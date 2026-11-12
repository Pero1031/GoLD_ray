/**
 * @file IO/PlyLoader.hpp
 * @brief Stanford PLY reader (ascii, binary little endian, binary big endian) into a TriangleMesh.
 *
 * This is the format that the mesh extraction of 3D Gaussian Splatting (2DGS, GOF, SuGaR via MeshLab) and most scanners write.
 * - The "vertex" element gives the positions (x y z); normals (nx ny nz) and texture coordinates (s t, u v, texture_u texture_v
 *   or texture_s texture_t) are read when present. Any other property (colors, confidence, ...) is skipped: the renderer
 *   takes the color from the material (docs/CAPTURE.md).
 * - The "face" element gives the polygons (a list named vertex_indices or vertex_index), fan-triangulated. A per-face
 *   "texcoord" list (2 numbers per corner, as MeshLab writes a textured mesh) splits the corners by their uv.
 * - Other elements (edges, materials, ...) are skipped.
 * - Texture coordinates are flipped to the renderer's convention (v grows downward, Core/Texture.hpp), as loadObj does.
 * - "comment TextureFile name" (MeshLab) is returned through outTexture, so a scene can find the photographed texture.
 */

#pragma once

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "Geometry/TriangleMesh.hpp"
#include "IO/ObjLoader.hpp"

namespace rayt::io {

    namespace plydetail {

        enum class Type { Char, UChar, Short, UShort, Int, UInt, Float, Double, None };

        inline Type typeOf(const std::string& s) {
            if (s == "char" || s == "int8") return Type::Char;
            if (s == "uchar" || s == "uint8") return Type::UChar;
            if (s == "short" || s == "int16") return Type::Short;
            if (s == "ushort" || s == "uint16") return Type::UShort;
            if (s == "int" || s == "int32") return Type::Int;
            if (s == "uint" || s == "uint32") return Type::UInt;
            if (s == "float" || s == "float32") return Type::Float;
            if (s == "double" || s == "float64") return Type::Double;
            return Type::None;
        }

        inline size_t sizeOf(Type t) {
            switch (t) {
            case Type::Char: case Type::UChar: return 1;
            case Type::Short: case Type::UShort: return 2;
            case Type::Int: case Type::UInt: case Type::Float: return 4;
            case Type::Double: return 8;
            default: return 0;
            }
        }

        struct Property {
            std::string name;
            Type type = Type::None;        // the item type (of the list items for a list)
            bool list = false;
            Type countType = Type::None;
        };

        struct Element {
            std::string name;
            size_t count = 0;
            std::vector<Property> props;
        };

        /// Reads one scalar after another, in the file's encoding.
        class Reader {
        public:
            Reader(std::istream& in, bool ascii, bool bigEndian) : m_in(in), m_ascii(ascii), m_big(bigEndian) {}

            bool read(Type t, double& v) {
                if (m_ascii) {
                    m_in >> v;
                    return static_cast<bool>(m_in);
                }
                unsigned char buf[8];
                const size_t n = sizeOf(t);
                m_in.read(reinterpret_cast<char*>(buf), static_cast<std::streamsize>(n));
                if (!m_in) return false;
                if (m_big) std::reverse(buf, buf + n);
                switch (t) {
                case Type::Char: { std::int8_t x; std::memcpy(&x, buf, 1); v = x; break; }
                case Type::UChar: { std::uint8_t x; std::memcpy(&x, buf, 1); v = x; break; }
                case Type::Short: { std::int16_t x; std::memcpy(&x, buf, 2); v = x; break; }
                case Type::UShort: { std::uint16_t x; std::memcpy(&x, buf, 2); v = x; break; }
                case Type::Int: { std::int32_t x; std::memcpy(&x, buf, 4); v = x; break; }
                case Type::UInt: { std::uint32_t x; std::memcpy(&x, buf, 4); v = x; break; }
                case Type::Float: { float x; std::memcpy(&x, buf, 4); v = x; break; }
                case Type::Double: { double x; std::memcpy(&x, buf, 8); v = x; break; }
                default: return false;
                }
                return true;
            }

        private:
            std::istream& m_in;
            bool m_ascii, m_big;
        };

    } // namespace plydetail

    /**
     * @brief Loads a PLY file. Returns false with a message in outError on failure.
     * @param outTexture  the file named by "comment TextureFile ..." (MeshLab), if any; untouched otherwise.
     */
    inline bool loadPly(const std::string& filename, TriangleMesh& out, std::string* outError = nullptr, std::string* outTexture = nullptr) {
        using namespace plydetail;
        auto fail = [&](const std::string& what) {
            if (outError) *outError = what + " (" + filename + ")";
            return false;
        };

        std::ifstream f(filename, std::ios::binary);
        if (!f) return fail("cannot open");

        // --- header ---
        std::string line;
        if (!std::getline(f, line)) return fail("empty file");
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line != "ply") return fail("not a PLY file (the first line is not \"ply\")");

        bool ascii = false, big = false, haveFormat = false, ended = false;
        std::vector<Element> elements;
        while (std::getline(f, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            std::istringstream ss(line);
            std::string tag;
            ss >> tag;
            if (tag == "format") {
                std::string fmt;
                ss >> fmt;
                if (fmt == "ascii") ascii = true;
                else if (fmt == "binary_little_endian") { ascii = false; big = false; }
                else if (fmt == "binary_big_endian") { ascii = false; big = true; }
                else return fail("unknown format \"" + fmt + "\"");
                haveFormat = true;
            }
            else if (tag == "comment" || tag == "obj_info") {
                std::string key;
                ss >> key;
                if (tag == "comment" && key == "TextureFile" && outTexture) {
                    std::string rest;
                    std::getline(ss, rest);
                    const size_t start = rest.find_first_not_of(" \t");
                    if (start != std::string::npos) *outTexture = rest.substr(start);
                }
            }
            else if (tag == "element") {
                Element e;
                ss >> e.name >> e.count;
                if (!ss) return fail("bad element line \"" + line + "\"");
                elements.push_back(e);
            }
            else if (tag == "property") {
                if (elements.empty()) return fail("a property before any element");
                Property p;
                std::string t;
                ss >> t;
                if (t == "list") {
                    std::string ct, it;
                    ss >> ct >> it >> p.name;
                    p.list = true;
                    p.countType = typeOf(ct);
                    p.type = typeOf(it);
                    if (p.countType == Type::None || p.type == Type::None) return fail("unknown list types in \"" + line + "\"");
                }
                else {
                    ss >> p.name;
                    p.type = typeOf(t);
                    if (p.type == Type::None) return fail("unknown property type \"" + t + "\"");
                }
                if (!ss && p.name.empty()) return fail("bad property line \"" + line + "\"");
                elements.back().props.push_back(p);
            }
            else if (tag == "end_header") { ended = true; break; }
            else if (tag.empty()) continue;
            else return fail("unexpected header line \"" + line + "\"");
        }
        if (!ended) return fail("no end_header");
        if (!haveFormat) return fail("no format line");

        // --- body ---
        Reader reader(f, ascii, big);
        std::vector<Point3> pos;
        std::vector<Vector3> nor;
        std::vector<Point2> tex;
        std::vector<std::vector<int>> faces;         // polygons, vertex indices
        std::vector<std::vector<Point2>> faceUv;     // per-corner uv of each polygon (MeshLab's "texcoord"), empty if none
        bool haveNormals = false, haveUv = false, haveFaceUv = false;

        for (const Element& e : elements) {
            if (e.name == "vertex") {
                int ix = -1, iy = -1, iz = -1, inx = -1, iny = -1, inz = -1, iu = -1, iv = -1;
                for (size_t k = 0; k < e.props.size(); ++k) {
                    const std::string& n = e.props[k].name;
                    const int i = static_cast<int>(k);
                    if (n == "x") ix = i; else if (n == "y") iy = i; else if (n == "z") iz = i;
                    else if (n == "nx") inx = i; else if (n == "ny") iny = i; else if (n == "nz") inz = i;
                    else if (n == "s" || n == "u" || n == "texture_u" || n == "texture_s") iu = i;
                    else if (n == "t" || n == "v" || n == "texture_v" || n == "texture_t") iv = i;
                }
                if (ix < 0 || iy < 0 || iz < 0) return fail("the vertex element has no x y z");
                haveNormals = inx >= 0 && iny >= 0 && inz >= 0;
                haveUv = iu >= 0 && iv >= 0;
                pos.reserve(e.count);
                if (haveNormals) nor.reserve(e.count);
                if (haveUv) tex.reserve(e.count);
                std::vector<double> vals(e.props.size(), 0.0);
                for (size_t i = 0; i < e.count; ++i) {
                    for (size_t k = 0; k < e.props.size(); ++k) {
                        const Property& p = e.props[k];
                        if (p.list) {   // a list in a vertex: skipped
                            double cnt = 0;
                            if (!reader.read(p.countType, cnt)) return fail("truncated vertex data");
                            double dummy;
                            for (int j = 0; j < static_cast<int>(cnt); ++j) if (!reader.read(p.type, dummy)) return fail("truncated vertex data");
                        }
                        else if (!reader.read(p.type, vals[k])) return fail("truncated vertex data");
                    }
                    pos.emplace_back(vals[static_cast<size_t>(ix)], vals[static_cast<size_t>(iy)], vals[static_cast<size_t>(iz)]);
                    if (haveNormals) nor.emplace_back(vals[static_cast<size_t>(inx)], vals[static_cast<size_t>(iny)], vals[static_cast<size_t>(inz)]);
                    if (haveUv) tex.emplace_back(vals[static_cast<size_t>(iu)], Real(1) - vals[static_cast<size_t>(iv)]);   // v grows downward in the renderer
                }
            }
            else if (e.name == "face") {
                faces.reserve(e.count);
                std::vector<int> poly;
                std::vector<Point2> uvs;
                for (size_t i = 0; i < e.count; ++i) {
                    poly.clear();
                    uvs.clear();
                    for (const Property& p : e.props) {
                        if (p.list) {
                            double cnt = 0;
                            if (!reader.read(p.countType, cnt)) return fail("truncated face data");
                            const int n = static_cast<int>(cnt);
                            if (p.name == "vertex_indices" || p.name == "vertex_index") {
                                for (int j = 0; j < n; ++j) {
                                    double v;
                                    if (!reader.read(p.type, v)) return fail("truncated face data");
                                    poly.push_back(static_cast<int>(v));
                                }
                            }
                            else if (p.name == "texcoord") {
                                for (int j = 0; j + 1 < n; j += 2) {
                                    double u, v;
                                    if (!reader.read(p.type, u) || !reader.read(p.type, v)) return fail("truncated face data");
                                    uvs.emplace_back(u, Real(1) - v);
                                }
                                if (n % 2 == 1) { double dummy; if (!reader.read(p.type, dummy)) return fail("truncated face data"); }
                                haveFaceUv = true;
                            }
                            else {
                                double dummy;
                                for (int j = 0; j < n; ++j) if (!reader.read(p.type, dummy)) return fail("truncated face data");
                            }
                        }
                        else {
                            double dummy;
                            if (!reader.read(p.type, dummy)) return fail("truncated face data");
                        }
                    }
                    if (poly.size() < 3) return fail("a face with fewer than 3 vertices");
                    faces.push_back(poly);
                    faceUv.push_back(uvs);
                }
            }
            else {   // another element: read past it
                for (size_t i = 0; i < e.count; ++i) {
                    for (const Property& p : e.props) {
                        double v;
                        if (p.list) {
                            if (!reader.read(p.countType, v)) return fail("truncated data in element " + e.name);
                            const int n = static_cast<int>(v);
                            for (int j = 0; j < n; ++j) if (!reader.read(p.type, v)) return fail("truncated data in element " + e.name);
                        }
                        else if (!reader.read(p.type, v)) return fail("truncated data in element " + e.name);
                    }
                }
            }
        }
        if (pos.empty()) return fail("no vertices");
        if (faces.empty()) return fail("no faces");

        // --- build the mesh ---
        out = TriangleMesh();
        const int nv = static_cast<int>(pos.size());
        auto emitFan = [&](const std::vector<int>& corners) {
            for (size_t k = 1; k + 1 < corners.size(); ++k) {
                out.indices.push_back(corners[0]);
                out.indices.push_back(corners[k]);
                out.indices.push_back(corners[k + 1]);
            }
        };

        if (haveFaceUv) {   // split the corners by (vertex, uv), as an OBJ with v/vt does
            std::map<std::tuple<int, float, float>, int> remap;
            std::vector<int> corners;
            for (size_t i = 0; i < faces.size(); ++i) {
                const std::vector<int>& poly = faces[i];
                const std::vector<Point2>& uvs = faceUv[i];
                corners.clear();
                for (size_t k = 0; k < poly.size(); ++k) {
                    const int vi = poly[k];
                    if (vi < 0 || vi >= nv) return fail("bad vertex index");
                    const Point2 uv = k < uvs.size() ? uvs[k] : Point2(0, 0);
                    const auto key = std::make_tuple(vi, static_cast<float>(uv.x), static_cast<float>(uv.y));
                    auto it = remap.find(key);
                    if (it == remap.end()) {
                        it = remap.emplace(key, static_cast<int>(out.p.size())).first;
                        out.p.push_back(pos[static_cast<size_t>(vi)]);
                        if (haveNormals) out.n.push_back(glm::normalize(nor[static_cast<size_t>(vi)]));
                        out.uv.push_back(uv);
                    }
                    corners.push_back(it->second);
                }
                emitFan(corners);
            }
        }
        else {
            out.p = std::move(pos);
            if (haveNormals) {
                out.n = std::move(nor);
                for (auto& n : out.n) n = glm::normalize(n);
            }
            if (haveUv) out.uv = std::move(tex);
            for (const std::vector<int>& poly : faces) {
                for (int vi : poly) if (vi < 0 || vi >= nv) return fail("bad vertex index");
                emitFan(poly);
            }
        }
        if (out.hasNormals()) {   // a zero normal (a degenerate vertex) must not poison the shading
            for (auto& n : out.n) if (!(glm::dot(n, n) > Real(0.5))) n = Vector3(0, 1, 0);
        }
        return true;
    }

    /// Loads an OBJ or a PLY by the file's extension (docs/SCENE_FILE.md, "mesh").
    inline bool loadMesh(const std::string& filename, TriangleMesh& out, std::string* outError = nullptr, std::string* outTexture = nullptr) {
        std::string ext;
        const size_t dot = filename.find_last_of('.');
        if (dot != std::string::npos) ext = filename.substr(dot + 1);
        std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (ext == "ply") return loadPly(filename, out, outError, outTexture);
        return loadObj(filename, out, outError);
    }

} // namespace rayt::io
