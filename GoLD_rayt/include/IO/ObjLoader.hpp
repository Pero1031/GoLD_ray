/**
 * @file IO/ObjLoader.hpp
 * @brief Minimal Wavefront OBJ reader (v, vt, vn, f) into a TriangleMesh.
 *
 * - Polygons are triangulated as a fan; "v", "v/vt", "v//vn" and "v/vt/vn" face vertices and
 *   negative (relative) indices are supported. Groups, materials and smoothing groups are ignored.
 * - Vertices are merged by their (position, uv, normal) index triple.
 * - If any face vertex has no normal, the whole mesh is flat-shaded (normals dropped).
 */

#pragma once

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "Geometry/TriangleMesh.hpp"

namespace rayt::io {

    inline bool loadObj(const std::string& filename, TriangleMesh& out, std::string* outError = nullptr) {
        std::ifstream f(filename);
        if (!f) {
            if (outError) *outError = "cannot open " + filename;
            return false;
        }

        std::vector<Point3> pos;
        std::vector<Vector3> nor;
        std::vector<Point2> tex;
        std::map<std::tuple<int, int, int>, int> remap;
        bool missingNormal = false;

        out = TriangleMesh();
        std::vector<Vector3> outN;
        std::vector<Point2> outUV;

        auto resolve = [](int idx, size_t count) { return idx > 0 ? idx - 1 : static_cast<int>(count) + idx; };

        std::string line;
        while (std::getline(f, line)) {
            std::istringstream ss(line);
            std::string tag;
            ss >> tag;
            if (tag == "v") { double x, y, z; if (ss >> x >> y >> z) pos.emplace_back(x, y, z); }
            else if (tag == "vn") { double x, y, z; if (ss >> x >> y >> z) nor.emplace_back(x, y, z); }
            else if (tag == "vt") { double u, v; if (ss >> u >> v) tex.emplace_back(u, v); }
            else if (tag == "f") {
                std::vector<int> face;
                std::string tok;
                while (ss >> tok) {
                    int vi = 0, ti = 0, ni = 0;
                    // "a", "a/b", "a//c", "a/b/c"
                    const size_t s1 = tok.find('/');
                    vi = std::stoi(tok.substr(0, s1));
                    if (s1 != std::string::npos) {
                        const size_t s2 = tok.find('/', s1 + 1);
                        const std::string mid = tok.substr(s1 + 1, s2 == std::string::npos ? std::string::npos : s2 - s1 - 1);
                        if (!mid.empty()) ti = std::stoi(mid);
                        if (s2 != std::string::npos && s2 + 1 < tok.size()) ni = std::stoi(tok.substr(s2 + 1));
                    }
                    const int rv = resolve(vi, pos.size());
                    const int rt = ti ? resolve(ti, tex.size()) : -1;
                    const int rn = ni ? resolve(ni, nor.size()) : -1;
                    if (rv < 0 || rv >= static_cast<int>(pos.size())) {
                        if (outError) *outError = "bad vertex index in " + filename;
                        return false;
                    }
                    if (rn < 0 || rn >= static_cast<int>(nor.size())) missingNormal = true;

                    const auto key = std::make_tuple(rv, rt, rn);
                    auto it = remap.find(key);
                    if (it == remap.end()) {
                        it = remap.emplace(key, static_cast<int>(out.p.size())).first;
                        out.p.push_back(pos[static_cast<size_t>(rv)]);
                        outN.push_back((rn >= 0 && rn < static_cast<int>(nor.size())) ? glm::normalize(nor[static_cast<size_t>(rn)]) : Vector3(0, 0, 0));
                        outUV.push_back((rt >= 0 && rt < static_cast<int>(tex.size())) ? tex[static_cast<size_t>(rt)] : Point2(0, 0));
                    }
                    face.push_back(it->second);
                }
                for (size_t k = 1; k + 1 < face.size(); ++k) {   // fan triangulation
                    out.indices.push_back(face[0]);
                    out.indices.push_back(face[k]);
                    out.indices.push_back(face[k + 1]);
                }
            }
        }

        if (out.indices.empty()) {
            if (outError) *outError = "no faces in " + filename;
            return false;
        }
        if (!missingNormal) out.n = std::move(outN);
        if (!tex.empty()) out.uv = std::move(outUV);
        return true;
    }

} // namespace rayt::io
