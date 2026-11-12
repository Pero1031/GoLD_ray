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
 * The reader is in src/IO/PlyLoader.cpp; this header holds the header parser's types and the scalar reader.
 */

#pragma once

#include <algorithm>
#include <istream>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "Geometry/TriangleMesh.hpp"

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

        /// Parses the header of an open (binary-mode) stream, leaving it at the first byte of the data. False with `err` set on failure.
        bool parseHeader(std::istream& f, bool& ascii, bool& big, std::vector<Element>& elements, std::string* outTexture, std::string& err);

    } // namespace plydetail

    /**
     * @brief Loads a PLY file. Returns false with a message in outError on failure.
     * @param outTexture  the file named by "comment TextureFile ..." (MeshLab), if any; untouched otherwise.
     */
    bool loadPly(const std::string& filename, TriangleMesh& out, std::string* outError = nullptr, std::string* outTexture = nullptr);

    /// Loads an OBJ or a PLY by the file's extension (docs/SCENE_FILE.md, "mesh").
    bool loadMesh(const std::string& filename, TriangleMesh& out, std::string* outError = nullptr, std::string* outTexture = nullptr);

} // namespace rayt::io
