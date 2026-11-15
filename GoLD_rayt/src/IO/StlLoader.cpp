// The STL reader (IO/StlLoader.hpp): binary and ascii, vertex welding, smooth normals from the winding.
#include "pch.h"
#include "IO/StlLoader.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace rayt::io {

    namespace {

        /// A vertex position as three float bit patterns (so the weld is exact, and -0 equals +0).
        using PosKey = std::array<std::uint32_t, 3>;

        struct PosKeyHash {
            size_t operator()(const PosKey& k) const {
                std::uint64_t h = 1469598103934665603ull;
                for (std::uint32_t w : k) { h ^= w; h *= 1099511628211ull; }
                return static_cast<size_t>(h ^ (h >> 29));
            }
        };

        PosKey keyOf(const float xyz[3]) {
            PosKey k;
            for (int i = 0; i < 3; ++i) {
                float v = xyz[i];
                if (v == 0.0f) v = 0.0f;   // -0 -> +0
                std::memcpy(&k[static_cast<size_t>(i)], &v, 4);
            }
            return k;
        }

        /// Collects the triangle soup, welds the corners and builds the smooth normals.
        struct Welder {
            std::unordered_map<PosKey, int, PosKeyHash> index;
            TriangleMesh* mesh;

            explicit Welder(TriangleMesh* m) : mesh(m) {}

            int corner(const float xyz[3]) {
                const auto it = index.find(keyOf(xyz));
                if (it != index.end()) return it->second;
                const int id = static_cast<int>(mesh->p.size());
                mesh->p.emplace_back(xyz[0], xyz[1], xyz[2]);
                index.emplace(keyOf(xyz), id);
                return id;
            }

            void triangle(const float v[9]) {
                const int a = corner(v), b = corner(v + 3), c = corner(v + 6);
                if (a == b || b == c || a == c) return;   // collapsed after welding
                mesh->indices.push_back(a);
                mesh->indices.push_back(b);
                mesh->indices.push_back(c);
            }

            /// Area-weighted vertex normals (the cross product's length is twice the area), from the counter-clockwise winding.
            void finish() {
                std::vector<Vector3> n(mesh->p.size(), Vector3(0.0));
                for (size_t i = 0; i < mesh->numTriangles(); ++i) {
                    const int a = mesh->indices[3 * i], b = mesh->indices[3 * i + 1], c = mesh->indices[3 * i + 2];
                    const Vector3 f = glm::cross(Vector3(mesh->p[static_cast<size_t>(b)] - mesh->p[static_cast<size_t>(a)]),
                                                 Vector3(mesh->p[static_cast<size_t>(c)] - mesh->p[static_cast<size_t>(a)]));
                    n[static_cast<size_t>(a)] += f; n[static_cast<size_t>(b)] += f; n[static_cast<size_t>(c)] += f;
                }
                for (Vector3& v : n) v = glm::dot(v, v) > 0.0 ? glm::normalize(v) : Vector3(0, 1, 0);   // an unused vertex gets a harmless normal
                mesh->n = std::move(n);
            }
        };

    } // namespace

    bool loadStl(const std::string& filename, TriangleMesh& out, std::string* outError) {
        auto fail = [&](const std::string& what) {
            if (outError) *outError = what + " (" + filename + ")";
            return false;
        };

        std::ifstream f(filename, std::ios::binary | std::ios::ate);
        if (!f) return fail("cannot open");
        const std::streamoff size = f.tellg();
        f.seekg(0);
        if (size < 15) return fail("too small to be an STL file");
        std::vector<char> data(static_cast<size_t>(size));
        if (!f.read(data.data(), size)) return fail("read error");

        out = TriangleMesh();
        Welder weld(&out);

        std::uint32_t count = 0;
        if (size >= 84) std::memcpy(&count, data.data() + 80, 4);   // little endian (all STL writers, and this PC)
        const bool binary = size >= 84 && static_cast<std::uint64_t>(size) == 84ull + 50ull * count;

        if (binary) {
            out.p.reserve(static_cast<size_t>(count) / 2 + 16);
            out.indices.reserve(static_cast<size_t>(count) * 3);
            const char* rec = data.data() + 84;
            for (std::uint32_t i = 0; i < count; ++i, rec += 50) {
                float v[9];
                std::memcpy(v, rec + 12, 36);   // skip the stored normal (12 bytes); the attribute (2 bytes) follows the vertices
                for (float x : v) if (!std::isfinite(x)) return fail("a vertex is not a finite number (triangle " + std::to_string(i) + ")");
                weld.triangle(v);
            }
        }
        else {
            if (std::strncmp(data.data(), "solid", 5) != 0) return fail("neither a binary STL (the size does not match the triangle count) nor an ascii one (no \"solid\")");
            std::istringstream in(std::string(data.begin(), data.end()));
            std::string word;
            float v[9];
            int k = 0;
            while (in >> word) {
                if (word != "vertex") continue;
                float x, y, z;
                if (!(in >> x >> y >> z)) return fail("a vertex line is incomplete");
                if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) return fail("a vertex is not a finite number");
                v[3 * k] = x; v[3 * k + 1] = y; v[3 * k + 2] = z;
                if (++k == 3) { weld.triangle(v); k = 0; }
            }
            if (k != 0) return fail("the last facet has fewer than 3 vertices");
        }

        if (out.indices.empty()) return fail("no triangles");
        weld.finish();
        return true;
    }

} // namespace rayt::io
