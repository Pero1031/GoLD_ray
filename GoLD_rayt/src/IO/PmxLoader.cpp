// The PMX reader (IO/PmxLoader.hpp): a bounds-checked little-endian reader over the file in memory, and the sections of the format.
#include "pch.h"
#include "IO/PmxLoader.hpp"

#include <cstring>
#include <fstream>
#include <iterator>

namespace rayt::io {

    namespace {

        void appendUtf8(std::string& out, uint32_t cp) {
            if (cp < 0x80) out += static_cast<char>(cp);
            else if (cp < 0x800) { out += static_cast<char>(0xC0 | (cp >> 6)); out += static_cast<char>(0x80 | (cp & 0x3F)); }
            else if (cp < 0x10000) { out += static_cast<char>(0xE0 | (cp >> 12)); out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F)); out += static_cast<char>(0x80 | (cp & 0x3F)); }
            else { out += static_cast<char>(0xF0 | (cp >> 18)); out += static_cast<char>(0x80 | ((cp >> 12) & 0x3F)); out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F)); out += static_cast<char>(0x80 | (cp & 0x3F)); }
        }

        /// Bounds-checked little-endian reader over an in-memory file.
        class PmxReader {
        public:
            PmxReader(const std::vector<unsigned char>& bytes) : m_data(bytes) {}

            bool ok() const { return m_ok; }
            size_t position() const { return m_pos; }

            template <class T>
            T read() {
                T v{};
                if (m_pos + sizeof(T) > m_data.size()) { m_ok = false; m_pos = m_data.size(); return v; }
                std::memcpy(&v, m_data.data() + m_pos, sizeof(T));
                m_pos += sizeof(T);
                return v;
            }

            void skip(size_t n) {
                if (m_pos + n > m_data.size()) { m_ok = false; m_pos = m_data.size(); return; }
                m_pos += n;
            }

            Point2 vec2() { const float x = read<float>(), y = read<float>(); return Point2(x, y); }
            Vector3 vec3() { const float x = read<float>(), y = read<float>(), z = read<float>(); return Vector3(x, y, z); }
            Vec4 vec4() { const float x = read<float>(), y = read<float>(), z = read<float>(), w = read<float>(); return Vec4(x, y, z, w); }

            /// A signed index of the given byte size (-1 means "none").
            int index(int size) {
                switch (size) {
                case 1: return read<int8_t>();
                case 2: return read<int16_t>();
                default: return read<int32_t>();
                }
            }

            /// A vertex index: unsigned for 1 and 2 bytes (so 65535 vertices fit in 2 bytes).
            int vertexIndex(int size) {
                switch (size) {
                case 1: return read<uint8_t>();
                case 2: return read<uint16_t>();
                default: return read<int32_t>();
                }
            }

            /// int32 byte length + text in the model's encoding, returned as UTF-8.
            std::string text(int encoding) {
                const int32_t len = read<int32_t>();
                if (len < 0 || m_pos + static_cast<size_t>(len) > m_data.size()) { m_ok = false; return {}; }
                std::string out;
                const unsigned char* p = m_data.data() + m_pos;
                if (encoding == 1) {
                    out.assign(reinterpret_cast<const char*>(p), static_cast<size_t>(len));
                }
                else {   // UTF-16LE with surrogate pairs
                    for (int i = 0; i + 1 < len; i += 2) {
                        uint32_t u = p[i] | (static_cast<uint32_t>(p[i + 1]) << 8);
                        if (u >= 0xD800 && u < 0xDC00 && i + 3 < len) {
                            const uint32_t lo = p[i + 2] | (static_cast<uint32_t>(p[i + 3]) << 8);
                            if (lo >= 0xDC00 && lo < 0xE000) { u = 0x10000 + ((u - 0xD800) << 10) + (lo - 0xDC00); i += 2; }
                        }
                        appendUtf8(out, u);
                    }
                }
                m_pos += static_cast<size_t>(len);
                return out;
            }

        private:
            const std::vector<unsigned char>& m_data;
            size_t m_pos = 0;
            bool m_ok = true;
        };

    } // namespace

    bool loadPmx(const std::filesystem::path& file, PmxModel& out, std::string* outError) {
        auto fail = [&](const std::string& msg) { if (outError) *outError = msg; return false; };

        std::ifstream f(file, std::ios::binary);
        if (!f) return fail("cannot open the PMX file");
        std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        if (bytes.size() < 16 || std::memcmp(bytes.data(), "PMX ", 4) != 0)
            return fail("not a PMX file (missing 'PMX ' signature; PMD is not supported)");

        PmxReader r(bytes);
        out = PmxModel();
        r.skip(4);
        out.version = r.read<float>();
        if (out.version < 2.0f || out.version > 2.15f) return fail("unsupported PMX version " + std::to_string(out.version));

        const int globalsCount = r.read<uint8_t>();
        if (globalsCount < 8) return fail("PMX header has fewer than 8 global bytes");
        int g[16] = {};
        for (int i = 0; i < globalsCount; ++i) { const int v = r.read<uint8_t>(); if (i < 16) g[i] = v; }
        const int encoding = g[0];
        out.additionalUv = g[1];
        const int vertexIndexSize = g[2], textureIndexSize = g[3], boneIndexSize = g[5];
        for (int s : { vertexIndexSize, textureIndexSize, g[4], boneIndexSize, g[6], g[7] })
            if (s != 1 && s != 2 && s != 4) return fail("invalid index size in the PMX header");

        out.name = r.text(encoding);
        out.nameEn = r.text(encoding);
        out.comment = r.text(encoding);
        out.commentEn = r.text(encoding);

        // ---- vertices ----
        const int32_t vertexCount = r.read<int32_t>();
        if (!r.ok() || vertexCount < 0 || vertexCount > 20000000) return fail("implausible vertex count");
        out.vertices.resize(static_cast<size_t>(vertexCount));
        for (PmxVertex& v : out.vertices) {
            v.position = Point3(r.vec3());
            v.normal = r.vec3();
            v.uv = r.vec2();
            for (int i = 0; i < out.additionalUv; ++i) r.vec4();
            v.weightType = r.read<uint8_t>();
            switch (v.weightType) {
            case 0: v.bones[0] = r.index(boneIndexSize); v.weights[0] = 1; break;                // BDEF1
            case 1: {                                                                            // BDEF2
                v.bones[0] = r.index(boneIndexSize); v.bones[1] = r.index(boneIndexSize);
                const Real w = r.read<float>();
                v.weights[0] = w; v.weights[1] = Real(1) - w;
                break;
            }
            case 2: case 4:                                                                      // BDEF4 / QDEF
                for (int i = 0; i < 4; ++i) v.bones[i] = r.index(boneIndexSize);
                for (int i = 0; i < 4; ++i) v.weights[i] = r.read<float>();
                break;
            case 3: {                                                                            // SDEF (treated as BDEF2)
                v.bones[0] = r.index(boneIndexSize); v.bones[1] = r.index(boneIndexSize);
                const Real w = r.read<float>();
                v.weights[0] = w; v.weights[1] = Real(1) - w;
                r.skip(9 * sizeof(float));                                                       // C, R0, R1
                break;
            }
            default: return fail("unknown skinning type " + std::to_string(v.weightType));
            }
            v.edgeScale = r.read<float>();
            if (!r.ok()) return fail("file ended inside the vertex list");
        }

        // ---- faces ----
        const int32_t indexCount = r.read<int32_t>();
        if (!r.ok() || indexCount < 0 || indexCount % 3 != 0 || indexCount > 60000000) return fail("implausible face index count");
        out.indices.resize(static_cast<size_t>(indexCount));
        for (int& idx : out.indices) {
            idx = r.vertexIndex(vertexIndexSize);
            if (idx < 0 || idx >= vertexCount) return fail("a face refers to a vertex that does not exist");
        }

        // ---- textures ----
        const int32_t textureCount = r.read<int32_t>();
        if (!r.ok() || textureCount < 0 || textureCount > 100000) return fail("implausible texture count");
        for (int i = 0; i < textureCount; ++i) out.textures.push_back(r.text(encoding));

        // ---- materials ----
        const int32_t materialCount = r.read<int32_t>();
        if (!r.ok() || materialCount < 0 || materialCount > 100000) return fail("implausible material count");
        int offset = 0;
        for (int i = 0; i < materialCount; ++i) {
            PmxMaterial m;
            m.name = r.text(encoding);
            m.nameEn = r.text(encoding);
            m.diffuse = r.vec4();
            m.specular = r.vec3();
            m.shininess = r.read<float>();
            m.ambient = r.vec3();
            m.flags = r.read<uint8_t>();
            m.edgeColor = r.vec4();
            m.edgeSize = r.read<float>();
            m.texture = r.index(textureIndexSize);
            m.sphereTexture = r.index(textureIndexSize);
            m.sphereMode = r.read<uint8_t>();
            m.sharedToon = r.read<uint8_t>() != 0;
            m.toon = m.sharedToon ? static_cast<int>(r.read<uint8_t>()) : r.index(textureIndexSize);
            m.memo = r.text(encoding);
            m.indexCount = r.read<int32_t>();
            m.firstIndex = offset;
            offset += m.indexCount;
            if (!r.ok()) return fail("file ended inside the material list");
            if (m.texture >= static_cast<int>(out.textures.size())) m.texture = -1;
            if (m.sphereTexture >= static_cast<int>(out.textures.size())) m.sphereTexture = -1;
            out.materials.push_back(std::move(m));
        }
        if (offset != indexCount) return fail("the materials do not cover exactly the face list (" + std::to_string(offset) + " vs " + std::to_string(indexCount) + " indices)");

        // ---- bones (the skeleton the vertices are skinned to) ----
        const int32_t boneCount = r.read<int32_t>();
        if (!r.ok() || boneCount < 0 || boneCount > 100000) return fail("implausible bone count");
        for (int i = 0; i < boneCount; ++i) {
            PmxBone b;
            b.name = r.text(encoding);
            b.nameEn = r.text(encoding);
            b.position = Point3(r.vec3());
            b.parent = r.index(boneIndexSize);
            b.layer = r.read<int32_t>();
            b.flags = r.read<uint16_t>();
            if (b.flags & 0x1u) b.tailBone = r.index(boneIndexSize); else b.tailOffset = r.vec3();
            if (b.flags & 0x300u) { b.grantParent = r.index(boneIndexSize); b.grantWeight = r.read<float>(); }
            if (b.flags & 0x400u) b.fixedAxis = r.vec3();
            if (b.flags & 0x800u) { b.localAxisX = r.vec3(); b.localAxisZ = r.vec3(); }
            if (b.flags & 0x2000u) r.read<int32_t>();   // external parent key
            if (b.flags & 0x20u) {
                b.ikTarget = r.index(boneIndexSize);
                b.ikLoops = r.read<int32_t>();
                b.ikLimit = r.read<float>();
                const int32_t links = r.read<int32_t>();
                if (links < 0 || links > 1000) return fail("implausible IK link count");
                for (int k = 0; k < links; ++k) {
                    PmxIkLink l;
                    l.bone = r.index(boneIndexSize);
                    l.hasLimit = r.read<uint8_t>() != 0;
                    if (l.hasLimit) { l.lower = r.vec3(); l.upper = r.vec3(); }
                    b.ikLinks.push_back(l);
                }
            }
            if (!r.ok()) return fail("file ended inside the bone list");
            out.bones.push_back(std::move(b));
        }
        for (const PmxVertex& v : out.vertices)
            for (int k = 0; k < 4; ++k)
                if (v.bones[k] >= boneCount) return fail("a vertex is skinned to a bone that does not exist");

        // ---- morphs ----
        const int materialIndexSize = g[4], morphIndexSize = g[6], rigidIndexSize = g[7];
        const int32_t morphCount = r.read<int32_t>();
        if (!r.ok() || morphCount < 0 || morphCount > 100000) return fail("implausible morph count");
        for (int i = 0; i < morphCount; ++i) {
            PmxMorph m;
            m.name = r.text(encoding);
            m.nameEn = r.text(encoding);
            m.panel = r.read<uint8_t>();
            m.type = r.read<uint8_t>();
            const int32_t n = r.read<int32_t>();
            if (!r.ok() || n < 0 || n > 20000000) return fail("implausible morph offset count");
            for (int k = 0; k < n; ++k) {
                switch (m.type) {
                case 0: case 9: {   // group / flip
                    PmxGroupOffset o; o.morph = r.index(morphIndexSize); o.rate = r.read<float>(); m.groupOffsets.push_back(o); break;
                }
                case 1: { PmxVertexOffset o; o.vertex = r.vertexIndex(vertexIndexSize); o.offset = r.vec3(); m.vertexOffsets.push_back(o); break; }
                case 2: { PmxBoneOffset o; o.bone = r.index(boneIndexSize); o.translation = r.vec3(); o.rotation = r.vec4(); m.boneOffsets.push_back(o); break; }
                case 3: case 4: case 5: case 6: case 7: {
                    PmxUvOffset o; o.vertex = r.vertexIndex(vertexIndexSize); o.offset = r.vec4(); m.uvOffsets.push_back(o); break;
                }
                case 8: {
                    PmxMaterialOffset o;
                    o.material = r.index(materialIndexSize);
                    o.op = r.read<uint8_t>();
                    o.diffuse = r.vec4(); o.specular = r.vec3(); o.shininess = r.read<float>(); o.ambient = r.vec3();
                    o.edgeColor = r.vec4(); o.edgeSize = r.read<float>();
                    o.textureTint = r.vec4(); o.sphereTint = r.vec4(); o.toonTint = r.vec4();
                    m.materialOffsets.push_back(o);
                    break;
                }
                case 10:            // impulse (rigid body, PMX 2.1): index, flag, velocity, torque
                    r.index(rigidIndexSize); r.read<uint8_t>(); r.skip(6 * sizeof(float));
                    break;
                default: return fail("unknown morph type " + std::to_string(m.type));
                }
                if (!r.ok()) return fail("file ended inside a morph");
            }
            out.morphs.push_back(std::move(m));
        }
        // display frames, rigid bodies and joints follow; not needed for a pose / expression.
        return true;
    }

} // namespace rayt::io
