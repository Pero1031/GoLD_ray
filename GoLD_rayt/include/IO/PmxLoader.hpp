/**
 * @file IO/PmxLoader.hpp
 * @brief Reader for MikuMikuDance PMX models (version 2.0 / 2.1): geometry, textures, materials.
 *
 * PMX is a little-endian binary format:
 *
 *   header      "PMX ", float version, byte count + that many "global" bytes:
 *                 [0] text encoding (0 = UTF-16LE, 1 = UTF-8)
 *                 [1] number of additional UV vec4 per vertex
 *                 [2] vertex index size   [3] texture index size   [4] material index size
 *                 [5] bone index size     [6] morph index size     [7] rigid body index size
 *               (index sizes are 1, 2 or 4 bytes; a smaller model uses smaller indices)
 *               then four texts: model name / English name / comment / English comment.
 *   text        int32 byte length + bytes in the declared encoding.
 *   vertices    int32 count; each: position, normal, uv, additional UVs, skinning
 *               (BDEF1 / BDEF2 / BDEF4 / SDEF / QDEF with bone indices and weights), edge scale.
 *   faces       int32 count of INDICES (3 per triangle), each a vertex index (unsigned for 1/2 bytes).
 *   textures    int32 count; each a text (relative path, usually with '\' separators).
 *   materials   int32 count; each: names, diffuse RGBA, specular RGB + shininess, ambient RGB,
 *               flags (1 double-sided, 2 ground shadow, 4 shadow-map caster, 8 shadow receiver, 16 edge),
 *               edge color + size, texture index, sphere-map texture index + mode, toon (shared 0..9
 *               or own texture), memo, and the number of indices (faces x 3) the material covers,
 *               consecutively in the face list.
 *   bones      int32 count; each: names, rest position, parent, deform level, flags (rotatable, translatable,
 *               IK, rotation/translation grant, fixed axis, local axis ...), tail, grant parent + weight, IK data.
 *               Vertices are skinned to these bones (indices + weights, stored in PmxVertex).
 *   then morphs, display frames, rigid bodies, joints (and soft bodies in 2.1): expressions / physics, not
 *   needed for a pose, so this loader stops after the bones.
 *
 * Coordinates are left-handed (DirectX style), +Y up; 1 unit is about 8 cm. Front faces wind
 * clockwise. Converting to the renderer's right-handed system is done by the caller.
 */

#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "Core/Types.hpp"

namespace rayt::io {

    using Vec4 = glm::vec<4, Real, glm::defaultp>;

    struct PmxVertex {
        Point3 position = Point3(0.0);
        Vector3 normal = Vector3(0, 1, 0);
        Point2 uv = Point2(0.0);
        int weightType = 0;        ///< 0 BDEF1, 1 BDEF2, 2 BDEF4, 3 SDEF, 4 QDEF
        Real edgeScale = 1;
        int bones[4] = { -1, -1, -1, -1 };        ///< skinning bones (-1 = unused)
        Real weights[4] = { 0, 0, 0, 0 };         ///< their weights (sum to 1). SDEF is treated as BDEF2, QDEF as BDEF4.
    };

    /// One link of an IK chain.
    struct PmxIkLink {
        int bone = -1;
        bool hasLimit = false;
        Vector3 lower = Vector3(0.0), upper = Vector3(0.0);   ///< angle limits (radians), only if hasLimit
    };

    struct PmxBone {
        std::string name, nameEn;
        Point3 position = Point3(0.0);          ///< rest position (model space)
        int parent = -1;
        int layer = 0;                          ///< deform order level
        unsigned flags = 0;
        int tailBone = -1;                      ///< connected bone (flag 0x1) or -1
        Vector3 tailOffset = Vector3(0.0);
        int grantParent = -1;                   ///< bone whose rotation / translation is inherited (flags 0x100 / 0x200)
        Real grantWeight = 0;
        Vector3 fixedAxis = Vector3(0.0);       ///< flag 0x400
        Vector3 localAxisX = Vector3(1, 0, 0), localAxisZ = Vector3(0, 0, 1);   ///< flag 0x800
        int ikTarget = -1;                      ///< IK bones (flag 0x20)
        int ikLoops = 0;
        Real ikLimit = 0;
        std::vector<PmxIkLink> ikLinks;

        bool rotatable() const { return (flags & 0x2u) != 0; }
        bool translatable() const { return (flags & 0x4u) != 0; }
        bool isIk() const { return (flags & 0x20u) != 0; }
        bool inheritsRotation() const { return (flags & 0x100u) != 0; }
        bool inheritsTranslation() const { return (flags & 0x200u) != 0; }
        bool afterPhysics() const { return (flags & 0x1000u) != 0; }
    };

    // ---- morphs (facial expressions, material switches, ...) ----
    struct PmxVertexOffset { int vertex = 0; Vector3 offset = Vector3(0.0); };
    struct PmxUvOffset { int vertex = 0; Vec4 offset = Vec4(0.0); };
    struct PmxBoneOffset { int bone = -1; Vector3 translation = Vector3(0.0); Vec4 rotation = Vec4(0, 0, 0, 1); };   ///< rotation = quaternion x, y, z, w
    struct PmxGroupOffset { int morph = -1; Real rate = 0; };
    struct PmxMaterialOffset {
        int material = -1;                      ///< -1 = every material
        int op = 0;                             ///< 0 multiply, 1 add
        Vec4 diffuse = Vec4(0.0);
        Vector3 specular = Vector3(0.0);
        Real shininess = 0;
        Vector3 ambient = Vector3(0.0);
        Vec4 edgeColor = Vec4(0.0);
        Real edgeSize = 0;
        Vec4 textureTint = Vec4(0.0), sphereTint = Vec4(0.0), toonTint = Vec4(0.0);
    };

    struct PmxMorph {
        std::string name, nameEn;
        int panel = 0;                          ///< 1 eyebrows, 2 eyes, 3 mouth, 4 other
        int type = 0;                           ///< 0 group, 1 vertex, 2 bone, 3-7 UV, 8 material, 9 flip, 10 impulse
        std::vector<PmxVertexOffset> vertexOffsets;
        std::vector<PmxUvOffset> uvOffsets;
        std::vector<PmxBoneOffset> boneOffsets;
        std::vector<PmxGroupOffset> groupOffsets;   ///< group (type 0) and flip (type 9)
        std::vector<PmxMaterialOffset> materialOffsets;
        size_t offsetCount() const { return vertexOffsets.size() + uvOffsets.size() + boneOffsets.size() + groupOffsets.size() + materialOffsets.size(); }
    };

    struct PmxMaterial {
        std::string name;                       ///< UTF-8
        std::string nameEn;
        Vec4 diffuse = Vec4(1.0);         ///< RGBA, display (gamma) space
        Vector3 specular = Vector3(0.0);
        Real shininess = 0;
        Vector3 ambient = Vector3(0.0);
        unsigned flags = 0;
        Vec4 edgeColor = Vec4(0, 0, 0, 1);
        Real edgeSize = 1;
        int texture = -1;                       ///< index into PmxModel::textures, -1 = none
        int sphereTexture = -1;
        int sphereMode = 0;                     ///< 0 off, 1 multiply, 2 add, 3 sub-texture
        bool sharedToon = true;
        int toon = -1;                          ///< shared toon number 0..9, or texture index
        std::string memo;
        bool alphaFromMorph = false;            ///< set when a morph changed the diffuse alpha (see applyMorphs)
        int indexCount = 0;                     ///< number of face indices of this material
        int firstIndex = 0;                     ///< offset into PmxModel::indices

        bool doubleSided() const { return (flags & 1u) != 0; }
        int triangleCount() const { return indexCount / 3; }
    };

    struct PmxModel {
        float version = 0;
        std::string name, nameEn, comment, commentEn;
        std::vector<PmxVertex> vertices;
        std::vector<int> indices;               ///< 3 per triangle, in file order (clockwise from the front)
        std::vector<std::string> textures;      ///< UTF-8 relative paths as stored in the file
        std::vector<PmxMaterial> materials;
        std::vector<PmxBone> bones;
        std::vector<PmxMorph> morphs;
        int additionalUv = 0;
    };

    namespace detail {

        inline void appendUtf8(std::string& out, uint32_t cp) {
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

    } // namespace detail

    /// Loads geometry, textures and materials of a PMX 2.0 / 2.1 file.
    inline bool loadPmx(const std::filesystem::path& file, PmxModel& out, std::string* outError = nullptr) {
        auto fail = [&](const std::string& msg) { if (outError) *outError = msg; return false; };

        std::ifstream f(file, std::ios::binary);
        if (!f) return fail("cannot open the PMX file");
        std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        if (bytes.size() < 16 || std::memcmp(bytes.data(), "PMX ", 4) != 0)
            return fail("not a PMX file (missing 'PMX ' signature; PMD is not supported)");

        detail::PmxReader r(bytes);
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
