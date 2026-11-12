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
 *
 * The reader itself is in src/IO/PmxLoader.cpp; this header holds the model's data structures.
 */

#pragma once

#include <cstdint>
#include <filesystem>
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

    /// Loads geometry, textures and materials of a PMX 2.0 / 2.1 file.
    bool loadPmx(const std::filesystem::path& file, PmxModel& out, std::string* outError = nullptr);

} // namespace rayt::io
