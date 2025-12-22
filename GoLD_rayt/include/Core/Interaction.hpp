#pragma once

#include "Core/Core.hpp"
#include "Core/Math.hpp"

namespace rayt {

    // Forward declaration to avoid circular dependency
    // note:declaration in Core.hpp
    class Material;

    // SurfaceInteraction stores all geometric information at the intersection point.
    struct SurfaceInteraction {
        Point3 p;          // Intersection point
        Vector3 n;         // Geometric normal
        Vector3 wo;        // Outgoing direction (vector towards camera)
        UV uv;             // Texture coordinates
        Real t;            // Ray parameter distance

        const Material* matPtr = nullptr; // Pointer to the material at this point

        // Differential Geometry (For Normal Mapping / Anisotropy)
        Vector3 dpdu;      // Tangent (u方向の接線)
        Vector3 dpdv;      // Bitangent (v方向の接線)
        Vector3 gn;        // Geometric Normal (本来の幾何形状の法線)

        // Methods
        
        // Initialize with basic data
        SurfaceInteraction() = default;

        // Helper to orient normal correctly
        // rayDir: Ray direction (Incident, pointing into surface)
        // geometricNormal: Raw surface normal
        void setFaceNormal(const Vector3& rayDir, const Vector3& geometricNormal) {
            // PBRT: "FaceForward" logic
            // レイと法線が逆向き（内積が負）なら表、同じ向きなら裏（内部からのヒット）
            bool frontFace = glm::dot(rayDir, geometricNormal) < 0;

            // シェーディング法線と幾何法線の両方を設定
            gn = frontFace ? geometricNormal : -geometricNormal;
            n = gn; // 初期状態ではシェーディング法線 = 幾何法線
        }

        // compute tangent を定義する必要
        // 接空間 for 金属反射や法線マップ（Normal Mapping）、異方性反射（Anisotropy）
    };

}