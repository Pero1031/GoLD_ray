#pragma once

#include "Core.hpp"

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

        // Helper to orient normal correctly (so it always points against the ray)
        void setFaceNormal(const Vector3& rayDir, const Vector3& geometricNormal) {
            bool frontFace = glm::dot(rayDir, geometricNormal) < 0;
            n = frontFace ? geometricNormal : -geometricNormal;
        }
    };

}