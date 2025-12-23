#pragma once

/**
* @file Interaction.hpp
* @brief Data structures for surface and volumetric interactions.
*/

#include "Core/Core.hpp"
#include "Core/Math.hpp"

namespace rayt {

    // Forward declaration to avoid circular dependency.
    // Note: Already declared in Core/Forward.hpp, but listed here for clarity.
    class Material;

    /**
     * @brief SurfaceInteraction stores all geometric and shading information at an intersection point.
     * It acts as a bridge between the geometry (Shapes) and the shading (Materials/BSDFs).
     */
    struct SurfaceInteraction {
        Point3 p;           // The hit point in world space.
        Vector3 n;          // Shading normal (may be modified by normal maps or smoothing).
        Vector3 wo;         // Outgoing direction (pointing away from the surface, toward the camera/ray origin).
        UV uv;              // 2D texture coordinates.
        Real t;             // Distance along the ray (parametric distance).

        const Material* matPtr = nullptr; // Pointer to the material property at hit point

        // ---------------------------------------------------------------------
        // Differential Geometry (For Normal Mapping / Anisotropy)
        // ---------------------------------------------------------------------
        // These vectors define the local tangent space (TBN frame).
        Vector3 dpdu;       // Tangent vector: partial derivative of position with respect to 'u'.
        Vector3 dpdv;       // Bitangent vector: partial derivative of position with respect to 'v'.
        Vector3 gn;         // Geometric normal: the true perpendicular vector of the underlying geometry.

        // ---------------------------------------------------------------------
        // Methods
        // ---------------------------------------------------------------------
        
        // Initialize with basic data
        SurfaceInteraction() = default;

        /**
         * @brief Orients the geometric and shading normals based on the incident ray.
         * Ensures that the normal always faces the side from which the ray arrived (FaceForward).
         * * @param rayDir The incident ray direction (pointing toward the surface).
         * @param geometricNormal The raw normal provided by the shape's geometry.
         */
        void setFaceNormal(const Vector3& rayDir, const Vector3& geometricNormal) {
            // Logic: If the dot product is negative, the ray and normal are in opposite directions (Front Face).
            bool frontFace = glm::dot(rayDir, geometricNormal) < 0;

            // Set both geometric and shading normals to face the ray.
            gn = frontFace ? geometricNormal : -geometricNormal;
            n = gn; // Initially, the shading normal matches the geometric normal.
        }

        // TODO: Implement coordinate system construction for tangent space.
        // Required for Normal Mapping, Anisotropy, and Microfacet distributions.
    };

}