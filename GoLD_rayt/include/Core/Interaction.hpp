/**
 * @file Core/Interaction.hpp
 * @brief Data structures for surface and volumetric interactions.
 *
 * This header defines interaction records used throughout the renderer.
 * SurfaceInteraction is the primary record for surface hits and is also
 * responsible for spawning numerically safe rays (shadow acne avoidance).
 */

#pragma once

#include <algorithm>

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Constants.hpp"
#include "Core/Math.hpp"
#include "Core/Ray.hpp"   

namespace rayt {

    // Forward declaration to avoid circular dependency.
    // Note: Already declared in Core/Forward.hpp, but listed here for clarity.
    class Material;

    /**
     * @brief Flip vector v so that it lies in the same hemisphere as ref.
     * Useful for keeping normals consistently oriented.
     */
    inline Vector3 FaceForward(const Vector3& v, const Vector3& ref) {
        return (glm::dot(v, ref) < 0) ? -v : v;
    }

    /**
     * @brief Minimal base interaction record (optional).
     * This can be expanded later (e.g., time/medium) or used as a base type
     * if you decide to model Surface/Medium interactions polymorphically.
     */
    struct Interaction {
        Point3 p{};
        Real t = 0;
        Vector3 wo{};
    };

    /**
     * @brief A simple shading frame (TBN). Currently unused but reserved for future use.
     */
    struct ShadingFrame { Vector3 n, t, b; };

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
        const Light* areaLight = nullptr; // Enter if the hit surface is an area light

        // ---------------------------------------------------------------------
        // Differential Geometry (For Normal Mapping / Anisotropy)
        // ---------------------------------------------------------------------
        
        // These vectors define the local tangent space (TBN frame).
        Vector3 dpdu;       // Tangent vector: partial derivative of position with respect to 'u'.
        Vector3 dpdv;       // Bitangent vector: partial derivative of position with respect to 'v'.
        Vector3 gn;         // Geometric normal: the true perpendicular vector of the underlying geometry.

        const Medium* medium = nullptr;   // The medium surrounding the interaction point.
        Real time = 0.0;                  // Ray time for motion blur / time-varying transforms.
        bool frontFace = true;            // True if the ray hit the "front" side of the surface.

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
            frontFace = glm::dot(rayDir, geometricNormal) < 0;

            // Set both geometric and shading normals to face the ray.
            gn = frontFace ? geometricNormal : -geometricNormal;
            n = gn; // Initially, the shading normal matches the geometric normal.
        }

        /**
         * @brief Offset the ray origin to avoid self-intersections (shadow acne).
         *
         * PBRT-style note:
         * A more robust version can include scale-dependent epsilon and nextafter(),
         * but this simple approach is often sufficient for small/medium scenes.
         *
         * @param w Ray direction (typically wi).
         * @return Offset origin point.
         */
        Point3 OffsetRayOrigin(const Vector3& w) const {
            // Use the GEOMETRIC normal (gn), not the shading normal (n),
            // to prevent normal-mapping artifacts in ray offsets.
            Vector3 offsetN = (glm::dot(gn, w) > 0.0f) ? gn : -gn;
            return p + offsetN * constants::RAY_EPSILON;
        }

        /**
         * @brief Spawn a new ray starting from this surface interaction.
         *
         * The ray origin is offset along the geometric normal to avoid self-intersections.
         * Medium and time are propagated to support volumetric rendering and motion blur.
         *
         * @param wi New ray direction.
         * @return A numerically safe ray starting at the interaction point.
         */
        Ray SpawnRay(const Vector3& wi) const {
            // When origin is already offset, tMin can often be 0 (or a very small value).
            // Keep it minimal to avoid missing near hits.
            constexpr Real kRayTMin = Real(0);

            Ray r(OffsetRayOrigin(wi), wi, kRayTMin, medium);
            r.time = time; // Propagate ray time (important for motion blur).
            return r;
        }

        /**
         * @brief Spawn a shadow ray toward a target with an upper bound tMax.
         *
         * This is used for visibility checks (occlusion). We subtract a small epsilon
         * from tMax to avoid falsely intersecting the light geometry itself.
         *
         * @param wi Ray direction toward the target.
         * @param tMax Maximum distance/parameter to test intersections.
         * @return A shadow ray with a finite intersection interval.
         */
        Ray SpawnRay(const Vector3& wi, Real tMax) const {
            Ray r = SpawnRay(wi);
            r.tMax = std::max(Real(0), tMax - constants::RAY_EPSILON);
            return r;
        }

        // TODO: Implement coordinate system construction for tangent space.
        // Required for Normal Mapping, Anisotropy, and Microfacet distributions.
    };

}