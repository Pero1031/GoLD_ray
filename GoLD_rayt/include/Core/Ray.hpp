#pragma once

/**
* @file Ray.hpp
* @brief Ray and RayDifferential structures.
* * Fundamental primitives for light transport. RayDifferential provides
* auxiliary information used for antialiasing and texture filtering (LOD).
*/

#include <cmath>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Core/Forward.hpp"
#include "Core/Math.hpp"

namespace rayt {

    // Ray class represents a semi-infinite line starting from an origin.
    // In PBRT, rays also carry information about time (for motion blur)
    // and the medium they are traveling through.
    /**
     * @brief Represents a semi-infinite line used for ray tracing.
     * Defined by an origin 'o' and a direction 'd'.
     * The ray is parameterized as: P(t) = o + t * d
     * * Designed as a 'struct' (public members) for direct access during
     * heavy mathematical computations, following PBRT standards.
     */
    struct Ray {
    public:
        Point3 o;       // Origin
        Vector3 d;      // Ray direction (not necessarily normalized, though common)

        /**
         * @brief Upper bound of the valid intersection interval.
         * Marked 'mutable' to allow updates within const methods (e.g., during traversal).
         * As the ray hits closer objects, this value decreases to prune the search.
         */
        mutable Real tMax;

        /**
         * @brief Lower bound of the valid intersection interval.
         * Used to prevent self-intersection artifacts (shadow acne).
         * Intersections closer than tMin are ignored.
         */
        Real tMin;

        // Participating Media
        // Enable this when rendering effects like "foggy environments," 
        // moving beyond simple surface-level metallic reflections.
        const Medium* medium = nullptr;

        /**
         * @brief Default constructor.
         * Initializes an invalid ray with infinite range.
         */
        Ray()
            : o(0.0), d(0.0, 0.0, 1.0)
            , tMin(constants::RAY_EPSILON), tMax(constants::INFINITY_VAL) {}

        /**
         * @brief Primary constructor.
         * @param o      The origin of the ray.
         * @param d      The direction of the ray.
         * @param tMin   Start distance to avoid self-intersection (default: EPSILON).
         * @param medium The medium containing the ray (default: nullptr).
         */
        Ray(const Point3& o, const Vector3& d, Real tMin = constants::RAY_EPSILON, const Medium* medium = nullptr)
            : o(o), d(d), tMin(tMin), tMax(constants::INFINITY_VAL), medium(medium) {
        }

        // Public Methods

        /**
         * @brief Calculates the point along the ray at parameter t.
         * @param t Distance from the origin.
         * @return Point3 Formula: P(t) = o + t * d
         */
        Point3 at(Real t) const {
            return o + d * t;
        }

        /**
         * @brief Debug utility to check for invalid numerical data.
         * @return True if any component is NaN or Infinity.
         */
        bool HasNaN() const {
            return math::hasNaNs(o) || math::hasNaNs(d) || std::isnan(tMin) || std::isnan(tMax);
        }
    };

    // -------------------------------------------------------------------------
    // RayDifferential Class (Critical for Texture Filtering)
    // -------------------------------------------------------------------------
    
    /**
     * @brief Extends the base Ray to include auxiliary rays for adjacent pixels.
     * Used for ray differentials, allowing the integrator to estimate the
     * footprint of the ray for high-quality texture filtering (LOD).
     */
    struct RayDifferential : public Ray {

        bool hasDifferentials = false;
        Point3 rxOrigin = Point3(0.0);
        Point3 ryOrigin = Point3(0.0);
        Vector3 rxDirection = Vector3(0.0);
        Vector3 ryDirection = Vector3(0.0);

        RayDifferential() = default;

        RayDifferential(const Point3& o, const Vector3& d, Real tMin = constants::RAY_EPSILON, const Medium* medium = nullptr)
            : Ray(o, d, tMin, medium) {
            hasDifferentials = false;
        }

        /**
         * @brief Implicit conversion from a standard Ray.
         * Useful for passing standard rays into functions expecting differentials.
         */
        RayDifferential(const Ray& ray) : Ray(ray) {
            hasDifferentials = false;
        }

        /**
         * @brief Scales the differential information.
         * Used to account for the spreading of rays during specular reflections or
         * transmissions, ensuring proper footprint estimation.
         * @param s Scaling factor.
         */
        void ScaleDifferentials(Real s) {
            if (!hasDifferentials) return;

            rxOrigin = o + (rxOrigin - o) * s;
            ryOrigin = o + (ryOrigin - o) * s;
            rxDirection = d + (rxDirection - d) * s;
            ryDirection = d + (ryDirection - d) * s;
        }
    };

    // -------------------------------------------------------------------------
    // Utility Functions
    // -------------------------------------------------------------------------

    /**
     * @brief Spawns a new ray from a surface point.
     * Correctly handles numerical offsets to avoid self-intersection (Shadow Acne).
     * * @param p   The surface hit point (origin for the new ray).
     * @param n   The surface normal at the hit point.
     * @param wi  The new ray direction (already calculated, e.g., via BRDF sampling).
     * @param med The medium containing the new ray.
     * @return Ray The newly constructed ray with proper safety offsets.
     */
    inline Ray SpawnRay(const Point3& p, const Vector3& n, const Vector3& wi, const Medium* med = nullptr) {

        // Offset the starting point along +n if wi is in the same hemisphere, or -n if refracted.
        // Without this, glass rendering will be plagued by self-intersection noise.
        Vector3 offset = glm::dot(n, wi) > 0.0 ? n : -n;

        // Offset the ray origin slightly from the hit point p along the ray direction.
        Point3 origin = p + offset * constants::RAY_EPSILON;

        return Ray(origin, wi, constants::RAY_EPSILON, med);
    }

} // namespace rayt