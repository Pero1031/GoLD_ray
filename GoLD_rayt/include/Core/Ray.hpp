/**
* @file Core/Ray.hpp
* @brief Ray and RayDifferential structures.
* * Fundamental primitives for light transport. RayDifferential provides
* auxiliary information used for antialiasing and texture filtering (LOD).
*/

#pragma once

#include <cmath>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Core/Forward.hpp"
#include "Core/Math.hpp"

namespace rayt {

    /**
     * @brief Represents a semi-infinite line used for ray tracing.
     * Defined by an origin 'o' and a direction 'd'.
     * The ray is parameterized as: P(t) = o + t * d
     * * Designed as a 'struct' (public members) for direct access during
     * heavy mathematical computations, following PBRT standards.
     * 
     * @Coution Upper bound of the valid interval. Used to limit intersection tests (e.g., shadow rays).
     * It is NOT modified by intersection routines; callers manage closest hits.
     * change mutable tMax → non mutable
     */
    struct Ray {
    public:
        Point3 o;       // Origin
        Vector3 d;      // Ray direction (not necessarily normalized, though common)

        /**
         * @brief Upper bound of the valid intersection interval.
         * As the ray hits closer objects, this value decreases to prune the search.
         */
        Real tMax;

        /**
         * @brief Lower bound of the valid intersection interval.
         * Used to prevent self-intersection artifacts (shadow acne).
         * Intersections closer than tMin are ignored.
         */
        Real tMin;

        /**
        * @brief The medium currently surrounding the ray.
        * If not null, the ray will account for volumetric effects (e.g., fog, smoke,
        * or subsurface scattering) as it prepagates through space.
        */
        const Medium* medium = nullptr;

        /**
        * @brief Time value associated with this ray.
        * Specifies the temporal position of the ray within the shutter interval,
         * enabling motion blur and time-varying geometry. This parameter does not
        * affect rendering unless explicitly used by time-dependent primitives
        * or transformations.
        */
        Real time = 0.0;

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
            return math::hasNonFinite(o) || math::hasNonFinite(d) || std::isnan(tMin) || std::isnan(tMax);
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
            // ★ここで事故防止のリセット！
            // ただし「tMaxを引き継ぎたいレアケース」がつぶれるトレードオフはある
            // 変なバグがあったらここを見ろ！
            // tMax = constants::INFINITY_VAL;
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

} // namespace rayt