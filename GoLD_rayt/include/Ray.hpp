#pragma once

#include "Core.hpp"

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
        Vector3 d;      // Direction

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

        /**
         * @brief Default constructor.
         * Initializes an invalid ray with infinite range.
         */
        Ray()
            : o(0.0), d(0.0, 0.0, 1.0)
            , tMin(Constants::RAY_EPSILON), tMax(Constants::INFINITY_VAL) {}

        /**
         * @brief Primary constructor.
         * @param o     The origin of the ray.
         * @param d     The direction of the ray.
         * @param tMin  Start distance to avoid self-intersection (default: EPSILON).
         */
        Ray(const Point3& o, const Vector3& d, Real tMin = Constants::RAY_EPSILON)
            : o(o), d(d), tMin(tMin), tMax(Constants::INFINITY_VAL) {
        }

        // Public Methods

        /**
         * @brief Calculates the point along the ray at parameter t.
         * Formula: P = o + t * d
         */
        Point3 at(Real t) const {
            return o + d * t;
        }
    };

} // namespace rayt