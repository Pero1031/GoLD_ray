#pragma once

/**
 * @file AABB.hpp
 * @brief Axis-Aligned Bounding Box (AABB) implementation.
 * * AABBs are the fundamental building blocks for acceleration structures (BVH).
 * They provide a fast way to cull groups of objects that a ray cannot possibly hit.
 */

#include <algorithm>
#include <cmath>

#include "Core/Types.hpp"
#include "Core/Ray.hpp"
#include "Core/Constants.hpp"

namespace rayt {

    /**
     * @brief Represents an Axis-Aligned Bounding Box.
     * * Defined by its minimum and maximum corners in world space.
     */
    class AABB {
    public:
        Vector3 min;   // Lower-left-front corner
        Vector3 max;   // Upper-right-back corner

        /**
         * @brief Default constructor creates an INVALID (empty) AABB.
         * Initialize min to +Infinity and max to -Infinity.
         * This ensures that any subsequent 'unite' operation will overwrite these values.
         */
        AABB() {
            const Real inf = constants::INFINITY_VAL;
            min = Vector3(inf, inf, inf);
            max = Vector3(-inf, -inf, -inf);
        }

        /**
         * @brief Constructs an AABB from two points.
         */
        AABB(const Vector3& pMin, const Vector3& pMax)
            : min(pMin), max(pMax) {}

        /**
         * @brief Ray–AABB intersection test using the Slab Method.
         * * This is a highly optimized conservative test to determine if a ray
         * potentially intersects any geometry within the box.
         * * @param r    The incident ray.
         * @param tMin The start of the ray's valid interval.
         * @param tMax The end of the ray's valid interval.
         * @return True if the ray's path overlaps with the AABB's volume.
         * * @note This implementation handles parallel rays and division-by-zero
         * through floating-point standard behavior (infinity).
         */
        bool intersect(const Ray& r, Real tMin, Real tMax) const {

            for (int axis = 0; axis < 3; ++axis) {
                const Real origin = r.o[axis];
                const Real direction = r.d[axis];

                if (std::abs(direction) < constants::INTERSECT_TOLERANCE) {
                    if (origin < min[axis] || origin > max[axis]) return false;
                    continue;
                }

                // Compute the distance to the two slabs on the current axis.
                const Real invD = Real(1) / direction;
                Real t0 = (min[axis] - origin) * invD;
                Real t1 = (max[axis] - origin) * invD;

                // If the ray direction is negative, swap the entry and exit points.
                if (invD < Real(0)) std::swap(t0, t1);

                // Tighten the intersection interval.
                tMin = (t0 > tMin) ? t0 : tMin;
                tMax = (t1 < tMax) ? t1 : tMax;

                // If the entry point is beyond the exit point, there is no intersection.
                if (tMax < tMin)
                    return false;
            }
            return true;
        }

        /**
         * @brief Merges two AABBs into a single bounding box that encloses both.
         * @return A new AABB representing the union of a and b.
         */
        static AABB unite(const AABB& a, const AABB& b) {
            return AABB(
                glm::min(a.min, b.min),
                glm::max(a.max, b.max)
            );
        }

        /**
         * @brief Returns the center point of the AABB.
         * Useful for BVH construction heuristics (e.g., sorting by centers).
         */
        Vector3 center() const {
            return (min + max) * Real(0.5);
        }

        /**
         * @brief Returns the diagonal extent of the AABB.
         */
        Vector3 extent() const {
            return max - min;
        }
    };

} // namespace rayt