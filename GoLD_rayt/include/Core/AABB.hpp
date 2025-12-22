#pragma once

#include <algorithm>
#include <limits>

#include "Core/Types.hpp"
#include "Core/Ray.hpp"

namespace rayt {

    class AABB {
    public:
        Vector3 min;   // lower corner
        Vector3 max;   // upper corner

        AABB() = default;

        AABB(const Vector3& pMin, const Vector3& pMax)
            : min(pMin), max(pMax) {}

        // Ray–AABB intersection.
        // Uses a slab test. The ray itself is NOT modified.
        // tMin/tMax are passed by value and locally updated.
        bool intersect(const Ray& r, Real tMin, Real tMax) const {
            for (int axis = 0; axis < 3; ++axis) {
                const Real invD = Real(1) / r.d[axis];
                Real t0 = (min[axis] - r.o[axis]) * invD;
                Real t1 = (max[axis] - r.o[axis]) * invD;

                if (invD < 0) std::swap(t0, t1);

                tMin = std::max(tMin, t0);
                tMax = std::min(tMax, t1);

                if (tMax <= tMin)
                    return false;
            }
            return true;
        }

        // Returns a new AABB that encloses both a and b.
        static AABB unite(const AABB& a, const AABB& b) {
            return AABB(
                glm::min(a.min, b.min),
                glm::max(a.max, b.max)
            );
        }

        // Optional helper: get center and extent (useful for BVH split heuristics)
        Vector3 center() const {
            return (min + max) * Real(0.5);
        }

        Vector3 extent() const {
            return max - min;
        }
    };

} // namespace rayt
