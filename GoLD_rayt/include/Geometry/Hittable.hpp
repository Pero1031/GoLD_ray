#pragma once

/**
 * @file Hittable.hpp
 * @brief Abstract interface for geometric objects.
 * * This header defines the 'Hittable' base class, which provides a unified
 * interface for any object that can be intersected by a ray. It serves as
 * the core abstraction for primitives, aggregates (like BVH), and instances.
 */

#include "Core/Forward.hpp"
#include "Core/Ray.hpp"
#include "Core/Interaction.hpp"
#include "Core/AABB.hpp"

namespace rayt {

    /**
     * @brief Abstract base class for any object that a ray can intersect.
     * * Following the PBRT-style architecture, the intersection range is strictly
     * managed within the Ray object (tMin and tMax). This allows for efficient
     * culling during acceleration structure traversal.
     */
    class Hittable {
    public:
        virtual ~Hittable() = default;

        /**
         * @brief Determines if a ray hits this object within its valid interval.
         * * If an intersection is found, the 'tMax' of the ray should be updated
         * to the hit distance to facilitate pruning of subsequent intersection tests.
         * * @param r   The incident ray. Its tMin/tMax define the valid search interval.
         * @param rec Output structure to store intersection details (normal, UV, etc.).
         * @return True if an intersection occurs within the ray's valid interval [tMin, tMax].
         */
        virtual bool hit(const Ray& r, SurfaceInteraction& rec) const = 0;

        /**
         * @brief Returns the world-space Axis-Aligned Bounding Box (AABB) of the object.
         * * This is a required interface for building acceleration structures like BVH.
         * The bounding box must encompass the entire geometric extent of the object.
         * @return AABB The bounding volume of the object.
         */
        virtual AABB bounds() const = 0;
    };

} // namespace rayt