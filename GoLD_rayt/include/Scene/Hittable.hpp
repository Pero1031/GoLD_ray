/**
 * @file Scene/Hittable.hpp
 * @brief Abstract interface for geometric objects.
 * * This header defines the 'Hittable' base class, which provides a unified
 * interface for any object that can be intersected by a ray. It serves as
 * the core abstraction for primitives, aggregates (like BVH), and instances.
 */

#pragma once

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
        * * The ray defines a valid search interval [tMin, tMax] (e.g., for shadow rays).
        * * This method MUST NOT modify the ray. If a hit is found, it writes the hit
        * distance to rec.t and fills other fields in rec.
        * * Selecting the closest hit and pruning farther tests is the caller's
        * responsibility (e.g., Aggregate/BVH manages a "closestSoFar" value and may
        * pass a reduced tMax via a copied Ray).
        * @param r   The incident ray. Its tMin/tMax define the valid search interval.
        * @param rec Output structure to store intersection details (normal, UV, etc.).
        * @return True if an intersection occurs within [r.tMin, r.tMax].
        */
        virtual bool hit(const Ray& r, SurfaceInteraction& rec) const = 0;

        /**
         * @brief Any-hit query for shadow rays: is there ANY intersection within [r.tMin, r.tMax]?
         * Implementations stop at the first one found and fill no hit record, which is much
         * cheaper than hit() for occlusion tests. Default: hit() with a scratch record.
         */
        virtual bool occluded(const Ray& r) const {
            SurfaceInteraction scratch;
            return hit(r, scratch);
        }

        /**
         * @brief Closest-hit distance only, without building a hit record.
         *
         * The BVH uses it while traversing: it only keeps (best primitive, t) and calls hit() once
         * on the winner, so the expensive record (normals, material, ...) is built once per ray
         * instead of once per candidate. Default: hit() with a scratch record.
         */
        virtual bool intersectT(const Ray& r, Real& t) const {
            SurfaceInteraction scratch;
            if (!hit(r, scratch)) return false;
            t = scratch.t;
            return true;
        }

        /**
         * @brief Returns the world-space Axis-Aligned Bounding Box (AABB) of the object.
         * * This is a required interface for building acceleration structures like BVH.
         * The bounding box must encompass the entire geometric extent of the object.
         * @return AABB The bounding volume of the object.
         */
        virtual AABB bounds() const = 0;
    };

} // namespace rayt