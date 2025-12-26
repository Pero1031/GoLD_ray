#pragma once

/**
 * @file BVH.hpp
 * @brief Bounding Volume Hierarchy (BVH) implementation.
 * * BVH is a spatial acceleration structure that organizes objects into a
 * hierarchical tree of AABBs. This allows the intersection complexity to be
 * reduced from O(N) to O(log N), making complex scene rendering feasible.
 */

#include <memory>
#include <vector>
#include <algorithm>
#include <random>
#include <limits>

#include "Core/Types.hpp"
#include "Core/Ray.hpp"
#include "Core/Interaction.hpp"
#include "Core/AABB.hpp"
#include "Geometry/Hittable.hpp"   

namespace rayt {

    /**
     * @brief A node in the Bounding Volume Hierarchy.
     * * BVHNode acts as both an internal node (containing two children) and
     * a potential leaf (containing a primitive). It recursively partitions
     * objects along the axis of maximum extent to maintain a balanced tree.
     */
    class BVHNode : public Hittable {
    public:
        std::shared_ptr<Hittable> left;
        std::shared_ptr<Hittable> right;
        AABB box;
        int splitAxis = 0;

        BVHNode() = default;

        /**
         * @brief Recursively constructs a BVH tree from a list of objects.
         * * This implementation uses a robust Median Split heuristic:
         * 1. Calculate the centroid bounds of all objects.
         * 2. Choose the axis with the maximum extent to minimize overlap.
         * 3. Sort and partition objects into two relatively equal sets.
         * * @param objects A vector of Hittable objects to partition.
         * @param start   The starting index in the objects vector.
         * @param end     The ending index in the objects vector.
         */
        BVHNode(std::vector<std::shared_ptr<Hittable>>& objects,
            size_t start, size_t end)
        {
            // ここは「とりあえず動く」median split（SAHは後でOK）
            const size_t span = end - start;

            // Compute the centroid bounding box to determine the optimal split axis.
            AABB centroidBox;
            {
                Vector3 cmin(std::numeric_limits<Real>::infinity());
                Vector3 cmax(-std::numeric_limits<Real>::infinity());
                for (size_t i = start; i < end; ++i) {
                    const AABB b = objects[i]->bounds();
                    const Vector3 c = (b.min + b.max) * Real(0.5);
                    cmin = glm::min(cmin, c);
                    cmax = glm::max(cmax, c);
                }
                centroidBox = AABB(cmin, cmax);
            }

            // Select the split axis based on the largest centroid extent (Heuristic).
            const Vector3 e = centroidBox.extent();
            int axis = 0;
            if (e.y > e.x) axis = 1;
            if (e.z > (axis == 0 ? e.x : e.y)) axis = 2;

            this->splitAxis = axis;

            // Comparator for sorting objects based on their centroids along the chosen axis.
            auto cmp = [axis](const std::shared_ptr<Hittable>& a,
                const std::shared_ptr<Hittable>& b) {
                    const AABB ba = a->bounds();
                    const AABB bb = b->bounds();
                    const Vector3 ca = (ba.min + ba.max) * Real(0.5);
                    const Vector3 cb = (bb.min + bb.max) * Real(0.5);
                    return ca[axis] < cb[axis];
                };

            if (span == 1) {
                left = objects[start];
                right = nullptr;
            }
            else if (span == 2) {
                // For exactly two objects, simply order them along the axis.
                if (cmp(objects[start], objects[start + 1])) {
                    left = objects[start];
                    right = objects[start + 1];
                }
                else {
                    left = objects[start + 1];
                    right = objects[start];
                }
            }
            else {
                // Median split: sort and divide the objects into two halves.
                std::sort(objects.begin() + start, objects.begin() + end, cmp);
                const size_t mid = start + span / 2;
                left = std::make_shared<BVHNode>(objects, start, mid);
                right = std::make_shared<BVHNode>(objects, mid, end);
            }

            // Consolidate the AABB to enclose all children.
            if (!right) {
                box = left->bounds();
            }
            else {
                box = AABB::unite(left->bounds(), right->bounds());
            }
        }

        /**
        * @brief Traverses the BVH tree to find the closest intersection.
        * * The ray is treated as read-only. This method keeps track of the closest
        * hit using a local variable (closestSoFar) and prunes farther traversal
        * by passing a reduced tMax via copied rays.
        * * Children are traversed in front-to-back order based on ray direction
        * to improve early pruning.
        */
        bool hit(const Ray& r, SurfaceInteraction& rec) const override {

            // Early exit if the ray doesn't hit this node's bounding box.
            if (!box.intersect(r, r.tMin, r.tMax))
                return false;

            // Leaf node case.
            if (!right) {
                return left->hit(r, rec);
            }

            // Internal node: Determine traversal order (Front-to-Back) based on ray direction.
            // This maximizes the chance of 'r.tMax' being shortened early.
            const bool dirNeg = r.d[splitAxis] < Real(0);

            // rightはnullptrではないことが保証されているので安全
            auto first = dirNeg ? right : left;
            auto second = dirNeg ? left : right;

            bool hitAnything = false;
            Real closestSoFar = r.tMax;

            SurfaceInteraction tempRec;

            Ray localRay = r;

            // ---- 1st child ----
            localRay.tMax = closestSoFar;
            if (first->hit(localRay, tempRec)) {
                hitAnything = true;
                closestSoFar = tempRec.t;
                rec = tempRec;
            }

            // ---- 2nd child ----
            localRay.tMax = closestSoFar;
            if (second->hit(localRay, tempRec)) {
                hitAnything = true;
                closestSoFar = tempRec.t;
                rec = tempRec;
            }

            /*
            // Test the closer child first.
            if (first->hit(r, rec)) hitAny = true;

            // Test the second child. 
            // Note: second->hit() will automatically use the updated (shortened) r.tMax.
            if (second->hit(r, rec)) hitAny = true;
            */

            return hitAnything;
        }

        /**
         * @brief Returns the bounding box for the entire BVH subtree.
         */
        AABB bounds() const override { return box; }
    };

} // namespace rayt
