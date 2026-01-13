/**
 * @file Accelerators/BVH.cpp
 * @brief Bounding Volume Hierarchy (BVH) implementation (definition).
 */

#include "pch.h"

#include <algorithm>
#include <limits>

#include "Core/Ray.hpp"
#include "Core/Interaction.hpp"
#include "Accelerators/BVH.hpp"

namespace rayt {

    BVHNode::BVHNode(std::vector<std::shared_ptr<Hittable>>& objects,
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

    bool BVHNode::hit(const Ray& r, SurfaceInteraction& rec) const {

        // Early exit if the ray doesn't hit this node's bounding box.
        if (!box.intersect(r, r.tMin, r.tMax))
            return false;

        // Leaf node case.
        if (!right) {
            return left->hit(r, rec);
        }

        // Internal node: Determine traversal order (Front-to-Back) based on ray direction.
        const bool dirNeg = r.d[splitAxis] < Real(0);

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

        return hitAnything;
    }

} // namespace rayt