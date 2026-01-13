/**
 * @file Accelerators/BVH.hpp
 * @brief Bounding Volume Hierarchy (BVH) implementation (declaration).
 */

#pragma once

#include <memory>
#include <vector>
#include <cstddef> // size_t

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/AABB.hpp"
#include "Scene/Hittable.hpp"

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
         * * Median split (SAH can be added later).
         * @param objects A vector of Hittable objects to partition.
         * @param start   The starting index in the objects vector.
         * @param end     The ending index in the objects vector.
         */
        BVHNode(std::vector<std::shared_ptr<Hittable>>& objects,
            size_t start, size_t end);

        /**
         * @brief Traverses the BVH tree to find the closest intersection.
         * * Children are traversed in front-to-back order based on ray direction.
         */
        bool hit(const Ray& r, SurfaceInteraction& rec) const override;

        /**
         * @brief Returns the bounding box for the entire BVH subtree.
         */
        AABB bounds() const override { return box; }
    };

} // namespace rayt