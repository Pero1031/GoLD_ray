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
     * @brief Simple recursive BVH (median split, one primitive per leaf). KEPT AS A REFERENCE for
     *        tests and benchmarks; use BVH below for rendering.
     *
     * A node in the Bounding Volume Hierarchy.
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

    /**
     * @brief Production BVH: binned SAH build, flat node array, iterative traversal.
     *
     *  - Build: surface area heuristic with 16 buckets on the widest centroid axis; leaves hold up
     *    to maxLeafSize primitives. Bounds are computed once per primitive (no virtual calls in the
     *    sort), partitioning is O(n) per level.
     *  - Layout: nodes in depth-first order in one contiguous array (the first child directly
     *    follows its parent, the second child is stored as an index); primitives are stored in leaf
     *    order as raw pointers.
     *  - hit(): stack-based, ray reciprocal computed once, slab test without divisions, the nearer
     *    child first, and the hit record is written in place (no scratch copies).
     *  - occluded(): same traversal that stops at the first hit (shadow rays).
     *
     * The BVH keeps its primitives alive through shared_ptr.
     *
     * recordInPlace (for the BVH over the objects of a whole scene): hit() calls each candidate's hit() directly and keeps
     * the closest record, as a plain list would. The default (distance first, one hit() on the winner) is right for
     * triangles, whose record is expensive; for an object that is itself a mesh BVH it would traverse that mesh twice.
     */
    class BVH : public Hittable {
    public:
        explicit BVH(std::vector<std::shared_ptr<Hittable>> objects, int maxLeafSize = 4, bool recordInPlace = false);

        bool hit(const Ray& r, SurfaceInteraction& rec) const override;
        bool occluded(const Ray& r) const override;
        AABB bounds() const override {
            if (m_nodes.empty()) return AABB(Point3(0), Point3(0));
            const Node& n = m_nodes[0];
            return AABB(Vector3(n.lo[0], n.lo[1], n.lo[2]), Vector3(n.hi[0], n.hi[1], n.hi[2]));
        }

        size_t numNodes() const { return m_nodes.size(); }
        size_t numPrimitives() const { return m_prims.size(); }
        int maxDepth() const { return m_maxDepth; }

    private:
        /// 56 bytes with plain arrays (no SIMD alignment padding): more nodes per cache line.
        struct Node {
            Real lo[3] = { 0, 0, 0 };
            Real hi[3] = { 0, 0, 0 };
            int offset = 0;              ///< leaf: index of the first primitive; interior: index of the second child
            unsigned short nPrims = 0;   ///< 0 for interior nodes
            unsigned char axis = 0;      ///< split axis of interior nodes
        };

        struct BuildPrim {
            AABB bounds;
            Vector3 centroid;
            int index = 0;
        };

        int build(std::vector<BuildPrim>& bp, size_t start, size_t end, int depth,
            std::vector<int>& order);

        std::vector<Node> m_nodes;
        std::vector<std::shared_ptr<Hittable>> m_owned;   ///< keeps the primitives alive
        std::vector<const Hittable*> m_prims;             ///< leaf order
        int m_maxLeaf = 4;
        int m_maxDepth = 0;
        bool m_recordInPlace = false;
    };

} // namespace rayt