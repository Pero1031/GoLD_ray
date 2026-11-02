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

    // =====================================================================================
    // BVH: binned SAH, flat array, iterative traversal
    // =====================================================================================

    namespace {

        constexpr int kBuckets = 16;
        constexpr int kMaxDepthLimit = 60;     // traversal stack holds 64 entries
        constexpr int kStackSize = 64;

        inline Real surfaceArea(const AABB& b) {
            const Vector3 d = b.max - b.min;
            return Real(2) * (d.x * d.y + d.y * d.z + d.z * d.x);
        }

        inline AABB unite(const AABB& a, const AABB& b) { return AABB::unite(a, b); }

        /// Slab test with a precomputed reciprocal direction (no divisions). NaNs from 0 * inf are ignored.
        inline bool intersectBox(const Real* lo, const Real* hi, const Vector3& o, const Vector3& invD,
            Real tMin, Real tMax)
        {
            for (int a = 0; a < 3; ++a) {
                const Real t0 = (lo[a] - o[a]) * invD[a];
                const Real t1 = (hi[a] - o[a]) * invD[a];
                tMin = std::max(tMin, std::min(t0, t1));
                tMax = std::min(tMax, std::max(t0, t1) * Real(1.000000000000002));   // round-off safety
            }
            return tMin <= tMax;
        }

    } // namespace

    BVH::BVH(std::vector<std::shared_ptr<Hittable>> objects, int maxLeafSize)
        : m_owned(std::move(objects)), m_maxLeaf(std::max(1, maxLeafSize))
    {
        const size_t n = m_owned.size();
        if (n == 0) return;

        // Bounds and centroids once per primitive.
        std::vector<BuildPrim> bp(n);
        for (size_t i = 0; i < n; ++i) {
            bp[i].bounds = m_owned[i]->bounds();
            bp[i].centroid = (bp[i].bounds.min + bp[i].bounds.max) * Real(0.5);
            bp[i].index = static_cast<int>(i);
        }

        m_nodes.reserve(2 * n / static_cast<size_t>(std::max(1, m_maxLeaf / 2)) + 1);
        std::vector<int> order;
        order.reserve(n);
        build(bp, 0, n, 0, order);

        m_prims.resize(n);
        for (size_t i = 0; i < n; ++i) m_prims[i] = m_owned[static_cast<size_t>(order[i])].get();
        m_nodes.shrink_to_fit();
    }

    int BVH::build(std::vector<BuildPrim>& bp, size_t start, size_t end, int depth,
        std::vector<int>& order)
    {
        const int nodeIndex = static_cast<int>(m_nodes.size());
        m_nodes.emplace_back();
        m_maxDepth = std::max(m_maxDepth, depth);

        const size_t count = end - start;

        AABB bounds;
        AABB centroidBounds(Vector3(std::numeric_limits<Real>::infinity()),
            Vector3(-std::numeric_limits<Real>::infinity()));
        for (size_t i = start; i < end; ++i) {
            bounds = unite(bounds, bp[i].bounds);
            centroidBounds.min = glm::min(centroidBounds.min, bp[i].centroid);
            centroidBounds.max = glm::max(centroidBounds.max, bp[i].centroid);
        }
        {
            Node& nb = m_nodes[static_cast<size_t>(nodeIndex)];
            for (int a = 0; a < 3; ++a) { nb.lo[a] = bounds.min[a]; nb.hi[a] = bounds.max[a]; }
        }

        auto makeLeaf = [&]() {
            Node& node = m_nodes[static_cast<size_t>(nodeIndex)];
            node.offset = static_cast<int>(order.size());
            node.nPrims = static_cast<unsigned short>(count);
            for (size_t i = start; i < end; ++i) order.push_back(bp[i].index);
            return nodeIndex;
        };

        if (count == 1 || depth >= kMaxDepthLimit || (count <= 2 && m_maxLeaf >= 2)) return makeLeaf();

        // Split axis: widest centroid extent.
        const Vector3 ext = centroidBounds.max - centroidBounds.min;
        int axis = 0;
        if (ext.y > ext.x) axis = 1;
        if (ext.z > ext[axis]) axis = 2;

        // All centroids coincide: cannot separate them spatially.
        if (ext[axis] <= Real(0)) {
            if (count <= static_cast<size_t>(std::numeric_limits<unsigned short>::max())) return makeLeaf();
            // Absurdly large degenerate set: split in the middle by index.
            const size_t mid = start + count / 2;
            build(bp, start, mid, depth + 1, order);
            m_nodes[static_cast<size_t>(nodeIndex)].offset = build(bp, mid, end, depth + 1, order);
            m_nodes[static_cast<size_t>(nodeIndex)].axis = static_cast<unsigned char>(axis);
            return nodeIndex;
        }

        size_t mid = start;
        if (count <= 2) {
            // Two primitives: order along the axis.
            if (bp[start].centroid[axis] > bp[start + 1].centroid[axis]) std::swap(bp[start], bp[start + 1]);
            mid = start + 1;
        }
        else {
            // Binned SAH.
            struct Bucket { size_t count = 0; AABB bounds; };
            Bucket buckets[kBuckets];
            const Real cmin = centroidBounds.min[axis];
            const Real scale = Real(kBuckets) / ext[axis];
            auto bucketOf = [&](const BuildPrim& p) {
                return std::min(static_cast<int>((p.centroid[axis] - cmin) * scale), kBuckets - 1);
            };
            for (size_t i = start; i < end; ++i) {
                Bucket& bk = buckets[bucketOf(bp[i])];
                ++bk.count;
                bk.bounds = unite(bk.bounds, bp[i].bounds);
            }

            // Cost of splitting after bucket i (prefix / suffix sweeps).
            Real cost[kBuckets - 1];
            {
                AABB b;
                size_t c = 0;
                for (int i = 0; i < kBuckets - 1; ++i) {
                    b = unite(b, buckets[i].bounds);
                    c += buckets[i].count;
                    cost[i] = (c > 0) ? static_cast<Real>(c) * surfaceArea(b) : Real(0);
                }
                b = AABB();
                c = 0;
                for (int i = kBuckets - 1; i >= 1; --i) {
                    b = unite(b, buckets[i].bounds);
                    c += buckets[i].count;
                    cost[i - 1] += (c > 0) ? static_cast<Real>(c) * surfaceArea(b) : Real(0);
                }
            }

            int best = 0;
            Real bestCost = std::numeric_limits<Real>::infinity();
            for (int i = 0; i < kBuckets - 1; ++i) {
                if (cost[i] < bestCost) { bestCost = cost[i]; best = i; }
            }
            // PBRT's constants: traversal 1/8, intersection 1.
            const Real total = surfaceArea(bounds);
            const Real splitCost = Real(0.125) + bestCost / std::max(total, Real(1e-300));
            const Real leafCost = static_cast<Real>(count);

            if (count <= static_cast<size_t>(m_maxLeaf) && splitCost >= leafCost) return makeLeaf();

            const auto it = std::partition(bp.begin() + static_cast<std::ptrdiff_t>(start),
                bp.begin() + static_cast<std::ptrdiff_t>(end),
                [&](const BuildPrim& p) { return bucketOf(p) <= best; });
            mid = static_cast<size_t>(it - bp.begin());
            if (mid == start || mid == end) {
                // SAH put everything on one side (identical bucket): split at the median instead.
                mid = start + count / 2;
                std::nth_element(bp.begin() + static_cast<std::ptrdiff_t>(start),
                    bp.begin() + static_cast<std::ptrdiff_t>(mid),
                    bp.begin() + static_cast<std::ptrdiff_t>(end),
                    [axis](const BuildPrim& a, const BuildPrim& b) { return a.centroid[axis] < b.centroid[axis]; });
            }
        }

        build(bp, start, mid, depth + 1, order);                     // first child = nodeIndex + 1
        const int second = build(bp, mid, end, depth + 1, order);
        Node& node = m_nodes[static_cast<size_t>(nodeIndex)];
        node.offset = second;
        node.axis = static_cast<unsigned char>(axis);
        return nodeIndex;
    }

    bool BVH::hit(const Ray& r, SurfaceInteraction& rec) const {
        if (m_nodes.empty()) return false;

        const Vector3 invD(Real(1) / r.d.x, Real(1) / r.d.y, Real(1) / r.d.z);
        const bool negDir[3] = { invD.x < Real(0), invD.y < Real(0), invD.z < Real(0) };

        Ray ray = r;
        const Hittable* best = nullptr;
        int stack[kStackSize];
        int sp = 0;
        int cur = 0;

        for (;;) {
            const Node& node = m_nodes[static_cast<size_t>(cur)];
            if (intersectBox(node.lo, node.hi, r.o, invD, r.tMin, ray.tMax)) {
                if (node.nPrims > 0) {
                    for (int i = 0; i < node.nPrims; ++i) {
                        const Hittable* prim = m_prims[static_cast<size_t>(node.offset + i)];
                        Real t;
                        if (prim->intersectT(ray, t)) {   // distance only; the record is built once, below
                            best = prim;
                            ray.tMax = t;                 // later tests only accept closer hits
                        }
                    }
                    if (sp == 0) break;
                    cur = stack[--sp];
                }
                else if (negDir[node.axis]) {   // visit the second child first
                    stack[sp++] = cur + 1;
                    cur = node.offset;
                }
                else {
                    stack[sp++] = node.offset;
                    cur = cur + 1;
                }
            }
            else {
                if (sp == 0) break;
                cur = stack[--sp];
            }
        }
        if (!best) return false;

        // Build the hit record for the winner only. The same arithmetic finds the same t, so a
        // range ending exactly at it accepts it; the slightly wider retry is only a round-off guard.
        Ray finalRay = r;
        finalRay.tMax = ray.tMax;
        if (!best->hit(finalRay, rec)) {
            finalRay.tMax = ray.tMax * Real(1.000000001) + Real(1e-12);
            if (!best->hit(finalRay, rec)) return false;
        }
        return true;
    }

    bool BVH::occluded(const Ray& r) const {
        if (m_nodes.empty()) return false;

        const Vector3 invD(Real(1) / r.d.x, Real(1) / r.d.y, Real(1) / r.d.z);
        const bool negDir[3] = { invD.x < Real(0), invD.y < Real(0), invD.z < Real(0) };

        int stack[kStackSize];
        int sp = 0;
        int cur = 0;

        for (;;) {
            const Node& node = m_nodes[static_cast<size_t>(cur)];
            if (intersectBox(node.lo, node.hi, r.o, invD, r.tMin, r.tMax)) {
                if (node.nPrims > 0) {
                    for (int i = 0; i < node.nPrims; ++i) {
                        if (m_prims[static_cast<size_t>(node.offset + i)]->occluded(r)) return true;
                    }
                    if (sp == 0) break;
                    cur = stack[--sp];
                }
                else if (negDir[node.axis]) {
                    stack[sp++] = cur + 1;
                    cur = node.offset;
                }
                else {
                    stack[sp++] = node.offset;
                    cur = cur + 1;
                }
            }
            else {
                if (sp == 0) break;
                cur = stack[--sp];
            }
        }
        return false;
    }

} // namespace rayt