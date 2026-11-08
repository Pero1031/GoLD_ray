/**
 * @file Geometry/HairCurves.hpp
 * @brief Many hair fibers as chains of round cylinders, with their own BVH.
 *
 * A strand is a polyline with a radius per vertex. Each segment is intersected as a finite round
 * cylinder (radius = the mean of its end radii): one quadratic per segment, no caps (a strand tip
 * tapers to nothing). Only the entry hit counts: a ray that starts inside a fiber (the continuation of a
 * transmitted ray, which the integrator offsets to the inside by SurfaceInteraction::rayOffset) does not hit
 * the fiber's own far wall; the Hair material has already accounted for the light inside.
 *
 * The hit record follows the convention of Materials/Hair.hpp: dpdu = the fiber tangent (interpolated
 * between the vertex tangents, so the shading is smooth along a strand), n = gn = the cylinder's outward
 * normal at the hit, made perpendicular to that tangent. uv = (position along the strand 0..1, 0.5).
 *
 * Why a dedicated BVH: a head of hair is 10^5 strands of ~20 segments; one Hittable per segment would cost
 * a shared_ptr, a vtable and a material pointer each. Here a segment is its first vertex index (4 bytes).
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <vector>

#include "Core/Types.hpp"
#include "Core/AABB.hpp"
#include "Core/Interaction.hpp"
#include "Geometry/Shape.hpp"
#include "Materials/Material.hpp"

namespace rayt {

    class HairCurves : public Shape {
    public:
        explicit HairCurves(std::shared_ptr<Material> material) : m_material(std::move(material)) {}

        /// Adds a strand (at least 2 points); the radius goes linearly from rootRadius to tipRadius. Call build() when done.
        void addStrand(const std::vector<Point3>& pts, Real rootRadius, Real tipRadius) {
            const size_t n = pts.size();
            if (n < 2) return;
            const uint32_t base = static_cast<uint32_t>(m_p.size());
            for (size_t i = 0; i < n; ++i) {
                m_p.push_back(pts[i]);
                const Real s = Real(i) / Real(n - 1);
                m_r.push_back(rootRadius + (tipRadius - rootRadius) * s);
                // central difference tangent (one-sided at the ends)
                const Point3& a = pts[i == 0 ? 0 : i - 1];
                const Point3& b = pts[i + 1 < n ? i + 1 : n - 1];
                Vector3 t = b - a;
                const Real len = glm::length(t);
                m_t.push_back(len > Real(0) ? t / len : Vector3(0, -1, 0));
            }
            for (size_t i = 0; i + 1 < n; ++i) m_segments.push_back(base + static_cast<uint32_t>(i));
            ++m_numStrands;
        }

        /// Builds the BVH over all segments added so far.
        void build() {
            const size_t n = m_segments.size();
            m_nodes.clear();
            if (n == 0) return;
            std::vector<Prim> prims(n);
            for (size_t k = 0; k < n; ++k) {
                const uint32_t i = m_segments[k];
                // pad by the radius plus the overlap at the joints (a quarter radius along the axis): sqrt(1 + 1/16) r < 1.1 r
                const Real r = std::max(m_r[i], m_r[i + 1]) * Real(1.1);
                Prim& pr = prims[k];
                pr.seg = i;
                for (int a = 0; a < 3; ++a) {
                    pr.lo[a] = std::min(m_p[i][a], m_p[i + 1][a]) - r;
                    pr.hi[a] = std::max(m_p[i][a], m_p[i + 1][a]) + r;
                    pr.c[a] = Real(0.5) * (pr.lo[a] + pr.hi[a]);
                }
            }
            m_nodes.reserve(n / 2 + 16);
            buildNode(prims, 0, n);
            for (size_t k = 0; k < n; ++k) m_segments[k] = prims[k].seg;   // leaf order
        }

        size_t numStrands() const { return m_numStrands; }
        size_t numSegments() const { return m_segments.size(); }

        // ---- Hittable ----

        bool hit(const Ray& r, SurfaceInteraction& rec) const override {
            Real t = 0, u = 0;
            uint32_t seg = 0;
            if (!closest(r, t, u, seg)) return false;
            fillRecord(r, t, u, seg, rec);
            return true;
        }

        bool occluded(const Ray& r) const override { Real t, u; uint32_t s; return traverse(r, true, t, u, s); }

        bool intersectT(const Ray& r, Real& tOut) const override {
            Real u; uint32_t s;
            return closest(r, tOut, u, s);
        }

        AABB bounds() const override {
            if (m_nodes.empty()) return AABB(Point3(0.0), Point3(0.0));
            const Node& n = m_nodes[0];
            return AABB(Point3(n.lo[0], n.lo[1], n.lo[2]), Point3(n.hi[0], n.hi[1], n.hi[2]));
        }

        // ---- Shape (a hair is not an emitter; sampling is not needed) ----
        Real area() const override {
            Real a = 0;
            for (uint32_t i : m_segments) a += Real(2) * constants::PI * Real(0.5) * (m_r[i] + m_r[i + 1]) * glm::length(m_p[i + 1] - m_p[i]);
            return a;
        }
        std::optional<ShapeSample> sampleSurface(const Point2&) const override { return std::nullopt; }
        Real pdfSurface(const SurfaceInteraction&, const Vector3&) const override { return Real(0); }

    private:
        struct Node {
            Real lo[3], hi[3];
            uint32_t first = 0;    // leaf: first segment (index into m_segments); inner: index of the second child
            uint32_t count = 0;    // leaf: number of segments; inner: 0
        };
        struct Prim { Real lo[3], hi[3], c[3]; uint32_t seg; };

        static constexpr size_t kLeafSize = 4;

        uint32_t buildNode(std::vector<Prim>& prims, size_t begin, size_t end) {
            const uint32_t idx = static_cast<uint32_t>(m_nodes.size());
            m_nodes.emplace_back();
            Real lo[3] = { INF, INF, INF }, hi[3] = { -INF, -INF, -INF };
            Real clo[3] = { INF, INF, INF }, chi[3] = { -INF, -INF, -INF };
            for (size_t k = begin; k < end; ++k)
                for (int a = 0; a < 3; ++a) {
                    lo[a] = std::min(lo[a], prims[k].lo[a]); hi[a] = std::max(hi[a], prims[k].hi[a]);
                    clo[a] = std::min(clo[a], prims[k].c[a]); chi[a] = std::max(chi[a], prims[k].c[a]);
                }
            for (int a = 0; a < 3; ++a) { m_nodes[idx].lo[a] = lo[a]; m_nodes[idx].hi[a] = hi[a]; }

            const size_t count = end - begin;
            int axis = 0;
            for (int a = 1; a < 3; ++a) if (chi[a] - clo[a] > chi[axis] - clo[axis]) axis = a;
            if (count <= kLeafSize || !(chi[axis] > clo[axis])) {
                m_nodes[idx].first = static_cast<uint32_t>(begin);
                m_nodes[idx].count = static_cast<uint32_t>(count);
                return idx;
            }
            const size_t mid = begin + count / 2;
            std::nth_element(prims.begin() + begin, prims.begin() + mid, prims.begin() + end,
                [axis](const Prim& x, const Prim& y) { return x.c[axis] < y.c[axis]; });
            buildNode(prims, begin, mid);
            const uint32_t second = buildNode(prims, mid, end);
            m_nodes[idx].first = second;
            m_nodes[idx].count = 0;
            return idx;
        }

        static constexpr Real INF = std::numeric_limits<Real>::infinity();

        /// Entry hit of a ray with the finite round cylinder of the segment starting at vertex i.
        bool intersectSegment(const Ray& r, Real tMax, uint32_t i, Real& tOut, Real& uOut) const {
            const Vector3& p0 = m_p[i];
            const Vector3 a = m_p[i + 1] - p0;
            const Real len2 = glm::dot(a, a);
            if (!(len2 > Real(0))) return false;
            const Real len = std::sqrt(len2);
            const Vector3 ax = a / len;
            const Real rad = Real(0.5) * (m_r[i] + m_r[i + 1]);

            const Vector3 oc = r.o - p0;
            const Real dA = glm::dot(r.d, ax), oA = glm::dot(oc, ax);
            const Vector3 dp = r.d - ax * dA, op = oc - ax * oA;
            const Real A = glm::dot(dp, dp);
            if (A < Real(1e-20) * glm::dot(r.d, r.d)) return false;   // parallel to the axis
            const Real B = glm::dot(op, dp);
            const Real C = glm::dot(op, op) - rad * rad;
            const Real disc = B * B - A * C;
            if (disc < Real(0)) return false;
            const Real t = (-B - std::sqrt(disc)) / A;   // the entry; C < 0 (origin inside) gives t < 0, rejected below
            if (t < r.tMin || t > tMax) return false;
            const Real s = oA + t * dA;
            const Real slack = Real(0.25) * rad;         // the segments overlap a little at the joints (no pinholes on the outside of a bend)
            if (s < -slack || s > len + slack) return false;
            tOut = t;
            uOut = std::clamp(s / len, Real(0), Real(1));
            return true;
        }

        bool boxHit(const Node& n, const Vector3& o, const Vector3& invD, Real tMin, Real tMax) const {
            for (int a = 0; a < 3; ++a) {
                Real t0 = (n.lo[a] - o[a]) * invD[a];
                Real t1 = (n.hi[a] - o[a]) * invD[a];
                if (t0 > t1) std::swap(t0, t1);
                tMin = std::max(tMin, t0);
                tMax = std::min(tMax, t1);
                if (tMax < tMin) return false;
            }
            return true;
        }

        /// Closest hit (anyHit = false) or any hit (true).
        bool traverse(const Ray& r, bool anyHit, Real& tBest, Real& uBest, uint32_t& segBest) const {
            if (m_nodes.empty()) return false;
            const Vector3 invD(Real(1) / r.d.x, Real(1) / r.d.y, Real(1) / r.d.z);
            Real tMax = r.tMax;
            bool found = false;
            uint32_t stack[64];
            int sp = 0;
            stack[sp++] = 0;
            while (sp > 0) {
                const uint32_t ni = stack[--sp];
                const Node& n = m_nodes[ni];
                if (!boxHit(n, r.o, invD, r.tMin, tMax)) continue;
                if (n.count > 0) {
                    for (uint32_t k = 0; k < n.count; ++k) {
                        const uint32_t seg = m_segments[n.first + k];
                        Real t, u;
                        if (intersectSegment(r, tMax, seg, t, u)) {
                            if (anyHit) return true;
                            tMax = t; tBest = t; uBest = u; segBest = seg; found = true;
                        }
                    }
                }
                else {
                    // visit the nearer child first
                    const Node& c0 = m_nodes[ni + 1];
                    const Node& c1 = m_nodes[n.first];
                    const Real d0 = (c0.lo[0] + c0.hi[0]) * r.d.x + (c0.lo[1] + c0.hi[1]) * r.d.y + (c0.lo[2] + c0.hi[2]) * r.d.z;
                    const Real d1 = (c1.lo[0] + c1.hi[0]) * r.d.x + (c1.lo[1] + c1.hi[1]) * r.d.y + (c1.lo[2] + c1.hi[2]) * r.d.z;
                    if (d0 <= d1) { stack[sp++] = n.first; stack[sp++] = ni + 1; }
                    else { stack[sp++] = ni + 1; stack[sp++] = n.first; }
                }
            }
            return found;
        }

        bool closest(const Ray& r, Real& t, Real& u, uint32_t& seg) const { return traverse(r, false, t, u, seg); }

        void fillRecord(const Ray& r, Real t, Real u, uint32_t i, SurfaceInteraction& rec) const {
            const Point3 p = r.at(t);
            const Vector3& p0 = m_p[i];
            const Vector3 a = m_p[i + 1] - p0;
            const Vector3 ax = a / glm::length(a);
            // smooth tangent along the strand
            Vector3 T = m_t[i] * (Real(1) - u) + m_t[i + 1] * u;
            const Real tl = glm::length(T);
            T = tl > Real(0) ? T / tl : ax;
            // outward normal: the radial direction, made perpendicular to the smooth tangent
            const Real s = glm::dot(p - p0, ax);
            Vector3 radial = p - (p0 + ax * s);
            radial -= T * glm::dot(radial, T);
            const Real rl = glm::length(radial);
            const Vector3 n = rl > Real(0) ? radial / rl : glm::normalize(glm::cross(T, Vector3(0, 0, 1)) + Vector3(1e-9));

            rec.t = t;
            rec.p = p;
            rec.n = n;
            rec.gn = n;
            rec.frontFace = true;
            rec.dpdu = T;
            rec.uv = UV(u, Real(0.5));
            rec.matPtr = m_material.get();
            rec.rayOffset = Real(1e-3) * Real(0.5) * (m_r[i] + m_r[i + 1]);
        }

        std::shared_ptr<Material> m_material;
        std::vector<Vector3> m_p, m_t;
        std::vector<Real> m_r;
        std::vector<uint32_t> m_segments;   ///< first vertex of each segment (in leaf order after build())
        std::vector<Node> m_nodes;
        size_t m_numStrands = 0;
    };

} // namespace rayt
