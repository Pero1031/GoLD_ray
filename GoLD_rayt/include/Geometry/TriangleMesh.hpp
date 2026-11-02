/**
 * @file Geometry/TriangleMesh.hpp
 * @brief Indexed triangle mesh, the Triangle shape, a procedural icosphere and a BVH builder.
 *
 * Winding: counter-clockwise seen from outside defines the geometric normal
 * cross(p1 - p0, p2 - p0). Per-vertex normals (optional) are interpolated for shading and also
 * decide which side is "outside" if they disagree with the winding (as in PBRT).
 *
 * A Triangle is a Shape, so an emissive mesh can be used as area lights (one AreaLight per triangle).
 */

#pragma once

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <vector>

#include "Core/Constants.hpp"
#include "Core/Math.hpp"
#include "Core/Interaction.hpp"
#include "Core/AABB.hpp"
#include "Accelerators/BVH.hpp"
#include "Core/Texture.hpp"
#include "Geometry/Shape.hpp"
#include "Scene/Hittable.hpp"

namespace rayt {

    struct TriangleMesh {
        std::vector<Point3> p;       ///< vertex positions
        std::vector<Vector3> n;      ///< per-vertex normals (empty = flat shading)
        std::vector<Point2> uv;      ///< per-vertex texture coordinates (may be empty)
        std::vector<int> indices;    ///< 3 vertex indices per triangle

        size_t numTriangles() const { return indices.size() / 3; }
        bool hasNormals() const { return n.size() == p.size() && !n.empty(); }

        /// Applies  x -> m * x + t  (m must be a rotation / uniform scale; normals are re-normalized).
        void transform(const glm::dmat3& m, const Vector3& t) {
            for (auto& v : p) v = Point3(m * Vector3(v) + t);
            for (auto& v : n) v = glm::normalize(m * Vector3(v));
        }

        /// Sum of triangle areas.
        Real surfaceArea() const {
            Real a = 0;
            for (size_t i = 0; i < numTriangles(); ++i) {
                const Vector3 e1 = p[indices[3 * i + 1]] - p[indices[3 * i]];
                const Vector3 e2 = p[indices[3 * i + 2]] - p[indices[3 * i]];
                a += Real(0.5) * glm::length(glm::cross(e1, e2));
            }
            return a;
        }

        /// Signed volume (positive for a closed mesh with outward counter-clockwise winding).
        Real signedVolume() const {
            Real v = 0;
            for (size_t i = 0; i < numTriangles(); ++i) {
                const Vector3& a = p[indices[3 * i]];
                const Vector3& b = p[indices[3 * i + 1]];
                const Vector3& c = p[indices[3 * i + 2]];
                v += glm::dot(a, glm::cross(b, c)) / Real(6);
            }
            return v;
        }

        /// Icosphere: icosahedron subdivided `subdiv` times (20 * 4^subdiv triangles).
        /// smooth = true stores per-vertex normals (the exact sphere normals), false = flat shading.
        static TriangleMesh Icosphere(int subdiv, Real radius, const Point3& center, bool smooth) {
            const Real t = (Real(1) + std::sqrt(Real(5))) / Real(2);
            std::vector<Vector3> v = {
                {-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0},
                {0, -1, t}, {0, 1, t}, {0, -1, -t}, {0, 1, -t},
                {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1} };
            for (auto& x : v) x = glm::normalize(x);
            std::vector<std::array<int, 3>> f = {
                {0,11,5},{0,5,1},{0,1,7},{0,7,10},{0,10,11},{1,5,9},{5,11,4},{11,10,2},{10,7,6},{7,1,8},
                {3,9,4},{3,4,2},{3,2,6},{3,6,8},{3,8,9},{4,9,5},{2,4,11},{6,2,10},{8,6,7},{9,8,1} };

            for (int s = 0; s < subdiv; ++s) {
                std::map<std::pair<int, int>, int> mid;
                auto midpoint = [&](int a, int b) {
                    const auto key = std::make_pair(std::min(a, b), std::max(a, b));
                    const auto it = mid.find(key);
                    if (it != mid.end()) return it->second;
                    v.push_back(glm::normalize(v[static_cast<size_t>(a)] + v[static_cast<size_t>(b)]));
                    return mid[key] = static_cast<int>(v.size()) - 1;
                };
                std::vector<std::array<int, 3>> next;
                for (const auto& tri : f) {
                    const int a = midpoint(tri[0], tri[1]), b = midpoint(tri[1], tri[2]), c = midpoint(tri[2], tri[0]);
                    next.push_back({ tri[0], a, c });
                    next.push_back({ tri[1], b, a });
                    next.push_back({ tri[2], c, b });
                    next.push_back({ a, b, c });
                }
                f.swap(next);
            }

            TriangleMesh m;
            for (const auto& x : v) {
                m.p.push_back(center + x * radius);
                if (smooth) m.n.push_back(x);
            }
            for (const auto& tri : f) for (int k = 0; k < 3; ++k) m.indices.push_back(tri[static_cast<size_t>(k)]);
            return m;
        }
    };

    /**
     * @brief One triangle of a shared mesh.
     */
    class Triangle : public Shape {
    public:
        Triangle(std::shared_ptr<const TriangleMesh> mesh, size_t index, std::shared_ptr<Material> mat)
            : m_mesh(std::move(mesh)), m_material(std::move(mat))
        {
            m_i0 = m_mesh->indices[3 * index];
            m_i1 = m_mesh->indices[3 * index + 1];
            m_i2 = m_mesh->indices[3 * index + 2];
            // Edges cached in the triangle: intersect() then needs no pointer chasing into the mesh.
            m_p0 = p0();
            m_e1 = p1() - p0();
            m_e2 = p2() - p0();
            const Vector3 e1 = m_e1, e2 = m_e2;
            const Vector3 c = glm::cross(e1, e2);
            m_area = Real(0.5) * glm::length(c);
            m_ng = m_area > Real(0) ? c / (Real(2) * m_area) : Vector3(0, 1, 0);
            // With vertex normals the side they point to is "outside" (PBRT flips ng toward ns).
            if (m_mesh->hasNormals()) {
                const Vector3 ns = m_mesh->n[static_cast<size_t>(m_i0)] + m_mesh->n[static_cast<size_t>(m_i1)] + m_mesh->n[static_cast<size_t>(m_i2)];
                if (glm::dot(m_ng, ns) < Real(0)) m_ng = -m_ng;
            }
        }

        bool hit(const Ray& r, SurfaceInteraction& rec) const override {
            Real t, u, v;
            if (!intersect(r, t, u, v)) return false;
            fill(r, t, u, v, rec);
            return true;
        }

        /**
         * @brief Builds the hit record for a hit that was found elsewhere (Embree, single precision).
         *
         * (t, u, v) are the external distance and barycentric weights of p1 and p2. They are
         * re-evaluated with this triangle's own (Real precision) test when that confirms the hit;
         * otherwise (an edge case where float and double disagree) the given values are used as is,
         * so an external hit is never lost.
         */
        void fillHit(const Ray& r, Real t, Real u, Real v, SurfaceInteraction& rec) const {
            Ray refine = r;
            refine.tMax = t * Real(1.0001) + Real(1e-6);
            Real t2, u2, v2;
            if (intersect(refine, t2, u2, v2)) { t = t2; u = u2; v = v2; }
            fill(r, t, u, v, rec);
        }

        bool occluded(const Ray& r) const override {
            if (!m_castsShadow) return false;
            Real t, u, v;
            return intersect(r, t, u, v);
        }

        /// false: shadow rays pass through (translucent overlay layers must not shade the surface under them).
        void setCastsShadow(bool casts) { m_castsShadow = casts; }

        bool intersectT(const Ray& r, Real& t) const override {
            Real u, v;
            return intersect(r, t, u, v);
        }

        AABB bounds() const override {
            const Vector3 lo = glm::min(p0(), glm::min(p1(), p2()));
            const Vector3 hi = glm::max(p0(), glm::max(p1(), p2()));
            const Vector3 pad(Real(1e-5));   // flat triangles have zero thickness
            return AABB(lo - pad, hi + pad);
        }

        Real area() const override { return m_area; }

        /// Makes texels of `tex` with alpha below `cutoff` transparent for every ray (hit, any-hit, closest).
        void setAlphaTexture(std::shared_ptr<const Texture2D> tex, Real cutoff = Real(0.5)) {
            m_alphaTex = std::move(tex);
            m_alphaCutoff = cutoff;
        }

        void normalBounds(Vector3& w, Real& cosTheta) const override {
            w = m_ng;
            cosTheta = Real(1);
        }

        std::optional<ShapeSample> sampleSurface(const Point2& u) const override {
            if (m_area <= Real(0)) return std::nullopt;
            const Real su = std::sqrt(u.x);
            const Real b0 = Real(1) - su, b1 = u.y * su;   // uniform barycentric coordinates
            ShapeSample s;
            s.si.p = b0 * p0() + b1 * p1() + (Real(1) - b0 - b1) * p2();
            s.si.gn = m_ng;
            s.si.n = m_ng;
            s.si.frontFace = true;
            s.pdf = Real(1) / m_area;
            return s;
        }

        Real pdfSurface(const SurfaceInteraction& ref, const Vector3& wi) const override {
            Ray ray(ref.OffsetRayOrigin(wi), wi);
            SurfaceInteraction isect;
            if (!hit(ray, isect)) return Real(0);
            const Real dist2 = glm::dot(isect.p - ref.p, isect.p - ref.p);
            const Real cosTheta = std::abs(glm::dot(m_ng, wi));
            if (cosTheta <= Real(1e-9)) return Real(0);
            return dist2 / (cosTheta * m_area);
        }

    private:
        /// Fills the hit record for a hit at distance t with barycentric weights (u, v).
        void fill(const Ray& r, Real t, Real u, Real v, SurfaceInteraction& rec) const {
            rec.t = t;
            rec.p = r.at(t);
            rec.setFaceNormal(r.d, m_ng);
            if (m_mesh->hasNormals()) {
                const Real w = Real(1) - u - v;
                const Vector3 ns = glm::normalize(w * m_mesh->n[static_cast<size_t>(m_i0)]
                    + u * m_mesh->n[static_cast<size_t>(m_i1)] + v * m_mesh->n[static_cast<size_t>(m_i2)]);
                rec.n = rec.frontFace ? ns : -ns;
            }
            if (m_mesh->uv.size() == m_mesh->p.size()) rec.uv = uvAt(u, v);
            rec.matPtr = m_material.get();
        }

        /// Interpolated texture coordinates for barycentric weights (u, v) of p1, p2.
        UV uvAt(Real u, Real v) const {
            const Real w = Real(1) - u - v;
            return UV(w * Real(m_mesh->uv[static_cast<size_t>(m_i0)].x) + u * Real(m_mesh->uv[static_cast<size_t>(m_i1)].x) + v * Real(m_mesh->uv[static_cast<size_t>(m_i2)].x),
                w * Real(m_mesh->uv[static_cast<size_t>(m_i0)].y) + u * Real(m_mesh->uv[static_cast<size_t>(m_i1)].y) + v * Real(m_mesh->uv[static_cast<size_t>(m_i2)].y));
        }

        /// Moller-Trumbore. On success t is the distance and (u, v) the barycentric weights of p1, p2.
        bool intersect(const Ray& r, Real& t, Real& u, Real& v) const {
            const Point3& a = m_p0;
            const Vector3& e1 = m_e1;
            const Vector3& e2 = m_e2;
            const Vector3 pvec = glm::cross(r.d, e2);
            const Real det = glm::dot(e1, pvec);
            if (std::abs(det) < Real(1e-14)) return false;   // parallel
            const Real invDet = Real(1) / det;

            const Vector3 tvec = r.o - a;
            u = glm::dot(tvec, pvec) * invDet;
            if (u < Real(0) || u > Real(1)) return false;

            const Vector3 qvec = glm::cross(tvec, e1);
            v = glm::dot(r.d, qvec) * invDet;
            if (v < Real(0) || u + v > Real(1)) return false;

            t = glm::dot(e2, qvec) * invDet;
            if (t < r.tMin || t > r.tMax) return false;

            // Alpha test: texels below the cutoff are holes (hair strands, lashes, lace ...).
            if (m_alphaTex && m_mesh->uv.size() == m_mesh->p.size()) {
                const UV uv = uvAt(u, v);
                if (m_alphaTex->alpha(uv.x, uv.y) < m_alphaCutoff) return false;
            }
            return true;
        }

        const Point3& p0() const { return m_mesh->p[static_cast<size_t>(m_i0)]; }
        const Point3& p1() const { return m_mesh->p[static_cast<size_t>(m_i1)]; }
        const Point3& p2() const { return m_mesh->p[static_cast<size_t>(m_i2)]; }

        std::shared_ptr<const TriangleMesh> m_mesh;
        std::shared_ptr<Material> m_material;
        int m_i0 = 0, m_i1 = 0, m_i2 = 0;
        bool m_castsShadow = true;
        std::shared_ptr<const Texture2D> m_alphaTex;   // optional alpha mask
        Real m_alphaCutoff = Real(0.5);
        Point3 m_p0 = Point3(0.0);
        Vector3 m_e1 = Vector3(0.0), m_e2 = Vector3(0.0);   // p1 - p0, p2 - p0
        Vector3 m_ng = Vector3(0, 1, 0);
        Real m_area = 0;
    };

    /**
     * @brief Wraps every triangle of the mesh in a BVH. `lights` (optional) receives the triangles
     *        as shapes, e.g. to turn an emissive mesh into area lights.
     */
    inline std::shared_ptr<Hittable> BuildMeshBVH(const std::shared_ptr<const TriangleMesh>& mesh,
        const std::shared_ptr<Material>& material,
        std::vector<std::shared_ptr<Triangle>>* triangles = nullptr)
    {
        std::vector<std::shared_ptr<Hittable>> objects;
        objects.reserve(mesh->numTriangles());
        for (size_t i = 0; i < mesh->numTriangles(); ++i) {
            auto t = std::make_shared<Triangle>(mesh, i, material);
            if (triangles) triangles->push_back(t);
            objects.push_back(t);
        }
        if (objects.empty()) return nullptr;
        return std::make_shared<BVH>(std::move(objects));
    }

} // namespace rayt
