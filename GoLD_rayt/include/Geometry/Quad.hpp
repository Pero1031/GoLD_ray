/**
 * @file Geometry/Quad.hpp
 * @brief Parallelogram (rectangle) shape: corner + two edge vectors.
 *
 * The surface is { corner + a * edge1 + b * edge2 : a, b in [0, 1] }.
 * The shape normal is normalize(cross(edge1, edge2)); a one-sided emitter (DiffuseLight) emits
 * toward that side only. Swap the edges to flip it.
 *
 * Sampling: uniform in area (pdf 1 / area), converted to solid angle by Shape::sampleFrom().
 */

#pragma once

#include <cmath>
#include <memory>

#include "Core/Constants.hpp"
#include "Core/Math.hpp"
#include "Core/Interaction.hpp"
#include "Core/AABB.hpp"
#include "Geometry/Shape.hpp"
#include "Scene/Hittable.hpp"

namespace rayt {

    class Quad : public Shape {
    public:
        Quad(const Point3& corner, const Vector3& edge1, const Vector3& edge2,
            std::shared_ptr<Material> mat)
            : m_corner(corner), m_e1(edge1), m_e2(edge2), m_material(std::move(mat))
        {
            m_nn = glm::cross(m_e1, m_e2);
            m_area = glm::length(m_nn);
            m_nnLen2 = glm::dot(m_nn, m_nn);
            m_n = m_area > Real(0) ? m_nn / m_area : Vector3(0, 1, 0);
        }

        bool hit(const Ray& r, SurfaceInteraction& rec) const override {
            const Real denom = glm::dot(m_n, r.d);
            if (std::abs(denom) < Real(1e-12)) return false;   // parallel

            const Real t = glm::dot(m_n, m_corner - r.o) / denom;
            if (t < r.tMin || t > r.tMax) return false;

            const Point3 p = r.at(t);
            const Vector3 rel = p - m_corner;
            const Real a = glm::dot(glm::cross(rel, m_e2), m_nn) / m_nnLen2;
            const Real b = glm::dot(glm::cross(m_e1, rel), m_nn) / m_nnLen2;
            if (a < Real(0) || a > Real(1) || b < Real(0) || b > Real(1)) return false;

            rec.t = t;
            rec.p = p;
            rec.setFaceNormal(r.d, m_n);
            rec.matPtr = m_material.get();
            return true;
        }

        AABB bounds() const override {
            const Point3 c[4] = { m_corner, m_corner + m_e1, m_corner + m_e2, m_corner + m_e1 + m_e2 };
            Vector3 lo = c[0], hi = c[0];
            for (int i = 1; i < 4; ++i) { lo = glm::min(lo, c[i]); hi = glm::max(hi, c[i]); }
            const Vector3 pad(Real(1e-4));   // a flat quad has zero thickness along its normal
            return AABB(lo - pad, hi + pad);
        }

        Real area() const override { return m_area; }

        void normalBounds(Vector3& w, Real& cosTheta) const override {
            w = m_n;
            cosTheta = Real(1);   // flat: one normal
        }

        std::optional<ShapeSample> sampleSurface(const Point2& u) const override {
            if (m_area <= Real(0)) return std::nullopt;
            ShapeSample s;
            s.si.p = m_corner + m_e1 * u.x + m_e2 * u.y;
            s.si.gn = m_n;
            s.si.n = m_n;
            s.si.frontFace = true;
            s.pdf = Real(1) / m_area;
            return s;
        }

        Real pdfSurface(const SurfaceInteraction& ref, const Vector3& wi) const override {
            Ray ray(ref.OffsetRayOrigin(wi), wi);
            SurfaceInteraction isect;
            if (!hit(ray, isect)) return Real(0);

            const Real dist2 = glm::dot(isect.p - ref.p, isect.p - ref.p);
            const Real cosTheta = std::abs(glm::dot(m_n, wi));
            if (cosTheta <= Real(1e-9)) return Real(0);
            return dist2 / (cosTheta * m_area);
        }

    private:
        Point3 m_corner;
        Vector3 m_e1, m_e2;
        Vector3 m_nn;        // cross(e1, e2), not normalized
        Vector3 m_n;         // unit normal
        Real m_area = 0;
        Real m_nnLen2 = 0;
        std::shared_ptr<Material> m_material;
    };

} // namespace rayt
