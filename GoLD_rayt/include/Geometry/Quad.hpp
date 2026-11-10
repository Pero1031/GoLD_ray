/**
 * @file Geometry/Quad.hpp
 * @brief Parallelogram (rectangle) shape: corner + two edge vectors.
 *
 * The surface is { corner + a * edge1 + b * edge2 : a, b in [0, 1] }.
 * The shape normal is normalize(cross(edge1, edge2)); a one-sided emitter (DiffuseLight) emits
 * toward that side only. Swap the edges to flip it.
 *
 * Sampling: a rectangle is sampled uniformly in solid angle as seen from the shading point
 * (spherical rectangle sampling, pdf 1 / S); a general parallelogram, a tiny solid angle or a
 * point in the quad's plane uses uniform area sampling (pdf 1 / area) converted to solid angle
 * by Shape::sampleFrom().
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>

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

            // Spherical rectangle sampling needs a true rectangle (perpendicular edges).
            m_lx = glm::length(m_e1);
            m_ly = glm::length(m_e2);
            if (m_lx > Real(0) && m_ly > Real(0)) {
                m_x = m_e1 / m_lx;
                m_y = m_e2 / m_ly;
                m_isRect = std::abs(glm::dot(m_x, m_y)) < Real(1e-6);
            }
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
            rec.dpdu = m_lx > Real(0) ? m_x : Vector3(0.0);   // tangent for anisotropic materials: along the first edge
            rec.uv = UV(a, b);                                 // texture coordinates: a along edge1, b along edge2
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

        /**
         * Solid-angle sampling of the rectangle as seen from `ref` (Urena et al., "An Area-Preserving
         * Parametrization for Spherical Rectangles", EGSR 2013). Directions are uniform over the
         * solid angle S the rectangle covers, so pdf = 1 / S and the estimator weight
         * (Le * f * cos) / pdf = Le * f * cos * S stays bounded. Uniform area sampling instead
         * carries a 1 / distance^2 factor that explodes for a point close to a bright light.
         * Falls back to area sampling for a tiny solid angle, a non-rectangular quad, or a
         * reference point in the plane of the quad. pdfSurface() makes the same choice.
         */
        std::optional<ShapeSample> sampleFrom(const Point3& ref, const Point2& u) const override {
            SphericalRect r;
            if (!sphericalRect(ref, r)) return Shape::sampleFrom(ref, u);

            // Pick x on the rectangle (cu = cos of the angle that fixes x), then y given x.
            const Real au = u.x * r.S + r.k;
            const Real fu = (std::cos(au) * r.b0 - r.b1) / std::sin(au);
            Real cu = Real(1) / std::sqrt(fu * fu + r.b0 * r.b0) * (fu > Real(0) ? Real(1) : Real(-1));
            cu = std::clamp(cu, Real(-1), Real(1));
            Real xu = -(cu * r.z0) / std::sqrt(std::max(Real(1) - cu * cu, Real(0)));
            xu = std::clamp(xu, r.x0, r.x1);

            const Real d2 = xu * xu + r.z0 * r.z0;
            const Real h0 = r.y0 / std::sqrt(d2 + r.y0 * r.y0);
            const Real h1 = r.y1 / std::sqrt(d2 + r.y1 * r.y1);
            const Real hv = h0 + u.y * (h1 - h0);
            const Real hv2 = hv * hv;
            Real yv = (hv2 < Real(1) - Real(1e-9)) ? (hv * std::sqrt(d2)) / std::sqrt(Real(1) - hv2) : r.y1;
            yv = std::clamp(yv, r.y0, r.y1);

            ShapeSample s;
            s.si.p = ref + m_x * xu + m_y * yv + r.ez * r.z0;
            s.si.gn = m_n;
            s.si.n = m_n;

            const Vector3 d = s.si.p - ref;
            const Real dist2 = glm::dot(d, d);
            if (dist2 <= Real(0)) return std::nullopt;
            const Real cosL = glm::dot(m_n, -d / std::sqrt(dist2));
            if (std::abs(cosL) <= Real(1e-9)) return std::nullopt;

            s.si.frontFace = cosL > Real(0);
            s.pdf = Real(1) / r.S;
            return s;
        }

        Real pdfSurface(const SurfaceInteraction& ref, const Vector3& wi) const override {
            Ray ray(ref.OffsetRayOrigin(wi), wi);
            SurfaceInteraction isect;
            if (!hit(ray, isect)) return Real(0);

            SphericalRect r;
            if (sphericalRect(ref.p, r)) return Real(1) / r.S;

            const Real dist2 = glm::dot(isect.p - ref.p, isect.p - ref.p);
            const Real cosTheta = std::abs(glm::dot(m_n, wi));
            if (cosTheta <= Real(1e-9)) return Real(0);
            return dist2 / (cosTheta * m_area);
        }

    private:
        /// The rectangle as seen from a reference point, in a local frame centered on that point.
        struct SphericalRect {
            Real x0 = 0, y0 = 0, z0 = 0, x1 = 0, y1 = 0;   ///< corner (x0, y0) and far corner (x1, y1); plane at z0 < 0
            Real b0 = 0, b1 = 0;                           ///< z of two edge-plane normals (see sampleFrom)
            Real k = 0;                                    ///< 2 pi - g2 - g3
            Real S = 0;                                    ///< solid angle of the rectangle (sr)
            Vector3 ez;                                    ///< local z axis in world space, flipped so that z0 < 0
        };

        // Below this solid angle uniform area sampling is just as good, and S = (sum of angles) - 2 pi
        // starts to lose precision (the threshold PBRT uses for spherical sampling).
        static constexpr Real kMinSolidAngle = Real(3e-4);

        /// Fills `r` and returns true when spherical rectangle sampling applies from `ref`.
        bool sphericalRect(const Point3& ref, SphericalRect& r) const {
            if (!m_isRect) return false;

            const Vector3 d = m_corner - ref;
            r.x0 = glm::dot(d, m_x);
            r.y0 = glm::dot(d, m_y);
            r.z0 = glm::dot(d, m_n);
            r.ez = m_n;
            if (r.z0 > Real(0)) { r.z0 = -r.z0; r.ez = -r.ez; }   // look at the rectangle from its -z side
            if (r.z0 > -Real(1e-9)) return false;                  // ref lies in the plane of the quad
            r.x1 = r.x0 + m_lx;
            r.y1 = r.y0 + m_ly;

            // Unit normals of the four planes through ref and an edge; their dihedral angles give S.
            const Vector3 v00(r.x0, r.y0, r.z0), v01(r.x0, r.y1, r.z0), v10(r.x1, r.y0, r.z0), v11(r.x1, r.y1, r.z0);
            const Vector3 n0 = glm::normalize(glm::cross(v00, v10));
            const Vector3 n1 = glm::normalize(glm::cross(v10, v11));
            const Vector3 n2 = glm::normalize(glm::cross(v11, v01));
            const Vector3 n3 = glm::normalize(glm::cross(v01, v00));
            const auto angle = [](const Vector3& a, const Vector3& b) {
                return std::acos(std::clamp(-glm::dot(a, b), Real(-1), Real(1)));
            };
            const Real g0 = angle(n0, n1), g1 = angle(n1, n2), g2 = angle(n2, n3), g3 = angle(n3, n0);

            r.b0 = n0.z;
            r.b1 = n2.z;
            r.k = Real(2) * constants::PI - g2 - g3;
            r.S = g0 + g1 - r.k;
            return r.S >= kMinSolidAngle;
        }

        Point3 m_corner;
        Vector3 m_e1, m_e2;
        Vector3 m_nn;        // cross(e1, e2), not normalized
        Vector3 m_n;         // unit normal
        Real m_area = 0;
        Real m_nnLen2 = 0;
        Vector3 m_x = Vector3(1, 0, 0), m_y = Vector3(0, 1, 0);   // unit edge directions
        Real m_lx = 0, m_ly = 0;                                  // edge lengths
        bool m_isRect = false;
        std::shared_ptr<Material> m_material;
    };

} // namespace rayt
