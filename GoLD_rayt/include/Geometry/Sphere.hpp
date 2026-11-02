/**
 * @file Geometry/Sphere.hpp
 * @brief Sphere primitive implementation.
 * * This header defines a standard geometric sphere. It uses an algebraic
 * solution for the ray-sphere intersection test, which is both precise
 * and efficient for primary and shadow rays.
 */

#pragma once

#include <cmath>
#include <memory>
#include <vector>

#include "Core/Constants.hpp"
#include "Core/Math.hpp"
#include "Core/Interaction.hpp"
#include "Core/AABB.hpp"
#include "Geometry/Shape.hpp"
#include "Scene/Hittable.hpp"

namespace rayt {

    /**
     * @brief Sphere primitive implementing the Hittable interface.
     * * The sphere is defined by its center and radius. It supports material
     * assignment and provides an AABB for acceleration structures.
     */
    class Sphere : public Shape {
    public:
        /**
         * @brief Constructs a sphere.
         * @param center   Center point of the sphere in world space.
         * @param radius   Radius of the sphere.
         * @param mat      Shared pointer to the material applied to the sphere's surface.
         */
        Sphere(Point3 center, Real radius, std::shared_ptr<Material> mat)
            : m_center(center), m_radius(radius), m_material(mat) {}

        /**
         * @brief Performs a ray-sphere intersection test.
         * * Solves the quadratic equation: |(o + td) - c|^2 = R^2.
         * Uses the optimized "half-b" form of the quadratic formula to reduce
         * floating-point operations.
         * * @param r   The incident ray.
         * @param rec Output structure to store intersection details if a hit occurs.
         * @return True if the ray hits the sphere within the valid interval [tMin, tMax].
         */
        bool hit(const Ray& r, SurfaceInteraction& rec) const override {
            // Vector from sphere center to ray origin
            Vector3 oc = r.o - m_center;

            // Quadratic coefficients (using half-b form)
            // a*t^2 + b*t + c = 0
            Real a = glm::dot(r.d, r.d); 
            if (a == Real(0)) return false;  // Prevent division by zero

            Real half_b = glm::dot(oc, r.d);
            Real c = glm::dot(oc, oc) - m_radius * m_radius;
            
            // Discriminant check (b^2 - 4ac, but scaled for half_b)
            Real discriminant = half_b * half_b - a * c;
            if (discriminant < 0) return false;
            Real sqrtd = std::sqrt(discriminant);

            // Find the nearest root that lies in the acceptable range [tMin, tMax].
            // We check the smaller root (-half_b - sqrtd) first as it represents the closer hit.
            Real root = (-half_b - sqrtd) / a;
            if (root < r.tMin || r.tMax < root) {
                root = (-half_b + sqrtd) / a;
                if (root < r.tMin || r.tMax < root)
                    return false;
            }

            // Populate Interaction record
            rec.t = root;
            rec.p = r.at(rec.t);

            // Calculate the geometric outward normal. 
            // Normalizing by dividing by radius is efficient for unit/uniform spheres.
            Vector3 outwardNormal = (rec.p - m_center) / m_radius;

            // Orient the normal based on whether the ray hit from the outside or inside.
            rec.setFaceNormal(r.d, outwardNormal);

            // Assign the material property
            rec.matPtr = m_material.get();

            // TODO: Implement spherical UV mapping for texture lookup (e.g., lat/long).

            return true;
        }

        /**
         * @brief Returns the Axis-Aligned Bounding Box (AABB) of the sphere.
         * @return AABB spanning from (center - radius) to (center + radius).
         */
        AABB bounds() const override {
            Vector3 rad(m_radius);
            return AABB(m_center - rad, m_center + rad);
        }

        Real area() const override {
            return Real(4) * constants::PI * m_radius * m_radius;
        }

        std::optional<ShapeSample> sampleSurface(const Point2& u) const override {
            Real z = Real(1) - Real(2) * u.x;
            Real r = math::safe_sqrt(std::max(Real(0), Real(1) - z * z));
            Real phi = Real(2) * constants::PI * u.y;

            Vector3 dir(r * std::cos(phi), r * std::sin(phi), z);
            Point3 p = m_center + dir * m_radius;

            SurfaceInteraction si;
            si.p = p;
            si.gn = dir;
            si.n = dir;
            si.frontFace = true;

            ShapeSample sample;
            sample.si = si;
            sample.pdf = math::safe_recip(area());
            return sample;
        }

        /**
         * @brief Cone sampling: a direction is drawn uniformly inside the cone the sphere subtends
         *        at `ref`, so every sample lands on the visible part of the sphere (area sampling
         *        wastes the far half). pdf = 1 / (2 pi (1 - cos(thetaMax))) w.r.t. solid angle.
         *        Falls back to area sampling when `ref` is inside the sphere.
         */
        std::optional<ShapeSample> sampleFrom(const Point3& ref, const Point2& u) const override {
            const Vector3 toCenter = m_center - ref;
            const Real dc2 = glm::dot(toCenter, toCenter);
            const Real r2 = m_radius * m_radius;
            if (dc2 <= r2 * (Real(1) + Real(1e-6))) return Shape::sampleFrom(ref, u);

            const Real dc = std::sqrt(dc2);
            const Vector3 wc = toCenter / dc;
            const Real sin2Max = r2 / dc2;
            const Real cosMax = math::safe_sqrt(Real(1) - sin2Max);
            const Real oneMinusCosMax = sin2Max / (Real(1) + cosMax);   // stable for small spheres

            // Uniform in cos(theta) inside the cone.
            const Real cosT = Real(1) - u.x * oneMinusCosMax;
            const Real sin2T = u.x * oneMinusCosMax * (Real(1) + cosT);
            const Real sinMax = std::sqrt(sin2Max);

            // Angle alpha at the sphere center between wc and the sampled point (law of sines).
            const Real cosA = sin2T / sinMax + cosT * math::safe_sqrt(Real(1) - sin2T / sin2Max);
            const Real sinA = math::safe_sqrt(Real(1) - cosA * cosA);
            const Real phi = Real(2) * constants::PI * u.y;

            // Orthonormal basis around wc.
            const Real sign = std::copysign(Real(1), wc.z);
            const Real a = -Real(1) / (sign + wc.z);
            const Real b = wc.x * wc.y * a;
            const Vector3 bx(Real(1) + sign * wc.x * wc.x * a, sign * b, -sign * wc.x);
            const Vector3 by(b, sign + wc.y * wc.y * a, -wc.y);

            const Vector3 n = -(bx * (sinA * std::cos(phi)) + by * (sinA * std::sin(phi)) + wc * cosA);

            ShapeSample sample;
            sample.si.p = m_center + n * m_radius;
            sample.si.gn = n;
            sample.si.n = n;
            sample.si.frontFace = true;   // the visible side faces `ref`
            sample.pdf = Real(1) / (Real(2) * constants::PI * oneMinusCosMax);
            return sample;
        }

        Real pdfSurface(const SurfaceInteraction& ref, const Vector3& wi) const override {
            // Same cone as sampleFrom(): constant density inside it.
            const Vector3 toCenter = m_center - ref.p;
            const Real dc2 = glm::dot(toCenter, toCenter);
            const Real r2 = m_radius * m_radius;
            if (dc2 > r2 * (Real(1) + Real(1e-6))) {
                const Real sin2Max = r2 / dc2;
                const Real cosMax = math::safe_sqrt(Real(1) - sin2Max);
                const Real cosT = glm::dot(wi, toCenter) / std::sqrt(dc2);
                if (cosT < cosMax) return Real(0);
                return (Real(1) + cosMax) / (Real(2) * constants::PI * sin2Max);
            }

            Ray ray(ref.OffsetRayOrigin(wi), wi);
            SurfaceInteraction isect;
            if (!hit(ray, isect)) {
                return Real(0);
            }

            Real dist2 = glm::dot(isect.p - ref.p, isect.p - ref.p);
            Real cosTheta = std::abs(glm::dot(isect.gn, -wi));
            if (cosTheta <= Real(0)) {
                return Real(0);
            }

            return dist2 / (cosTheta * area());
        }

    private:
        Point3 m_center;
        Real m_radius;
        std::shared_ptr<Material> m_material;
    };

} // namespace rayt