#pragma once

/**
 * @file Sphere.hpp
 * @brief Sphere primitive implementation.
 * * This header defines a standard geometric sphere. It uses an algebraic
 * solution for the ray-sphere intersection test, which is both precise
 * and efficient for primary and shadow rays.
 */

#include "Geometry/Hittable.hpp"
#include "Core/Interaction.hpp"
#include "Core/AABB.hpp"

#include <memory>
#include <vector>

namespace rayt {

    /**
     * @brief Sphere primitive implementing the Hittable interface.
     * * The sphere is defined by its center and radius. It supports material
     * assignment and provides an AABB for acceleration structures.
     */
    class Sphere : public Hittable {
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
        virtual bool hit(const Ray& r, SurfaceInteraction& rec) const override {
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

            // Critical: Update the ray's maximum valid distance. 
            // This ensures subsequent intersection tests in a list or BVH prune farther objects.
            r.tMax = rec.t;

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

    private:
        Point3 m_center;
        Real m_radius;
        std::shared_ptr<Material> m_material;
    };

} // namespace rayt