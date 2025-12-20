#pragma once

#include "Hittable.hpp"
#include "Interaction.hpp"

#include <memory>
#include <vector>

namespace rayt {

    // Sphere primitive implementing ray–sphere intersection
    class Sphere : public Hittable {
    public:
        // center : sphere center
        // radius : sphere radius
        // mat    : surface material
        Sphere(Point3 center, Real radius, std::shared_ptr<Material> mat)
            : m_center(center), m_radius(radius), m_material(mat) {}

        // Ray–sphere intersection test
        // If hit, fills SurfaceInteraction and returns true
        virtual bool hit(const Ray& r, SurfaceInteraction& rec) const override {
            // Vector from sphere center to ray origin
            Vector3 oc = r.o - m_center;

            // Quadratic coefficients (using half-b form)
            Real a = glm::length2(r.d); // requires <glm/gtx/norm.hpp> or dot(d, d)
            Real half_b = glm::dot(oc, r.d);
            Real c = glm::length2(oc) - m_radius * m_radius;
            
            // Discriminant check
            Real discriminant = half_b * half_b - a * c;
            if (discriminant < 0) return false;
            Real sqrtd = std::sqrt(discriminant);

            // Find the nearest root that lies in the acceptable range.
            Real root = (-half_b - sqrtd) / a;
            if (root < r.tMin || r.tMax < root) {
                root = (-half_b + sqrtd) / a;
                if (root < r.tMin || r.tMax < root)
                    return false;
            }

            // Populate Interaction record
            rec.t = root;
            rec.p = r.at(rec.t);

            // Geometric outward normal
            Vector3 outwardNormal = (rec.p - m_center) / m_radius;

            // Orient normal against ray direction
            rec.setFaceNormal(r.d, outwardNormal);

            // Assign material
            rec.matPtr = m_material.get();

            // TODO: Calculate UV coordinates for sphere mapping here

            // Update ray's range of validity
            //r.tMax = rec.t;

            return true;
        }

    private:
        Point3 m_center;
        Real m_radius;
        std::shared_ptr<Material> m_material;
    };

}