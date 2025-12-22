#pragma once

#include "Core/Core.hpp"
#include "Core/Ray.hpp"
#include "Core/Interaction.hpp"
#include "Core/AABB.hpp"

namespace rayt {

    // Abstract base class for any object that a ray can intersect.
    class Hittable {
    public:
        virtual ~Hittable() = default;

        // Determines if a ray hits this object.
        // The valid interval is now strictly controlled by r.tMin and r.tMax.
        // rec: Output structure to store intersection details.
        virtual bool hit(const Ray& r, SurfaceInteraction& rec) const = 0;

        // Abstract base class for BVH nodes/objects.
        // Returns the Axis-Aligned Bounding Box (AABB) of the object.
        virtual AABB bounds() const = 0;
    };

}