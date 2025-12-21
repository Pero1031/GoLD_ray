#pragma once

#include "Geometry/Hittable.hpp"

#include <memory>
#include <vector>

namespace rayt {

    // Scene class holds all information about the environment being rendered.
    // Currently, it just wraps the geometry (aggregate), but in the future,
    // it will also hold the list of Light sources.
    class Scene {
    public:
        // Constructor takes the root object of the scene (usually a HittableList)
        Scene(std::shared_ptr<Hittable> aggregate)
            : m_aggregate(aggregate) {}

        // Forward the intersection request to the aggregate geometry
        bool hit(const Ray& r, SurfaceInteraction& rec) const {
            return m_aggregate->hit(r, rec);
        }

        // Future extension:
        // const std::vector<std::shared_ptr<Light>>& lights() const { return m_lights; }

    private:
        // The root of the geometry hierarchy (BVH or List)
        std::shared_ptr<Hittable> m_aggregate;

        // std::vector<std::shared_ptr<Light>> m_lights;
    };

}