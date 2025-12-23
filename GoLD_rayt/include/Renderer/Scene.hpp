#pragma once

/**
 * @file Scene.hpp
 * @brief Scene class representing the entire environment to be rendered.
 * * The Scene acts as a high-level container that orchestrates geometry,
 * light sources, and environmental properties. It provides a simplified
 * interface for the Integrator to query the world state.
 */

#include "Geometry/Hittable.hpp"

#include <memory>
#include <vector>

namespace rayt {

    /**
     * @brief The Scene class holds all information about the virtual environment.
     * * Currently, it manages the geometric aggregate (e.g., a BVH). In the future,
     * it will be expanded to manage a collection of Light sources for importance
     * sampling and direct lighting calculations.
     */
    class Scene {
    public:
        /**
         * @brief Constructs a scene with a geometric root.
         * @param aggregate The root of the geometry hierarchy (usually a BVHNode or HittableList).
         */
        Scene(std::shared_ptr<Hittable> aggregate)
            : m_aggregate(aggregate) {}

        /**
         * @brief Queries the scene for the closest ray-geometry intersection.
         * * This method delegates the intersection test to the internal aggregate
         * acceleration structure, ensuring efficient spatial queries.
         * * @param r   The incident ray.
         * @param rec Output structure to store the intersection details.
         * @return True if any geometry in the scene was hit within the ray's interval.
         */
        bool hit(const Ray& r, SurfaceInteraction& rec) const {
            return m_aggregate->hit(r, rec);
        }

        // ---------------------------------------------------------------------
        // Future Extensions
        // ---------------------------------------------------------------------
        // TODO: Implement light source management for direct illumination.
        // const std::vector<std::shared_ptr<Light>>& lights() const { return m_lights; }

    private:
        /**
         * @brief The root of the spatial acceleration structure.
         * Typically a BVHNode that provides logarithmic-time intersection tests.
         */
        std::shared_ptr<Hittable> m_aggregate;

        // Future member for light sources
        // std::vector<std::shared_ptr<Light>> m_lights;
    };

} // namespace rayt