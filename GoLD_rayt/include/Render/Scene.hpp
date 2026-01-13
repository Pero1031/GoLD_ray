/**
 * @file Renderer/Scene.hpp
 * @brief Scene class representing the entire environment to be rendered.
 * * The Scene acts as a high-level container that orchestrates geometry,
 * light sources, and environmental properties. It provides a simplified
 * interface for the Integrator to query the world state.
 */

#pragma once

#include <memory>
#include <vector>
#include <string>

#include "Geometry/Hittable.hpp"
#include "Lights/Light.hpp"
#include "Lights/LightSampler.hpp"
#include "Lights/EnvLight.hpp"

namespace rayt {

    /**
     * @brief The Scene class holds all information about the virtual environment.
     *
     * The Scene manages:
     * - Geometric aggregate (BVH or similar acceleration structure)
     * - Collection of light sources
     * - Light sampling strategy via LightSampler
     * - Optional infinite/environment light
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

        /**
         * @brief Add a light source to the scene.
         * @param light The light to add
         *
         * @note After adding all lights, you must call finalize() to build
         *       the light sampler.
         */
        void addLight(std::shared_ptr<Light> l) { m_lights.push_back(l); }

        /**
         * @brief Set the infinite/environment light.
         * @param envLight The environment light (can be nullptr)
         */
        void setEnvLight(std::shared_ptr<EnvLight> envLight) {
            m_envLight = envLight;
        }

        /**
         * @brief Finalize the scene after all lights have been added.
         *
         * This builds the light sampler and prepares the scene for rendering.
         * Must be called before rendering begins.
         *
         * @param samplerType Type of light sampler to use (default: uniform)
         */
        void finalize(const std::string& samplerType = "uniform") {
            if (samplerType == "uniform") {
                m_lightSampler = std::make_unique<UniformLightSampler>(m_lights);
            }
            else if (samplerType == "power") {
                // Future: power-based sampling
                m_lightSampler = std::make_unique<PowerLightSampler>(m_lights);
            }
            else {
                // Default to uniform
                m_lightSampler = std::make_unique<UniformLightSampler>(m_lights);
            }
        }

        // ========== Light Query Interface ==========

        /**
         * @brief Check if the scene has any lights (excluding environment).
         */
        bool hasLights() const {
            return m_lightSampler && m_lightSampler->hasLights();
        }

        /**
         * @brief Get the number of lights (excluding environment).
         */
        size_t numLights() const {
            return m_lightSampler ? m_lightSampler->numLights() : 0;
        }

        /**
         * @brief Check if an environment light is present.
         */
        bool hasEnvLight() const {
            return m_envLight != nullptr;
        }

        /**
         * @brief Get the environment light.
         */
        const EnvLight* envLight() const {
            return m_envLight.get();
        }

        /**
         * @brief Access the light sampler directly.
         */
        const LightSampler* lightSampler() const {
            return m_lightSampler.get();
        }

        // ========== Direct Light Access (for debugging) ==========

        /**
         * @brief Get all lights in the scene.
         * @note Prefer using the LightSampler interface for rendering.
         */
        const std::vector<std::shared_ptr<Light>>& lights() const {
            return m_lights;
        }

    private:
        /**
         * @brief The root of the spatial acceleration structure.
         * Typically a BVHNode that provides logarithmic-time intersection tests.
         */
        std::shared_ptr<Hittable> m_aggregate;

        // member for light sources
        std::vector<std::shared_ptr<Light>> m_lights;

        // Optional environment/infinite light
        std::shared_ptr<EnvLight> m_envLight;

        // Light sampling strategy
        std::unique_ptr<LightSampler> m_lightSampler;
    };

} // namespace rayt