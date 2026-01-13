/**
 * @file Lights/LightSampler.hpp
 * @brief Light sampling strategies for importance sampling in rendering.
 *
 * The LightSampler encapsulates the logic for selecting which light(s) to sample
 * during direct illumination calculations. This separation allows for different
 * sampling strategies (uniform, power-based, spatial, etc.) without modifying
 * the core Scene or Integrator classes.
 */

#pragma once

#include <memory>
#include <vector>
#include <optional>

#include "Core/Types.hpp"
#include "Core/Sampling.hpp"
#include "Lights/Light.hpp"

namespace rayt {

    /**
     * @brief Result of sampling a light source.
     */
    struct LightSampleResult {
        const Light* light = nullptr;   ///< Pointer to the sampled light
        Real pdf = 0;                   ///< PDF of selecting this light
        int index = -1;                 ///< Index in the light array (useful for debugging)
    };

    /**
     * @brief Abstract base class for light sampling strategies.
     *
     * Different implementations can provide various sampling strategies:
     * - UniformLightSampler: Sample all lights with equal probability
     * - PowerLightSampler: Sample based on emitted power (future)
     * - BVHLightSampler: Spatial light sampling using BVH (future)
     */
    class LightSampler {
    public:
        virtual ~LightSampler() = default;

        /**
         * @brief Sample a single light source.
         * @param u Random sample in [0,1)
         * @return LightSampleResult containing the sampled light and its PDF
         */
        virtual std::optional<LightSampleResult> sample(Real u) const = 0;

        /**
         * @brief Compute the PDF of sampling a specific light.
         * @param lightIndex Index of the light in question
         * @return Probability of selecting this light
         */
        virtual Real pdf(int lightIndex) const = 0;

        /**
         * @brief Get total number of lights managed by this sampler.
         */
        virtual size_t numLights() const = 0;

        /**
         * @brief Check if any lights are available.
         */
        virtual bool hasLights() const = 0;
    };

    /**
     * @brief Uniform light sampler - selects lights with equal probability.
     *
     * This is the simplest and most commonly used strategy. Each light has
     * probability 1/N where N is the total number of lights.
     */
    class UniformLightSampler : public LightSampler {
    public:
        /**
         * @brief Construct from a collection of lights.
         * @param lights Vector of light sources to sample from
         */
        explicit UniformLightSampler(const std::vector<std::shared_ptr<Light>>& lights)
            : m_lights(lights) {}

        std::optional<LightSampleResult> sample(Real u) const override {
            if (m_lights.empty()) {
                return std::nullopt;
            }

            // Uniform sampling: select light index uniformly
            int index = std::min(
                static_cast<int>(u * m_lights.size()),
                static_cast<int>(m_lights.size()) - 1
            );

            LightSampleResult result;
            result.light = m_lights[index].get();
            result.pdf = Real(1) / Real(m_lights.size());
            result.index = index;

            return result;
        }

        Real pdf(int lightIndex) const override {
            if (m_lights.empty() || lightIndex < 0 ||
                lightIndex >= static_cast<int>(m_lights.size())) {
                return 0;
            }
            return Real(1) / Real(m_lights.size());
        }

        size_t numLights() const override {
            return m_lights.size();
        }

        bool hasLights() const override {
            return !m_lights.empty();
        }

    private:
        std::vector<std::shared_ptr<Light>> m_lights;
    };

    /**
     * @brief Power-based light sampler (future implementation).
     *
     * Samples lights proportional to their emitted power, which can reduce
     * variance when lights have very different intensities.
     *
     * @note This is a placeholder for future implementation.
     */
    class PowerLightSampler : public LightSampler {
    public:
        explicit PowerLightSampler(const std::vector<std::shared_ptr<Light>>& lights) {
            // TODO: Compute power for each light and build discrete distribution
            // m_distribution = DiscreteDistribution1D(powers);
        }

        std::optional<LightSampleResult> sample(Real u) const override {
            // TODO: Sample from power distribution
            return std::nullopt;
        }

        Real pdf(int lightIndex) const override {
            // TODO: Return power-based PDF
            return 0;
        }

        size_t numLights() const override { return 0; }
        bool hasLights() const override { return false; }

    private:
        // TODO: Store discrete distribution
        // DiscreteDistribution1D m_distribution;
        // std::vector<std::shared_ptr<Light>> m_lights;
    };

} // namespace rayt