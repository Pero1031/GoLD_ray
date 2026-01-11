/**
 * @file EnvLight.hpp
 * @brief Infinite/Environment light source based on an environment map.
 *
 * This light wraps an EnvMap (from IO/EnvMap.hpp) and exposes it through
 * the Light interface, enabling importance sampling for Next Event Estimation.
 */

#pragma once

#include <memory>
#include <optional>

#include "Core/Types.hpp"
#include "Core/Ray.hpp"
#include "Lights/Light.hpp"
#include "IO/EnvMap.hpp"

namespace rayt {

    /**
     * @brief Environment light source for Image-Based Lighting (IBL).
     *
     * The EnvLight represents an infinitely distant light source that
     * surrounds the entire scene. It uses an equirectangular environment
     * map for both evaluation and importance sampling.
     *
     * Features:
     * - Importance sampling based on environment map luminance
     * - Support for arbitrary HDR environment maps
     * - Proper solid angle PDF computation
     * - Integration with the Light interface
     */
    class EnvLight final : public Light {
    public:
        /**
         * @brief Construct an environment light from an EnvMap.
         * @param envMap The environment map (equirectangular)
         * @param scale Global intensity multiplier (default: 1.0)
         */
        explicit EnvLight(std::shared_ptr<EnvMap> env, Real scale = 1.0)
            : m_envMap(std::move(env)), m_scale(scale) {}

        /**
         * @brief Get the light type.
         */
        LightType type() const override { return LightType::Infinite; }

        /**
         * @brief Sample the environment light for Next Event Estimation.
         *
         * This method importance samples a direction based on the environment
         * map's luminance distribution and returns the radiance arriving from
         * that direction.
         *
         * @param ctx Reference point context (position, normal, etc.)
         * @param u Random sample in [0,1)^2
         * @param lambda Wavelength samples (unused in RGB mode)
         * @param allowIncompletePDF Allow incomplete PDF (unused)
         * @return Light sample containing direction, radiance, and PDF
         */
        std::optional<LightLiSample>
        sampleLi(const LightSampleContext& ctx,
                 const Point2& u,
                 const SampledWavelengths& /*lambda*/,
                 bool /*allowIncompletePDF*/ = false) const override
        {
            if (!m_envMap) return std::nullopt;

            // Importance sample direction from environment map
            Vector3 wi;
            Real pdfW = 0;
            Vector3 rgb = m_envMap->sample(u, wi, pdfW); 
            if (pdfW <= Real(0)) return std::nullopt;

            LightLiSample ls;
            ls.wi = wi;            // Sampled direction (world space)
            ls.pdf = pdfW;         // PDF w.r.t. solid angle
            ls.isDelta = false;    // Environment lights are not delta

            // For infinite lights, distance is infinite
            ls.tMax = std::numeric_limits<Real>::infinity();

            // Dummy position (not used for infinite lights)
            ls.pLight = ctx.p + wi * ls.tMax; 

            // Radiance from sampled direction
            ls.Li = Spectrum(rgb.x, rgb.y, rgb.z) * m_scale;

            return ls;
        }

        /**
         * @brief Compute the PDF of sampling a specific direction.
         *
         * @param ctx Reference point context
         * @param wi Direction to evaluate (should be normalized)
         * @param allowIncompletePDF Allow incomplete PDF (unused)
         * @return PDF w.r.t. solid angle
         */
        Real pdfLi(const LightSampleContext& /*ctx*/,
            const Vector3& wi,
            bool /*allowIncompletePDF*/ = false) const override
        {
            if (!m_envMap) return Real(0);
            return m_envMap->pdf(wi); // solid angle pdf
        }

        /**
         * @brief Evaluate radiance from a given ray direction.
         *
         * This is called when a camera ray or path misses all geometry
         * and hits the background.
         *
         * @param ray The ray that missed geometry
         * @param lambda Wavelength samples (unused in RGB mode)
         * @return Background radiance in the ray direction
         */
        Spectrum Le(const Ray& ray, const SampledWavelengths& /*lambda*/) const override
        {
            if (!m_envMap) return Spectrum(0.0);
            Vector3 rgb = m_envMap->eval(ray.d);
            return Spectrum(rgb.x, rgb.y, rgb.z) * m_scale;
        }

        /**
         * @brief Sample a ray emitted from the light (for bidirectional methods).
         *
         * @note Not yet implemented for environment lights.
         * This would be used in bidirectional path tracing or light tracing.
         */
        std::optional<LightLeSample>
            sampleLe(const Point2& /*u1*/,
                const Point2& /*u2*/,
                SampledWavelengths& /*lambda*/,
                Real /*time*/) const override
        {
            // TODO: For bidirectional path tracing
            // Would need to:
            // 1. Sample a point on the scene bounding sphere
            // 2. Sample a direction from the environment map
            // 3. Return ray pointing inward
            return std::nullopt;
        }

        /**
         * @brief Compute PDFs for emitted rays (for bidirectional methods).
         *
         * @note Not yet implemented for environment lights.
         */
        void pdfLe(const Ray& /*ray*/, Real* pdfPos, Real* pdfDir) const override
        {
            if (pdfPos) *pdfPos = 0;
            if (pdfDir) *pdfDir = 0;
        }

        // ========== Additional Utilities ==========

        /**
         * @brief Get the underlying environment map.
         */
        const EnvMap* envMap() const {
            return m_envMap.get();
        }

        /**
         * @brief Get the intensity scale factor.
         */
        Real scale() const {
            return m_scale;
        }

        /**
         * @brief Set a new intensity scale factor.
         */
        void setScale(Real scale) {
            m_scale = scale;
        }

    private:
        std::shared_ptr<EnvMap> m_envMap;
        Real m_scale = Real(1);
    };

} // namespace rayt