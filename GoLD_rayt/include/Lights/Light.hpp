/**
* @file Lights/Light.hpp
* Common light interface for sampling emitted radiance in a renderer.
* Designed primarily for Next Event Estimation (NEE), with future support
* for spectral rendering and bidirectional techniques.
*/

#pragma once

#include <optional>

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Interaction.hpp"
#include "Core/Ray.hpp"

namespace rayt {

    class SampledWavelengths;

    /**
     * @brief Classification of light sources.
     *
     * The distinction between delta and non-delta lights is essential
     * for correct PDF handling and Multiple Importance Sampling (MIS).
     */
    enum class LightType {
        DeltaPosition,   ///< Point-like light with delta position (e.g. point light)
        DeltaDirection,  ///< Directional light with delta direction (e.g. sun)
        Area,            ///< Finite area light emitting from a surface
        Infinite         ///< Environment / infinite distant light
    };

    /**
     * @brief Minimal reference context required for light sampling.
     *
     * This structure encapsulates all information needed by a light source
     * to evaluate its contribution toward a given reference point.
     */
    struct LightSampleContext {
        Point3 p;                        ///< Reference position (usually a surface intersection)
        Normal3 gn;                      ///< Geometry normal at the reference point
        const Medium* medium = nullptr;  ///< Medium containing the reference point
        Real time = 0;                   ///< Time for motion blur
    };

    /**
     * @brief Result of sampling incident radiance from a light source.
     *
     * Returned by Light::sampleLi() and consumed by the integrator
     * during Next Event Estimation.
     */
    struct LightLiSample {
        Vector3 wi;        ///< Direction from reference point toward the light (normalized)
        Spectrum Li;       ///< Incident radiance arriving along wi
        Real pdf = 0;      ///< PDF with respect to solid angle at the reference point
        Point3 pLight;     ///< Sampled point on the light (area lights only)
        Real tMax = 0;     ///< Distance to the light sample (for shadow ray max t)
        bool isDelta = false;  ///< True if the light is delta-distributed
    }; 

    /**
     * @brief Result of sampling emitted radiance from a light source.
     *
     * Used by bidirectional algorithms such as BDPT or MLT.
     */
    struct LightLeSample {
        Ray ray;           ///< Emitted ray from the light
        Spectrum Le;       ///< Radiance carried by the emitted ray
        Real pdfPos = 0;   ///< PDF of sampling the emission position
        Real pdfDir = 0;   ///< PDF of sampling the emission direction
        bool isDelta = false; ///< True if the light has delta position
    };

    /**
     * @brief Abstract base class for all light sources.
     *
     * Lights are responsible for sampling incident radiance toward
     * a reference point and providing correct probability densities
     * for Monte Carlo integration.
     */
    class Light {
    public:
        virtual ~Light() = default;

        virtual LightType type() const = 0;

        /**
         * @brief Sample incident radiance arriving at a reference point.
         *
         * This function is the core interface for Next Event Estimation (NEE).
         * @param ctx Reference point context (position, normal, medium, time)
         * @param u 2D uniform random sample
         * @param lambda Sampled wavelengths (spectral rendering)
         * @param allowIncompletePDF Allow approximate or incomplete PDFs
         * @return Optional LightLiSample; std::nullopt if sampling fails
         */
        virtual std::optional<LightLiSample>
            sampleLi(const LightSampleContext& ctx, const Point2& u,
                const SampledWavelengths& lambda,
                bool allowIncompletePDF = false) const = 0;

        /**
         * @brief Probability density of sampling direction wi toward this light.
         *
         * The returned PDF must be expressed with respect to solid angle
         * at the reference point.
         * @param ctx Reference point context
         * @param wi Direction toward the light
         * @param allowIncompletePDF Allow approximate or incomplete PDFs
         */
        virtual Real pdfLi(const LightSampleContext& ctx, const Vector3& wi,
            bool allowIncompletePDF = false) const = 0;

        /**
         * @brief Radiance emitted along a ray that misses all geometry.
         * This is only meaningful for infinite lights (environment maps).
         * For other lights, the default implementation returns zero.
         */
        virtual Spectrum Le(const Ray& ray, const SampledWavelengths& lambda) const {
            return Spectrum(0.0);
        }

        /**
         * @brief Sample a ray emitted from the light source.
         * This interface is intended for bidirectional rendering algorithms.
         */
        virtual std::optional<LightLeSample>
            sampleLe(const Point2& u1, const Point2& u2,
                SampledWavelengths& lambda, Real time) const {
            return std::nullopt;
        }

        /**
         * @brief Compute PDFs for a previously sampled emitted ray.
         *
         * @param ray Emitted ray
         * @param pdfPos Output PDF for position sampling
         * @param pdfDir Output PDF for direction sampling
         */
        virtual void pdfLe(const Ray& ray, Real* pdfPos, Real* pdfDir) const {
            if (pdfPos) *pdfPos = 0;
            if (pdfDir) *pdfDir = 0;
        }

        /**
         * @brief Emitted radiance from a point on an area light surface.
         *
         * This is equivalent to Material::emitted() and is only meaningful
         * for area lights. Other lights return zero by default.
         * @param pLight Surface interaction on the light
         * @param w Outgoing direction from the light surface
         * @param lambda Sampled wavelengths
         */
        virtual Spectrum L(const SurfaceInteraction& pLight,
            const Vector3& w,
            const SampledWavelengths& lambda) const {
            return Spectrum(0.0);
        }
    };

} // namespace rayt