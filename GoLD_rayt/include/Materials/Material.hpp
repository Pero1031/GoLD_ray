#pragma once

/**
 * @file Material.hpp
 * @brief Abstract interface for physical materials and BxDF functions.
 * * This header defines the Material base class, which encapsulates the
 * light-scattering properties of a surface (BSDF). It supports evaluation,
 * importance sampling, and PDF calculation required for robust Monte Carlo
 * path tracing.
 */

#include "Core/Core.hpp"
#include "Core/Interaction.hpp"
#include "Core/Ray.hpp"

#include <optional>

// Include your spectral data handler
// #include "IORInterpolator.hpp" 

namespace rayt {

    /**
     * @brief Flags representing the properties of a BxDF (BSDF/BRDF/BTDF).
     * * Used to categorize light scattering behavior, which is essential for
     * Multiple Importance Sampling (MIS) and handling delta distributions
     * (perfect mirrors/glass).
     */
    enum BxDFType {
        BSDF_REFLECTION = 1 << 0,
        BSDF_TRANSMISSION = 1 << 1,
        BSDF_DIFFUSE = 1 << 2,
        BSDF_GLOSSY = 1 << 3,
        BSDF_SPECULAR = 1 << 4,
        BSDF_ALL = 0xFF
    };

    /**
     * @brief Light transport mode for asymmetric scattering.
     * * Distinguishes between light moving from sources (Radiance) or from
     * the camera (Importance). Necessary for handling refraction correctly
     * in bidirectional methods.
     */
    enum class TransportMode {
        Radiance,
        Importance
    };

    /**
     * @brief Result of a BSDF sampling operation.
     */
    struct BSDFSample {
        Spectrum f;           // The evaluated BSDF value (throughput)
        Vector3 wi;           // The sampled incident direction (World Space).
        Real pdf = 0;         // The probability density for the sampled direction.
        BxDFType sampledType; // The type of lobe that was sampled.

        /**
         * @brief Checks if the sampled interaction is a perfectly specular (delta) scattering.
         */
        bool isSpecular() const { return sampledType & BSDF_SPECULAR; }
    };

    /**
     * @brief Abstract base class for all materials.
     * * Materials define how light interacts with a surface by providing
     * a Bidirectional Scattering Distribution Function (BSDF).
     */
    class Material {
    public:
        virtual ~Material() = default;

        /**
         * @brief Evaluates the BSDF value f(wo, wi) for a given pair of directions.
         * * @param rec  The surface interaction details (normal, tangents, UVs).
         * @param wo   The outgoing (view) direction in World Space.
         * @param wi   The incident (light) direction in World Space.
         * @param mode The transport mode (usually Radiance).
         * @return Spectrum The BSDF value [1/sr].
         * * @note Following standard PBR conventions, this returns the pure BSDF value.
         * The cosine term (n·wi) is typically applied by the Integrator.
         * For delta distributions (perfect Specular), this returns 0.
         */
        virtual Spectrum eval(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            TransportMode mode = TransportMode::Radiance) const = 0;

        /**
         * @brief Importance samples a new incident direction wi based on the BSDF.
         * * This method is crucial for variance reduction in path tracing.
         * @param rec    The surface interaction details.
         * @param wo     The outgoing (view) direction in World Space.
         * @param u      A 2D random sample from the RNG [0, 1)^2.
         * @param mode   The transport mode.
         * @return std::optional<BSDFSample> The result of the sampling,
         * or nullopt if sampling fails (e.g., Total Internal Reflection).
         */
        virtual std::optional<BSDFSample> sample(const SurfaceInteraction& rec,
            const Vector3& wo,
            const Point2& u, // ランダムシード
            TransportMode mode = TransportMode::Radiance) const = 0;

        /**
         * @brief Evaluates the Probability Density Function (PDF) for a given direction.
         * * Essential for Multiple Importance Sampling (MIS).
         * @param rec The surface interaction details.
         * @param wo  The outgoing direction.
         * @param wi  The incident direction to evaluate.
         * @return Real The PDF value with respect to solid angle.
         */
        virtual Real pdf(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi) const = 0;

        // -----------------------------------------------------------
        // 以下のヘルパーや、以前の簡易的な関数は必要に応じて残すか削除します
        // -----------------------------------------------------------

        /**
         * @brief Returns the emitted radiance for area lights.
         * @param rec The surface interaction details.
         * @param wo  The outgoing direction.
         * @return Spectrum The self-emitted radiance (default: black).
         */
        virtual Spectrum emitted(const SurfaceInteraction& rec,
            const Vector3& wo) const {
            return Spectrum(0.0);
        }

        /**
         * @brief Optimization hint to check if the material is perfectly specular.
         */
        virtual bool isSpecular() const { return false; }
    };

} // namespace rayt