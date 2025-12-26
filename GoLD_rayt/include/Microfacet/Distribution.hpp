#pragma once

/**
 * @file Distribution.hpp
 * @brief Abstract interface for Microfacet Distributions (NDF).
 * * Defines the statistical distribution of micro-geometry on a rough surface.
 * Key components:
 * - D(): Distribution of Normals (NDF)
 * - G(): Geometric Shadowing/Masking (Smith Model)
 * - sample_wh(): Importance sampling of the micro-normal
 */

#include "Core/Core.hpp"
#include "Core/Math.hpp"

namespace rayt {

    /**
     * @brief Abstract base class for microfacet distributions.
     * All vectors (wh, wo, wi) are assumed to be in **Local Tangent Space**,
     * where the macroscopic surface normal is (0, 0, 1).
     */
    class MicrofacetDistribution {
    protected:
        // Roughness parameters (linear scale typically mapped to [0,1])
        // alpha = roughness^2 is a common convention (perceptually linear).
        Real alpha_x;
        Real alpha_y;

    public:
        /**
         * @brief Constructs a distribution with anisotropic roughness.
         * @param ax Roughness in the tangent direction (alpha_x).
         * @param ay Roughness in the bitangent direction (alpha_y).
         */
        MicrofacetDistribution(Real ax, Real ay)
            : alpha_x(std::max(Real(1e-4), ax)),
            alpha_y(std::max(Real(1e-4), ay)) {}

        virtual ~MicrofacetDistribution() = default;

        /**
         * @brief The Normal Distribution Function (NDF), D(wh).
         * Describes the statistical concentration of micro-normals 'wh'.
         * @param wh The micro-surface normal (half vector).
         * @return Probability density of micro-facets aligned with wh.
         */
        virtual Real D(const Vector3& wh) const = 0;

        /**
         * @brief The Smith Lambda function (Auxiliary function for G).
         * Represents the area of the micro-surface that is blocked/shadowed per unit projected area.
         * Used to compute G1 and G2.
         * @param w The direction vector (either wo or wi).
         */
        virtual Real lambda(const Vector3& w) const = 0;

        /**
         * @brief The Geometric Shadowing-Masking function, G(wo, wi).
         * Uses the Smith model separation: G(wo, wi) = 1 / (1 + lambda(wo) + lambda(wi)).
         * @param wo Outgoing (view) direction.
         * @param wi Incident (light) direction.
         * @return Probability that a micro-facet is visible from both directions.
         */
        Real G(const Vector3& wo, const Vector3& wi) const {
            return Real(1.0) / (Real(1.0) + lambda(wo) + lambda(wi));
        }

        /**
         * @brief The Geometric attenuation for a single direction, G1(w).
         * G1(w) = 1 / (1 + lambda(w)).
         */
        Real G1(const Vector3& w) const {
            return Real(1.0) / (Real(1.0) + lambda(w));
        }

        /**
         * @brief Importance samples a micro-normal (wh).
         * Recommended to use Visible Normal Distribution Function (VNDF) sampling.
         * @param wo The outgoing (view) direction.
         * @param u  Random numbers [0,1)^2.
         * @return Sampled micro-normal 'wh'.
         */
        virtual Vector3 sample_wh(const Vector3& wo, const Point2& u) const = 0;

        /**
         * @brief Computes the PDF of sampling 'wh' given 'wo'.
         * For VNDF: pdf = G1(wo) * max(0, wo.wh) * D(wh) / wo.z
         * @param wo View direction.
         * @param wh Sampled micro-normal.
         */
        virtual Real pdf(const Vector3& wo, const Vector3& wh) const = 0;

        /**
         * @brief Helper to convert user-friendly "roughness" [0,1] to alpha value.
         * Typically: alpha = roughness * roughness.
         */
        static Real roughnessToAlpha(Real roughness) {
            roughness = std::clamp(roughness, Real(1e-3), Real(1.0));
            return roughness * roughness;
        }
    };

} // namespace rayt