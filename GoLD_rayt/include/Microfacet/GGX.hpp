/**
 * @file Microfacet/GGX.hpp
 * @brief GGX (Trowbridge-Reitz) Microfacet Distribution declaration.
 * * References:
 * - "Microfacet Models for Refraction through Rough Surfaces" (Walter et al. 2007)
 * - "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs" (Heitz 2014)
 * - "Sampling the GGX Distribution of Visible Normals" (Heitz 2018)
 */

#pragma once

#include "Microfacet/Distribution.hpp"

namespace rayt {

    /**
     * @brief GGX (Trowbridge-Reitz) microfacet distribution.
     * * Implements NDF D(wh), Smith lambda(w), VNDF sampling sample_wh(), and its pdf().
     * * All vectors are assumed to be in local tangent space (normal = (0,0,1)).
     */
    class GGXDistribution : public MicrofacetDistribution {
    public:
        /**
         * @brief Constructs an anisotropic GGX distribution.
         * @param ax Roughness in X direction (alpha_x).
         * @param ay Roughness in Y direction (alpha_y).
         */
        GGXDistribution(Real ax, Real ay);

        /**
         * @brief Evaluates the GGX Normal Distribution Function (NDF).
         * Formula: D(m) = 1 / ( PI * ax * ay * cos^4(theta) * (1 + tan^2(theta) * (cos^2(phi)/ax^2 + sin^2(phi)/ay^2))^2 )
         * Optimized algebraic form is used in the implementation.
         */
        Real D(const Vector3& wh) const override;

        /**
         * @brief Evaluates the Smith Lambda function for GGX.
         * Approximation: lambda(w) = (-1 + sqrt(1 + alpha^2 * tan^2(theta))) / 2
         */
        Real lambda(const Vector3& w) const override;

        /**
         * @brief Samples the Visible Normal Distribution Function (VNDF).
         * Uses the exact sampling method by Eric Heitz (2018).
         * This technique accounts for masking, significantly reducing variance at grazing angles.
         */
        Vector3 sample_wh(const Vector3& wo, const Point2& u) const override;

        /**
         * @brief Computes the PDF of the sampled micro-normal 'wh'.
         * pdf = G1(wo) * abs(wo·wh) * D(wh) / wo.z
         * Note: This is the PDF wrt solid angle of wh, NOT wi.
         */
        Real pdf(const Vector3& wo, const Vector3& wh) const override;
    };

} // namespace rayt