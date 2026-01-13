/**
 * @file Microfacet/GGX.hpp
 * @brief GGX (Trowbridge-Reitz) Microfacet Distribution implementation.
 * * References:
 * - "Microfacet Models for Refraction through Rough Surfaces" (Walter et al. 2007)
 * - "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs" (Heitz 2014)
 * - "Sampling the GGX Distribution of Visible Normals" (Heitz 2018)
 */

#pragma once

#include <cmath>

#include "Core/Math.hpp"
#include "Microfacet/Distribution.hpp"

namespace rayt {

    class GGXDistribution : public MicrofacetDistribution {
    public:
        /**
         * @brief Constructs an anisotropic GGX distribution.
         * @param ax Roughness in X direction (alpha_x).
         * @param ay Roughness in Y direction (alpha_y).
         */
        GGXDistribution(Real ax, Real ay) : MicrofacetDistribution(ax, ay) {}

        /**
         * @brief Evaluates the GGX Normal Distribution Function (NDF).
         * Formula: D(m) = 1 / ( PI * ax * ay * cos^4(theta) * (1 + tan^2(theta) * (cos^2(phi)/ax^2 + sin^2(phi)/ay^2))^2 )
         * Optimized algebraic form is used below.
         */
        Real D(const Vector3& wh) const override {

            if (wh.z <= Real(0.0)) return Real(0.0);

            Real ax2 = math::sqr(alpha_x);
            Real ay2 = math::sqr(alpha_y);

            Real e =
                math::sqr(wh.x) / ax2 +
                math::sqr(wh.y) / ay2 +
                math::sqr(wh.z);

            return Real(1.0) / (constants::PI * alpha_x * alpha_y * math::sqr(e));
        }

        /**
         * @brief Evaluates the Smith Lambda function for GGX.
         * Approximation: lambda(w) = (-1 + sqrt(1 + alpha^2 * tan^2(theta))) / 2
         */
        Real lambda(const Vector3& w) const override {
            Real absCosTheta = std::abs(w.z);
            if (absCosTheta <= Real(0.0))
                return std::numeric_limits<Real>::infinity();

            Real sin2Theta = std::max(Real(0.0), Real(1.0) - math::sqr(absCosTheta));
            if (sin2Theta <= Real(0.0))
                return Real(0.0);

            Real tan2Theta = sin2Theta * math::safe_recip(math::sqr(absCosTheta));

            Real wxy2 = math::sqr(w.x) + math::sqr(w.y);
            Real ax2 = math::sqr(alpha_x);
            Real ay2 = math::sqr(alpha_y);

            // Compute effective alpha for the azimuthal direction of w
            Real alpha2 = (wxy2 > Real(0.0))
                ? (math::sqr(w.x) * ax2 + math::sqr(w.y) * ay2) / wxy2
                : ax2;

            return Real(0.5) * (math::safe_sqrt(Real(1.0) + alpha2 * tan2Theta) - Real(1.0));
        }

        /**
         * @brief Samples the Visible Normal Distribution Function (VNDF).
         * Uses the exact sampling method by Eric Heitz (2018).
         * This technique accounts for masking, significantly reducing variance at grazing angles.
         */
        Vector3 sample_wh(const Vector3& wo, const Point2& u) const override {
            if (wo.z <= Real(0.0)) return Vector3(0, 0, 1);

            // 1. Transform the view vector to the hemisphere configuration
            //    (stretching by roughness)
            Vector3 Vh = glm::normalize(Vector3(alpha_x * wo.x, alpha_y * wo.y, wo.z));

            // 2. Orthonormal basis (with special handling for Vh.z)
            //    We construct a basis (T1, T2, Vh) without using square roots if possible.
            Real lenSq = math::sqr(Vh.x) + math::sqr(Vh.y);
            Vector3 T1 = lenSq > 0 ? Vector3(-Vh.y, Vh.x, 0) * math::safe_recip(math::safe_sqrt(lenSq))
                : Vector3(1, 0, 0);
            Vector3 T2 = glm::cross(Vh, T1);

            // 3. Sample the unit disk (r, phi)
            Real r = math::safe_sqrt(u.x);
            Real phi = Real(2.0) * constants::PI * u.y;
            Real t1 = r * std::cos(phi);
            Real t2 = r * std::sin(phi);

            // 4. Warping to fit the projection of the visible hemisphere
            //    s is a parameter mixing the disk sample and the vertical component.
            Real s = Real(0.5) * (Real(1.0) + Vh.z);
            t2 = (Real(1.0) - s) * math::safe_sqrt(Real(1.0) - math::sqr(t1)) + s * t2;

            // 5. Compute the sampled micro-normal in the stretched space (Nh)
            Vector3 Nh = t1 * T1 + t2 * T2 +
                math::safe_sqrt(std::max(Real(0.0), Real(1.0) - math::sqr(t1) - math::sqr(t2))) * Vh;

            // 6. Transform back to the original space (unstretching)
            //    m = normalize( (x*ax, y*ay, z) ) -- careful with normalization logic
            Vector3 ne = Vector3(
                alpha_x * Nh.x,
                alpha_y * Nh.y,
                std::max(Real(0.0), Nh.z) // Ensure strictly positive Z
            );

            return glm::normalize(ne);
        }

        /**
         * @brief Computes the PDF of the sampled micro-normal 'wh'.
         * pdf = G1(wo) * abs(wo·wh) * D(wh) / wo.z
         * Note: This is the PDF wrt solid angle of wh, NOT wi.
         */
        Real pdf(const Vector3& wo, const Vector3& wh) const override {
            if (wo.z <= Real(0.0) || wh.z <= Real(0.0))
                return Real(0.0);

            return G1(wo) * std::abs(glm::dot(wo, wh)) * D(wh) * math::safe_recip(wo.z);
        }
    };

} // namespace rayt