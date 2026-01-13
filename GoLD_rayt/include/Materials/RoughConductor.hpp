/**
 * @file Materials/RoughConductor.hpp
 * @brief Rough conductor material based on Microfacet theory (Cook-Torrance).
 * Implements GGX distribution with VNDF sampling and complex Fresnel terms.
 */

#pragma once

#include <complex>
#include <optional>

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Constants.hpp"
#include "Core/Interaction.hpp"
#include "Core/Math.hpp"
#include "Core/Fresnel.hpp"

#include "Geometry/Frame.hpp"
#include "Materials/Material.hpp"
#include "Microfacet/GGX.hpp" 

namespace rayt {

    /**
     * @class RoughConductor
     * @brief Models a rough conductive surface using the Cook-Torrance BRDF.
     * * The BRDF is defined as: fr = (D * G * F) / (4 * (n.wi) * (n.wo)).
     * Features:
     * - Trowbridge-Reitz (GGX) Normal Distribution Function.
     * - Smith Shadowing-Masking function.
     * - Fresnel equations for conductors using complex refractive indices.
     * - Visible Normal Distribution Function (VNDF) importance sampling.
     */
    class RoughConductor : public Material {
        Spectrum eta;   // Index of Refraction (Real)
        Spectrum k;     // Extinction Coefficient
        Real alpha_x;   // Roughness X
        Real alpha_y;   // Roughness Y (same as X if isotropic)

    public:
        /**
         * @brief Constructs a rough conductor.
         * @param eta Real part of IOR.
         * @param k Imaginary part of IOR.
         * @param roughness Roughness value [0, 1].
         * @param anisotropy Anisotropy factor [-1, 1] (0 for isotropic).
         */
        RoughConductor(const Spectrum& eta, const Spectrum& k, Real roughness, Real anisotropy = 0.0)
            : eta(eta), k(k) {

            // Convert perceptual roughness to linear alpha using an anisotropic mapping.
            Real aspect = math::safe_sqrt(1.0 - anisotropy * 0.9);
            alpha_x = MicrofacetDistribution::roughnessToAlpha(roughness / aspect);
            alpha_y = MicrofacetDistribution::roughnessToAlpha(roughness * aspect);
        }

        /**
         * @brief Evaluates the BRDF for a given pair of directions.
         * @param rec Surface interaction details (position, normal, etc.).
         * @param wo Outgoing direction (towards camera).
         * @param wi Incident direction (towards light).
         * @param mode Transport mode (Radiance or Importance).
         * @return The evaluated BRDF spectrum.
         */
        Spectrum eval(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            TransportMode mode) const override {

            // Ensure both directions are in the same hemisphere as the geometric normal.
            if (glm::dot(rec.gn, wo) <= 0 || glm::dot(rec.gn, wi) <= 0) return Spectrum(0.0);

            Real cosThetaO = std::max(Real(0), glm::dot(rec.n, wo));
            Real cosThetaI = std::max(Real(0), glm::dot(rec.n, wi));

            if (cosThetaI <= 0 || cosThetaO <= 0) return Spectrum(0.0);

            // 1. Calculate the half-vector and ensure it's valid.
            Vector3 wh = wi + wo;
            Real whLen2 = glm::dot(wh, wh);
            if (whLen2 <= Real(0)) return Spectrum(0.0);
            wh *= math::safe_recip(math::safe_sqrt(whLen2));

            // 2. Transform world vectors to the local shading frame.
            // The Distribution class assumes operations in Tangent Space.
            // We need to transform World Space vectors to Local Space for D() and G().
            frame::Frame frame(rec.n);
            Vector3 wo_local = frame.worldToLocal(wo);
            Vector3 wi_local = frame.worldToLocal(wi);
            Vector3 wh_local = frame.worldToLocal(wh);

            // 3. Create Distribution (On-the-fly is fine as it's lightweight)
            GGXDistribution dist(alpha_x, alpha_y);

            // 4. Evaluate Terms
            Real D = dist.D(wh_local);
            Real G = dist.G(wo_local, wi_local);

            Real cosThetaD = std::abs(glm::dot(wi, wh));
            Spectrum F = fresnel::fresnelConductor(cosThetaD, eta, k); 

            // 5. Cook-Torrance Formula
            Real denominator = 4.0 * cosThetaI * cosThetaO;
            if (denominator < 1e-6) return Spectrum(0.0);

            return (D * G * F) / denominator;
        }

        /**
         * @brief Importance samples the BRDF using the Visible Normal Distribution Function (VNDF).
         * @param rec Surface interaction details.
         * @param wo Outgoing direction.
         * @param u Uniform random sample in [0, 1]^2.
         * @return A BSDFSample containing the sampled direction, weight, and PDF, or nullopt on failure.
         */
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec,
            const Vector3& wo,
            const Point2& u,
            TransportMode mode) const override {

            // Reflection only
            if (glm::dot(rec.gn, wo) <= 0)
                return std::nullopt;

            BSDFSample bsdfSample;

            // 1. Setup Frame & Distribution
            rayt::frame::Frame frame(rec.n);
            Vector3 wo_local = frame.worldToLocal(wo);
            GGXDistribution dist(alpha_x, alpha_y);

            // 2. Sample a visible microfacet normal (wh) and transform to world space.
            Vector3 wh_local = dist.sample_wh(wo_local, u);
            Vector3 wh = frame.localToWorld(wh_local);

            // 3. Reflect wo about wh to determine the incident direction wi.
            bsdfSample.wi = math::reflectOutward(wo, wh);

            if (glm::dot(rec.gn, bsdfSample.wi) <= 0)
                return std::nullopt;

            Vector3 wi_local = frame.worldToLocal(bsdfSample.wi);

            // 4. Validate vectors against the local shading frame.
            if (wo_local.z == 0 || wi_local.z == 0) return std::nullopt;

            // 5. Compute the PDF of the sampled direction.
            // dwh/dwi Jacobian: 1 / (4 * (wo . wh))
            Real dot_wo_wh = std::abs(glm::dot(wo_local, wh_local));
            if (dot_wo_wh <= 0) return std::nullopt;

            Real pdf_wh = dist.pdf(wo_local, wh_local);
            bsdfSample.pdf = pdf_wh / (Real(4) * dot_wo_wh);

            // 6. Evaluate the sample weight (throughput).
            bsdfSample.f = eval(rec, wo, bsdfSample.wi, mode);

            if (bsdfSample.pdf <= 1e-6 || math::hasNonFinite(bsdfSample.f)) return std::nullopt;

            // Flags 
            bsdfSample.flags =
                BxDFFlags::Reflection |
                BxDFFlags::Glossy;


            return bsdfSample;
        }

        /**
         * @brief Computes the probability density function (PDF) for a given pair of directions.
         * @param rec Surface interaction details.
         * @param wo Outgoing direction.
         * @param wi Incident direction.
         * @return The PDF value with respect to solid angle.
         */
        Real pdf(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi) const override {

            if (glm::dot(rec.gn, wi) <= 0 || glm::dot(rec.gn, wo) <= 0) return 0.0;

            Vector3 wh = wo + wi;
            Real whLen2 = glm::dot(wh, wh);
            if (whLen2 <= Real(0.0)) return Real(0.0);
            wh = wh * math::safe_recip(math::safe_sqrt(whLen2));

            rayt::frame::Frame frame(rec.n);
            Vector3 wo_local = frame.worldToLocal(wo);
            Vector3 wh_local = frame.worldToLocal(wh);

            GGXDistribution dist(alpha_x, alpha_y);
            Real pdf_wh = dist.pdf(wo_local, wh_local);

            // Transform PDF from half-vector to solid angle
            Real dot_wo_wh = std::abs(glm::dot(wo_local, wh_local));
            if (dot_wo_wh <= Real(0.0)) return Real(0.0);
            return pdf_wh / (Real(4) * dot_wo_wh);
        }

        /// @brief Returns false as this is a glossy material, not a perfect specular mirror.
        bool isSpecular() const override { return false; } // It is Glossy, not delta-Specular

    };

} // namespace rayt