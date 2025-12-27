#pragma once

/**
 * @file RoughConductor.hpp
 * @brief Rough conductor material using the Microfacet theory.
 * * Implements the Cook-Torrance BRDF model:
 * fr = (D * G * F) / (4 * (n.wi) * (n.wo))
 */

#include <complex>

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

            // Convert roughness to alpha (perceptually linear mapping)
            Real aspect = std::sqrt(1.0 - anisotropy * 0.9);
            alpha_x = MicrofacetDistribution::roughnessToAlpha(roughness / aspect);
            alpha_y = MicrofacetDistribution::roughnessToAlpha(roughness * aspect);
        }

        /**
         * @brief Evaluates the Cook-Torrance BRDF.
         */
        Spectrum eval(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            TransportMode mode) const override {

            Real cosThetaO = std::abs(glm::dot(rec.n, wo)); // n dot wo
            Real cosThetaI = std::abs(glm::dot(rec.n, wi)); // n dot wi

            // Ignore grazing angles or back-facing
            if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0.0);
            if (glm::dot(rec.gn, wi) <= 0 || glm::dot(rec.gn, wo) <= 0) return Spectrum(0.0); // Geometric normal check

            // 1. Half vector
            Vector3 wh = wi + wo;
            if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.0);
            wh = glm::normalize(wh);

            // 2. Local Frame conversion
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
            Spectrum F = fresnel::fresnelConductor(glm::dot(wh, wi), eta, k); // F using wh

            // 5. Cook-Torrance Formula
            return (D * G * F) / (4.0f * cosThetaI * cosThetaO);
        }

        /**
         * @brief Importance samples the GGX distribution (VNDF).
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

            // 2. Sample micro-normal (wh) using VNDF
            Vector3 wh_local = dist.sample_wh(wo_local, u);
            Vector3 wh = frame.localToWorld(wh_local);

            // 3. Reflect wo about wh to get wi
            bsdfSample.wi = glm::reflect(-wo, wh);

            if (glm::dot(rec.gn, bsdfSample.wi) <= 0)
                return std::nullopt;

            Vector3 wi_local = frame.worldToLocal(bsdfSample.wi);

            // 4. Sanity checks (geometric & shading normals)
            if (glm::dot(rec.gn, bsdfSample.wi) <= 0) return std::nullopt; // Below surface
            if (wo_local.z == 0 || wi_local.z == 0) return std::nullopt;

            // 5. Compute PDF
            // Jacobian transformation: dwh / dwi = 1 / (4 * (wo . wh))
            Real dot_wo_wh = glm::dot(wo, wh);
            if (dot_wo_wh <= 0) return std::nullopt;

            Real pdf_wh = dist.pdf(wo_local, wh_local);
            bsdfSample.pdf = pdf_wh / (4.0f * dot_wo_wh);

            // 6. Evaluate Weight (f / pdf) directly to reduce variance
            // f = D * G * F / (4 * cosI * cosO)
            // pdf = D * G1 * (wo.wh) / (wo.z * 4 * (wo.wh)) ... (VNDF pdf simplified)
            //
            // However, to be safe and verify correctness first, we can just call eval() / pdf.
            // But let's calculate F and G explicitly for clarity.

            Spectrum F = fresnel::fresnelConductor(dot_wo_wh, eta, k);
            Real G = dist.G(wo_local, wi_local);
            Real D = dist.D(wh_local); // Usually cancels out in weight, but computed for completeness

            // Correct Cook-Torrance weight with VNDF usually simplifies to: F * G2 / G1
            // Let's stick to the standard definition for now:

            bsdfSample.f = eval(rec, wo, bsdfSample.wi, mode);
            // bsdfSample.sampledType = BxDFType(BSDF_GLOSSY | BSDF_REFLECTION);

            if (bsdfSample.pdf <= 1e-6f || math::hasNaNs(bsdfSample.f)) return std::nullopt;

            // Flags (重要)
            bsdfSample.flags =
                BxDFFlags::Reflection |
                BxDFFlags::Glossy;


            return bsdfSample;
        }

        /**
         * @brief PDF evaluation.
         */
        Real pdf(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi) const override {

            if (glm::dot(rec.gn, wi) <= 0 || glm::dot(rec.gn, wo) <= 0) return 0.0;

            Vector3 wh = glm::normalize(wo + wi);
            if (glm::dot(wh, wh) == 0)
                return 0.0;
            wh = glm::normalize(wh);

            rayt::frame::Frame frame(rec.n);
            Vector3 wo_local = frame.worldToLocal(wo);
            Vector3 wh_local = frame.worldToLocal(wh);

            GGXDistribution dist(alpha_x, alpha_y);
            Real pdf_wh = dist.pdf(wo_local, wh_local);

            // Transform PDF from half-vector to solid angle
            return pdf_wh / (4.0f * std::abs(glm::dot(wo, wh)));
        }

        bool isSpecular() const override { return false; } // It is Glossy, not delta-Specular

    };

} // namespace rayt