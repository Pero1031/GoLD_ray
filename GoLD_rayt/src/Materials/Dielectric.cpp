// Member functions of Dielectric that are defined here (not in Materials/Dielectric.hpp) so that changing them recompiles one file.
#include "pch.h"

#include "Materials/Dielectric.hpp"
#include "Core/Sampling.hpp"

namespace rayt {

    Spectrum Dielectric::eval(const SurfaceInteraction& rec,
        const Vector3& wo, const Vector3& wi,
        const SampledWavelengths& lambda,
        TransportMode mode) const {

        // Delta lobes: eval() returns 0 (Dirac delta)
        if (smoothAt(rec)) return Spectrum(0.0);

        const Real ior = iorAt(lambda);

        Real cosThetaO = glm::dot(rec.n, wo);
        Real cosThetaI = glm::dot(rec.n, wi);

        // Determine reflection vs transmission in shading-normal hemisphere sense
        bool isReflection = cosThetaO * cosThetaI > 0;

        // For reflection, both must be on the geometric front side
        // For transmission, allow opposite sides (don't kill it here)
        if (isReflection) {
            if (glm::dot(rec.gn, wo) <= 0 || glm::dot(rec.gn, wi) <= 0)
                return Spectrum(0.0);
        }
        else {
            if (glm::dot(rec.gn, wo) * glm::dot(rec.gn, wi) > 0) {
                return Spectrum(0.0);
            }
        }

        // PBRT logic for eta/etap
        bool entering = rec.frontFace;
        Real eta = entering ? (1.0 / ior) : ior;          // etaI / etaT  (for refractOutward)
        Real etap = entering ? ior : (1.0 / ior);          // etaT / etaI  (for half-vector reconstruction)

        // Compute half-vector (wh)
        Vector3 wh;
        if (isReflection) {
            wh = wo + wi;
        }
        else {
            wh = wo + wi * etap;
        }
        if (glm::length2(wh) == 0) return Spectrum(0.0);
        wh = glm::normalize(wh);

        if (glm::dot(wh, wo) < 0) wh = -wh;

        // Local frame conversion for GGX
        rayt::frame::Frame frame(rec.n);
        Vector3 wo_local = frame.worldToLocal(wo);
        Vector3 wi_local = frame.worldToLocal(wi);
        Vector3 wh_local = frame.worldToLocal(wh);

        Real ax, ay;
        alphasAt(rec, ax, ay);
        GGXDistribution dist(ax, ay);
        Real D = dist.D(wh_local);
        Real G = dist.G(wo_local, wi_local);
        Real F = entering
            ? rayt::fresnel::fresnelDielectric(std::abs(glm::dot(wo, wh)), 1.0, ior)
            : rayt::fresnel::fresnelDielectric(std::abs(glm::dot(wo, wh)), ior, 1.0);

        if (isReflection) {
            Real denom = std::abs(4.0 * cosThetaI * cosThetaO);
            if (denom < 1e-8) return Spectrum(0.0);
            return Spectrum(D * G * F / denom);
        }
        else {
            // Refraction BSDF (Walter 2007 / PBRT v4)
            Real dot_wi_wh = glm::dot(wi, wh);
            Real dot_wo_wh = glm::dot(wo, wh);

            Real sqrtDenom = dot_wi_wh * etap + dot_wo_wh;
            Real denom = rayt::math::sqr(sqrtDenom) * cosThetaI * cosThetaO; // signed

            if (std::abs(denom) < 1e-8) return Spectrum(0.0);

            // NOTE: `denom` already carries an implicit etap^2 factor because
            // wh is reconstructed as (wo + wi * etap) instead of PBRT's
            // (wo / etap + wi). An extra `/ etap^2` here would double-apply
            // the radiance non-symmetry correction, so none is needed.
            Real val = D * G * (1.0 - F) *
                std::abs(dot_wi_wh * dot_wo_wh / denom);

            return Spectrum(val);
        }
    }

    std::optional<BSDFSample> Dielectric::sample(const SurfaceInteraction& rec,
        const Vector3& wo,
        const Point2& u,
        const SampledWavelengths& lambda,
        TransportMode mode) const {

        BSDFSample bsdfSample;

        const Real ior = iorAt(lambda);
        Real cosThetaO = glm::dot(rec.n, wo);
        bool entering = rec.frontFace;

        Real eta = entering ? (1.0 / ior) : ior;          // for refractOutward
        Real etap = entering ? ior : (1.0 / ior);          // for scaling in BTDF
        // refractOutward() needs a normal on the same side as wo (dot(wo, n) > 0). rec.n already
        // faces the incoming ray for entering AND exiting rays, so it must NOT be flipped here:
        // flipping it for exiting rays bent the light back into the glass (it never came out).
        Vector3 n_eff = rec.n;

        // --- Case A: Smooth (delta) ---
        if (smoothAt(rec)) {
            Real F = entering
                ? rayt::fresnel::fresnelDielectric(std::abs(cosThetaO), 1.0, ior)
                : rayt::fresnel::fresnelDielectric(std::abs(cosThetaO), ior, 1.0);

            if (u.x < F) {
                // Specular reflection
                bsdfSample.wi = rayt::math::reflectOutward(wo, n_eff);
                bsdfSample.pdf = F;

                // IMPORTANT:
                // Integrator's specular branch uses: beta *= f;
                // With Fresnel-selection sampling, throughput multiplier is 1 (F cancels).
                bsdfSample.f = Spectrum(1.0);

                bsdfSample.flags = BxDFFlags::Specular | BxDFFlags::Reflection;
            }
            else {
                // Specular transmission
                Vector3 wi;
                bool tir = !rayt::math::refractOutward(wo, n_eff, eta, wi);
                if (tir) return std::nullopt;

                bsdfSample.wi = wi;
                bsdfSample.pdf = 1.0 - F;

                Spectrum ft(1.0);
                if (mode == TransportMode::Radiance)
                    ft /= rayt::math::sqr(etap);
                bsdfSample.f = ft;

                bsdfSample.flags = BxDFFlags::Specular | BxDFFlags::Transmission;
            }
            return bsdfSample;
        }

        // --- Case B: Rough (microfacet) ---
        rayt::frame::Frame frame(rec.n);
        Vector3 wo_local = frame.worldToLocal(wo);
        Real ax, ay;
        alphasAt(rec, ax, ay);
        GGXDistribution dist(ax, ay);

        Vector3 wo_sampling = wo_local;
        if (wo_sampling.z < 0) wo_sampling = -wo_sampling;

        // 1) sample wh (uses u)
        Vector3 wh_local = dist.sample_wh(wo_sampling, u);
        Vector3 wh = frame.localToWorld(wh_local);
        if (glm::dot(wh, wo) < 0) {
            wh = -wh;
            wh_local = -wh_local;
        }

        // 2) Fresnel
        Real dot_wo_wh = glm::dot(wo, wh);
        if (dot_wo_wh == 0) return std::nullopt;

        Real F = entering
            ? rayt::fresnel::fresnelDielectric(std::abs(dot_wo_wh), 1.0, ior)
            : rayt::fresnel::fresnelDielectric(std::abs(dot_wo_wh), ior, 1.0);

        // IMPORTANT:
        // Don't reuse u.x (already used in VNDF sampling) for lobe selection.
        // Ideally sample should take Point3. Minimal fix: draw one extra RNG here.
        Real uLobe = rayt::sampling::Random();

        // 3) choose reflection vs transmission
        if (uLobe < F) {
            // --- Reflection ---
            bsdfSample.wi = rayt::math::reflectIncident(-wo, wh);

            // For reflection, must stay on geometric front side
            if (glm::dot(rec.gn, wo) <= 0 || glm::dot(rec.gn, bsdfSample.wi) <= 0)
                return std::nullopt;

            bsdfSample.flags = BxDFFlags::Glossy | BxDFFlags::Reflection;

            // f via eval() to keep eval/pdf/sample consistent
            bsdfSample.f = eval(rec, wo, bsdfSample.wi, lambda, mode);

            // PDF
            Real pdf_wh = dist.pdf(wo_sampling, wh_local);
            bsdfSample.pdf = (pdf_wh / (4.0 * std::abs(dot_wo_wh))) * F;
        }
        else {
            // --- Transmission ---
            Vector3 wi;
            bool tir = !rayt::math::refractOutward(wo, wh, eta, wi);
            if (tir) return std::nullopt;

            // transmission must go to opposite side in shading-normal sense
            if (glm::dot(rec.n, wi) * cosThetaO > 0) return std::nullopt;

            bsdfSample.wi = wi;
            bsdfSample.flags = BxDFFlags::Glossy | BxDFFlags::Transmission;

            // f via eval() to keep consistent
            bsdfSample.f = eval(rec, wo, wi, lambda, mode);

            // PDF (PBRT v4 Jacobian)
            Real dot_wi_wh = glm::dot(wi, wh);
            Real sqrtDenom = dot_wi_wh * etap + dot_wo_wh;
            if (sqrtDenom == 0) return std::nullopt;

            Real dwh_dwi = std::abs(dot_wi_wh) * rayt::math::sqr(etap) / rayt::math::sqr(sqrtDenom);
            Real pdf_wh = dist.pdf(wo_sampling, wh_local);
            bsdfSample.pdf = pdf_wh * dwh_dwi * (1.0 - F);
        }

        if (bsdfSample.pdf < 1e-8 || rayt::math::hasNonFinite(bsdfSample.f))
            return std::nullopt;

        return bsdfSample;
    }

    Real Dielectric::pdf(const SurfaceInteraction& rec,
        const Vector3& wo, const Vector3& wi,
        const SampledWavelengths& lambda) const {

        if (smoothAt(rec)) return 0.0;

        const Real ior = iorAt(lambda);

        Real cosThetaO = glm::dot(rec.n, wo);
        Real cosThetaI = glm::dot(rec.n, wi);

        bool isReflection = cosThetaO * cosThetaI > 0;

        bool entering = rec.frontFace;
        Real etap = entering ? ior : (1.0 / ior);

        // Reconstruct wh correctly depending on reflection/transmission
        Vector3 wh;
        if (isReflection) {
            wh = wo + wi;
        }
        else {
            wh = wo + wi * etap;
        }
        if (glm::length2(wh) == 0) return 0.0;
        wh = glm::normalize(wh);
        if (glm::dot(wh, wo) < 0) wh = -wh;

        Real dot_wo_wh = glm::dot(wo, wh);
        if (dot_wo_wh == 0) return 0.0;

        Real F = entering
            ? rayt::fresnel::fresnelDielectric(std::abs(dot_wo_wh), 1.0, ior)
            : rayt::fresnel::fresnelDielectric(std::abs(dot_wo_wh), ior, 1.0);

        rayt::frame::Frame frame(rec.n);
        Real ax, ay;
        alphasAt(rec, ax, ay);
        GGXDistribution dist(ax, ay);

        Vector3 wo_local = frame.worldToLocal(wo);
        if (wo_local.z < 0) wo_local = -wo_local;

        Real pdf_wh = dist.pdf(wo_local, frame.worldToLocal(wh));

        if (isReflection) {
            return (pdf_wh / (4.0 * std::abs(dot_wo_wh))) * F;
        }
        else {
            Real dot_wi_wh = glm::dot(wi, wh);
            Real sqrtDenom = dot_wi_wh * etap + dot_wo_wh;
            if (sqrtDenom == 0) return 0.0;

            Real dwh_dwi = std::abs(dot_wi_wh) * rayt::math::sqr(etap) / rayt::math::sqr(sqrtDenom);
            return (pdf_wh * dwh_dwi) * (1.0 - F);
        }
    }

} // namespace rayt
