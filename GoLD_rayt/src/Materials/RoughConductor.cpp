// Member functions of RoughConductor that are defined here (not in Materials/RoughConductor.hpp) so that changing them recompiles one file.
#include "pch.h"

#include "Materials/RoughConductor.hpp"
#include "Core/Sampling.hpp"

namespace rayt {

    Spectrum RoughConductor::eval(const SurfaceInteraction& rec,
        const Vector3& wo, const Vector3& wi,
        const SampledWavelengths& lambda,
        TransportMode mode) const {

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
        Spectrum F = fresnelAt(cosThetaD, lambda);

        // 5. Cook-Torrance Formula
        Real denominator = 4.0 * cosThetaI * cosThetaO;
        if (denominator < 1e-6) return Spectrum(0.0);

        Spectrum f = (D * G * F) / denominator;
        if (m_msActive) f += msBrdf(cosThetaO, cosThetaI, lambda);
        return f;
    }

    std::optional<BSDFSample> RoughConductor::sample(const SurfaceInteraction& rec,
        const Vector3& wo,
        const Point2& u,
        const SampledWavelengths& lambda,
        TransportMode mode) const {

        // Reflection only
        if (glm::dot(rec.gn, wo) <= 0)
            return std::nullopt;

        BSDFSample bsdfSample;

        // 1. Setup Frame & Distribution
        rayt::frame::Frame frame(rec.n);
        Vector3 wo_local = frame.worldToLocal(wo);
        GGXDistribution dist(alpha_x, alpha_y);

        // Choose the lobe: GGX single scattering, or the multiple-scattering lobe (cosine-weighted)
        // with probability 1 - E(mu_o), the energy the GGX lobe cannot account for.
        const Real pMs = m_msActive ? (Real(1) - eAt(wo_local.z)) : Real(0);
        Point2 uLobe = u;
        bool useMs = false;
        if (pMs > Real(0)) {
            if (u.x < pMs) { useMs = true; uLobe.x = u.x / pMs; }
            else { uLobe.x = (u.x - pMs) / (Real(1) - pMs); }
        }

        if (useMs) {
            const Real r = math::safe_sqrt(uLobe.x);
            const Real phi = Real(2) * constants::PI * uLobe.y;
            const Vector3 wi_ms(r * std::cos(phi), r * std::sin(phi), math::safe_sqrt(Real(1) - uLobe.x));
            bsdfSample.wi = frame.localToWorld(wi_ms);
        }
        else {
            // 2. Sample a visible microfacet normal (wh) and transform to world space.
            Vector3 wh_local = dist.sample_wh(wo_local, uLobe);
            Vector3 wh = frame.localToWorld(wh_local);

            // 3. Reflect wo about wh to determine the incident direction wi.
            bsdfSample.wi = math::reflectOutward(wo, wh);
        }

        if (glm::dot(rec.gn, bsdfSample.wi) <= 0)
            return std::nullopt;

        Vector3 wi_local = frame.worldToLocal(bsdfSample.wi);

        // 4. Validate vectors against the local shading frame.
        if (wo_local.z == 0 || wi_local.z == 0) return std::nullopt;

        // 5. The PDF of the sampled direction: the mixture of both lobes (see pdf()).
        bsdfSample.pdf = pdf(rec, wo, bsdfSample.wi, lambda);

        // 6. Evaluate the sample weight (throughput).
        bsdfSample.f = eval(rec, wo, bsdfSample.wi, lambda, mode);

        if (bsdfSample.pdf <= 1e-6 || math::hasNonFinite(bsdfSample.f)) return std::nullopt;

        // Flags 
        bsdfSample.flags =
            BxDFFlags::Reflection |
            BxDFFlags::Glossy;


        return bsdfSample;
    }

    Real RoughConductor::pdf(const SurfaceInteraction& rec,
        const Vector3& wo, const Vector3& wi,
        const SampledWavelengths& lambda) const {

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
        const Real pdfGgx = (dot_wo_wh > Real(0.0)) ? pdf_wh / (Real(4) * dot_wo_wh) : Real(0.0);

        // Mixture with the multiple-scattering lobe (cosine-weighted), chosen with probability 1 - E(mu_o).
        if (!m_msActive) return pdfGgx;
        const Real pMs = Real(1) - eAt(wo_local.z);
        const Real pdfMs = std::max(glm::dot(rec.n, wi), Real(0)) / constants::PI;
        return (Real(1) - pMs) * pdfGgx + pMs * pdfMs;
    }

} // namespace rayt
