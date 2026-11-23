// Member functions of RoughConductor that are defined here (not in Materials/RoughConductor.hpp) so that changing them recompiles one file.
#include "pch.h"

#include "Materials/RoughConductor.hpp"
#include "Core/Sampling.hpp"

namespace rayt {

    frame::Frame RoughConductor::shadingFrame(const SurfaceInteraction& rec) const {
        frame::Frame f;
        const bool aniso = alpha_x != alpha_y;
        if (aniso && glm::dot(rec.dpdu, rec.dpdu) > Real(0)) f.buildFromNormalAndTangent(rec.n, rec.dpdu);
        else f.buildFromW(rec.n);
        if (aniso && m_rotation != Real(0)) {
            const Real c = std::cos(m_rotation), s = std::sin(m_rotation);
            const Vector3 s0 = f.s, t0 = f.t;
            f.s = c * s0 + s * t0;
            f.t = -s * s0 + c * t0;
        }
        return f;
    }

    bool RoughConductor::glintPlanes(const SurfaceInteraction& rec, const Vector3& wo, GlintPlanes& planes) const {
        if (!m_glint || !(rec.pixelSpacing > Real(0))) return false;
        // The pixel's filter on the surface: a disc across the ray (radius = the filter's standard deviation), projected along the ray
        // onto the tangent plane: an ellipse with its long axis in the plane of incidence.
        const Vector3 n = rec.gn;
        // the pixel's filter, widened by what the path's earlier bounces blur
        const Real r = std::sqrt(math::sqr(m_glintFilter * rec.pixelSpacing) + math::sqr(rec.pathBlur));
        Vector3 t1 = glm::cross(n, wo);
        const Real l1 = glm::length(t1);
        if (l1 > Real(1e-6)) t1 /= l1;
        else t1 = glm::normalize(glm::cross(n, std::abs(n.x) < Real(0.9) ? Vector3(1, 0, 0) : Vector3(0, 1, 0)));
        const Vector3 t2 = glm::cross(t1, n);
        const Vector3 a = t1 * r, b = t2 * (r / std::max(std::abs(glm::dot(n, wo)), Real(0.02)));
        // Three arrangements of facets, one on each coordinate plane of the world, mixed by a roulette with the weights |n_k| / sum
        // (the paper's extension to surfaces without texture coordinates). A plane the surface is nearly perpendicular to is left out.
        Vector3 w(std::abs(n.x), std::abs(n.y), std::abs(n.z));
        for (int k = 0; k < 3; ++k) if (w[k] < Real(0.02)) w[k] = 0;
        w /= w.x + w.y + w.z;
        Real selectSum = 0;
        for (int k = 0; k < 3; ++k) {
            if (!(w[k] > Real(0))) continue;
            const int i = (k + 1) % 3, j = (k + 2) % 3, c = planes.count++;
            planes.stream[c] = k + 1;
            planes.weight[c] = w[k];
            planes.uv[c] = Point2(rec.p[i], rec.p[j]);
            planes.fp[c] = GlintFootprint{ a[i] * a[i] + b[i] * b[i], a[i] * a[j] + b[i] * b[j], a[j] * a[j] + b[j] * b[j] };
            planes.select[c] = m_glint->samplable(planes.fp[c]) ? w[k] : Real(0);
            selectSum += planes.select[c];
        }
        planes.samplable = selectSum > Real(0);
        if (planes.samplable) for (int c = 0; c < planes.count; ++c) planes.select[c] /= selectSum;
        return true;
    }

    Real RoughConductor::glintRatio(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& whLocal) const {
        GlintPlanes planes;
        if (!glintPlanes(rec, wo, planes)) return 1;
        Real ratio = 0;
        for (int c = 0; c < planes.count; ++c)
            ratio += m_glint->ratio(whLocal, planes.uv[c], planes.fp[c], planes.weight[c], static_cast<std::uint32_t>(planes.stream[c]));
        return ratio;
    }

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
        const frame::Frame frame = shadingFrame(rec);
        Vector3 wo_local = frame.worldToLocal(wo);
        Vector3 wi_local = frame.worldToLocal(wi);
        Vector3 wh_local = frame.worldToLocal(wh);

        // 3. Create Distribution (On-the-fly is fine as it's lightweight)
        GGXDistribution dist(alpha_x, alpha_y);

        // 4. Evaluate Terms
        Real D = dist.D(wh_local);
        if (m_glint) D *= glintRatio(rec, wo, wh_local);
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
        const frame::Frame frame = shadingFrame(rec);
        Vector3 wo_local = frame.worldToLocal(wo);
        GGXDistribution dist(alpha_x, alpha_y);

        // Choose the lobe: GGX single scattering, or the multiple-scattering lobe (cosine-weighted)
        // with probability 1 - E(mu_o), the energy the GGX lobe cannot account for.
        const Real pMs = m_msActive ? (Real(1) - eAt(wo_local.z)) : Real(0);
        Point2 uLobe = u;
        // Glints: a share of the samples from the lobes of the facets under the pixel (its own random numbers: more than two).
        bool useGlint = false;
        GlintPlanes planes;
        if (m_glint && glintPlanes(rec, wo, planes)) {
            const Real pGlint = glintSampleShare(planes);
            if (pGlint > Real(0)) {
                if (uLobe.x < pGlint) useGlint = true;
                else uLobe.x = (uLobe.x - pGlint) / (Real(1) - pGlint);
            }
        }
        bool useMs = false;
        if (!useGlint && pMs > Real(0)) {
            if (uLobe.x < pMs) { useMs = true; uLobe.x = uLobe.x / pMs; }
            else { uLobe.x = (uLobe.x - pMs) / (Real(1) - pMs); }
        }

        if (useGlint) {
            int c = 0;
            Real pick = sampling::Random();
            while (c + 1 < planes.count && pick >= planes.select[c]) pick -= planes.select[c++];
            Vector3 wh_local;
            if (!(planes.select[c] > Real(0))
                || !m_glint->sample(planes.uv[c], planes.fp[c], planes.weight[c], static_cast<std::uint32_t>(planes.stream[c]), [] { return sampling::Random(); }, wh_local))
                return std::nullopt;
            bsdfSample.wi = math::reflectOutward(wo, frame.localToWorld(wh_local));
        }
        else if (useMs) {
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

    bool RoughConductor::scatterMueller(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
        const SampledWavelengths& lambda, pol::ScatterEvent& out) const {

        if (glm::dot(rec.gn, wo) <= 0 || glm::dot(rec.gn, wi) <= 0) return false;
        Vector3 wh = wi + wo;
        const Real whLen2 = glm::dot(wh, wh);
        if (whLen2 <= Real(0)) return false;
        wh *= math::safe_recip(math::safe_sqrt(whLen2));

        const Real cosD = std::abs(glm::dot(wi, wh));
        out.s = pol::planeNormalS(wo, wh);
        for (int i = 0; i < pol::kN; ++i) {
#if defined(RAYT_SPECTRAL_MODE)
            const pol::Complex n(lutAt(etaLut, lambda[i]), lutAt(kLut, lambda[i]));
#else
            const pol::Complex n(eta[i], k[i]);
#endif
            const pol::Amplitudes a = pol::fresnelAmplitudes(cosD, n);
            pol::setInterface(out.m, i, a.rs, a.rp);
        }

        if (m_msActive) {
            const Spectrum fms = msBrdf(std::max(Real(0), glm::dot(rec.n, wo)), std::max(Real(0), glm::dot(rec.n, wi)), lambda);
            const Spectrum f = eval(rec, wo, wi, lambda, TransportMode::Radiance);
            for (int i = 0; i < pol::kN; ++i) {
                if (!(f[i] > Real(0))) continue;
                const Real w = std::clamp((f[i] - fms[i]) / f[i], Real(0), Real(1));   // the polarizing share of the lobe
                for (int r = 0; r < 4; ++r)
                    for (int c = 0; c < 4; ++c) out.m.e[r][c][i] = w * out.m.e[r][c][i] + (r == 0 && c == 0 ? Real(1) - w : Real(0));
            }
        }
        return true;
    }

    Real RoughConductor::pdf(const SurfaceInteraction& rec,
        const Vector3& wo, const Vector3& wi,
        const SampledWavelengths& lambda) const {

        if (glm::dot(rec.gn, wi) <= 0 || glm::dot(rec.gn, wo) <= 0) return 0.0;

        Vector3 wh = wo + wi;
        Real whLen2 = glm::dot(wh, wh);
        if (whLen2 <= Real(0.0)) return Real(0.0);
        wh = wh * math::safe_recip(math::safe_sqrt(whLen2));

        const frame::Frame frame = shadingFrame(rec);
        Vector3 wo_local = frame.worldToLocal(wo);
        Vector3 wh_local = frame.worldToLocal(wh);

        GGXDistribution dist(alpha_x, alpha_y);
        Real pdf_wh = dist.pdf(wo_local, wh_local);

        // Transform PDF from half-vector to solid angle
        Real dot_wo_wh = std::abs(glm::dot(wo_local, wh_local));
        const Real pdfGgx = (dot_wo_wh > Real(0.0)) ? pdf_wh / (Real(4) * dot_wo_wh) : Real(0.0);

        // Mixture with the multiple-scattering lobe (cosine-weighted), chosen with probability 1 - E(mu_o).
        Real pdfSmooth = pdfGgx;
        if (m_msActive) {
            const Real pMs = Real(1) - eAt(wo_local.z);
            const Real pdfMs = std::max(glm::dot(rec.n, wi), Real(0)) / constants::PI;
            pdfSmooth = (Real(1) - pMs) * pdfGgx + pMs * pdfMs;
        }
        // Glints: and with the samples drawn from the facets' lobes (see sample()).
        GlintPlanes planes;
        if (!m_glint || !glintPlanes(rec, wo, planes)) return pdfSmooth;
        const Real pGlint = glintSampleShare(planes);
        if (!(pGlint > Real(0))) return pdfSmooth;
        Real pdfWh = 0;
        for (int c = 0; c < planes.count; ++c)
            if (planes.select[c] > Real(0))
                pdfWh += planes.select[c] * m_glint->pdf(wh_local, planes.uv[c], planes.fp[c], planes.weight[c], static_cast<std::uint32_t>(planes.stream[c]));
        const Real pdfGlint = (dot_wo_wh > Real(0.0)) ? pdfWh / (Real(4) * dot_wo_wh) : Real(0.0);
        return pGlint * pdfGlint + (Real(1) - pGlint) * pdfSmooth;
    }

} // namespace rayt
