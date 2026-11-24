// A thin leaf with PROSPECT's diffuse reflection and transmission and a GGX cuticle on each face (Materials/Leaf.hpp).
#include "pch.h"

#include "Materials/Leaf.hpp"

#include <algorithm>
#include "Core/Constants.hpp"
#include "Core/Fresnel.hpp"
#include "Core/Math.hpp"
#include "Core/Sampling.hpp"
#include "Color/SpectralConversion.hpp"
#include "Color/SpectralCurve.hpp"
#include "Geometry/Frame.hpp"
#include "Microfacet/GGX.hpp"

namespace rayt {

    namespace {
        Real meanOf(const Spectrum& s) {
#if defined(RAYT_SPECTRAL_MODE)
            Real m = 0;
            for (int i = 0; i < kSpectralSamples; ++i) m += s[i];
            return m / Real(kSpectralSamples);
#else
            return (s[0] + s[1] + s[2]) / Real(3);
#endif
        }
        constexpr int kChannels =
#if defined(RAYT_SPECTRAL_MODE)
            kSpectralSamples;
#else
            3;
#endif
    }

    Leaf::Leaf(std::shared_ptr<const leaf::ProspectTable> table, const leaf::Pigments& pigments, Real roughnessTop, Real roughnessBottom)
        : m_table(std::move(table)), m_p(pigments)
        , m_alphaTop(MicrofacetDistribution::roughnessToAlpha(roughnessTop))
        , m_alphaBottom(MicrofacetDistribution::roughnessToAlpha(roughnessBottom))
    {
        color::SpectralCurve R, T, Rd, Td;
        for (int nm = 360; nm <= 830; ++nm) {
            const leaf::ReflTrans rt = leaf::prospect(*m_table, m_p, Real(nm));
            R.add(Real(nm), rt.R);
            T.add(Real(nm), rt.T);
            Rd.add(Real(nm), std::max(rt.R - (Real(1) - rt.t12), Real(0)));
            Td.add(Real(nm), rt.T);
        }
        for (auto* c : { &R, &T, &Rd, &Td }) c->finalize();
        m_rgbR = color::reflectanceCurveToRgb(R);
        m_rgbT = color::reflectanceCurveToRgb(T);
        m_rgbRd = glm::max(color::reflectanceCurveToRgb(Rd), Vector3(0.0));
        m_rgbTd = glm::max(color::reflectanceCurveToRgb(Td), Vector3(0.0));
        m_n550 = m_table->at(Real(550)).n;
        m_esTop = cuticleTable(m_alphaTop, m_n550);
        m_esBottom = cuticleTable(m_alphaBottom, m_n550);
    }

    Leaf::CuticleTable Leaf::cuticleTable(Real alpha, Real n) {
        // E(mu) = Integral F D G / (4 mu mu_i) mu_i dwi, by VNDF samples on a fixed 96 x 96 grid: weight F G2 / G1(wo)
        CuticleTable t;
        const GGXDistribution dist(alpha, alpha);
        constexpr int g = 96;
        for (int k = 0; k < CuticleTable::kSize; ++k) {
            const Real mu = std::max(Real(k) / Real(CuticleTable::kSize - 1), Real(1e-4));
            const Vector3 wo(math::safe_sqrt(Real(1) - mu * mu), Real(0), mu);
            Real sum = 0;
            for (int i = 0; i < g; ++i)
                for (int j = 0; j < g; ++j) {
                    const Vector3 wh = dist.sample_wh(wo, Point2((Real(i) + Real(0.5)) / g, (Real(j) + Real(0.5)) / g));
                    const Real oh = glm::dot(wo, wh);
                    const Vector3 wi = Real(2) * oh * wh - wo;
                    if (wi.z <= 0 || oh <= 0) continue;
                    sum += fresnel::fresnelDielectric(oh, Real(1), n) * dist.G(wo, wi) / std::max(dist.G1(wo), Real(1e-12));
                }
            t.e[k] = sum / Real(g * g);
        }
        Real avg = 0;   // 2 Integral E(mu) mu dmu, midpoint rule on a fine grid
        constexpr int m = 512;
        for (int i = 0; i < m; ++i) {
            const Real mu = (Real(i) + Real(0.5)) / m;
            avg += t.at(mu) * mu;
        }
        t.avg = Real(2) * avg / m;
        return t;
    }

    Real Leaf::CuticleTable::at(Real mu) const {
        const Real x = std::clamp(mu, Real(0), Real(1)) * Real(kSize - 1);
        const int i = std::min(int(x), kSize - 2);
        const Real f = x - Real(i);
        return e[i] + (e[i + 1] - e[i]) * f;
    }

    Leaf::Optics Leaf::opticsAt(const SurfaceInteraction& /*rec*/, const SampledWavelengths& lambda) const {
        Optics o;
#if defined(RAYT_SPECTRAL_MODE)
        for (int i = 0; i < kSpectralSamples; ++i) {
            const leaf::ReflTrans rt = leaf::prospect(*m_table, m_p, lambda[i]);
            o.Rd[i] = std::max(rt.R - (Real(1) - rt.t12), Real(0));
            o.T[i] = rt.T;
            o.t12[i] = rt.t12;
            o.n[i] = rt.n;
        }
#else
        o.Rd = RgbColor(m_rgbRd).reflectance(lambda);
        o.T = RgbColor(m_rgbTd).reflectance(lambda);
        o.t12 = Spectrum(leaf::tav(Real(90), m_n550));   // one index for the three channels: t12 must belong to the index of F
        o.n = Spectrum(m_n550);
#endif
        return o;
    }

    Vector3 Leaf::aovAlbedo(const SurfaceInteraction& /*rec*/) const {
        return glm::clamp(m_rgbR + m_rgbT, Vector3(0.0), Vector3(1.0));
    }

    Real Leaf::lobeProbability(const SurfaceInteraction& rec, const Optics& o, Real cosO, Real& pReflect) const {
        const CuticleTable& es = rec.frontFace ? m_esTop : m_esBottom;
        const Real F = es.at(cosO);
        const Real rd = meanOf(o.Rd), t = meanOf(o.T);
        const Real d = (Real(1) - F) * (rd + t) / (Real(1) - es.avg);
        pReflect = (rd + t) > Real(1e-8) ? std::clamp(rd / (rd + t), Real(0.05), Real(0.95)) : Real(0.5);
        return std::clamp(F / std::max(F + d, Real(1e-6)), Real(0.05), Real(0.95));
    }

    Spectrum Leaf::specular(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const Optics& o) const {
        const Real cosO = glm::dot(rec.n, wo), cosI = glm::dot(rec.n, wi);
        if (cosO <= 0 || cosI <= 0) return Spectrum(0.0);
        Vector3 wh = wo + wi;
        const Real len2 = glm::dot(wh, wh);
        if (len2 <= Real(0)) return Spectrum(0.0);
        wh /= std::sqrt(len2);
        const frame::Frame frame(rec.n);
        const GGXDistribution dist(alphaFor(rec), alphaFor(rec));
        const Vector3 o_l = frame.worldToLocal(wo), i_l = frame.worldToLocal(wi), h_l = frame.worldToLocal(wh);
        const Real DG = dist.D(h_l) * dist.G(o_l, i_l) / (Real(4) * cosO * cosI);
        const Real cosD = std::abs(glm::dot(wi, wh));
        Spectrum f(0.0);
        for (int c = 0; c < kChannels; ++c) f[c] = fresnel::fresnelDielectric(cosD, Real(1), o.n[c]) * DG;
        return f;
    }

    Spectrum Leaf::diffuseLobe(const SurfaceInteraction& rec, const Optics& o, Real cosO, Real cosI, bool reflect) const {
        // the face the light leaves by (the viewer's) and the face it arrives at (the same one, or the other for a transmission)
        const CuticleTable& out = rec.frontFace ? m_esTop : m_esBottom;
        const CuticleTable& in = (rec.frontFace == reflect) ? m_esTop : m_esBottom;
        const Real k = (Real(1) - out.at(cosO)) * (Real(1) - in.at(cosI)) / (constants::PI * (Real(1) - out.avg) * (Real(1) - in.avg));
        Spectrum f(0.0);
        for (int c = 0; c < kChannels; ++c) f[c] = (reflect ? o.Rd[c] : o.T[c]) * k;
        return f;
    }

    Spectrum Leaf::evalDiffuse(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const SampledWavelengths& lambda) const {
        const Real go = glm::dot(rec.gn, wo), gi = glm::dot(rec.gn, wi);
        const Real cosO = glm::dot(rec.n, wo), cosI = glm::dot(rec.n, wi);
        if (go <= 0 || cosO <= 0) return Spectrum(0.0);
        const bool reflect = gi > 0;
        if (reflect ? cosI <= 0 : (gi >= 0 || cosI >= 0)) return Spectrum(0.0);
        const Optics o = opticsAt(rec, lambda);
        return diffuseLobe(rec, o, cosO, std::abs(cosI), reflect);
    }

    Spectrum Leaf::eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const SampledWavelengths& lambda,
        TransportMode /*mode*/) const
    {
        const Real go = glm::dot(rec.gn, wo), gi = glm::dot(rec.gn, wi);
        const Real cosO = glm::dot(rec.n, wo), cosI = glm::dot(rec.n, wi);
        if (go <= 0 || cosO <= 0) return Spectrum(0.0);
        const bool reflect = gi > 0;
        if (reflect ? cosI <= 0 : (gi >= 0 || cosI >= 0)) return Spectrum(0.0);

        const Optics o = opticsAt(rec, lambda);
        Spectrum f = diffuseLobe(rec, o, cosO, std::abs(cosI), reflect);
        if (reflect) f = f + specular(rec, wo, wi, o);
        return f;
    }

    Real Leaf::pdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const SampledWavelengths& lambda) const {
        const Real go = glm::dot(rec.gn, wo), gi = glm::dot(rec.gn, wi);
        const Real cosO = glm::dot(rec.n, wo), cosI = glm::dot(rec.n, wi);
        if (go <= 0 || cosO <= 0) return Real(0);
        const Optics o = opticsAt(rec, lambda);
        Real pR = 0;
        const Real pS = lobeProbability(rec, o, cosO, pR);
        if (gi > 0) {
            if (cosI <= 0) return Real(0);
            Real pdfSpec = 0;
            Vector3 wh = wo + wi;
            const Real len2 = glm::dot(wh, wh);
            if (len2 > Real(0)) {
                wh /= std::sqrt(len2);
                const frame::Frame frame(rec.n);
                const GGXDistribution dist(alphaFor(rec), alphaFor(rec));
                const Vector3 o_l = frame.worldToLocal(wo), h_l = frame.worldToLocal(wh);
                const Real dotOH = std::abs(glm::dot(o_l, h_l));
                if (dotOH > 0) pdfSpec = dist.pdf(o_l, h_l) / (Real(4) * dotOH);
            }
            return pS * pdfSpec + (Real(1) - pS) * pR * cosI / constants::PI;
        }
        if (gi >= 0 || cosI >= 0) return Real(0);
        return (Real(1) - pS) * (Real(1) - pR) * (-cosI) / constants::PI;
    }

    std::optional<BSDFSample> Leaf::sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u, const SampledWavelengths& lambda,
        TransportMode mode) const
    {
        const Real cosO = glm::dot(rec.n, wo);
        if (glm::dot(rec.gn, wo) <= 0 || cosO <= 0) return std::nullopt;
        const Optics o = opticsAt(rec, lambda);
        Real pR = 0;
        const Real pS = lobeProbability(rec, o, cosO, pR);
        const frame::Frame frame(rec.n);

        BSDFSample s;
        if (u.x < pS) {
            const GGXDistribution dist(alphaFor(rec), alphaFor(rec));
            const Vector3 wh = frame.localToWorld(dist.sample_wh(frame.worldToLocal(wo), Point2(u.x / pS, u.y)));
            s.wi = math::reflectOutward(wo, wh);
            // a reflection that points into the leaf is lost (the lobe's own shadowing): it must not pass for a transmission
            if (glm::dot(rec.gn, s.wi) <= 0 || glm::dot(rec.n, s.wi) <= 0) return std::nullopt;
            s.flags = BxDFFlags::Reflection | BxDFFlags::Glossy;
        }
        else {
            Real ux = (u.x - pS) / (Real(1) - pS);
            const bool reflect = ux < pR;
            ux = reflect ? ux / pR : (ux - pR) / (Real(1) - pR);
            Vector3 d = sampling::CosineSampleHemisphere(Point2(std::min(ux, Real(1) - Real(1e-9)), u.y));
            if (!reflect) d.z = -d.z;
            s.wi = frame.localToWorld(d);
            s.flags = BxDFFlags::Diffuse | (reflect ? BxDFFlags::Reflection : BxDFFlags::Transmission);
        }
        s.pdf = pdf(rec, wo, s.wi, lambda);
        s.f = eval(rec, wo, s.wi, lambda, mode);
        if (s.pdf <= Real(1e-8) || math::hasNonFinite(s.f)) return std::nullopt;
        return s;
    }

} // namespace rayt
