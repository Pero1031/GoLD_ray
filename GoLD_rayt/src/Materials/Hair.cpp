// The hair fiber material (Materials/Hair.hpp): set-up of the lobes, evaluation, sampling and the pdf.
#include "pch.h"
#include "Materials/Hair.hpp"

namespace rayt {

    using namespace hairmath;

    Vector3 Hair::sigmaFromColor(const Vector3& c, Real betaN) {
        const Real b = betaN;
        const Real denom = Real(5.969) - Real(0.215) * b + Real(2.532) * b * b - Real(10.73) * b * b * b + Real(5.574) * b * b * b * b + Real(0.245) * b * b * b * b * b;
        Vector3 s;
        for (int i = 0; i < 3; ++i) s[i] = sqr(std::log(std::clamp(c[i], Real(1e-4), Real(1))) / denom);
        return s;
    }

    Hair::Hair(const Params& p) : m_p(p) {
        // Longitudinal variances of the lobes R, TT, TRT, rest (Chiang et al. 2016, fitted to Marschner's measurements).
        const Real bm = std::clamp(p.betaM, Real(0.01), Real(1));
        m_v[0] = sqr(Real(0.726) * bm + Real(0.812) * bm * bm + Real(3.7) * std::pow(bm, Real(20)));
        m_v[1] = Real(0.25) * m_v[0];
        m_v[2] = Real(4) * m_v[0];
        m_v[3] = m_v[2];
        // Azimuthal logistic scale.
        const Real bn = std::clamp(p.betaN, Real(0.01), Real(1));
        m_s = Real(0.626657069) * (Real(0.265) * bn + Real(1.194) * bn * bn + Real(5.372) * std::pow(bn, Real(22)));
        // Cuticle scale tilt: the lobes R, TT, TRT are rotated by -2a, a, 4a around the fiber's axis.
        const Real a = p.alphaDeg * constants::PI / Real(180);
        m_sin2kAlpha[0] = std::sin(a);
        m_cos2kAlpha[0] = safeSqrt(Real(1) - sqr(m_sin2kAlpha[0]));
        for (int i = 1; i < 3; ++i) {
            m_sin2kAlpha[i] = Real(2) * m_cos2kAlpha[i - 1] * m_sin2kAlpha[i - 1];
            m_cos2kAlpha[i] = sqr(m_cos2kAlpha[i - 1]) - sqr(m_sin2kAlpha[i - 1]);
        }
#if defined(RAYT_SPECTRAL_MODE)
        buildSigmaLut();
#endif
        // Denoiser albedo: the color of a fiber seen through ~1.5 radii of pigment.
        for (int i = 0; i < 3; ++i) m_albedoRgb[i] = std::clamp(std::exp(-p.sigmaA[i] * Real(1.5)), Real(0), Real(1));
    }

    // ------------------------------------------------------------------------------------
    // Material interface
    // ------------------------------------------------------------------------------------

    Spectrum Hair::eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
        const SampledWavelengths& lambda, TransportMode) const
    {
        Frame f; makeFrame(rec, f);
        const Vector3 woL = f.toLocal(wo), wiL = f.toLocal(wi);
        Common c; if (!setup(woL, lambda, c)) return Spectrum(0.0);
        return evalLocal(c, wiL);
    }

    std::optional<BSDFSample> Hair::sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u,
        const SampledWavelengths& lambda, TransportMode) const
    {
        Frame fr; makeFrame(rec, fr);
        const Vector3 woL = fr.toLocal(wo);
        Common c; if (!setup(woL, lambda, c)) return std::nullopt;

        // Choose the lobe p in proportion to the (wavelength-averaged) energy of A_p.
        int p = kMaxP;
        {
            const Real up = sampling::Random();
            Real acc = 0;
            for (int i = 0; i <= kMaxP; ++i) {
                acc += c.apPdf[i];
                if (up < acc) { p = i; break; }
            }
        }
        Real sinThetaOp, cosThetaOp;
        tiltedTheta(p, c.sinThetaO, c.cosThetaO, sinThetaOp, cosThetaOp);
        cosThetaOp = std::abs(cosThetaOp);

        // theta_i from M_p.
        const Real cosTheta = Real(1) + m_v[p] * std::log(std::max(u.x, Real(1e-5)) + (Real(1) - u.x) * std::exp(Real(-2) / m_v[p]));
        const Real sinTheta = safeSqrt(Real(1) - sqr(cosTheta));
        const Real cosPhiM = std::cos(Real(2) * constants::PI * u.y);
        const Real sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhiM * cosThetaOp;
        const Real cosThetaI = safeSqrt(Real(1) - sqr(sinThetaI));

        // phi_i from N_p.
        Real dphi;
        if (p < kMaxP) dphi = Phi(p, c.gammaO, c.gammaT) + SampleTrimmedLogistic(sampling::Random(), m_s, -constants::PI, constants::PI);
        else dphi = Real(2) * constants::PI * sampling::Random();
        const Real phiI = c.phiO + dphi;
        const Vector3 wiL(sinThetaI, cosThetaI * std::cos(phiI), cosThetaI * std::sin(phiI));

        BSDFSample bs;
        bs.wi = fr.toWorld(wiL);
        bs.pdf = pdfLocal(c, wiL);
        if (!(bs.pdf > Real(0))) return std::nullopt;
        bs.f = evalLocal(c, wiL);
        bs.flags = BxDFFlags::Glossy | BxDFFlags::Reflection | BxDFFlags::Transmission;
        return bs;
    }

    Real Hair::pdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
        const SampledWavelengths& lambda) const
    {
        Frame f; makeFrame(rec, f);
        const Vector3 woL = f.toLocal(wo), wiL = f.toLocal(wi);
        Common c; if (!setup(woL, lambda, c)) return Real(0);
        return pdfLocal(c, wiL);
    }

    std::array<Spectrum, 4> Hair::attenuation(const Vector3& woLocal, const SampledWavelengths& lambda) const {
        Common c; std::array<Spectrum, 4> out{};
        if (setup(woLocal, lambda, c)) out = c.ap;
        return out;
    }

    // ------------------------------------------------------------------------------------
    // the lobes
    // ------------------------------------------------------------------------------------

    void Hair::makeFrame(const SurfaceInteraction& rec, Frame& f) {
        f.z = rec.n;
        Vector3 t = rec.dpdu - f.z * glm::dot(rec.dpdu, f.z);
        Real len2 = glm::dot(t, t);
        if (len2 < Real(1e-20)) {   // no tangent given: any direction perpendicular to the normal (an isotropic, arbitrary fiber)
            t = std::abs(f.z.x) < Real(0.9) ? glm::cross(f.z, Vector3(1, 0, 0)) : glm::cross(f.z, Vector3(0, 1, 0));
            len2 = glm::dot(t, t);
        }
        f.x = t / std::sqrt(len2);
        f.y = glm::cross(f.z, f.x);
    }

    bool Hair::setup(const Vector3& woL, const SampledWavelengths& lambda, Common& c) const {
        c.sinThetaO = std::clamp(woL.x, Real(-1), Real(1));
        c.cosThetaO = safeSqrt(Real(1) - sqr(c.sinThetaO));
        if (c.cosThetaO < Real(1e-6)) return false;   // along the fiber: no cross-section
        c.phiO = std::atan2(woL.z, woL.y);
        c.h = std::clamp(-woL.y / c.cosThetaO, Real(-1), Real(1));   // = sin(gamma_o), see the file comment
        c.gammaO = safeAsin(c.h);

        const Real eta = m_p.eta;
        const Real sinThetaT = c.sinThetaO / eta;
        const Real cosThetaT = safeSqrt(Real(1) - sqr(sinThetaT));
        const Real etap = safeSqrt(sqr(eta) - sqr(c.sinThetaO)) / c.cosThetaO;   // the index seen in the cross-section
        const Real sinGammaT = c.h / etap;
        const Real cosGammaT = safeSqrt(Real(1) - sqr(sinGammaT));
        c.gammaT = safeAsin(sinGammaT);

        // Transmittance of one pass through the fiber (path length in radii: 2 cos(gamma_t) / cos(theta_t)).
        const Real pathLen = Real(2) * cosGammaT / cosThetaT;
        Spec T;
        for (int i = 0; i < kN; ++i) T[i] = std::exp(-sigmaA(lambda, i) * pathLen);

        // Attenuation of each lobe: Fresnel at the surface, and T per pass.
        const Real cosGammaO = safeSqrt(Real(1) - sqr(c.h));
        const Real fr = fresnel::fresnelDielectric(c.cosThetaO * cosGammaO, Real(1), eta);
        c.ap[0] = Spec(fr);
        c.ap[1] = T * sqr(Real(1) - fr);
        c.ap[2] = c.ap[1] * T * fr;
        // lobes p >= 3 form a geometric series: sum_{p>=3} A_2 (f T)^{p-2} = A_2 f T / (1 - f T)
        for (int i = 0; i < kN; ++i) c.ap[3][i] = c.ap[2][i] * fr * T[i] / (Real(1) - T[i] * fr);

        for (int k = 0; k <= kMaxP; ++k)   // switched-off lobes: no f, and no sampling of them either
            if (!((m_p.lobeMask >> k) & 1u)) c.ap[k] = Spec(0.0);

        Real sum = 0;
        std::array<Real, 4> avg{};
        for (int k = 0; k <= kMaxP; ++k) {
            Real a = 0;
            for (int i = 0; i < kN; ++i) a += c.ap[k][i];
            avg[k] = a / Real(kN);
            sum += avg[k];
        }
        for (int k = 0; k <= kMaxP; ++k) c.apPdf[k] = sum > Real(0) ? avg[k] / sum : Real(0.25);
        return true;
    }

    /// sin / cos of theta_o rotated by the scale tilt of lobe p.
    void Hair::tiltedTheta(int p, Real sinThetaO, Real cosThetaO, Real& sinOut, Real& cosOut) const {
        if (p == 0) {
            sinOut = sinThetaO * m_cos2kAlpha[1] - cosThetaO * m_sin2kAlpha[1];
            cosOut = cosThetaO * m_cos2kAlpha[1] + sinThetaO * m_sin2kAlpha[1];
        }
        else if (p == 1) {
            sinOut = sinThetaO * m_cos2kAlpha[0] + cosThetaO * m_sin2kAlpha[0];
            cosOut = cosThetaO * m_cos2kAlpha[0] - sinThetaO * m_sin2kAlpha[0];
        }
        else if (p == 2) {
            sinOut = sinThetaO * m_cos2kAlpha[2] + cosThetaO * m_sin2kAlpha[2];
            cosOut = cosThetaO * m_cos2kAlpha[2] - sinThetaO * m_sin2kAlpha[2];
        }
        else { sinOut = sinThetaO; cosOut = cosThetaO; }
    }

    Spectrum Hair::evalLocal(const Common& c, const Vector3& wiL) const {
        const Real sinThetaI = std::clamp(wiL.x, Real(-1), Real(1));
        const Real cosThetaI = safeSqrt(Real(1) - sqr(sinThetaI));
        const Real phiI = std::atan2(wiL.z, wiL.y);
        const Real phi = phiI - c.phiO;

        Spec fsum(0.0);
        for (int p = 0; p < kMaxP; ++p) {
            Real sinThetaOp, cosThetaOp;
            tiltedTheta(p, c.sinThetaO, c.cosThetaO, sinThetaOp, cosThetaOp);
            cosThetaOp = std::abs(cosThetaOp);
            fsum += c.ap[p] * (Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, m_v[p]) * Np(phi, p, m_s, c.gammaO, c.gammaT));
        }
        fsum += c.ap[kMaxP] * (Mp(cosThetaI, c.cosThetaO, sinThetaI, c.sinThetaO, m_v[kMaxP]) / (Real(2) * constants::PI));

        const Real absCos = std::abs(wiL.z);
        if (absCos > Real(0)) fsum = fsum / absCos;
        return Sanitize(fsum);
    }

    Real Hair::pdfLocal(const Common& c, const Vector3& wiL) const {
        const Real sinThetaI = std::clamp(wiL.x, Real(-1), Real(1));
        const Real cosThetaI = safeSqrt(Real(1) - sqr(sinThetaI));
        const Real phiI = std::atan2(wiL.z, wiL.y);
        const Real dphi = phiI - c.phiO;

        Real pdf = 0;
        for (int p = 0; p < kMaxP; ++p) {
            Real sinThetaOp, cosThetaOp;
            tiltedTheta(p, c.sinThetaO, c.cosThetaO, sinThetaOp, cosThetaOp);
            cosThetaOp = std::abs(cosThetaOp);
            pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, m_v[p]) * c.apPdf[p] * Np(dphi, p, m_s, c.gammaO, c.gammaT);
        }
        pdf += Mp(cosThetaI, c.cosThetaO, sinThetaI, c.sinThetaO, m_v[kMaxP]) * c.apPdf[kMaxP] / (Real(2) * constants::PI);
        return pdf;
    }

    // ------------------------------------------------------------------------------------
    // Absorption per wavelength
    // ------------------------------------------------------------------------------------
#if defined(RAYT_SPECTRAL_MODE)
    void Hair::buildSigmaLut() {
        // log-log interpolation through (610, R), (545, G), (465, B); the end segments are extended.
        const Real nm[3] = { Real(610), Real(545), Real(465) };
        Real ls[3];
        for (int i = 0; i < 3; ++i) ls[i] = std::log(std::max(m_p.sigmaA[i], Real(1e-9)));
        m_sigmaLut.resize(size_t(kLutLastNm - kLutFirstNm + 1));
        for (int l = kLutFirstNm; l <= kLutLastNm; ++l) {
            const Real x = std::log(Real(l));
            int a = (l >= 545) ? 0 : 1;   // segment (610, 545) or (545, 465)
            const Real x0 = std::log(nm[a]), x1 = std::log(nm[a + 1]);
            const Real t = (x - x0) / (x1 - x0);
            const Real v = std::exp(ls[a] + (ls[a + 1] - ls[a]) * t);
            m_sigmaLut[size_t(l - kLutFirstNm)] = (v < Real(2e-9)) ? Real(0) : v;
        }
    }
#endif

} // namespace rayt
