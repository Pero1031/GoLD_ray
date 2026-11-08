/**
 * @file Materials/Hair.hpp
 * @brief Scattering from a single hair fiber: the Marschner et al. 2003 lobes with the
 *        energy-conserving, importance-sampled form of Chiang et al. 2016 (the one pbrt-v4 uses).
 *
 * A hair is a dielectric cylinder (index eta ~ 1.55) that contains pigment (melanin). Light that hits it
 * leaves through a small number of paths p:
 *   p = 0  R    reflection off the surface
 *   p = 1  TT   transmission in, transmission out (the glow of a backlit hair)
 *   p = 2  TRT  in, one internal reflection, out (the secondary, colored highlight)
 *   p >= 3      everything with more internal bounces, summed in closed form and treated as diffuse in azimuth
 * The scattering function factorizes per lobe:
 *   f(wo, wi) = sum_p  M_p(theta_i, theta_o) * A_p(h, ...) * N_p(phi) / |cos theta_i|
 * M_p is the longitudinal distribution (angle along the fiber, a Gaussian-like Bessel form, variance v_p),
 * A_p the attenuation (Fresnel at every interface and absorption inside, per wavelength), N_p the azimuthal
 * distribution (a logistic around the geometric exit angle Phi(p, h)).
 *
 * Geometry convention (the shape must provide it, see Geometry/HairCurves.hpp):
 *   rec.dpdu = fiber tangent, rec.n = the cylinder's outward normal at the hit (perpendicular to the tangent).
 * The offset h in [-1, 1] across the fiber's width is not passed in: for a hit on a cylinder it follows from
 * the outgoing direction, h = sin(gamma_o) = -wo_y / cos(theta_o) in the frame (x = tangent, z = normal, y = z x x).
 *
 * The division by |cos theta_i| (= |wi . n|) cancels the cosine that the integrator multiplies in: the hair
 * model is defined for the projected width of the fiber, which does not depend on the azimuth.
 *
 * Absorption sigma_a is given per hair radius (the path through the fiber is measured in radii), as in pbrt.
 * In the spectral build the RGB values are the values at 610 / 545 / 465 nm, interpolated in log-log between
 * them (a simplification: real melanin spectra are smooth, monotone and close to a power law).
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <vector>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Core/Interaction.hpp"
#include "Core/Math.hpp"
#include "Core/Fresnel.hpp"
#include "Core/Sampling.hpp"
#include "Core/SpectrumUtils.hpp"
#include "Materials/Material.hpp"

namespace rayt {

    namespace hairmath {

        constexpr int kMaxP = 3;   ///< explicit lobes R, TT, TRT; lobe kMaxP is the sum of all the rest

        inline Real sqr(Real x) { return x * x; }
        inline Real safeSqrt(Real x) { return std::sqrt(std::max(x, Real(0))); }
        inline Real safeAsin(Real x) { return std::asin(std::clamp(x, Real(-1), Real(1))); }

        /// Modified Bessel function of the first kind, order 0 (power series; accurate for x up to ~12).
        inline Real I0(Real x) {
            Real val = 0, x2i = 1;
            long long ifact = 1;
            long long i4 = 1;
            for (int i = 0; i < 10; ++i) {
                if (i > 1) ifact *= i;
                val += x2i / (Real(i4) * sqr(Real(ifact)));
                x2i *= x * x;
                i4 *= 4;
            }
            return val;
        }

        /// log I0(x), with the asymptotic expansion for large x.
        inline Real LogI0(Real x) {
            if (x > Real(12))
                return x + Real(0.5) * (-std::log(Real(2) * constants::PI) + std::log(Real(1) / x) + Real(1) / (Real(8) * x));
            return std::log(I0(x));
        }

        /// Longitudinal scattering M_p (normalized: its integral over theta_i with the measure cos(theta_i) is 1).
        inline Real Mp(Real cosThetaI, Real cosThetaO, Real sinThetaI, Real sinThetaO, Real v) {
            const Real a = cosThetaI * cosThetaO / v;
            const Real b = sinThetaI * sinThetaO / v;
            if (v <= Real(0.1))
                return std::exp(LogI0(a) - b - Real(1) / v + Real(0.6931) + std::log(Real(1) / (Real(2) * v)));
            return (std::exp(-b) * I0(a)) / (std::sinh(Real(1) / v) * Real(2) * v);
        }

        inline Real Logistic(Real x, Real s) {
            x = std::abs(x);
            const Real e = std::exp(-x / s);
            return e / (s * sqr(Real(1) + e));
        }
        inline Real LogisticCdf(Real x, Real s) { return Real(1) / (Real(1) + std::exp(-x / s)); }
        inline Real TrimmedLogistic(Real x, Real s, Real a, Real b) {
            return Logistic(x, s) / (LogisticCdf(b, s) - LogisticCdf(a, s));
        }
        inline Real SampleTrimmedLogistic(Real u, Real s, Real a, Real b) {
            const Real k = LogisticCdf(b, s) - LogisticCdf(a, s);
            const Real p = std::clamp(u * k + LogisticCdf(a, s), Real(1e-12), Real(1) - Real(1e-12));
            const Real x = -s * std::log(Real(1) / p - Real(1));
            return std::clamp(x, a, b);
        }

        /// Geometric exit azimuth of the path with p internal segments, relative to the incoming azimuth.
        inline Real Phi(int p, Real gammaO, Real gammaT) {
            return Real(2) * Real(p) * gammaT - Real(2) * gammaO + Real(p) * constants::PI;
        }

        /// Azimuthal scattering N_p: a logistic around Phi, wrapped into [-pi, pi].
        inline Real Np(Real phi, int p, Real s, Real gammaO, Real gammaT) {
            Real dphi = phi - Phi(p, gammaO, gammaT);
            while (dphi > constants::PI) dphi -= Real(2) * constants::PI;
            while (dphi < -constants::PI) dphi += Real(2) * constants::PI;
            return TrimmedLogistic(dphi, s, -constants::PI, constants::PI);
        }

    } // namespace hairmath

    class Hair : public Material {
    public:
        struct Params {
            Real eta = Real(1.55);                    ///< index of refraction of the fiber
            Vector3 sigmaA = Vector3(0.4, 0.7, 1.4);  ///< absorption per hair radius at R / G / B (see sigmaFromMelanin)
            Real betaM = Real(0.3);                   ///< longitudinal roughness in [0, 1] (the width of the highlight along the fiber)
            Real betaN = Real(0.3);                   ///< azimuthal roughness in [0, 1] (the width around the fiber)
            Real alphaDeg = Real(2);                  ///< tilt of the cuticle scales in degrees (shifts R, TT and TRT apart)
            unsigned lobeMask = 0xF;                  ///< which lobes are on: bit 0 R, 1 TT, 2 TRT, 3 the rest (an illustration aid; the default is all)
        };

        // --- ways to choose the absorption ---

        /// Eumelanin / pheomelanin concentrations (about 1.3 / 0 for brown, 0.3 / 0 light brown, 8 / 0 black, 0 / 1 red, 0.05 / 0.2 blond;
        /// 0 / 0 is white). The two spectra are the RGB values of pbrt-v4's HairMaterial.
        static Vector3 sigmaFromMelanin(Real eumelanin, Real pheomelanin) {
            return eumelanin * Vector3(0.419, 0.697, 1.37) + pheomelanin * Vector3(0.187, 0.4, 1.05);
        }

        /// The absorption that gives a fiber the apparent color c after all the internal bounces (d'Eon et al. 2011),
        /// for the azimuthal roughness betaN. c in (0, 1] per channel.
        static Vector3 sigmaFromColor(const Vector3& c, Real betaN) {
            const Real b = betaN;
            const Real denom = Real(5.969) - Real(0.215) * b + Real(2.532) * b * b - Real(10.73) * b * b * b + Real(5.574) * b * b * b * b + Real(0.245) * b * b * b * b * b;
            Vector3 s;
            for (int i = 0; i < 3; ++i) s[i] = hairmath::sqr(std::log(std::clamp(c[i], Real(1e-4), Real(1))) / denom);
            return s;
        }

        explicit Hair(const Params& p) : m_p(p) {
            using namespace hairmath;
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

        const Params& params() const { return m_p; }

        // ------------------------------------------------------------------------------------
        // Material interface
        // ------------------------------------------------------------------------------------

        Spectrum eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda, TransportMode = TransportMode::Radiance) const override
        {
            Frame f; makeFrame(rec, f);
            const Vector3 woL = f.toLocal(wo), wiL = f.toLocal(wi);
            Common c; if (!setup(woL, lambda, c)) return Spectrum(0.0);
            return evalLocal(c, wiL);
        }

        std::optional<BSDFSample> sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u,
            const SampledWavelengths& lambda, TransportMode = TransportMode::Radiance) const override
        {
            using namespace hairmath;
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

        Real pdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda) const override
        {
            Frame f; makeFrame(rec, f);
            const Vector3 woL = f.toLocal(wo), wiL = f.toLocal(wi);
            Common c; if (!setup(woL, lambda, c)) return Real(0);
            return pdfLocal(c, wiL);
        }

        Vector3 aovAlbedo(const SurfaceInteraction&) const override { return m_albedoRgb; }

        /// Test access: the A_p of the four lobes at an outgoing direction in the hair frame (wo.x along the fiber, h = -wo.y / cos).
        std::array<Spectrum, 4> attenuation(const Vector3& woLocal, const SampledWavelengths& lambda) const {
            Common c; std::array<Spectrum, 4> out{};
            if (setup(woLocal, lambda, c)) out = c.ap;
            return out;
        }

    private:
        using Spec = Spectrum;

        /// Orthonormal frame at the hit: x = fiber tangent, z = surface normal, y = z x x.
        struct Frame {
            Vector3 x, y, z;
            Vector3 toLocal(const Vector3& w) const { return Vector3(glm::dot(w, x), glm::dot(w, y), glm::dot(w, z)); }
            Vector3 toWorld(const Vector3& w) const { return x * w.x + y * w.y + z * w.z; }
        };

        static void makeFrame(const SurfaceInteraction& rec, Frame& f) {
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

        /// Everything that depends on wo only (and the wavelengths): shared by f, pdf and sampling.
        struct Common {
            Real sinThetaO = 0, cosThetaO = 1, phiO = 0, h = 0, gammaO = 0, gammaT = 0;
            std::array<Spec, 4> ap{};
            std::array<Real, 4> apPdf{};
        };

        bool setup(const Vector3& woL, const SampledWavelengths& lambda, Common& c) const {
            using namespace hairmath;
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
        void tiltedTheta(int p, Real sinThetaO, Real cosThetaO, Real& sinOut, Real& cosOut) const {
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

        Spectrum evalLocal(const Common& c, const Vector3& wiL) const {
            using namespace hairmath;
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

        Real pdfLocal(const Common& c, const Vector3& wiL) const {
            using namespace hairmath;
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
        static constexpr int kN = kSpectralSamples;
        static constexpr int kLutFirstNm = 360, kLutLastNm = 830;
        std::vector<Real> m_sigmaLut;

        void buildSigmaLut() {
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
        Real sigmaA(const SampledWavelengths& lambda, int i) const {
            const Real x = std::clamp(lambda[i], Real(kLutFirstNm), Real(kLutLastNm)) - Real(kLutFirstNm);
            const int k = std::min(int(x), int(m_sigmaLut.size()) - 2);
            const Real f = x - Real(k);
            return m_sigmaLut[size_t(k)] * (Real(1) - f) + m_sigmaLut[size_t(k) + 1] * f;
        }
#else
        static constexpr int kN = 3;
        Real sigmaA(const SampledWavelengths&, int i) const { return m_p.sigmaA[i]; }
#endif

        Params m_p;
        std::array<Real, 4> m_v{};
        Real m_s = Real(0.3);
        std::array<Real, 3> m_sin2kAlpha{}, m_cos2kAlpha{};
        Vector3 m_albedoRgb = Vector3(0.5);
    };

} // namespace rayt
