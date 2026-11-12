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
 *
 * The lobe mathematics (hairmath) is in this header; the material's own code is in src/Materials/Hair.cpp.
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
        static Vector3 sigmaFromColor(const Vector3& c, Real betaN);

        explicit Hair(const Params& p);

        const Params& params() const { return m_p; }

        // ------------------------------------------------------------------------------------
        // Material interface
        // ------------------------------------------------------------------------------------

        Spectrum eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda, TransportMode = TransportMode::Radiance) const override;

        std::optional<BSDFSample> sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u,
            const SampledWavelengths& lambda, TransportMode = TransportMode::Radiance) const override;

        Real pdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda) const override;

        Vector3 aovAlbedo(const SurfaceInteraction&) const override { return m_albedoRgb; }

        /// Test access: the A_p of the four lobes at an outgoing direction in the hair frame (wo.x along the fiber, h = -wo.y / cos).
        std::array<Spectrum, 4> attenuation(const Vector3& woLocal, const SampledWavelengths& lambda) const;

    private:
        using Spec = Spectrum;

        /// Orthonormal frame at the hit: x = fiber tangent, z = surface normal, y = z x x.
        struct Frame {
            Vector3 x, y, z;
            Vector3 toLocal(const Vector3& w) const { return Vector3(glm::dot(w, x), glm::dot(w, y), glm::dot(w, z)); }
            Vector3 toWorld(const Vector3& w) const { return x * w.x + y * w.y + z * w.z; }
        };

        static void makeFrame(const SurfaceInteraction& rec, Frame& f);

        /// Everything that depends on wo only (and the wavelengths): shared by f, pdf and sampling.
        struct Common {
            Real sinThetaO = 0, cosThetaO = 1, phiO = 0, h = 0, gammaO = 0, gammaT = 0;
            std::array<Spec, 4> ap{};
            std::array<Real, 4> apPdf{};
        };

        bool setup(const Vector3& woL, const SampledWavelengths& lambda, Common& c) const;
        void tiltedTheta(int p, Real sinThetaO, Real cosThetaO, Real& sinOut, Real& cosOut) const;
        Spectrum evalLocal(const Common& c, const Vector3& wiL) const;
        Real pdfLocal(const Common& c, const Vector3& wiL) const;

        // ------------------------------------------------------------------------------------
        // Absorption per wavelength
        // ------------------------------------------------------------------------------------
#if defined(RAYT_SPECTRAL_MODE)
        static constexpr int kN = kSpectralSamples;
        static constexpr int kLutFirstNm = 360, kLutLastNm = 830;
        std::vector<Real> m_sigmaLut;

        void buildSigmaLut();
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
