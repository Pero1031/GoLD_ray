/**
 * @file Materials/Medium.hpp
 * @brief Participating media for volumetric path tracing (a random walk inside closed surfaces).
 *
 * A medium is described by its extinction coefficient sigma_t (1/m), its single-scattering albedo
 * alpha = sigma_s / sigma_t and a phase function (Henyey-Greenstein, anisotropy g).
 *
 * Free-flight sampling. The probability of travelling a distance t without a collision is the
 * transmittance T(t) = exp(-sigma_t t) (Beer-Lambert), so the distance to the next collision has
 * the density  p(t) = sigma_t exp(-sigma_t t)  and is sampled by  t = -ln(1 - u) / sigma_t.
 * At the collision the path scatters with probability alpha and is absorbed otherwise; a scattered
 * direction is drawn from the phase function with exactly its density. In a homogeneous medium all
 * these factors cancel, so the path weight stays 1: the estimator is unbiased without any
 * bookkeeping (Novak et al. 2018, "Monte Carlo methods for volumetric light transport simulation").
 *
 * Gray scattering media (HomogeneousMedium): sigma_t and alpha do not depend on the wavelength.
 * Chromatic absorption without scattering needs no sampling at all: see AbsorbingMedium (a transmittance weight).
 * Chromatic scattering media (ChromaticMedium): the free-flight distance is sampled with the extinction of ONE wavelength
 * of the path (the hero; one of the three channels in the RGB build), and every wavelength gets the weight
 * p_i / p_hero, where p_i is the density of the same event had wavelength i been the hero. The integrator keeps the
 * running ratios normalized so that dividing by their mean (the balance heuristic over the hero choices, as in
 * PBRT v4 "SampledSpectrum r_u") is already done: see PathIntegrator::Li.
 *
 * The medium that a ray travels in is Ray::medium. A closed boundary material (Dielectric) knows
 * the medium inside it (Material::mediumInside()); the integrator switches the ray's medium when
 * the ray crosses the boundary.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Core/Forward.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Color/SpectralCurve.hpp"
#include "Color/SpectralConversion.hpp"
#include "Geometry/Frame.hpp"

namespace rayt {

    /// A new direction from the Henyey-Greenstein phase function for a ray travelling along `dir` (g > 0: forward).
    inline Vector3 sampleHenyeyGreenstein(const Vector3& dir, Real g, const Point2& u) {
        // cos(theta) between the old and the new propagation direction, by inverting the
        // cumulative distribution of Henyey-Greenstein: g > 0 favors cos = 1 (forward).
        Real cosTheta;
        if (std::abs(g) < Real(1e-3)) {
            cosTheta = Real(1) - Real(2) * u.x;
        }
        else {
            const Real s = (Real(1) - g * g) / (Real(1) - g + Real(2) * g * u.x);
            cosTheta = (Real(1) + g * g - s * s) / (Real(2) * g);
        }
        cosTheta = std::clamp(cosTheta, Real(-1), Real(1));
        const Real sinTheta = std::sqrt(std::max(Real(0), Real(1) - cosTheta * cosTheta));
        const Real phi = Real(2) * constants::PI * u.y;
        const frame::Frame frame(dir);
        return frame.localToWorld(Vector3(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta));
    }

    /**
     * @brief Phase function: one Henyey-Greenstein lobe, or a mixture of two (Gkioulekas et al. 2013 use such mixtures:
     *        a forward lobe and a backward lobe).
     *
     *   p(cos) = w HG(cos; g1) + (1 - w) HG(cos; g2)
     *
     * w = 1 is a single lobe. The mean cosine (the "g" of the similarity relation sigma_s' = sigma_s (1 - g)) is
     * w g1 + (1 - w) g2: two phase functions with the same mean cosine look alike in thick media, and differ where few
     * collisions happen (thin parts, edges, a light close to the surface).
     */
    struct PhaseFunction {
        Real g1 = Real(0), g2 = Real(0), w = Real(1);
        PhaseFunction() = default;
        PhaseFunction(Real g) : g1(g), g2(g), w(Real(1)) {}   // implicit: a plain g is a single lobe
        PhaseFunction(Real forward, Real backward, Real weight) : g1(forward), g2(backward), w(weight) {}

        Real meanCosine() const { return w * g1 + (Real(1) - w) * g2; }

        /// A new direction for a ray travelling along `dir`, drawn with exactly the phase function's density.
        Vector3 sample(const Vector3& dir, const Point2& u) const {
            // The first number picks the lobe and, rescaled, is reused for the angle (it stays uniform in [0, 1)).
            if (u.x < w) return sampleHenyeyGreenstein(dir, g1, Point2(w > Real(0) ? u.x / w : Real(0), u.y));
            return sampleHenyeyGreenstein(dir, g2, Point2((u.x - w) / (Real(1) - w), u.y));
        }
    };

    class Medium {
    public:
        virtual ~Medium() = default;

        virtual Real sigmaT() const = 0;      ///< extinction coefficient (1/m)
        virtual Real albedo() const = 0;      ///< probability that a collision scatters (else: absorbs)

        /// Distance to the next collision for the uniform number u in [0, 1).
        virtual Real sampleDistance(Real u) const = 0;

        /// New propagation direction after a scattering collision of a ray travelling along `dir`.
        virtual Vector3 samplePhase(const Vector3& dir, const Point2& u) const = 0;

        /// True for a medium whose coefficients depend on the wavelength and that scatters: the integrator then samples with
        /// the hero wavelength and weights the others (sigmaTSpectrum / sigmaSSpectrum). False: sigmaT() and albedo() are all.
        virtual bool isChromatic() const { return false; }

        /// Extinction and scattering coefficients (1/m) at the path's wavelengths (RGB build: per channel).
        virtual Spectrum sigmaTSpectrum(const SampledWavelengths& /*lambda*/) const { return Spectrum(sigmaT()); }
        virtual Spectrum sigmaSSpectrum(const SampledWavelengths& /*lambda*/) const { return Spectrum(sigmaT() * albedo()); }

        /// True for a medium that never scatters and whose absorption is applied as a weight along the whole
        /// segment instead of by sampling collisions (see AbsorbingMedium). The integrator then multiplies the
        /// path throughput by transmittance() and ignores sigmaT() / albedo() / sampleDistance() / samplePhase().
        virtual bool attenuatesContinuously() const { return false; }

        /// Beer-Lambert transmittance exp(-sigma_a(lambda) d) over `distance` metres, per wavelength.
        virtual Spectrum transmittance(Real /*distance*/, const SampledWavelengths& /*lambda*/) const { return Spectrum(1.0); }
    };

    /**
     * @brief Purely absorbing, wavelength-dependent medium (clear water, tinted glass, ...).
     *
     * Without scattering the light that survives a segment of length d is exactly exp(-sigma_a(lambda) d), so
     * there is nothing to sample: the path weight is multiplied by that factor for each carried wavelength.
     * This is deterministic (no variance, no collision sampling) and needs no spectral MIS, unlike a chromatic
     * medium with scattering. Spectral mode evaluates the coefficient table at the path's wavelengths; RGB mode
     * uses the mean coefficient over three bands (R 580-680, G 500-580, B 420-500 nm).
     */
    class AbsorbingMedium final : public Medium {
    public:
        /// @param sigmaA  absorption coefficient in 1/m as a function of wavelength (nm); clamped outside its range
        /// @param scale   multiplies the coefficient (an artistic control: 1 = the physical medium)
        AbsorbingMedium(color::SpectralCurve sigmaA, Real scale = Real(1))
            : m_sigmaA(std::move(sigmaA)), m_scale(scale)
        {
            auto bandMean = [&](Real lo, Real hi) {
                Real sum = 0; int n = 0;
                for (Real l = lo; l <= hi; l += Real(5)) { sum += m_sigmaA.evaluate(l, color::Extrapolate::Clamp); ++n; }
                return sum / Real(n);
            };
            m_sigmaRgb = Vector3(bandMean(580, 680), bandMean(500, 580), bandMean(420, 500)) * m_scale;
        }

        Real sigmaT() const override { return m_sigmaRgb.y; }   // not used: attenuatesContinuously()
        Real albedo() const override { return Real(0); }
        Real sampleDistance(Real) const override { return constants::INFINITY_VAL; }
        Vector3 samplePhase(const Vector3& dir, const Point2&) const override { return dir; }

        bool attenuatesContinuously() const override { return true; }

        /// Absorption coefficient (1/m) at one wavelength, scale included.
        Real sigmaA(Real lambda_nm) const { return m_scale * m_sigmaA.evaluate(lambda_nm, color::Extrapolate::Clamp); }
        const Vector3& sigmaRgb() const { return m_sigmaRgb; }

        Spectrum transmittance(Real distance, const SampledWavelengths& lambda) const override {
            Spectrum T(1.0);
#if defined(RAYT_SPECTRAL_MODE)
            for (int i = 0; i < kSpectralSamples; ++i) T[i] = std::exp(-sigmaA(lambda[i]) * distance);
#else
            (void)lambda;
            for (int c = 0; c < 3; ++c) T[c] = std::exp(-m_sigmaRgb[c] * distance);
#endif
            return T;
        }

    private:
        color::SpectralCurve m_sigmaA;
        Real m_scale;
        Vector3 m_sigmaRgb;
    };

    /// Homogeneous gray medium with a Henyey-Greenstein phase function.
    class HomogeneousMedium final : public Medium {
    public:
        /// @param sigmaT extinction coefficient (1/m)
        /// @param albedo single-scattering albedo in [0, 1]
        /// @param phase  Henyey-Greenstein anisotropy g in (-1, 1) (0 = isotropic), or a two-lobe PhaseFunction
        HomogeneousMedium(Real sigmaT, Real albedo, PhaseFunction phase = PhaseFunction())
            : m_sigmaT(sigmaT), m_albedo(albedo), m_phase(phase) {}

        Real sigmaT() const override { return m_sigmaT; }
        Real albedo() const override { return m_albedo; }
        Real anisotropy() const { return m_phase.meanCosine(); }
        const PhaseFunction& phase() const { return m_phase; }

        Real sampleDistance(Real u) const override {
            return -std::log(Real(1) - u) / m_sigmaT;
        }

        Vector3 samplePhase(const Vector3& dir, const Point2& u) const override {
            return m_phase.sample(dir, u);
        }

    private:
        Real m_sigmaT, m_albedo;
        PhaseFunction m_phase;
    };

    /**
     * @brief Homogeneous scattering medium whose extinction and single-scattering albedo depend on the wavelength
     *        (jade, wax, milk with a tint, ruby glass with scattering, ...). Gray phase function (Henyey-Greenstein).
     *
     * Given as RGB triplets (extinction in 1/m, albedo in [0, 1]). The RGB build uses them per channel; the spectral build
     * turns each into a smooth spectrum (the same basis as reflectance colors), so the medium has a continuous spectrum
     * and every wavelength of the path sees its own coefficients.
     */
    class ChromaticMedium final : public Medium {
    public:
        ChromaticMedium(const Vector3& sigmaTRgb, const Vector3& albedoRgb, PhaseFunction phase = PhaseFunction())
            : m_sigmaT(glm::max(sigmaTRgb, Vector3(1e-3))), m_albedo(glm::clamp(albedoRgb, Vector3(0.0), Vector3(1.0))), m_phase(phase) {}

        bool isChromatic() const override { return true; }
        Real anisotropy() const { return m_phase.meanCosine(); }
        const PhaseFunction& phase() const { return m_phase; }
        const Vector3& sigmaTRgb() const { return m_sigmaT; }
        const Vector3& albedoRgb() const { return m_albedo; }

        // Scalars (the mean over the channels) for code that does not know about chromatic media.
        Real sigmaT() const override { return (m_sigmaT.x + m_sigmaT.y + m_sigmaT.z) / Real(3); }
        Real albedo() const override { return (m_albedo.x + m_albedo.y + m_albedo.z) / Real(3); }
        Real sampleDistance(Real u) const override { return -std::log(Real(1) - u) / sigmaT(); }
        Vector3 samplePhase(const Vector3& dir, const Point2& u) const override { return m_phase.sample(dir, u); }

        Spectrum sigmaTSpectrum(const SampledWavelengths& lambda) const override {
#if defined(RAYT_SPECTRAL_MODE)
            Spectrum s;
            for (int i = 0; i < kSpectralSamples; ++i) s[i] = std::max(color::rgbToReflectance(m_sigmaT, lambda[i]), Real(1e-3));
            return s;
#else
            (void)lambda;
            return m_sigmaT;
#endif
        }
        Spectrum sigmaSSpectrum(const SampledWavelengths& lambda) const override {
#if defined(RAYT_SPECTRAL_MODE)
            Spectrum s;
            for (int i = 0; i < kSpectralSamples; ++i)
                s[i] = sigmaTSpectrum(lambda)[i] * std::min(color::rgbToReflectance(m_albedo, lambda[i]), Real(1));
            return s;
#else
            (void)lambda;
            return m_sigmaT * m_albedo;
#endif
        }

    private:
        Vector3 m_sigmaT, m_albedo;
        PhaseFunction m_phase;
    };

} // namespace rayt
