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
 * Kinds of media:
 *   - HomogeneousMedium: sigma_t and alpha do not depend on the wavelength (the cheapest walk: no weights).
 *   - ChromaticMedium: sigma_t and alpha depend on the wavelength, optionally on the position (a DensityField).
 *   - AbsorbingMedium: absorption only; needs no sampling at all (a transmittance weight).
 *
 * All scattering media are walked by the same code (PathIntegrator::advanceInMedium). The free-flight distance is sampled
 * with the extinction of ONE wavelength of the path (the hero; one of the three channels in the RGB build) at the densest
 * point, null collisions make up for a lower local density, and every wavelength gets the weight p_i / p_hero, where p_i is
 * the density of the same event had wavelength i been the hero. The integrator keeps the running ratios normalized so that
 * the division by their mean (the balance heuristic over the hero choices, as in PBRT v4 "SampledSpectrum r_u") is done.
 *
 * The medium that a ray travels in is Ray::medium. A closed boundary material (Dielectric) knows
 * the medium inside it (Material::mediumInside()); the integrator switches the ray's medium when
 * the ray crosses the boundary.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <memory>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Core/Forward.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Color/SpectralCurve.hpp"
#include "Color/SpectralConversion.hpp"
#include "Materials/DensityField.hpp"
#include "Materials/TabulatedPhase.hpp"
#include "Geometry/Frame.hpp"

namespace rayt {

    /// The Henyey-Greenstein density 1 / (4 pi) (1 - g^2) / (1 + g^2 - 2 g cos)^(3/2), cos between the old and the new propagation direction.
    inline Real henyeyGreenstein(Real cosTheta, Real g) {
        const Real d = Real(1) + g * g - Real(2) * g * cosTheta;
        return (Real(1) - g * g) / (Real(4) * constants::PI * d * std::sqrt(std::max(d, Real(1e-12))));
    }

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
        std::shared_ptr<const TabulatedPhase> table;          // if set, it replaces the Henyey-Greenstein lobes (Mie scattering, ...)
        PhaseFunction() = default;
        PhaseFunction(Real g) : g1(g), g2(g), w(Real(1)) {}   // implicit: a plain g is a single lobe
        PhaseFunction(Real forward, Real backward, Real weight) : g1(forward), g2(backward), w(weight) {}
        PhaseFunction(std::shared_ptr<const TabulatedPhase> t) : g1(t->meanCosine()), g2(t->meanCosine()), w(Real(1)), table(std::move(t)) {}

        Real meanCosine() const { return table ? table->meanCosine() : w * g1 + (Real(1) - w) * g2; }

        /// The density (per steradian) of scattering from the propagation direction d_in into d_out, with cosTheta = dot(d_in, d_out)
        /// (cosTheta = 1 is straight on: forward scattering for g > 0). Sample() draws from exactly this density.
        Real eval(Real cosTheta) const {
            if (table) return table->eval(cosTheta);
            return w * henyeyGreenstein(cosTheta, g1) + (Real(1) - w) * henyeyGreenstein(cosTheta, g2);
        }

        /// A new direction for a ray travelling along `dir`, drawn with exactly the phase function's density.
        Vector3 sample(const Vector3& dir, const Point2& u) const {
            if (table) {
                const Real cosTheta = table->sampleCosine(u.x);
                const Real sinTheta = std::sqrt(std::max(Real(0), Real(1) - cosTheta * cosTheta));
                const Real phi = Real(2) * constants::PI * u.y;
                const frame::Frame frame(dir);
                return frame.localToWorld(Vector3(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta));
            }
            // The first number picks the lobe and, rescaled, is reused for the angle (it stays uniform in [0, 1)).
            if (u.x < w) return sampleHenyeyGreenstein(dir, g1, Point2(w > Real(0) ? u.x / w : Real(0), u.y));
            return sampleHenyeyGreenstein(dir, g2, Point2((u.x - w) / (Real(1) - w), u.y));
        }
    };

    class Medium {
    public:
        virtual ~Medium() = default;

        /// Extinction and scattering coefficients (1/m) at the path's wavelengths (RGB build: per channel).
        /// For a heterogeneous medium: at the densest point (q = 1), the majorant of the walk.
        virtual Spectrum sigmaTSpectrum(const SampledWavelengths& lambda) const = 0;
        virtual Spectrum sigmaSSpectrum(const SampledWavelengths& lambda) const = 0;

        /// New propagation direction after a scattering collision of a ray travelling along `dir`.
        virtual Vector3 samplePhase(const Vector3& dir, const Point2& u) const = 0;

        /// Density (per steradian) of the phase function, cosTheta = dot(incoming propagation direction, outgoing propagation direction).
        /// The density that samplePhase() draws from. (Isotropic unless a medium says otherwise.)
        virtual Real phase(Real /*cosTheta*/) const { return Real(1) / (Real(4) * constants::PI); }

        /// True if the density is the same everywhere (q = 1): the transmittance along a segment is then exp(-sigma_t d), without sampling.
        virtual bool isHomogeneous() const { return false; }

        /// Relative density q(x) in [0, 1] at a world-space point: sigma_t(x) = q(x) sigmaTSpectrum(). Homogeneous media: 1.
        /// With q < 1 null collisions make up the difference to the majorant.
        virtual Real density(const Point3& /*p*/) const { return Real(1); }

        /// True if sigma_t and the albedo are the same for every wavelength: the weights of the walk are then all 1 and the
        /// integrator skips them (a gray medium, with or without a density field).
        virtual bool isSpectrallyConstant() const { return false; }

        /// True for a medium that never scatters and whose absorption is applied as a weight along the whole
        /// segment instead of by sampling collisions (see AbsorbingMedium). The integrator then multiplies the
        /// path throughput by transmittance() and does not walk.
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

        // Not used by the walk (attenuatesContinuously()); the coefficients are those of the physics: all absorption.
        Spectrum sigmaTSpectrum(const SampledWavelengths& lambda) const override {
            Spectrum s(0.0);
#if defined(RAYT_SPECTRAL_MODE)
            for (int i = 0; i < kSpectralSamples; ++i) s[i] = sigmaA(lambda[i]);
#else
            (void)lambda;
            s = m_sigmaRgb;
#endif
            return s;
        }
        Spectrum sigmaSSpectrum(const SampledWavelengths&) const override { return Spectrum(0.0); }
        Vector3 samplePhase(const Vector3& dir, const Point2&) const override { return dir; }

        bool attenuatesContinuously() const override { return true; }
        bool isHomogeneous() const override { return true; }

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

        Real sigmaT() const { return m_sigmaT; }
        Real albedo() const { return m_albedo; }
        Real anisotropy() const { return m_phase.meanCosine(); }
        const PhaseFunction& phase() const { return m_phase; }

        bool isSpectrallyConstant() const override { return true; }
        bool isHomogeneous() const override { return true; }
        Real phase(Real cosTheta) const override { return m_phase.eval(cosTheta); }
        Spectrum sigmaTSpectrum(const SampledWavelengths&) const override { return Spectrum(m_sigmaT); }
        Spectrum sigmaSSpectrum(const SampledWavelengths&) const override { return Spectrum(m_sigmaT * m_albedo); }

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
        /// @param density  optional density field q(x) in [0, 1] (null: homogeneous); sigmaTRgb is then the extinction where q = 1
        ChromaticMedium(const Vector3& sigmaTRgb, const Vector3& albedoRgb, PhaseFunction phase = PhaseFunction(),
            std::shared_ptr<const DensityField> density = nullptr)
            : m_sigmaT(glm::max(sigmaTRgb, Vector3(1e-3))), m_albedo(glm::clamp(albedoRgb, Vector3(0.0), Vector3(1.0))), m_phase(phase), m_density(std::move(density))
        {
            m_constant = m_sigmaT.x == m_sigmaT.y && m_sigmaT.y == m_sigmaT.z && m_albedo.x == m_albedo.y && m_albedo.y == m_albedo.z;
        }

        bool isSpectrallyConstant() const override { return m_constant; }
        bool isHomogeneous() const override { return m_density == nullptr; }
        Real phase(Real cosTheta) const override { return m_phase.eval(cosTheta); }
        Real density(const Point3& p) const override { return m_density ? std::clamp((*m_density)(p), Real(0), Real(1)) : Real(1); }
        bool isHeterogeneous() const { return m_density != nullptr; }
        Real anisotropy() const { return m_phase.meanCosine(); }
        const PhaseFunction& phase() const { return m_phase; }
        const Vector3& sigmaTRgb() const { return m_sigmaT; }
        const Vector3& albedoRgb() const { return m_albedo; }

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
        std::shared_ptr<const DensityField> m_density;
        bool m_constant = false;   // equal channels: see isSpectrallyConstant()
    };

} // namespace rayt
