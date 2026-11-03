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
 * Stage 1: gray scattering media (sigma_t and alpha do not depend on the wavelength). Chromatic scattering
 * media need one wavelength to sample the distance and a pdf-ratio weight for the others (spectral MIS).
 * Chromatic absorption without scattering needs none of that: see AbsorbingMedium (a transmittance weight).
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
#include "Geometry/Frame.hpp"

namespace rayt {

    class Medium {
    public:
        virtual ~Medium() = default;

        virtual Real sigmaT() const = 0;      ///< extinction coefficient (1/m)
        virtual Real albedo() const = 0;      ///< probability that a collision scatters (else: absorbs)

        /// Distance to the next collision for the uniform number u in [0, 1).
        virtual Real sampleDistance(Real u) const = 0;

        /// New propagation direction after a scattering collision of a ray travelling along `dir`.
        virtual Vector3 samplePhase(const Vector3& dir, const Point2& u) const = 0;

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
        /// @param g      Henyey-Greenstein anisotropy in (-1, 1); 0 = isotropic
        HomogeneousMedium(Real sigmaT, Real albedo, Real g = Real(0))
            : m_sigmaT(sigmaT), m_albedo(albedo), m_g(g) {}

        Real sigmaT() const override { return m_sigmaT; }
        Real albedo() const override { return m_albedo; }
        Real anisotropy() const { return m_g; }

        Real sampleDistance(Real u) const override {
            return -std::log(Real(1) - u) / m_sigmaT;
        }

        Vector3 samplePhase(const Vector3& dir, const Point2& u) const override {
            // cos(theta) between the old and the new propagation direction, by inverting the
            // cumulative distribution of Henyey-Greenstein: g > 0 favors cos = 1 (forward).
            Real cosTheta;
            if (std::abs(m_g) < Real(1e-3)) {
                cosTheta = Real(1) - Real(2) * u.x;
            }
            else {
                const Real s = (Real(1) - m_g * m_g) / (Real(1) - m_g + Real(2) * m_g * u.x);
                cosTheta = (Real(1) + m_g * m_g - s * s) / (Real(2) * m_g);
            }
            cosTheta = std::clamp(cosTheta, Real(-1), Real(1));
            const Real sinTheta = std::sqrt(std::max(Real(0), Real(1) - cosTheta * cosTheta));
            const Real phi = Real(2) * constants::PI * u.y;
            const frame::Frame frame(dir);
            return frame.localToWorld(Vector3(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta));
        }

    private:
        Real m_sigmaT, m_albedo, m_g;
    };

} // namespace rayt
