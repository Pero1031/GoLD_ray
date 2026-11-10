/**
 * @file Materials/SubsurfaceWalk.hpp
 * @brief Subsurface scattering by a random walk inside the object, with the same two controls as the BSSRDF
 *        (Materials/Subsurface.hpp): the diffuse color A and the scatter distance l per channel.
 *
 * The BSSRDF jumps from the entry point to an exit point on the surface; it knows nothing of how thick the object is.
 * Here the path really walks through the inside (PathIntegrator::advanceInMedium), so thin and curved parts let through
 * what they should. This file holds the two pieces that the walk needs:
 *
 *  1. The conversion (A, l, g) -> (sigma_t, sigma_s): the coefficients of a medium that, as a thick slab under a diffuse
 *     interface, returns the diffuse color A and spreads the light over the distance l.
 *       - s = sqrt(sigma_a / (sigma_a + sigma_s (1 - g))) from A by van de Hulst's formula for the albedo of a half space
 *         under diffuse light,  A = (1 - s)(1 - 0.139 s) / (1 + 1.17 s)   (solved for s: a quadratic). Cycles uses the same
 *         formula; Chiang et al. 2016 and Wrenninge et al. 2017 fit the same relation numerically.
 *       - 1 / l is the reduced extinction sigma_a + sigma_s (1 - g): the mean free path of the equivalent isotropic medium,
 *         which is what the BSSRDF's profile calls the scatter distance. Changing g keeps s and l, so the color and the
 *         spread of a thick part stay; what changes is the light after few collisions (thin parts, near the entry).
 *  2. SubsurfaceWalkMedium: that medium, per color channel (RGB build) or per wavelength (spectral build; the
 *     conversion is done at each wavelength, so the color is that of the BSSRDF there too).
 *
 * The boundary is a "diffuse interface" (the convention of production renderers): the light enters with a cosine
 * distribution and every walk that reaches the surface leaves through the BSSRDF's exit lobe. The Fresnel factors are the
 * BSSRDF's: (1 - F(wo)) going in, the normalized SubsurfaceExit lobe coming out; nothing is reflected back inside.
 * CoatedDiffuse::setSubsurfaceWalk() switches a skin material to this mode.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <memory>

#include "Core/Types.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Color/SpectralConversion.hpp"
#include "Materials/Medium.hpp"
#include "Materials/Subsurface.hpp"

namespace rayt {

    namespace sss {

        /// van de Hulst: the albedo of a half space under diffuse light, s = sqrt(sigma_a / (sigma_a + sigma_s (1 - g))) in [0, 1].
        inline Real walkAlbedo(Real s) {
            return (Real(1) - s) * (Real(1) - Real(0.139) * s) / (Real(1) + Real(1.17) * s);
        }

        /// The inverse of walkAlbedo(): s for the diffuse color A in [0, 1] (A = 1: s = 0, no absorption; A = 0: s = 1, no scattering).
        inline Real walkSimilarity(Real A) {
            A = std::clamp(A, Real(0), Real(1));
            // 0.139 s^2 - (1.139 + 1.17 A) s + (1 - A) = 0, the root in [0, 1]
            const Real b = Real(1.139) + Real(1.17) * A;
            const Real disc = std::max(b * b - Real(4) * Real(0.139) * (Real(1) - A), Real(0));
            return std::clamp((b - std::sqrt(disc)) / (Real(2) * Real(0.139)), Real(0), Real(1));
        }

        struct WalkCoefficients {
            Real sigmaT = 0;   ///< extinction (1/m)
            Real sigmaS = 0;   ///< scattering (1/m)
        };

        /// The medium for the diffuse color A, the scatter distance l (m) and the mean cosine g of the phase function.
        inline WalkCoefficients walkCoefficients(Real A, Real scatterDistance, Real g) {
            const Real s = walkSimilarity(A);
            const Real reduced = Real(1) / std::max(scatterDistance, Real(1e-9));   // sigma_a + sigma_s (1 - g)
            const Real sigmaA = s * s * reduced;
            const Real sigmaS = (Real(1) - s * s) * reduced / std::max(Real(1) - g, Real(1e-3));
            return { sigmaA + sigmaS, sigmaS };
        }

    } // namespace sss

    /// The inside of a skin-like object for the random walk: homogeneous, coefficients from (diffuse color, scatter distance, g).
    class SubsurfaceWalkMedium final : public Medium {
    public:
        SubsurfaceWalkMedium(const RgbColor& albedo, std::shared_ptr<const SubsurfaceProfile> profile, Real g)
            : m_albedo(albedo), m_profile(std::move(profile)), m_g(std::clamp(g, Real(-0.95), Real(0.95))), m_phase(m_g)
        {
            const Vector3 A = m_albedo.rgb(), l = m_profile->scatterDistanceRgb();
            for (int c = 0; c < 3; ++c) {
                const sss::WalkCoefficients k = sss::walkCoefficients(A[c], l[c], m_g);
                m_sigmaTRgb[c] = k.sigmaT;
                m_sigmaSRgb[c] = k.sigmaS;
            }
            m_constant = A.x == A.y && A.y == A.z && l.x == l.y && l.y == l.z;
        }

        Real anisotropy() const { return m_g; }
        const Vector3& sigmaTRgb() const { return m_sigmaTRgb; }
        const Vector3& sigmaSRgb() const { return m_sigmaSRgb; }

        bool isHomogeneous() const override { return true; }
        bool isSpectrallyConstant() const override { return m_constant; }
        Real phase(Real cosTheta) const override { return m_phase.eval(cosTheta); }
        Vector3 samplePhase(const Vector3& dir, const Point2& u) const override { return m_phase.sample(dir, u); }

        Spectrum sigmaTSpectrum(const SampledWavelengths& lambda) const override { return coefficients(lambda, false); }
        Spectrum sigmaSSpectrum(const SampledWavelengths& lambda) const override { return coefficients(lambda, true); }

    private:
        Spectrum coefficients(const SampledWavelengths& lambda, bool scattering) const {
#if defined(RAYT_SPECTRAL_MODE)
            if (m_constant) return Spectrum(scattering ? m_sigmaSRgb.x : m_sigmaTRgb.x);
            // The conversion is not linear in A: it is done at each wavelength, on the spectra the BSSRDF uses.
            const Spectrum A = m_albedo.reflectance(lambda), l = m_profile->scatterDistance(lambda);
            Spectrum out;
            for (int i = 0; i < kSpectralSamples; ++i) {
                const sss::WalkCoefficients k = sss::walkCoefficients(A[i], l[i], m_g);
                out[i] = scattering ? k.sigmaS : k.sigmaT;
            }
            return out;
#else
            (void)lambda;
            return Spectrum(scattering ? m_sigmaSRgb : m_sigmaTRgb);
#endif
        }

        RgbColor m_albedo;
        std::shared_ptr<const SubsurfaceProfile> m_profile;
        Real m_g;
        PhaseFunction m_phase;
        Vector3 m_sigmaTRgb{ 0.0 }, m_sigmaSRgb{ 0.0 };
        bool m_constant = false;
    };

} // namespace rayt
