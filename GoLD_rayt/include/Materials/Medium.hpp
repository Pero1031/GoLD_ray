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
 * Stage 1: gray media (sigma_t and alpha do not depend on the wavelength). Chromatic media need
 * one wavelength to sample the distance and a pdf-ratio weight for the others (spectral MIS).
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
