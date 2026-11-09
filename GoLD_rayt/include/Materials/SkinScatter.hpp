/**
 * @file Materials/SkinScatter.hpp
 * @brief The light that dry skin scatters at its surface: loose flakes of the stratum corneum, with air between their layers (docs/SKIN_MOISTURE.md, 2.3).
 *
 * A local, white, diffuse reflection of reflectance
 *     rho(p) = rho0 * sigma(p) * f_match * f_hyd.
 *  - rho0   : the reflectance of a flaking patch of dry skin (not measured: a free parameter).
 *  - sigma  : 0..1, where the flakes are (patches of a low-frequency noise; a higher threshold = fewer flakes).
 *  - f_match: a liquid that fills the gaps between the flakes takes away most of what scatters there. The flake / air interface reflects
 *             R(n_s | 1); with a liquid of index n_f it reflects R(n_s | n_f), so the scatter drops by the ratio of the two Fresnel reflectances,
 *             f = [ ((n_s - n_f) / (n_s + n_f)) / ((n_s - 1) / (n_s + 1)) ]^2  (water 0.125, glycerin 0.015, air 1).
 *             How much of the gaps the liquid fills is W = 1 - exp(-thickness / gap): 1 for a film many gaps thick.
 *  - f_hyd  : a hydrated stratum corneum scatters less inside (in vivo, the scattering falls when skin is hydrated; the size is not known):
 *             1 - kappaH * hydration.
 * The reflection is white, so it lowers the chroma of the skin; removing it makes the skin darker and more saturated, the direction of wet surfaces.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include "Core/Types.hpp"
#include "Materials/DensityField.hpp"

namespace rayt {

    class SkinScatter {
    public:
        struct Params {
            Real rho0 = Real(0.08);          ///< reflectance of the flakes where they are
            Real frequency = Real(300);      ///< patches of flakes per metre (a few millimetres)
            Real threshold = Real(0.52);     ///< noise level above which there are flakes (the noise is bell shaped around 0.5)
            Real softness = Real(0.08);      ///< width of the edge of a patch, in noise units
            Real hydration = Real(0);        ///< 0..1, how hydrated the stratum corneum is
            Real kappaH = Real(0.3);         ///< fraction of the scatter that full hydration takes away (unknown: to be varied)
            Real gap = Real(1e-6);           ///< the liquid fills the gaps over a depth of about this (m)
            Real skinIor = Real(1.55);       ///< refractive index of the flakes
            std::uint32_t seed = 1;
        };

        explicit SkinScatter(const Params& p)
            : m_p(p), m_noise(p.frequency, 3, p.threshold - p.softness, p.threshold + p.softness, p.seed + 211u) {}

        const Params& params() const { return m_p; }

        /// Where the flakes are, 0..1.
        Real flakes(const Point3& p) const { return m_noise(p); }

        /// The factor by which a liquid of index `liquidIor` that fills the gaps cuts the scatter: 1 for air (<= 1), 0 for the index of the flakes.
        static Real match(Real skinIor, Real liquidIor) {
            if (liquidIor <= Real(1)) return Real(1);
            const Real a = (skinIor - liquidIor) / (skinIor + liquidIor), b = (skinIor - Real(1)) / (skinIor + Real(1));
            return (a / b) * (a / b);
        }

        /// The reflectance of the scatter at p, with a liquid film of this mean thickness and index on the skin (thickness 0 = dry).
        Real reflectance(const Point3& p, Real liquidIor, Real filmThickness) const {
            const Real sigma = m_noise(p);
            if (sigma <= Real(0)) return Real(0);
            const Real wet = filmThickness > Real(0) ? Real(1) - std::exp(-filmThickness / m_p.gap) : Real(0);
            const Real fMatch = Real(1) - wet * (Real(1) - match(m_p.skinIor, liquidIor));
            const Real fHyd = std::max(Real(1) - m_p.kappaH * m_p.hydration, Real(0));
            return m_p.rho0 * sigma * fMatch * fHyd;
        }

    private:
        Params m_p;
        NoiseDensity m_noise;
    };

} // namespace rayt
