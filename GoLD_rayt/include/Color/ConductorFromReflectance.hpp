/**
 * @file Color/ConductorFromReflectance.hpp
 * @brief Builds a conductor's complex IOR table (eta(l), k(l)) from a measured normal-incidence
 *        reflectance R(l), using Gulbrandsen's reflectivity / edge-tint parameterization.
 *
 * What is exact and what is assumed:
 *  - EXACT: the reflectance at normal incidence. The returned (eta, k) satisfy
 *      F0(l) = ((eta-1)^2 + k^2) / ((eta+1)^2 + k^2) = R(l)
 *    for every wavelength, so the color seen head-on matches the measurement.
 *  - ASSUMED: how the reflectance changes with angle. One number alone cannot fix it, because
 *    many (eta, k) pairs give the same R. The edge tint g in [0,1] picks one of them:
 *      eta = g * (1-r)/(1+r) + (1-g) * (1+sqrt(r))/(1-sqrt(r)),   r = R(l)
 *      k   = sqrt( (r (eta+1)^2 - (eta-1)^2) / (1-r) )
 *    The default g = 0.9 was calibrated against the Johnson & Christy tables shipped in
 *    assets/Metal_data (Au, Ag, Cu), 400..700 nm. Maximum |dF| versus the true (eta, k):
 *    <= 0.013 up to 45 deg, <= 0.04 at 60 deg, <= 0.07 at 70 deg, 0.15..0.2 at 80..85 deg.
 *    The resulting eta and k are NOT the metal's physical constants, only a consistent pair.
 *
 * If reflectance was measured at several angles, fit (eta, k) instead of using this.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include "Core/Precision.hpp"
#include "Color/SpectralCurve.hpp"
#include "Color/ComplexIorTable.hpp"

namespace rayt::color {

    inline constexpr Real kDefaultEdgeTint = Real(0.9);

    /// (eta, k) with normal-incidence reflectance r and edge tint g (clamped to valid ranges).
    inline void complexIorFromReflectance(Real r, Real edgeTint, Real& eta, Real& k) {
        r = std::clamp(r, Real(0), Real(0.9999));
        const Real g = std::clamp(edgeTint, Real(0), Real(1));
        const Real sr = std::sqrt(r);
        eta = g * (Real(1) - r) / (Real(1) + r) + (Real(1) - g) * (Real(1) + sr) / (Real(1) - sr);
        const Real k2 = (r * (eta + Real(1)) * (eta + Real(1)) - (eta - Real(1)) * (eta - Real(1))) / (Real(1) - r);
        k = (k2 > Real(0)) ? std::sqrt(k2) : Real(0);
    }

    /**
     * @brief Converts a reflectance curve R(l) (fractions 0..1) to a complex IOR table.
     *
     * The table is sampled every stepNm over the curve's own wavelength range; outside that range
     * the table (like every ComplexIorTable) holds its end values, i.e. the reflectance at the
     * first / last measured wavelength.
     */
    inline ComplexIorTable complexIorTableFromReflectance(const SpectralCurve& reflectance,
        Real edgeTint = kDefaultEdgeTint, Real stepNm = Real(5))
    {
        ComplexIorTable table;
        if (reflectance.size() < 2) return table;

        const Real first = reflectance.firstNm(), last = reflectance.lastNm();
        for (Real l = first; ; l += stepNm) {
            const Real wl = std::min(l, last);
            Real eta, k;
            complexIorFromReflectance(reflectance.evaluate(wl, Extrapolate::Clamp), edgeTint, eta, k);
            table.add(wl, eta, k);
            if (l >= last) break;
        }
        table.finalize();
        return table;
    }

} // namespace rayt::color
