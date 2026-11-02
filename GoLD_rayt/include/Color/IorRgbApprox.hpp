/*
* @file Color/IorRgbApprox.hpp
*/

#pragma once

#include <algorithm>
#include <cmath>

#include "Core/Precision.hpp"
#include "Core/Types.hpp"
#include "Color/ComplexIorTable.hpp"
#include "Color/ColorSpace.hpp"

namespace rayt::color {

    /**
     * @brief Unpolarized Fresnel reflectance of a conductor at normal incidence (F0).
     *
     *   F0 = ((eta - 1)^2 + k^2) / ((eta + 1)^2 + k^2)
     */
    inline Real conductorNormalIncidenceReflectance(Real eta, Real k) {
        const Real a = eta - Real(1);
        const Real b = eta + Real(1);
        return (a * a + k * k) / (b * b + k * k);
    }

    /**
     * @brief Derives RGB eta/k from the FULL tabulated spectrum (CMF integration under D65).
     *
     * Why not simply integrate eta(l) and k(l) separately?  The visible color of a metal is
     * determined by its *reflectance* R(l), which is a non-linear function of (eta, k). Averaging
     * eta and k first and then computing R gives a different color than averaging R. So:
     *
     *  1. F0_rgb  = spectrumToRGB_D65( R(l) )            -- the color the metal really has
     *  2. eta_rgb = spectrumToRGB_D65( eta(l) )          -- smooth quantity, averaged like a color
     *  3. k_rgb   : solved per channel (closed form) so that F0(eta_rgb, k_rgb) == F0_rgb
     *               k^2 = (F0 (eta+1)^2 - (eta-1)^2) / (1 - F0)
     *
     * The angular behaviour is then still given by the per-channel Fresnel equations
     * (fresnelConductor), but F0 -- the dominant term for the look -- is spectrally correct.
     *
     * Clamping: the sRGB matrix has negative lobes, so a saturated spectrum can give F0 > 1 in a
     * channel (gold's red). F0 is clamped to [0, max_l R(l)] (max over 380..780 nm) so that the
     * solved k stays physical. If the target is below the Fresnel floor for the chosen eta
     * (k^2 < 0), k is set to 0.
     *
     * @param table       Complex IOR table (wavelengths in nm).
     * @param outEtaRgb   Output eta  as (R, G, B).
     * @param outKRgb     Output k    as (R, G, B).
     */
    inline void integrateIorToRgb_CMF(
        const ComplexIorTable& table,
        Vector3& outEtaRgb,
        Vector3& outKRgb)
    {
        if (table.empty()) {
            outEtaRgb = Vector3(1);
            outKRgb = Vector3(0);
            return;
        }

        constexpr Real kLambdaMin = Real(380);
        constexpr Real kLambdaMax = Real(780);

        const auto etaAt = [&table](Real lambda_nm) -> Real {
            return table.evaluateNK(lambda_nm).first;
        };
        const auto reflectanceAt = [&table](Real lambda_nm) -> Real {
            const auto [eta, k] = table.evaluateNK(lambda_nm);
            return conductorNormalIncidenceReflectance(eta, k);
        };

        // Upper bound for F0: no channel can reflect more than the best spectral reflectance.
        Real rMax = Real(0);
        for (Real l = kLambdaMin; l <= kLambdaMax; l += Real(1)) {
            rMax = std::max(rMax, reflectanceAt(l));
        }
        rMax = std::min(rMax, Real(0.9999));

        const Vector3 f0Rgb = spectrumToRGB_D65(reflectanceAt, kLambdaMin, kLambdaMax);
        const Vector3 etaRgb = spectrumToRGB_D65(etaAt, kLambdaMin, kLambdaMax);

        for (int c = 0; c < 3; ++c) {
            const Real eta = std::max(etaRgb[c], Real(0.01));
            const Real f0 = std::clamp(f0Rgb[c], Real(0), rMax);

            const Real k2 = (f0 * (eta + Real(1)) * (eta + Real(1)) - (eta - Real(1)) * (eta - Real(1)))
                / (Real(1) - f0);

            outEtaRgb[c] = eta;
            outKRgb[c] = (k2 > Real(0)) ? std::sqrt(k2) : Real(0);
        }
    }

    /**
     * @brief Extracts RGB-like eta/k by sampling the IOR table at 3 representative wavelengths.
     *
     * Quick approximation (kept for comparison). Prefer integrateIorToRgb_CMF(), which uses the
     * whole spectrum instead of three point samples.
     */
    inline void approxIorToRgb_3Wavelengths(
        const ComplexIorTable& table,
        Vector3& outEtaRgb,
        Vector3& outKRgb,
        Real lambdaB_nm = Real(450),
        Real lambdaG_nm = Real(550),
        Real lambdaR_nm = Real(650))
    {
        auto [etaB, kB] = table.evaluateNK(lambdaB_nm);
        auto [etaG, kG] = table.evaluateNK(lambdaG_nm);
        auto [etaR, kR] = table.evaluateNK(lambdaR_nm);

        // Note: We pack as (R,G,B).
        outEtaRgb = Vector3(etaR, etaG, etaB);
        outKRgb = Vector3(kR, kG, kB);
    }

} // namespace rayt::color