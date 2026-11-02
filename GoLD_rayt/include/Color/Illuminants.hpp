/**
 * @file Color/Illuminants.hpp
 * @brief Standard illuminant spectral power distributions (SPDs).
 *
 * Currently provides CIE standard illuminant D65 (the sRGB / Rec.709 white point).
 *
 * Why D65 is needed:
 * - A reflectance spectrum R(lambda) only becomes a color once it is lit by something.
 *   sRGB assumes D65 lighting, so integrating R(lambda) * D65(lambda) * CMF(lambda)
 *   gives the color the material has under "white" light in an sRGB pipeline.
 * - See Color/ColorSpace.hpp : spectrumToRGB_D65().
 *
 * Data:
 * - Relative SPD, normalized to 100 at 560 nm.
 * - 380..780 nm in 5 nm steps (81 samples), linearly interpolated, 0 outside the range.
 * - Sanity check (Python, Wyman 2013 CMF fit): chromaticity of the table is
 *   (x, y) = (0.3126, 0.3292) vs the official D65 (0.3127, 0.3290).
 */

#pragma once

#include <array>
#include <algorithm>
#include <cmath>

#include "Core/Precision.hpp"

namespace rayt::color {

    inline constexpr Real kD65_FirstNm = Real(380);
    inline constexpr Real kD65_LastNm = Real(780);
    inline constexpr Real kD65_StepNm = Real(5);

    inline constexpr std::array<Real, 81> kD65_SPD = {
        49.9755, 52.3118, 54.6482, 68.7015, 82.7549, 87.1204, 91.4860, 92.4589, 93.4318, 90.0570,   // 380-425
        86.6823, 95.7736, 104.865, 110.936, 117.008, 117.410, 117.812, 116.336, 114.861, 115.392,   // 430-475
        115.923, 112.367, 108.811, 109.082, 109.354, 108.578, 107.802, 106.296, 104.790, 106.239,   // 480-525
        107.689, 106.047, 104.405, 104.225, 104.046, 102.023, 100.000, 98.1671, 96.3342, 96.0611,   // 530-575
        95.7880, 92.2368, 88.6856, 89.3459, 90.0062, 89.8026, 89.5991, 88.6489, 87.6987, 85.4936,   // 580-625
        83.2886, 83.4939, 83.6992, 81.8630, 80.0268, 80.1207, 80.2146, 81.2462, 82.2778, 80.2810,   // 630-675
        78.2842, 74.0027, 69.7213, 70.6652, 71.6091, 72.9790, 74.3490, 67.9765, 61.6040, 65.7448,   // 680-725
        69.8856, 72.4863, 75.0870, 69.3398, 63.5927, 55.0054, 46.4182, 56.6118, 66.8054, 65.0941,   // 730-775
        63.3828                                                                                    // 780
    };

    /**
     * @brief Relative spectral power of CIE illuminant D65 at a wavelength (nm).
     * @return Linearly interpolated value (560 nm == 100), or 0 outside [380, 780] nm.
     */
    inline Real evaluateD65(Real lambda_nm) {
        if (lambda_nm < kD65_FirstNm || lambda_nm > kD65_LastNm) return Real(0);

        const Real f = (lambda_nm - kD65_FirstNm) / kD65_StepNm;
        const int last = static_cast<int>(kD65_SPD.size()) - 2;
        const int i = std::min(static_cast<int>(f), last);
        const Real t = f - Real(i);
        return kD65_SPD[i] + t * (kD65_SPD[i + 1] - kD65_SPD[i]);
    }

    /**
     * @brief Relative spectral radiance of a blackbody (Planck's law) at a wavelength (nm).
     *
     * Only the shape matters (the caller normalizes), so the absolute prefactor is dropped.
     * Illuminant A is a blackbody of 2856 K (xy = 0.4476, 0.4074).
     * @param kelvin Temperature in K (> 0).
     */
    inline Real evaluateBlackbody(Real lambda_nm, Real kelvin) {
        constexpr Real kC2 = Real(1.438777e7);          // hc/k in nm*K
        const Real l = lambda_nm;
        const Real l5 = l * l * l * l * l;
        return Real(1e15) / (l5 * std::expm1(kC2 / (l * kelvin)));   // 1e15 keeps values near 1..100
    }

} // namespace rayt::color
