/**
 * @file Color/ColorSpace.hpp
 * @brief CIE color matching and RGB conversion utilities.
 *
 * Notes / Assumptions:
 * - CIE 1931 2° CMFs: the 1 nm table of the standard (Color/CIE1931Table.hpp), linearly interpolated, 0 outside [360, 830] nm.
 * - spectralToRGB() integrates S(λ) * CMF(λ) over wavelength.
 *   What S(λ) represents matters:
 *     * If S(λ) is an SPD/radiance already including illuminant, XYZ will be meaningful up to a scale.
 *     * If S(λ) is *reflectance*, you should multiply by an illuminant SPD (e.g., D65) before integration.
 * - XYZ_to_sRGB() uses the standard sRGB matrix for D65 white point.
 *   If your XYZ is not D65-referenced (or you want colorimetry-grade results), you may need white-point adaptation / normalization.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Color/CIE1931Table.hpp"
#include "Color/Illuminants.hpp"

namespace rayt::color {

    // ========================================================================
    // CIE 1931 Standard Observer (2-degree)
    // ========================================================================

    /// CIE 1931 XYZ tristimulus values (or the three colour-matching functions at one wavelength).
    struct CIE_XYZ {
        Real X, Y, Z;

        CIE_XYZ(Real x = 0, Real y = 0, Real z = 0) : X(x), Y(y), Z(z) {}

        CIE_XYZ& operator+=(const CIE_XYZ& o) {
            X += o.X; Y += o.Y; Z += o.Z;
            return *this;
        }

        friend CIE_XYZ operator*(const CIE_XYZ& xyz, Real s) {
            return CIE_XYZ(xyz.X * s, xyz.Y * s, xyz.Z * s);
        }
    };

    /**
     * @brief CIE 1931 2° colour-matching functions at a wavelength (nm): the 1 nm table (Color/CIE1931Table.hpp), linearly
     *        interpolated, 0 outside 360..830 nm. At whole nanometres it returns the table's values exactly.
     */
    inline CIE_XYZ evaluateCIE_XYZ_1931(Real lambda_nm) {
        if (!(lambda_nm >= Real(kCIE1931_FirstNm) && lambda_nm <= Real(kCIE1931_LastNm))) return CIE_XYZ(0, 0, 0);
        const Real x = lambda_nm - Real(kCIE1931_FirstNm);
        const int i = std::min(static_cast<int>(x), kCIE1931_LastNm - kCIE1931_FirstNm - 1);
        const Real t = x - Real(i);
        const auto& a = kCIE1931_2deg[static_cast<size_t>(i)];
        const auto& b = kCIE1931_2deg[static_cast<size_t>(i) + 1];
        return CIE_XYZ(Real(a[0] + t * (b[0] - a[0])), Real(a[1] + t * (b[1] - a[1])), Real(a[2] + t * (b[2] - a[2])));
    }

    // ========================================================================
    // RGB Color Spaces
    // ========================================================================

    /**
     * @brief Convert CIE XYZ to linear sRGB (D65).
     *
     * This is the standard matrix for sRGB/BT.709 primaries with D65 white point.
     * Input XYZ should be in the same reference white (D65) to be strictly consistent.
     */
    inline Vector3 XYZ_to_sRGB(const CIE_XYZ& xyz) {
        Vector3 rgb;
        rgb.x = Real(3.2406) * xyz.X - Real(1.5372) * xyz.Y - Real(0.4986) * xyz.Z;
        rgb.y = -Real(0.9689) * xyz.X + Real(1.8758) * xyz.Y + Real(0.0415) * xyz.Z;
        rgb.z = Real(0.0557) * xyz.X - Real(0.2040) * xyz.Y + Real(1.0570) * xyz.Z;
        return rgb;
    }

    /**
    * @brief sRGB OETF (linear -> sRGB encoded) for a single channel.
    * @note Input is assumed linear-light. Negative values are not meaningful for display; caller should clamp if needed.
    */
    inline Real linearToSRGB(Real x) {
        if (x <= Real(0.0031308))
            return Real(12.92) * x;
        return Real(1.055) * std::pow(x, Real(1.0) / Real(2.4)) - Real(0.055);
    }

    /**
     * @brief sRGB OETF (linear -> sRGB encoded).
     *
     * Note: This is defined for non-negative linear values in typical display pipelines.
     * If linear values can go negative (e.g., out-of-gamut conversions), handle/clamp as needed
     * before encoding for display.
     */
    inline Vector3 linearToSRGB(const Vector3& linear) {
        Vector3 srgb;
        srgb[0] = linearToSRGB(linear[0]);
        srgb[1] = linearToSRGB(linear[1]);
        srgb[2] = linearToSRGB(linear[2]);
        return srgb;
    }

    /**
     * @brief Apply inverse sRGB gamma (EOTF).
     */
    inline Vector3 sRGBToLinear(const Vector3& srgb) {
        constexpr Real a = Real(0.04045);
        Vector3 linear;
        for (int i = 0; i < 3; ++i) {
            const Real c = srgb[i];
            if (c <= a) {
                linear[i] = c / Real(12.92);
            }
            else {
                linear[i] = std::pow((c + Real(0.055)) / Real(1.055), Real(2.4));
            }
        }
        return linear;
    }

    // ========================================================================
    // Spectral to RGB Conversion
    // ========================================================================

    /**
     * @brief Convert spectral samples to CIE XYZ using CIE 1931 2° CMFs.
     *
     * Integration:
     * - Uses trapezoidal rule over the provided wavelength samples.
     * - Assumes wavelengths are in nanometers and sorted ascending.
     *
     * @param wavelengths Wavelengths in nm (must be ascending)
     * @param values     Spectral values S(λ) at each wavelength
     * @return CIE XYZ tristimulus values (scale depends on S(λ) and any normalization you apply)
     */
    template<int N>
    inline CIE_XYZ spectralToXYZ_1931(const std::array<Real, N>& wavelengths,
        const std::array<Real, N>& values) {
        CIE_XYZ xyz(0, 0, 0);
        if constexpr (N < 2) return xyz;

        // Trapezoidal integration: ∫ f(λ) dλ ≈ Σ 0.5 (f_i + f_{i+1}) Δλ_i
        for (int i = 0; i < N - 1; ++i) {
            const Real l0 = wavelengths[i];
            const Real l1 = wavelengths[i + 1];
            const Real dl = (l1 - l0);

            const CIE_XYZ c0 = evaluateCIE_XYZ_1931(l0);
            const CIE_XYZ c1 = evaluateCIE_XYZ_1931(l1);

            const Real v0 = values[i];
            const Real v1 = values[i + 1];

            // f(λ) = S(λ) * CMF(λ)
            xyz.X += Real(0.5) * (v0 * c0.X + v1 * c1.X) * dl;
            xyz.Y += Real(0.5) * (v0 * c0.Y + v1 * c1.Y) * dl;
            xyz.Z += Real(0.5) * (v0 * c0.Z + v1 * c1.Z) * dl;
        }
        return xyz;
    }


    /**
     * @brief Convert spectral samples to linear sRGB using CIE 1931 2° CMFs.
     *
     * Returns *linear* sRGB (no gamma encoding). Use linearToSRGB() for display encoding.
     */
    template<int N>
    inline Vector3 spectralToRGB(const std::array<Real, N>& wavelengths,
        const std::array<Real, N>& values) {
        const CIE_XYZ xyz = spectralToXYZ_1931<N>(wavelengths, values);
        return XYZ_to_sRGB(xyz);
    }

    /**
     * @brief Convert a reflectance-like spectrum S(lambda) to linear sRGB under D65 light.
     *
     * Computes  RGB = M * Integral( S(l) * D65(l) * CMF(l) dl )  and divides it per channel by the
     * same integral for S(l) == 1, so that:
     *   - a perfectly white spectrum (S == 1) maps to exactly (1, 1, 1),
     *   - the result does not depend on the D65 table's normalization.
     *
     * Intended for smooth, "0..1-ish" spectra: reflectance, and also eta(lambda) (an RGB value
     * that stays close to the spectrum's average, like a color).
     *
     * @note The sRGB matrix has negative lobes, so saturated spectra (e.g. gold's red plateau)
     *       can land slightly outside [0, 1] (out of sRGB gamut). Callers decide how to clamp.
     * @note Uses 1 nm trapezoid integration over [lambdaMin, lambdaMax] (D65 table: 380..780 nm).
     *
     * @tparam SpectrumFn Callable Real(Real lambda_nm) returning S(lambda).
     */
    template<class SpectrumFn>
    inline Vector3 spectrumToRGB_D65(SpectrumFn&& spectrum,
        Real lambdaMin_nm = Real(380),
        Real lambdaMax_nm = Real(780),
        Real stepNm = Real(1))
    {
        CIE_XYZ xyzS(0, 0, 0);   // with the spectrum
        CIE_XYZ xyzW(0, 0, 0);   // reference white (S == 1)

        const int n = static_cast<int>(std::round((lambdaMax_nm - lambdaMin_nm) / stepNm));
        for (int i = 0; i <= n; ++i) {
            const Real lambda = lambdaMin_nm + stepNm * Real(i);
            const Real trapezoid = (i == 0 || i == n) ? Real(0.5) : Real(1);
            const Real weight = evaluateD65(lambda) * trapezoid * stepNm;

            const CIE_XYZ cmf = evaluateCIE_XYZ_1931(lambda);
            const Real s = static_cast<Real>(spectrum(lambda));

            xyzS += cmf * (weight * s);
            xyzW += cmf * weight;
        }

        const Vector3 rgbS = XYZ_to_sRGB(xyzS);
        const Vector3 rgbW = XYZ_to_sRGB(xyzW);
        return Vector3(rgbS.x / rgbW.x, rgbS.y / rgbW.y, rgbS.z / rgbW.z);
    }

} // namespace rayt::color