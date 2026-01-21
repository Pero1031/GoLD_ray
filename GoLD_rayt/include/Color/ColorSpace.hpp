/**
 * @file Color/ColorSpace.hpp
 * @brief CIE color matching and RGB conversion utilities.
 *
 * Notes / Assumptions:
 * - CIE 1931 2° CMFs are approximated using Wyman/Sloan/Shirley (2013) "multi-lobe" analytic fits.
 * - The CMF approximation is intended for wavelengths roughly in [360, 830] nm (same range used in the paper's viewer).
 * - spectralToRGB() integrates S(λ) * CMF(λ) over wavelength.
 *   What S(λ) represents matters:
 *     * If S(λ) is an SPD/radiance already including illuminant, XYZ will be meaningful up to a scale.
 *     * If S(λ) is *reflectance*, you should multiply by an illuminant SPD (e.g., D65) before integration.
 * - XYZ_to_sRGB() uses the standard sRGB matrix for D65 white point.
 *   If your XYZ is not D65-referenced (or you want colorimetry-grade results), you may need white-point adaptation / normalization.
 */

#pragma once

#include <array>
#include <cmath>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"

namespace rayt::color {

    // ========================================================================
    // CIE 1931 Standard Observer (2-degree)
    // ========================================================================

    /**
     * @brief CIE 1931 XYZ color matching functions.
     *
     * These are analytic approximations to the tabulated data.
     * More accurate than lookup tables for arbitrary wavelengths.
     *
     * Reference: "Simple Analytic Approximations to the CIE XYZ Color Matching Functions"
     * by Chris Wyman, Peter-Pike Sloan, and Peter Shirley (2013)
     */
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
     * @brief Asymmetric Gaussian used by Wyman et al. (2013).
     *
     * The paper implements terms in the form:
     *   d = (λ - μ) * (λ < μ ? k_left : k_right)
     *   value = A * exp(-0.5 * d^2)
     * Note: k is effectively (1/σ) and differs on left/right sides of the peak.
     */
    inline Real gaussAsym(Real lambda, Real mu, Real kLeft, Real kRight, Real A) {
        const Real d = (lambda - mu) * (lambda < mu ? kLeft : kRight);
        return A * std::exp(Real(-0.5) * d * d);
    }

    /**
     * @brief Wyman/Sloan/Shirley (2013) multi-lobe analytic fit for CIE 1931 2° CMFs.
     *
     * Validated in the paper against the tabulated CIE 1931 1nm data.
     * Intended usage range: approximately 360–830 nm.
     */
    inline CIE_XYZ evaluateCIE_XYZ_1931_MultiLobe(Real lambda_nm) {
        CIE_XYZ xyz;

        // X
        xyz.X =
            gaussAsym(lambda_nm, Real(442.0), Real(0.0624), Real(0.0374), Real(0.362)) +
            gaussAsym(lambda_nm, Real(599.8), Real(0.0264), Real(0.0323), Real(1.056)) -
            gaussAsym(lambda_nm, Real(501.1), Real(0.0490), Real(0.0382), Real(0.065));

        // Y
        xyz.Y =
            gaussAsym(lambda_nm, Real(568.8), Real(0.0213), Real(0.0247), Real(0.821)) +
            gaussAsym(lambda_nm, Real(530.9), Real(0.0613), Real(0.0322), Real(0.286));

        // Z
        xyz.Z =
            gaussAsym(lambda_nm, Real(437.0), Real(0.0845), Real(0.0278), Real(1.217)) +
            gaussAsym(lambda_nm, Real(459.0), Real(0.0385), Real(0.0725), Real(0.681));

        return xyz;
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
     * @brief Convert spectral samples to CIE XYZ using CIE 1931 2° CMFs (Wyman 2013 multi-lobe).
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

            const CIE_XYZ c0 = evaluateCIE_XYZ_1931_MultiLobe(l0);
            const CIE_XYZ c1 = evaluateCIE_XYZ_1931_MultiLobe(l1);

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
     * @brief Convert spectral samples to linear sRGB using CIE 1931 2° CMFs (Wyman 2013 multi-lobe).
     *
     * Returns *linear* sRGB (no gamma encoding). Use linearToSRGB() for display encoding.
     */
    template<int N>
    inline Vector3 spectralToRGB(const std::array<Real, N>& wavelengths,
        const std::array<Real, N>& values) {
        const CIE_XYZ xyz = spectralToXYZ_1931<N>(wavelengths, values);
        return XYZ_to_sRGB(xyz);
    }

} // namespace rayt::color