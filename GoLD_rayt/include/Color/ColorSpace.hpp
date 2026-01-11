/**
 * @file ColorSpace.hpp
 * @brief CIE color matching and RGB conversion utilities.
 *
 * Provides:
 * - CIE XYZ color matching functions
 * - Spectral to RGB conversion
 * - RGB color space transformations (sRGB, Rec.2020, etc.)
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
     * @brief Evaluate CIE XYZ color matching function at wavelength λ (nm).
     */
    inline CIE_XYZ evaluateCIE_XYZ(Real lambda) {
        // Gaussian approximation (simpler, faster)
        // For production, use piecewise Gaussian fits (Wyman et al. 2013)

        CIE_XYZ xyz;

        // X component (peaks around 600nm - red/orange)
        {
            Real t1 = (lambda - Real(442)) / Real(36);
            Real t2 = (lambda - Real(599.8)) / Real(37.9);
            Real t3 = (lambda - Real(501.1)) / Real(20.4);
            xyz.X = Real(0.362) * std::exp(Real(-0.5) * t1 * t1)
                + Real(1.056) * std::exp(Real(-0.5) * t2 * t2)
                - Real(0.065) * std::exp(Real(-0.5) * t3 * t3);
        }

        // Y component (peaks around 555nm - green, luminance)
        {
            Real t1 = (lambda - Real(568.8)) / Real(46.9);
            Real t2 = (lambda - Real(530.9)) / Real(16.3);
            xyz.Y = Real(0.821) * std::exp(Real(-0.5) * t1 * t1)
                + Real(0.286) * std::exp(Real(-0.5) * t2 * t2);
        }

        // Z component (peaks around 445nm - blue)
        {
            Real t1 = (lambda - Real(437)) / Real(11.8);
            Real t2 = (lambda - Real(459)) / Real(26);
            xyz.Z = Real(1.217) * std::exp(Real(-0.5) * t1 * t1)
                + Real(0.681) * std::exp(Real(-0.5) * t2 * t2);
        }

        return xyz;
    }

    // ========================================================================
    // RGB Color Spaces
    // ========================================================================

    /**
     * @brief XYZ to sRGB conversion matrix (D65 white point).
     *
     * Standard sRGB primaries with D65 illuminant.
     */
    inline glm::vec3 XYZ_to_sRGB(const CIE_XYZ& xyz) {
        // ITU-R BT.709 / sRGB conversion matrix
        glm::vec3 rgb;
        rgb.x = static_cast<float>(3.2406 * xyz.X - 1.5372 * xyz.Y - 0.4986 * xyz.Z);
        rgb.y = static_cast<float>(-0.9689 * xyz.X + 1.8758 * xyz.Y + 0.0415 * xyz.Z);
        rgb.z = static_cast<float>(0.0557 * xyz.X - 0.2040 * xyz.Y + 1.0570 * xyz.Z);
        return rgb;
    }

    /**
     * @brief Apply sRGB gamma correction (inverse EOTF).
     */
    inline glm::vec3 linearToSRGB(const glm::vec3& linear) {
        glm::vec3 srgb;
        for (int i = 0; i < 3; ++i) {
            float c = linear[i];
            if (c <= 0.0031308f) {
                srgb[i] = 12.92f * c;
            }
            else {
                srgb[i] = 1.055f * std::pow(c, 1.0f / 2.4f) - 0.055f;
            }
        }
        return srgb;
    }

    /**
     * @brief Apply inverse sRGB gamma (EOTF).
     */
    inline glm::vec3 sRGBToLinear(const glm::vec3& srgb) {
        glm::vec3 linear;
        for (int i = 0; i < 3; ++i) {
            float c = srgb[i];
            if (c <= 0.04045f) {
                linear[i] = c / 12.92f;
            }
            else {
                linear[i] = std::pow((c + 0.055f) / 1.055f, 2.4f);
            }
        }
        return linear;
    }

    // ========================================================================
    // Spectral to RGB Conversion
    // ========================================================================

    /**
     * @brief Convert spectral samples to RGB using CIE color matching.
     *
     * @param wavelengths Array of wavelengths in nm
     * @param values Spectral values at each wavelength
     * @param n Number of samples
     * @return Linear RGB color (not gamma corrected)
     */
    template<int N>
    glm::vec3 spectralToRGB(const std::array<Real, N>& wavelengths,
        const std::array<Real, N>& values) {
        CIE_XYZ xyz(0, 0, 0);

        // Riemann sum approximation of the integral
        // ∫ S(λ) * CMF(λ) dλ
        for (int i = 0; i < N; ++i) {
            Real lambda = wavelengths[i];
            Real value = values[i];

            CIE_XYZ cmf = evaluateCIE_XYZ(lambda);
            xyz += cmf * value;
        }

        // Normalize by number of samples (approximate integration)
        Real scale = (constants::LAMBDA_MAX - constants::LAMBDA_MIN) / Real(N);
        xyz.X *= scale;
        xyz.Y *= scale;
        xyz.Z *= scale;

        // Convert XYZ to linear RGB
        return XYZ_to_sRGB(xyz);
    }

} // namespace rayt::color