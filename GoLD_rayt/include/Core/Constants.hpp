#pragma once

#include <limits>
#include <numbers>

#include "Core/Types.hpp"  // to use Real = double

/**
 * @file constants.hpp
 * @brief Global constants for mathematical, computational, and physical calculations.
 * * Consistent with the project-wide 'Real' type definition.
 * Mathematical constants leverage C++20 standard library templates for maximum precision.
 */

namespace rayt::constants {

    // -------------------------------------------------------------------------
    // 1. Mathematical Constants
    // -------------------------------------------------------------------------

    // C++20: Use standard library constants for maximum precision and portability.
    constexpr Real PI = std::numbers::pi_v<Real>;
    constexpr Real INV_PI = std::numbers::inv_pi_v<Real>;          // 1/π
    constexpr Real INV_SQRT_PI = std::numbers::inv_sqrtpi_v<Real>; // 1/√π

    // Derived constants frequently used in ray tracing equations.
    constexpr Real TWO_PI = Real(2) * PI;             // Full circle (360 deg) in radians
    constexpr Real FOUR_PI = Real(4) * PI;            // Solid angle of a sphere
    constexpr Real INV_TWO_PI = Real(1) / TWO_PI;
    constexpr Real INV_FOUR_PI = Real(1) / FOUR_PI;   // Normalization factor for isotropic phase functions 1/4π

    // Unit Conversion: Angle degree/radian transformations
    constexpr Real DEG_TO_RAD = PI / Real(180.0);
    constexpr Real RAD_TO_DEG = Real(180.0) / PI;


    // -------------------------------------------------------------------------
    // 2. Computational Constants
    // -------------------------------------------------------------------------

    // Representation of infinity for initial ray distances and bounding box bounds.
    constexpr Real INFINITY_VAL = std::numeric_limits<Real>::infinity();

    // Numerical epsilon to prevent self-intersection artifacts (Shadow Acne).
    // Note: This value is scale-dependent. 1e-5 is typically suitable for meter-scale scenes.
    // For microscopic or millimeter-scale rendering, consider 1e-7 or smaller.
    constexpr Real RAY_EPSILON = 1e-5;

    // Numerical tolerance for floating-point comparisons in intersection tests.
    constexpr Real INTERSECT_TOLERANCE = 1e-8;


    // -------------------------------------------------------------------------
    // 3. Physical Constants for Spectral Rendering
    // -------------------------------------------------------------------------
    // These constants are required for defining light sources (e.g., Blackbody radiation)
    // when performing full spectral rendering.
    // Standard: SI Units (m, kg, s, J, K) based on the 2019 redefinition.

    // Speed of light in vacuum c [m/s]
    constexpr Real SPEED_OF_LIGHT = 299792458.0;

    // Planck constant h [J*s]
    constexpr Real PLANCK_CONSTANT = 6.62607015e-34;

    // Boltzmann constant k_B [J/K]
    constexpr Real BOLTZMANN_CONSTANT = 1.380649e-23;

    // Unit conversions
    constexpr Real NM_TO_M = 1e-9; // Nanometers to Meters
    constexpr Real M_TO_NM = 1e9;  // Meters to Nanometers

    // -------------------------------------------------------------------------
    // 4. Spectral Constants
    // -------------------------------------------------------------------------

    // Visible light spectrum range (Unit: nm).
    // Aligned with the CIE 1931 Standard Observer range (approx. 360nm to 830nm).
    constexpr Real LAMBDA_MIN = 360.0;
    constexpr Real LAMBDA_MAX = 830.0;

} // namespace constants