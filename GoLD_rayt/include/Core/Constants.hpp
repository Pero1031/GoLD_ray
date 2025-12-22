#pragma once

#include <limits>
#include <numbers>

#include "Core/Types.hpp"  // to use Real = double

/*
* @brief Constants for calculations
* @note Define Real as double consistently throughout the project.
*/
namespace rayt::constants {

    // -------------------------------------------------------------------------
    // 1. Mathematical Constants
    // -------------------------------------------------------------------------

    // Using alias 'Real' allows easy switching between double and float precision.
    // For scientific research and spectral rendering, 'double' is recommended.
    using Real = rayt::Real;

#if defined(__cpp_lib_math_constants)
    // C++20: Use standard library constants for maximum precision and portability.
    constexpr Real PI = std::numbers::pi_v<Real>;
    constexpr Real INV_PI = std::numbers::inv_pi_v<Real>;          // 1/π
    constexpr Real INV_SQRT_PI = std::numbers::inv_sqrtpi_v<Real>; // 1/√π
#else
    // Pre-C++20 fallback: Define constants manually with high precision.
    constexpr Real PI = 3.14159265358979323846;
    constexpr Real INV_PI = 1.0 / PI;
    constexpr Real INV_SQRT_PI = 0.564189583547756286948;
#endif

    // Derived constants frequently used in ray tracing equations.
    constexpr Real TWO_PI = Real(2) * PI;             // Full circle (360 deg) in radians
    constexpr Real FOUR_PI = Real(4) * PI;            // Solid angle of a sphere
    constexpr Real INV_TWO_PI = Real(1) / TWO_PI;
    constexpr Real INV_FOUR_PI = Real(1) / FOUR_PI;   // Normalization factor for isotropic phase functions 1/4π

    // Angle conversions
    constexpr Real DEG_TO_RAD = PI / Real(180.0);
    constexpr Real RAD_TO_DEG = Real(180.0) / PI;


    // -------------------------------------------------------------------------
    // 2. Computational Constants
    // -------------------------------------------------------------------------

    // Representation of infinity for initial ray distances, etc.
    constexpr Real INFINITY_VAL = std::numeric_limits<Real>::infinity();

    // Small epsilon to prevent self-intersection artifacts (Shadow Acne).
    // CAUTION: Depends on scene scale. 1e-5 is good for meter-scale scenes.
    // For microscopic/millimeter scale, consider 1e-7.
    constexpr Real RAY_EPSILON = 1e-5;

    // Tolerance for intersection tests.
    constexpr Real INTERSECT_TOLERANCE = 1e-6;


    // -------------------------------------------------------------------------
    // 3. Physical Constants for Spectral Rendering
    // -------------------------------------------------------------------------
    // Johnsonのデータを使う＝スペクトルレンダリングを行う場合、
    // 光源（黒体放射など）の定義に以下の物理定数が必要になります。
    // ※ 単位は SI系 (m, kg, s, J, K) に準拠
    // Essential for spectral rendering and blackbody radiation calculations.
    // Units: meter (m), kilogram (kg), second (s), Joule (J), Kelvin (K)
    // Values based on 2019 SI redefinition.

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

    // Visible light range definition (Unit: nm)
    // Covers the full CIE 1931 standard observer range.
    constexpr Real LAMBDA_MIN = 360.0;
    constexpr Real LAMBDA_MAX = 830.0;

} // namespace constants