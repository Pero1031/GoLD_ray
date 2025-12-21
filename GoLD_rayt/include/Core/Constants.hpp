#pragma once

#include <limits>
#include <numbers>

/*
* @brief 計算で使う定数
*/
namespace Constants {

    // -------------------------------------------------------------------------
    // 1. Mathematical Constants
    // -------------------------------------------------------------------------

    // Using alias 'Real' allows easy switching between double and float precision.
    // For scientific research and spectral rendering, 'double' is recommended.
    using Real = double;

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
    constexpr Real TWO_PI = 2.0 * PI;             // Full circle (360 deg) in radians
    constexpr Real FOUR_PI = 4.0 * PI;            // Solid angle of a sphere
    constexpr Real INV_TWO_PI = 1.0 / TWO_PI;
    constexpr Real INV_FOUR_PI = 1.0 / FOUR_PI;   // Normalization factor for isotropic phase functions 1/4π

    // Angle conversions
    constexpr Real DEG_TO_RAD = PI / 180.0;
    constexpr Real RAD_TO_DEG = 180.0 / PI;


    // -------------------------------------------------------------------------
    // 2. Computational Constants
    // -------------------------------------------------------------------------

    // Representation of infinity for initial ray distances, etc.
    constexpr Real INFINITY_VAL = std::numeric_limits<Real>::infinity();

    // Small epsilon to prevent self-intersection artifacts (Shadow Acne).
    // The optimal value depends on the scene scale; 1e-4 to 1e-6 is typical.
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

    // Speed of light in vacuum c [m/s]
    constexpr Real SPEED_OF_LIGHT = 299792458.0;

    // Planck constant h [J*s]
    constexpr Real PLANCK_CONSTANT = 6.62607015e-34;

    // Boltzmann constant k_B [J/K]
    constexpr Real BOLTZMANN_CONSTANT = 1.380649e-23;

    // -------------------------------------------------------------------------
    // 4. Spectral Constants
    // -------------------------------------------------------------------------

    // Visible light range definition (Unit: nm)
    constexpr Real LAMBDA_MIN = 360.0;
    constexpr Real LAMBDA_MAX = 830.0;

    // -------------------------------------------------------------------------
    // 5. Utility Functions (Compile-time)
    // -------------------------------------------------------------------------
    // 便利なヘルパー関数（度数法 <-> 弧度法）
    // Converts degrees to radians.
    /*constexpr Real toRadians(Real degrees) {
        return degrees * DEG_TO_RAD;
    }

    // Converts radians to degrees.
    constexpr Real toDegrees(Real radians) {
        return radians * RAD_TO_DEG;
    }*/
}