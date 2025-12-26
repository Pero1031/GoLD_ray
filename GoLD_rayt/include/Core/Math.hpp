#pragma once

/**
 * @file Math.hpp
 * @brief Core mathematical utilities and geometric routines for the rayt engine.
 * * This header provides a comprehensive suite of math helpers designed for
 * robustness in path tracing. Key features include:
 * * - Numerical Stability: "Safe" wrappers for sqrt and trigonometric functions to
 * handle floating-point drift.
 * - Multi-Style Convention: Support for both incident-based (GLSL) and
 * outward-based (PBRT/Mitsuba) reflection/refraction routines.
 * - Optimized Primitives: Efficient implementations of common power functions
 * (sqr, pow5) and linear interpolation.
 * - Validation: Comprehensive NaN and Infinity checks for both scalars and
 * GLM vector types.
 * * @note All functions are templated or inlined to ensure zero-overhead
 * integration with the rendering pipeline.
 */

#include <algorithm>
#include <cmath>
#include <limits>

#include "Core/Types.hpp"   // for type def
#include "Core/Assert.hpp"  // assert
#include "Core/Constants.hpp"

#include <glm/glm.hpp>

namespace rayt::math {

    // -------------------------------------------------------------------------
    // Math Helpers
    // -------------------------------------------------------------------------

    // Square function (x^2)
    // Significantly faster than std::pow(x, 2).
    template <typename T>
    inline T sqr(T x) {
        return x * x;
    }

    // Safe square root.
    // Clamps negative inputs caused by numerical error to zero.
    template <typename T>
    inline T safe_sqrt(T x) {
        return std::sqrt(std::max(x, T(0)));
    }

    // 5th power function (x^5)
    // Frequently used in Fresnel equations (Schlick's approximation).
    // More efficient than std::pow(x, 5.0f).
    template <typename T>
    inline T pow5(T x) {
        T x2 = x * x;
        return x2 * x2 * x;
    }

    // Clamp to [0, 1] range
    // Essential for keeping color values within valid range.
    template <typename T>
    inline T saturate(T x) {
        return glm::clamp(x, static_cast<T>(0), static_cast<T>(1));
    }

    // Check if a value is close to zero (for floating-point safety).
    // Useful for avoiding direct comparisons like (x == 0).
    template <typename T>
    inline bool isZero(T x, T eps = T(1e-6)) {
        return std::abs(x) < eps;
    }

    // Check if two floating-point values are nearly equal.
    // Essential for robust geometric and shading computations.
    template <typename T>
    inline bool nearlyEqual(T a, T b, T eps = T(1e-6)) {
        return std::abs(a - b) < eps;
    }

    // Linear Interpolation
    // Interpolates between 'a' and 'b' based on weight 't'.
    // (Note: GLM also provides glm::mix, but this is useful for scalar types)
    template <typename T, typename U>
    inline T lerp(const T& a, const T& b, U t) {
        return (static_cast<U>(1) - t) * a + t * b;
    }

    // Safe Reciprocal
    // Returns 1/x but avoids division by zero (Infinity/NaN results).
    // Useful for preventing artifacts when the dot product (N dot V) is zero.
    inline Real safe_recip(Real x) {
        constexpr Real EPS = Real(1e-6);
        if (std::abs(x) < EPS) {
            // Divide by epsilon preserving the sign
            return Real(1.0) / (x >= Real(0) ? EPS : -EPS);
        }
        return 1.0 / x;
    }

    // -------------------------------------------------------------------------
    // Angle Conversions
    // -------------------------------------------------------------------------
    // Convert degrees to radians
    template <typename T>
    inline T toRadians(T degrees) {
        return degrees * (constants::PI / static_cast<T>(180.0));
    }

    // Convert radians to degrees
    template <typename T>
    inline T toDegrees(T radians) {
        return radians * (static_cast<T>(180.0) / constants::PI);
    }

    //--------------------------------------------------------------------------
    // vector operation
    // -------------------------------------------------------------------------
    // Reflect 
    // GLSL-style: v points TOWARD the surface (Incident)
    // Result points AWAY from surface
    inline Vector3 reflectIncident(const Vector3& v, const Vector3& n) {
        return v - Real(2) * glm::dot(v, n) * n;
    }

    // pbrt-style: wo points AWAY from the surface (Outgoing)
    // Result points AWAY from surface (Incoming / Light direction)
    inline Vector3 reflectOutward(const Vector3& wo, const Vector3& n) {
        // Invert the input (-w) to get the "piercing direction".
        // Reflecting this direction naturally yields the "vector pointing toward the light source (Away)".
        return reflectIncident(-wo, n);
    }

    // [GLSL / Classic Style]
    // Input:
    //   incident: Points TOWARD the surface 
    //   n:        Surface Normal (Normalized)
    //   eta:      ratio of indices (eta_i / eta_t)
    // Output:
    //   refracted: Points AWAY from the surface (into the volume)
    // Returns: false if Total Internal Reflection (TIR) occurs.
    inline bool refractIncident(
        const Vector3& v,
        const Vector3& n,
        Real eta,                 
        Vector3& refracted)       
    {
        Real cosi = glm::dot(-v, n);  // cos > 0
        Real sin2_t = sqr(eta) * (Real(1) - sqr(cosi));

        if (sin2_t > Real(1)) {
            return false; // Total internal reflection
        }

        Real cost = safe_sqrt(Real(1) - sin2_t);
        refracted = eta * v + (eta * cosi - cost) * n;

        return true;
    }

    // [PBRT / Research Style]
    // Input:
    //   wo: direction (outward convention). In pbrt/Mitsuba, both wo/wi are directions 
    //       pointing away from the shading point.
    //   n:  geometric normal (normalized). Must be oriented such that dot(n, wo) >= 0
    //      (use faceforward if needed).
    //   eta: ratio of indices (eta_i / eta_t)
    // Output:
    //   wi: Incoming direction (Points AWAY from surface, into the volume)
    // Returns: false if Total Internal Reflection (TIR) occurs.
    inline bool refractOutward(
        const Vector3& wo,
        const Vector3& n,
        Real eta,
        Vector3& wi)
    {
        // Convert outward to incident-toward and reuse the robust core routine.
        return refractIncident(-wo, n, eta, wi);
    }

    // Safe Inverse Trigonometric Functions
    // Prevents NaN generation when inputs slightly exceed [-1, 1] due to float precision errors.
    template <typename T>
    inline T safe_asin(T x) {
        // Use asserts to catch logical errors (clear out-of-bounds),
        // while clamping to allow for minor floating-point inaccuracies.
        Assert(x >= T(-1.0001) && x <= T(1.0001));
        return std::asin(std::clamp(x, T(-1), T(1)));
    }

    template <typename T>
    inline T safe_acos(T x) {
        Assert(x >= T(-1.0001) && x <= T(1.0001));
        return std::acos(std::clamp(x, T(-1), T(1)));
    }

    // -------------------------------------------------------------------------
    // NaN / Inf Checks
    // -------------------------------------------------------------------------
    // 1. スカラ値用 (Real / float / double)
    template <typename T>
    inline bool hasNaNs(T x) {
        return std::isnan(x) || std::isinf(x);
    }

    // 2. ベクトル用 (Vector3, Point3, Normal3, Spectrum 兼用)
    // GLMの型(vec3など)なら何でも受け取れるテンプレート
    template <int L, typename T, glm::qualifier Q>
    inline bool hasNaNs(const glm::vec<L, T, Q>& v) {
        // glm::isnan(v) は boolのベクトル(bvec)を返すので、
        // glm::any() を使って「どれか一つでもNaNならtrue」に
        return glm::any(glm::isnan(v)) || glm::any(glm::isinf(v));
    }

    //-----------------------------------------------------------------------------
    // ベクトル(Spectrum)の中で最大の成分を返す関数
    template <typename T>
    inline T maxComponent(const glm::vec<3, T, glm::defaultp>& v) {
        return std::max({ v.x, v.y, v.z });
    }

} // namespace rayt::math