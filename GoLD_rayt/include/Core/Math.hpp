/**
 * @file Core/Math.hpp
 * @brief Core mathematical utilities and geometric routines for the rayt engine.
 *
 * This header provides a comprehensive suite of math helpers designed for
 * robustness in path tracing. Key features include:
 * - Numerical Stability: "Safe" wrappers for sqrt and trigonometric functions.
 * - Multi-Style Convention: Support for both incident-based and outward-based routines.
 * - Optimized Primitives: Efficient implementations of common power functions.
 * - Validation: Comprehensive NaN and Infinity checks.
 *
 * @note All functions are templated or inlined to ensure zero-overhead
 * integration with the rendering pipeline.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include "Core/Types.hpp"   
#include "Core/Assert.hpp"  
#include "Core/Constants.hpp"

namespace rayt::math {

    // -------------------------------------------------------------------------
    // Basic Helpers
    // -------------------------------------------------------------------------

    /// @brief Calculates the square of a value (x * x).
    /// @details Significantly faster than std::pow(x, 2) due to avoiding function call overhead and branching.
    template <typename T>
    inline T sqr(T x) {
        return x * x;
    }

    /// @brief Calculates the 5th power of a value (x^5).
    /// @details Efficient implementation (x^2 * x^2 * x).
    /// Frequently used in Fresnel equations (Schlick's approximation) for conductor/dielectric blending.
    template <typename T>
    inline T pow5(T x) {
        T x2 = x * x;
        return x2 * x2 * x;
    }

    /// @brief Linearly interpolates between 'a' and 'b' based on weight 't'.
    /// @details Formula: (1 - t) * a + t * b
    /// Useful for scalar types where glm::mix might be overkill or ambiguous.
    template <typename T, typename U>
    inline T lerp(const T& a, const T& b, U t) {
        return (static_cast<U>(1) - t) * a + t * b;
    }

    /// @brief Clamps a value to the [0, 1] range.
    /// @details Essential for keeping energy/color values within valid physical ranges or HDR tone mapping.
    template <typename T>
    inline T saturate(T x) {
        return glm::clamp(x, static_cast<T>(0), static_cast<T>(1));
    }

    /// @brief Calculates a safe reciprocal (1/x).
    /// @details Avoids division by zero (Infinity/NaN) by returning a large value with the correct sign.
    /// Critical for robust shading when the dot product (N dot V) approaches zero.
    inline Real safe_recip(Real x) {
        constexpr Real EPS = Real(1e-6);
        if (std::abs(x) < EPS) {
            // Divide by epsilon preserving the sign to avoid +/- Inf
            return Real(1.0) / (x >= Real(0) ? EPS : -EPS);
        }
        return 1.0 / x;
    }

    /// @brief Checks if a value is effectively zero within a given epsilon.
    /// @details Prefer this over (x == 0) for floating-point safety.
    template <typename T>
    inline bool isZero(T x, T eps = T(1e-6)) {
        return std::abs(x) < eps;
    }

    /// @brief Checks if two floating-point values are nearly equal.
    /// @details Essential for geometry intersections and convergence checks to handle floating-point drift.
    template <typename T>
    inline bool nearlyEqual(T a, T b, T eps = T(1e-6)) {
        return std::abs(a - b) < eps;
    }

    // -------------------------------------------------------------------------
    // Safe Elementary Functions
    // -------------------------------------------------------------------------

    /// @brief Calculates the square root safely by clamping negative inputs to zero.
    /// @details Prevents NaN generation when inputs are theoretically non-negative (e.g., 1 - dot^2) 
    /// but drift to slightly negative values (e.g., -1e-18) due to floating-point errors.
    /// @param x Input value.
    /// @return std::sqrt(max(0, x))
    template <typename T>
    inline T safe_sqrt(T x) {
        return std::sqrt(std::max(x, T(0)));
    }

    /// @brief Component-wise safe square root for Vector3.
    /// @details Applies safe_sqrt logic to each component independently. 
    /// Useful for manipulating spectral data or color vectors where components must be non-negative.
    /// @param v Input vector.
    /// @return Vector with component-wise square roots.
    inline Vector3 safe_sqrt(const Vector3& v) {
        // glm::max compares component-wise, glm::sqrt computes component-wise.
        return glm::sqrt(glm::max(v, Vector3(0.0)));
    }

    // -------------------------------------------------------------------------
    // Safe Inverse Trigonometric Functions
    // -------------------------------------------------------------------------

    /// @brief Safe arc sine (inverse sine). Clamps input to the [-1, 1] range.
    /// @details Prevents NaN when inputs slightly exceed range (e.g., 1.0000001) due to precision errors.
    /// Includes an assertion to catch significant logical errors (values far outside valid range).
    /// @param x Input value.
    /// @return std::asin(clamp(x, -1, 1))
    template <typename T>
    inline T safe_asin(T x) {
        // Use asserts to catch logical errors (clear out-of-bounds),
        // while clamping to allow for minor floating-point inaccuracies.
        Assert(x >= T(-1.0001) && x <= T(1.0001));
        return std::asin(std::clamp(x, T(-1), T(1)));
    }

    /// @brief Safe arc cosine (inverse cosine). Clamps input to the [-1, 1] range.
    /// @details Essential for spherical coordinate conversions (e.g., retrieving theta from z-coordinate).
    /// Prevents NaN when input is slightly > 1.0 or < -1.0.
    /// @param x Input value.
    /// @return std::acos(clamp(x, -1, 1))
    template <typename T>
    inline T safe_acos(T x) {
        Assert(x >= T(-1.0001) && x <= T(1.0001));
        return std::acos(std::clamp(x, T(-1), T(1)));
    }

    // -------------------------------------------------------------------------
    // Angle Conversions
    // -------------------------------------------------------------------------
    
    /// @brief Converts degrees to radians.
    /// @param degrees Angle in degrees (e.g., 90.0).
    /// @return Angle in radians (e.g., PI/2).
    template <typename T>
    inline T toRadians(T degrees) {
        return degrees * (constants::PI / static_cast<T>(180.0));
    }

    /// @brief Converts radians to degrees.
    /// @param radians Angle in radians.
    /// @return Angle in degrees.
    template <typename T>
    inline T toDegrees(T radians) {
        return radians * (static_cast<T>(180.0) / constants::PI);
    }

    // -------------------------------------------------------------------------
    // Reflection / Refraction (Convention Helpers)
    // -------------------------------------------------------------------------
    
    /// @brief Reflects a vector based on the "Incident" convention (GLSL style).
    /// @details 
    /// Input 'v' points TOWARD the surface.
    /// Result points AWAY from the surface.
    /// Formula: v - 2 * dot(v, n) * n
    /// @param v Incident vector (Points TOWARD the surface).
    /// @param n Surface Normal (Normalized).
    /// @return Reflected vector.
    inline Vector3 reflectIncident(const Vector3& v, const Vector3& n) {
        return v - Real(2) * glm::dot(v, n) * n;
    }

    /// @brief Reflects a vector based on the "Outward" convention (PBRT/Mitsuba style).
    /// @details 
    /// Input 'wo' points AWAY from the surface (e.g., towards the camera).
    /// Result points AWAY from the surface (e.g., towards the light/environment).
    /// Useful for BSDF evaluation where all vectors are in the local hemisphere.
    /// @param wo Outgoing direction (Points AWAY from surface).
    /// @param n Surface Normal (Normalized).
    /// @return Reflected direction (Incoming light direction, pointing AWAY).
    inline Vector3 reflectOutward(const Vector3& wo, const Vector3& n) {
        // Invert the input (-w) to get the "piercing direction".
        // Reflecting this direction naturally yields the "vector pointing toward the light source (Away)".
        return reflectIncident(-wo, n);
    }

    /// @brief Refracts an incident vector using Snell's Law (GLSL style).
    /// @details 
    /// Input 'v' points TOWARD the surface.
    /// Output 'refracted' points AWAY from the surface (into the volume).
    /// Handles Total Internal Reflection (TIR).
    /// @param v Incident vector (Points TOWARD the surface).
    /// @param n Surface Normal (Normalized).
    /// @param eta Ratio of indices (eta_i / eta_t). Example: 1.0/1.5 for Air->Glass.
    /// @param[out] refracted The computed refracted vector.
    /// @return false if Total Internal Reflection (TIR) occurs, true otherwise.
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

    /// @brief Refracts an outward vector (PBRT/Research style).
    /// @details
    /// Both 'wo' and 'wi' point AWAY from the shading point.
    /// This is a wrapper around refractIncident that handles direction flipping.
    /// @param wo Outgoing direction (Points AWAY from surface).
    /// @param n Geometric normal (Normalized).
    /// @param eta Ratio of indices (eta_i / eta_t).
    /// @param[out] wi Incoming direction (Points AWAY from surface, into the volume).
    /// @return false if Total Internal Reflection (TIR) occurs, true otherwise.
    inline bool refractOutward(
        const Vector3& wo,
        const Vector3& n,
        Real eta,
        Vector3& wi)
    {
        // Convert outward to incident-toward and reuse the robust core routine.
        return refractIncident(-wo, n, eta, wi);
    }

    // -------------------------------------------------------------------------
    // Non-finite Checks (NaN/Inf)
    // -------------------------------------------------------------------------
    
    /// @brief Checks if a scalar value is NOT finite (i.e., is NaN or Infinity).
    /// @details Returns true if the value is corrupted. 
    /// Essential for debugging "fireflies" (bright spots) or black NaN-propagation in the image.
    /// @param x The value to check.
    /// @return true if x is NaN or Infinity (Invalid), false if finite (Valid).
    template <typename T>
    inline bool hasNonFinite(T x) {
        return std::isnan(x) || std::isinf(x);
    }

    /// @brief Checks if any component of a GLM vector is NOT finite.
    /// @details Works for any GLM vector type (vec2, vec3, vec4, etc.).
    /// Uses glm::any to detect if at least one component is corrupted.
    /// @param v Input vector (e.g., Radiance, Position, Normal).
    /// @return true if the vector contains any NaN or Infinity components.
    template <int L, typename T, glm::qualifier Q>
    inline bool hasNonFinite(const glm::vec<L, T, Q>& v) {
        // glm::isnan returns a boolean vector (bvec).
        // glm::any checks if any component in the bvec is true.
        return glm::any(glm::isnan(v)) || glm::any(glm::isinf(v));
    }

    /// @brief Returns the maximum component of a 3D vector.
    /// @details Frequently used for:
    /// - Calculating Russian Roulette termination probabilities (max component of throughput).
    /// - Normalizing spectral power distributions.
    /// @param v Input vector (typically Spectrum or Throughput).
    /// @return The value of the largest component (std::max({x, y, z})).
    template <typename T>
    inline T maxComponent(const glm::vec<3, T, glm::defaultp>& v) {
        return std::max({ v.x, v.y, v.z });
    }

    // -------------------------------------------------------------------------
    // Multiple Importance Sampling (MIS)
    // -------------------------------------------------------------------------

    /// @brief Power heuristic for Multiple Importance Sampling (MIS) with exponent = 2.
    /// @details
    /// Balances contributions from two sampling strategies (e.g., Light sampling vs. BSDF sampling).
    /// Formula: f^2 / (f^2 + g^2)
    /// @param pf Probability density (PDF) of the strategy used to generate the sample.
    /// @param pg Probability density (PDF) of the alternative strategy.
    /// @return The weight [0, 1] to apply to the sample contribution.
    inline Real powerHeuristic(Real pf, Real pg) {
        Real f2 = pf * pf;
        Real g2 = pg * pg;

        // Prevent Divide-By-Zero if both PDFs are zero (e.g., singular delta distributions).
        if (f2 == Real(0) && g2 == Real(0)) return Real(0);

        return f2 / (f2 + g2);
    }

} // namespace rayt::math