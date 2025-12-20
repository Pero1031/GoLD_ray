#pragma once

#include <random>
#include <numbers>
#include <algorithm>
#include <cmath>

// --- Math / Geometry ---
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

#include "Constants.hpp"   // Include the constants header if available
#include "Core.hpp"

// template <typename>は型に依存しない関数を作るという宣言
// Use std::sqrt for scalar, glm::sqrt for vector.

/**
* @brief 汎用関数
* 
*/
namespace Utils {

    using Real = Constants::Real;

    // -------------------------------------------------------------------------
    // 1. Math Helpers
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

    // Linear Interpolation
    // Interpolates between 'a' and 'b' based on weight 't'.
    // (Note: GLM also provides glm::mix, but this is useful for scalar types)
    template <typename T>
    inline T lerp(T a, T b, Real t) {
        return a + t * (b - a);
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

    // -------------------------------------------------------------------------
    // 2. Random Number Generation (RNG)
    // -------------------------------------------------------------------------
    // Replaces legacy functions like drand48().
    // Uses modern C++ <random> which provides better distribution properties.
    // 'thread_local' ensures thread safety during parallel rendering.

    inline Real Random() {
        // Initialize generator once per thread
        static thread_local std::mt19937 generator(std::random_device{}());
        static thread_local std::uniform_real_distribution<Real> distribution(0.0, 1.0);
        return distribution(generator);
    }

    // Utility for Scattering
    // Returns a random float in range [min, max)
    inline Real Random(Real min, Real max) {
        return min + (max - min) * Random();
    }

    // Returns a random point inside a unit sphere (radius 1).
    // Used for diffuse scattering and rough reflection.
    inline rayt::Vector3 randomInUnitSphere() {
        while (true) {
            // Generate a random vector in [-1, 1]^3 cube
            rayt::Vector3 p(Random(Real(-1), Real(1)),
                Random(Real(-1), Real(1)),
                Random(Real(-1), Real(1)));

            // If inside sphere, return it
            if (glm::dot(p, p) < Real(1)) {
                return p;
            }
        }
    }

    /**
     * @brief Generates a random point inside a unit disk (radius 1) on the XY plane.
     * Z-coordinate is always 0.
     * Used for camera aperture sampling (Depth of Field).
     */
    inline rayt::Vector3 randomInUnitDisk() {
        while (true) {
            // Generate x, y in [-1, 1]
            rayt::Vector3 p(Utils::Random(Real(- 1.0), Real(1.0)), Utils::Random(Real(- 1.0), Real(1.0)), Real(0.0));
            // Check if the point is within the circle (x^2 + y^2 < 1)
            if (glm::length2(p) >= Real(1.0)) continue; // Retry if outside the disk
            return p;
        }
    }

    // Returns a random unit vector (on the sphere surface).
    // Used for ideal diffuse reflection (Lambertian).
    inline rayt::Vector3 randomUnitVector() {
        while (true) {
            rayt::Vector3 v = randomInUnitSphere();
            Real len2 = glm::dot(v, v);
            if (len2 > Real(1e-12)) return v / std::sqrt(len2);
        }
    }


    // -------------------------------------------------------------------------
    // 3. Debug & Utility
    // -------------------------------------------------------------------------

    // Gamma Correction (Linear space -> sRGB space)
    // Standard gamma value is approximately 2.2.
    // Should be applied to the final pixel color before saving to an image.
    inline Real linearToGamma(Real linearComponent) {
        if (linearComponent > 0) {
            return std::pow(linearComponent, 1.0 / 2.2);
        }
        return 0.0;
    }

    // -------------------------------------------------------------------------
    // 4. Angle Conversions
    // -------------------------------------------------------------------------

    // Convert degrees to radians
    template <typename T>
    inline T toRadians(T degrees) {
        return degrees * (Constants::PI / static_cast<T>(180.0));
    }

    // Convert radians to degrees
    template <typename T>
    inline T toDegrees(T radians) {
        return radians * (static_cast<T>(180.0) / Constants::PI);
    }

    // -------------------------------------------------------------------------
    // 5. Fresnel Equations (Conductor / Metal)
    // -------------------------------------------------------------------------

    /**
     * @brief Calculates Fresnel reflectance for conductors (metals).
     * * Uses the exact solution based on complex refractive indices (eta + i*k).
     * This is physically more accurate than Schlick's approximation, especially
     * for reproducing the color shift observed in metals at grazing angles.
     * * @param cosThetaI Cosine of the incident angle (N dot V).
     * @param eta       Real part of the refractive index (n) - per RGB channel.
     * @param k         Imaginary part of the refractive index (Extinction coefficient, k) - per RGB channel.
     * @return rayt::Vector3 Fresnel reflectance (Reflectance).
     */
    inline rayt::Vector3 fresnelConductor(Real cosThetaI, const rayt::Vector3& eta, const rayt::Vector3& k) {
        // Clamp cosine to [0, 1] to handle numerical errors
        cosThetaI = Utils::saturate(cosThetaI);

        // Pre-calculate frequently used terms
        Real cosThetaI2 = cosThetaI * cosThetaI;
        Real sinThetaI2 = Real(1.0) - cosThetaI2;

        rayt::Vector3 eta2 = eta * eta;
        rayt::Vector3 k2 = k * k;

        // --- Preparation for Complex Arithmetic ---
        // Calculate coefficients 'a' and 'b' based on:
        // t0 = eta^2 - k^2 - sin^2(theta)
        // a^2 + b^2 = sqrt(t0^2 + 4 * eta^2 * k^2)

        rayt::Vector3 t0 = eta2 - k2 - rayt::Vector3(sinThetaI2);
        rayt::Vector3 a2plusb2 = glm::sqrt(t0 * t0 + Real(4.0) * eta2 * k2);

        // a = sqrt(0.5 * (a^2 + b^2 + t0))
        rayt::Vector3 t1 = a2plusb2 + rayt::Vector3(cosThetaI2);
        rayt::Vector3 a = glm::sqrt(Real(0.5) * (a2plusb2 + t0));

        // --- Calculate Rs (S-polarized reflectance) ---
        // Rs = ((a^2 + b^2) + cos^2 - 2a*cos) / ((a^2 + b^2) + cos^2 + 2a*cos)
        rayt::Vector3 t2 = Real(2.0) * cosThetaI * a;
        rayt::Vector3 Rs = (t1 - t2) / (t1 + t2);

        // --- Calculate Rp (P-polarized reflectance) ---
        // Rp = Rs * ( (a^2+b^2)*cos^2 + sin^4 - 2a*cos*sin^2 ) / ...
        // Using a simplified form derived via algebraic manipulation:
        // Rp = Rs * (t3 - t4) / (t3 + t4)

        rayt::Vector3 t3 = cosThetaI2 * a2plusb2 + rayt::Vector3(sinThetaI2 * sinThetaI2);
        rayt::Vector3 t4 = t2 * sinThetaI2;
        rayt::Vector3 Rp = Rs * (t3 - t4) / (t3 + t4);

        // Return the average for unpolarized light
        return Real(0.5) * (Rs + Rp);
    }

    // -------------------------------------------------------------------------
    // 6. Coordinate System (Orthonormal Basis)
    // -------------------------------------------------------------------------

    /**
     * @brief Builds an orthonormal basis (T, B, N) from a given normal N.
     * * A robust, branchless orthonormal basis construction (Frisvad-style).
     * Essential for transforming sampled vectors from local space to world space.
     * * @param N  The unit normal vector (z-axis of the local frame).
     * @param T  (Output) The computed tangent vector.
     * @param B  (Output) The computed bitangent vector.
     */
    inline void makeOrthonormalBasis(const rayt::Vector3& N, rayt::Vector3& T, rayt::Vector3& B) {
        Real sign = std::copysign(Real(1.0), N.z);
        Real a = -Real(1.0) / (sign + N.z);
        Real b = N.x * N.y * a;

        T = rayt::Vector3(Real(1.0) + sign * N.x * N.x * a, sign * b, -sign * N.x);
        B = rayt::Vector3(b, sign + N.y * N.y * a, -N.y);
    }

    // -------------------------------------------------------------------------
    // 7. Fresnel Equations (Dielectric)
    // -------------------------------------------------------------------------

    /**
     * @brief Calculates Fresnel reflectance for dielectrics (glass, water, coating).
     * * @param cosThetaI Cosine of the incident angle (positive).
     * @param etaI      Refractive index of the incident medium (usually 1.0 for air).
     * @param etaT      Refractive index of the transmission medium (e.g., 1.5 for glass).
     * @return Real     Fresnel reflectance (probability of reflection).
     */
    inline Real fresnelDielectric(Real cosThetaI, Real etaI, Real etaT) {
        cosThetaI = Utils::saturate(cosThetaI);

        // Check for total internal reflection
        Real sinThetaI = std::sqrt(std::max(Real(0.0), Real(1.0) - cosThetaI * cosThetaI));
        Real sinThetaT = (etaI / etaT) * sinThetaI;

        if (sinThetaT >= Real(1.0)) {
            return Real(1.0); // Total Internal Reflection
        }

        Real cosThetaT = std::sqrt(std::max(Real(0.0), Real(1.0) - sinThetaT * sinThetaT));

        Real rParl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
            ((etaT * cosThetaI) + (etaI * cosThetaT));

        Real rPerp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
            ((etaI * cosThetaI) + (etaT * cosThetaT));

        return (rParl * rParl + rPerp * rPerp) * Real(0.5);
    }

    // -------------------------------------------------------------------------
    // 8. Color Utilities
    // -------------------------------------------------------------------------

    /**
     * @brief Calculates the luminance of a linear RGB color.
     * * Uses standard coefficients (Rec. 709) for human vision sensitivity.
     * Green contributes the most, followed by Red, then Blue.
     */
    inline Real luminance(const rayt::Vector3& color) {
        return glm::dot(color, rayt::Vector3(Real(0.2126), Real(0.7152), Real(0.0722)));
    }

    //--------------------------------------------------------------------------
    // 9. vector operation
    // -------------------------------------------------------------------------
    // 
    // Reflect an incident vector around a surface normal.
    // v: incident direction (pointing *towards* the surface)
    // n: surface normal (must be normalized)
    inline rayt::Vector3 reflect(
        const rayt::Vector3& v,
        const rayt::Vector3& n)
    {
        return v - Real(2) * glm::dot(v, n) * n;
    }

    // Compute refraction using Snell's law.
    // Returns false if total internal reflection occurs.
    inline bool refract(
        const rayt::Vector3& v,
        const rayt::Vector3& n,
        Real eta,                 // relative index of refraction (eta_i / eta_t)
        rayt::Vector3& refracted) // output direction
    {
        Real cosi = glm::dot(-v, n);
        Real sin2_t = eta * eta * (Real(1) - cosi * cosi);

        if (sin2_t > Real(1)) {
            return false; // Total internal reflection
        }

        Real cost = safe_sqrt(Real(1) - sin2_t);
        refracted = eta * v + (eta * cosi - cost) * n;

        return true;
    }

} // namespace Utils