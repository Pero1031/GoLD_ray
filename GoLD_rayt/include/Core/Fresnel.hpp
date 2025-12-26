#pragma once

#include <numbers>
#include <algorithm>
#include <cmath>

// --- Math / Geometry ---
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Constants.hpp"   // Include the constants header if available
#include "Core/Math.hpp"

namespace rayt::fresnel {

	using Real = rayt::Real;

    // -------------------------------------------------------------------------
    // Fresnel Equations (Conductor / Metal)
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
    inline Vector3 fresnelConductor(Real cosThetaI, const Vector3& eta, const Vector3& k) {
        // Clamp cosine to [0, 1] to handle numerical errors
        cosThetaI = rayt::math::saturate(cosThetaI);

        // Pre-calculate frequently used terms
        Real cosThetaI2 = cosThetaI * cosThetaI;
        Real sinThetaI2 = Real(1.0) - cosThetaI2;

        Vector3 eta2 = eta * eta;
        Vector3 k2 = k * k;

        // --- Preparation for Complex Arithmetic ---
        // Calculate coefficients 'a' and 'b' based on:
        // t0 = eta^2 - k^2 - sin^2(theta)
        // a^2 + b^2 = sqrt(t0^2 + 4 * eta^2 * k^2)

        Vector3 t0 = eta2 - k2 - Vector3(sinThetaI2);
        Vector3 a2plusb2 = glm::sqrt(t0 * t0 + Real(4.0) * eta2 * k2);

        // a = sqrt(0.5 * (a^2 + b^2 + t0))
        Vector3 t1 = a2plusb2 + Vector3(cosThetaI2);
        Vector3 a = glm::sqrt(Real(0.5) * (a2plusb2 + t0));

        // --- Calculate Rs (S-polarized reflectance) ---
        // Rs = ((a^2 + b^2) + cos^2 - 2a*cos) / ((a^2 + b^2) + cos^2 + 2a*cos)
        Vector3 t2 = Real(2.0) * cosThetaI * a;
        Vector3 Rs = (t1 - t2) / (t1 + t2);

        // --- Calculate Rp (P-polarized reflectance) ---
        // Rp = Rs * ( (a^2+b^2)*cos^2 + sin^4 - 2a*cos*sin^2 ) / ...
        // Using a simplified form derived via algebraic manipulation:
        // Rp = Rs * (t3 - t4) / (t3 + t4)

        Vector3 t3 = cosThetaI2 * a2plusb2 + Vector3(sinThetaI2 * sinThetaI2);
        Vector3 t4 = t2 * sinThetaI2;
        Vector3 Rp = Rs * (t3 - t4) / (t3 + t4);

        // Return the average for unpolarized light
        return Real(0.5) * (Rs + Rp);
    }

    // -------------------------------------------------------------------------
    // Fresnel Equations (Dielectric)
    // -------------------------------------------------------------------------

    /**
     * @brief Calculates Fresnel reflectance for dielectrics (glass, water, coating).
     * * @param cosThetaI Cosine of the incident angle (positive).
     * @param etaI      Refractive index of the incident medium (usually 1.0 for air).
     * @param etaT      Refractive index of the transmission medium (e.g., 1.5 for glass).
     * @return Real     Fresnel reflectance (probability of reflection).
     */
    inline Real fresnelDielectric(Real cosThetaI, Real etaI, Real etaT) {
        cosThetaI = rayt::math::saturate(cosThetaI);

        // Check for total internal reflection
        Real sinThetaI = glm::sqrt(std::max(Real(0.0), Real(1.0) - cosThetaI * cosThetaI));
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
}