/**
 * @file Core/Fresnel.hpp
 * @brief Fresnel reflection utilities for conductors and dielectrics.
 *
 * Provides physically-based Fresnel reflectance evaluation:
 * - Conductor (metal): exact formulation using complex IOR (eta + i*k)
 * - Dielectric: exact unpolarized Fresnel reflectance with TIR handling
 *
 * @note This module is currently RGB-oriented. It will be extended to support
 *       spectral rendering by evaluating eta(λ), k(λ) over sampled wavelengths.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"   
#include "Core/Math.hpp"
#include "Kernel/Microfacet.hpp"

namespace rayt::fresnel {
    // -------------------------------------------------------------------------
    // Fresnel Equations (Conductor / Metal)
    // -------------------------------------------------------------------------

    // The Vector3 overload below is for RGB; the scalar overload (per wavelength) is what the spectral build calls.

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
        // Per channel, by the shared formula (Kernel/Microfacet.hpp, also on the GPU): |cos| clamped to [0, 1];
        // t0 = eta^2 - k^2 - sin^2, a^2 + b^2 = sqrt(t0^2 + 4 eta^2 k^2), a = sqrt((a^2 + b^2 + t0) / 2),
        // Rs = (a^2 + b^2 + cos^2 - 2 a cos) / (... + 2 a cos), Rp = Rs (t3 - t4) / (t3 + t4) with t3 = cos^2 (a^2 + b^2) + sin^4,
        // t4 = 2 a cos sin^2; the result is (Rs + Rp) / 2.
        return Vector3(kernel::fresnelConductor(cosThetaI, eta.x, k.x), kernel::fresnelConductor(cosThetaI, eta.y, k.y),
            kernel::fresnelConductor(cosThetaI, eta.z, k.z));
    }

    inline Real fresnelConductor(Real cosThetaI, Real eta, Real k) {
        return kernel::fresnelConductor(cosThetaI, eta, k);
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
        cosThetaI = math::saturate(cosThetaI);

        // Check for total internal reflection
        Real sinThetaI = math::safe_sqrt(Real(1.0) - math::sqr(cosThetaI));
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