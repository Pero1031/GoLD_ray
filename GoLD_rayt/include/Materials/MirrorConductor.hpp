#pragma once

/**
 * @file MirrorConductor.hpp
 * @brief Perfectly specular conductor (metal) material.
 * * This material simulates smooth metallic surfaces using complex Fresnel
 * equations. By providing the refractive index (eta) and extinction
 * coefficient (k), it can accurately render the characteristic colors of
 * metals like Gold, Copper, or Chromium.
 */

#include "Core/Core.hpp"
#include "Materials/Material.hpp"
#include "Core/Interaction.hpp"
#include <complex>
#include <algorithm> // std::clampに必要
#include <cmath>

namespace rayt {

    /**
     * @brief A mirror-like conductor material.
     * * Handles perfect specular reflection with wavelength-dependent (RGB)
     * Fresnel reflectance for metallic surfaces.
     */
    class MirrorConductor : public Material {
        Spectrum eta; // Real part of the refractive index 'n'
        Spectrum k;   // Imaginary part (extinction coefficient) 'k'

    public:
        /**
         * @brief Constructs a conductor with complex IOR values.
         * @param eta Real part of IOR (RGB).
         * @param k   Extinction coefficient (RGB).
         */
        MirrorConductor(const Spectrum& eta, const Spectrum& k) : eta(eta), k(k) {}

        /**
         * @brief Evaluation for specular materials returns zero.
         * * Specular reflections are Dirac delta distributions and cannot be
         * evaluated via standard sampling; they must be generated via sample().
         */
        Spectrum eval(const SurfaceInteraction&, const Vector3&, const Vector3&, TransportMode) const override {
            return Spectrum(0.0);
        }

        /**
         * @brief PDF for a delta distribution is 0.
         */
        Real pdf(const SurfaceInteraction&, const Vector3&, const Vector3&) const override {
            return 0.0;
        }

        /**
         * @brief Samples the perfect specular reflection direction.
         * * Computes the reflection vector and scales the throughput by the
         * complex Fresnel reflectance for each RGB channel.
         */
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec,
            const Vector3& wo,
            const Point2& u,
            TransportMode mode) const override {

            BSDFSample bsdfSample;

            // 1. Direction Calculation: Perfect Specular Reflection
            // Reflecting the outgoing vector 'wo' about the normal 'n'.
            // Formula: wi = -wo + 2(n · wo)n
            bsdfSample.wi = glm::reflect(-wo, rec.n);

            // Geometric sanity check (ensure reflection is above the surface)
            Real cosTheta = glm::dot(bsdfSample.wi, rec.n);
            if (cosTheta <= 0) return std::nullopt;

            // 2. Metadata Flags
            // Marked as SPECULAR and REFLECTION for the Integrator's path handling.
            bsdfSample.sampledType = BxDFType(BSDF_SPECULAR | BSDF_REFLECTION);
            bsdfSample.pdf = 1.0; // Dummy value for delta distribution

            // 3. Complex Fresnel Throughput (f)
            // Calculates the exact reflectance based on the angle of incidence.
            // Since this is a delta distribution, we return (f / cosTheta) if 
            // the integrator multiplies by cosTheta, but PBRT style usually 
            // returns just the reflectance factor here.
            bsdfSample.f = Spectrum(
                fresnelConductorExact(cosTheta, eta.x, k.x),
                fresnelConductorExact(cosTheta, eta.y, k.y),
                fresnelConductorExact(cosTheta, eta.z, k.z)
            );

            return bsdfSample;
        }

        /**
         * @brief Checks if the material is specular.
         */
        bool isSpecular() const override { return true; }

    private:
        /**
         * @brief Exact Fresnel reflectance for conductors using complex arithmetic.
         * * Handles the transition from the medium (vacuum/air) to the conductor.
         * @param cosThetaI Cosine of the angle of incidence (n · wi).
         * @param etaVal    The real part of the IOR.
         * @param kVal      The extinction coefficient.
         * @return Real Average of s-polarized and p-polarized reflectance.
         */
        Real fresnelConductorExact(Real cosThetaI, Real etaVal, Real kVal) const {
            // Clamping to avoid numerical artifacts from precision errors.
            cosThetaI = std::clamp(cosThetaI, Real(0.0), Real(1.0));

            // Complex Refractive Index: N = eta + ik
            std::complex<Real> N(etaVal, kVal);

            Real cosThetaI2 = cosThetaI * cosThetaI;
            Real sinThetaI2 = Real(1.0) - cosThetaI2;

            // Compute complex cosThetaT using Snell's Law (Complex form)
            std::complex<Real> sinThetaT2 = sinThetaI2 / (N * N);
            std::complex<Real> cosThetaT = std::sqrt(Real(1.0) - sinThetaT2);

            // Fresnel coefficients for parallel (p) and perpendicular (s) polarization
            std::complex<Real> r_s = (cosThetaI - N * cosThetaT) / (cosThetaI + N * cosThetaT);
            std::complex<Real> r_p = (N * cosThetaI - cosThetaT) / (N * cosThetaI + cosThetaT);

            // ※注: 一般的な教科書では r_s の分子が (ni cos - nt cos) だったり (cos - N cos) だったり定義揺れがありますが、
            // 最終的に絶対値の二乗を取るため、符号の違いは結果に影響しません。
            // ここでは PBRT v3/v4 の実装に近い形式を採用しています。

            // Reflectance for unpolarized light is the average of squared magnitudes:
            // R = (|r_s|^2 + |r_p|^2) / 2
            return (std::norm(r_s) + std::norm(r_p)) * Real(0.5);
        }
    };

} // namespace rayt