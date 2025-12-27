#pragma once

/**
 * @file MirrorConductor.hpp
 * @brief Perfectly specular conductor (metal) material.
 * * This material simulates smooth metallic surfaces using complex Fresnel
 * equations. By providing the refractive index (eta) and extinction
 * coefficient (k), it can accurately render the characteristic colors of
 * metals like Gold, Copper, or Chromium.
 */

#include <complex>

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Fresnel.hpp"
#include "Materials/Material.hpp"
#include "Core/Interaction.hpp"

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

            // Geometric sanity: must go above the geometric surface
            if (glm::dot(rec.gn, bsdfSample.wi) <= Real(0))
                return std::nullopt;

            // 2. Metadata Flags
            // Marked as SPECULAR and REFLECTION for the Integrator's path handling.
            // Flags (new design)
            bsdfSample.flags = BxDFFlags::Specular | BxDFFlags::Reflection;
            bsdfSample.pdf = 1.0; // Dummy value for delta distribution

            // Fresnel conductor reflectance
            Real cosTheta = math::saturate(glm::dot(rec.n, bsdfSample.wi));

            // 3. Complex Fresnel Throughput (f)
            // Calculates the exact reflectance based on the angle of incidence.
            // Since this is a delta distribution, we return (f / cosTheta) if 
            // the integrator multiplies by cosTheta, but PBRT style usually 
            // returns just the reflectance factor here.
            bsdfSample.f = fresnel::fresnelConductor(cosTheta, eta, k);

            return bsdfSample;
        }

        /**
         * @brief Checks if the material is specular.
         */
        bool isSpecular() const override { return true; }
    
    };

} // namespace rayt