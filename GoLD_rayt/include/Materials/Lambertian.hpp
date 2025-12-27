#pragma once

/**
 * @file Lambertian.hpp
 * @brief Ideal diffuse (Lambertian) material implementation.
 * * Represents a perfectly matte surface that scatters light uniformly in all
 * directions. This implementation uses cosine-weighted importance sampling
 * to significantly reduce noise in diffuse global illumination.
 */

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Math.hpp"
#include "Materials/Material.hpp"
#include "Geometry/Frame.hpp"
#include "Core/Sampling.hpp"

namespace rayt {

    /**
     * @brief Perfectly diffuse material following Lambert's Cosine Law.
     * * The BRDF is a constant value: f = albedo / PI.
     */
    class Lambertian : public Material {
        Spectrum albedo;

        static inline Real cosNg(const SurfaceInteraction& rec, const Vector3& w) {
            return glm::dot(rec.gn, w);
        }

        static inline Real cosNs(const SurfaceInteraction& rec, const Vector3& w) {
            return glm::dot(rec.n, w);
        }

    public:
        /**
         * @brief Constructs a Lambertian material.
         * @param a The reflectance (color) of the surface.
         */
        Lambertian(const Spectrum& a) : albedo(a) {}

        /**
         * @brief Evaluates the Lambertian BRDF.
         * * Formula: f(wo, wi) = albedo / PI
         * @return The constant spectral reflectance divided by PI.
         */
        Spectrum eval(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            TransportMode mode) const override {

            // Rejects light coming from behind the surface (back-face)
            if (cosNg(rec, wi) <= 0) {  // using geometry normal
                return Spectrum(0.0);
            }

            return albedo * (1.0 / constants::PI);
        }

        /**
         * @brief Evaluates the PDF for cosine-weighted sampling.
         * * Formula: p(wi) = cos(theta) / PI
         * @return The probability density of choosing direction wi.
         */
        Real pdf(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi) const override {
            if (cosNg(rec, wi) <= 0) {
                return 0.0;
            }

            Real cosTheta = cosNs(rec, wi);
            if (cosTheta <= 0) {
                return 0.0;
            }

            return cosTheta * (1.0 / constants::PI);
        }

        /**
         * @brief Importance samples the hemisphere using a cosine distribution.
         * * This method aligns the sample density with the cosine term of the
         * rendering equation, resulting in zero-variance for the cosine term.
         * * @param u Random samples in [0, 1)^2.
         * @return A BSDFSample containing the direction, value, and PDF.
         */
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec,
            const Vector3& wo,
            const Point2& u,
            TransportMode mode) const override {

            BSDFSample bsdfSample;

            // 1. Generate a direction in local space using cosine-weighted sampling.
            //    Z-axis in local space corresponds to the surface normal.
            rayt::frame::Frame onb(rec.n);
            Vector3 localDir = sampling::CosineSampleHemisphere(u); 

            // 2. Transform the sampled direction to World Space.
            bsdfSample.wi = onb.localToWorld(localDir);

            // 3. Geometric sanity check: Ensure the direction is in the upper hemisphere.
            if (cosNg(rec, bsdfSample.wi) <= 0) return std::nullopt;

            // 4. PDF calculation: p(wi) = cos(theta) / PI.
            //    Since localDir is on a unit hemisphere, localDir.z is exactly cos(theta).
            bsdfSample.pdf = localDir.z * (Real(1.0) / constants::PI);

            // 5. BSDF value: constant for Lambertian surfaces.
            bsdfSample.f = albedo *(Real(1.0) / constants::PI);

            // 6. Set flags indicating a diffuse reflection interaction.
            bsdfSample.flags =
                BxDFFlags::Diffuse |
                BxDFFlags::Reflection;

            return bsdfSample;
        }
    };

} // namespace rayt