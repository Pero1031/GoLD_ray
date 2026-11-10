/**
 * @file Materials/Lambertian.hpp
 * @brief Ideal diffuse (Lambertian) material implementation.
 * * Represents a perfectly matte surface that scatters light uniformly in all
 * directions. This implementation uses cosine-weighted importance sampling
 * to significantly reduce noise in diffuse global illumination.
 */

#pragma once

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Math.hpp"
#include "Core/Sampling.hpp"
#include "Materials/Material.hpp"
#include "Geometry/Frame.hpp"
#include "Color/SpectralConversion.hpp"
#include "Core/Texture.hpp"
#include "Core/SpectralTexture.hpp"

namespace rayt {

    /**
     * @brief Perfectly diffuse material following Lambert's Cosine Law.
     * * The BRDF is a constant value: f = albedo / PI.
     */
    class Lambertian : public Material {
        RgbColor albedo;
        std::shared_ptr<const Texture2D> m_tex;   // optional: albedo is multiplied by the texture color at rec.uv
        Real m_opacity = Real(1);                 // < 1: translucent layer (see Material::opacity)
        bool m_translucent = false;

        std::shared_ptr<const color::SpectralCurve> m_curve;   // optional: a measured reflectance R(lambda)
        std::shared_ptr<const SpectralTexture> m_stex;         // optional: a measured reflectance per texel R(u, v, lambda)

        RgbColor albedoAt(const SurfaceInteraction& rec) const {
            if (m_stex) return RgbColor(albedo.rgb() * m_stex->rgb(rec.uv.x, rec.uv.y));
            if (!m_tex) return albedo;
            return RgbColor(albedo.rgb() * m_tex->rgb(rec.uv.x, rec.uv.y));
        }

        /// The reflectance at the path's wavelengths: the measured curve or spectral texture if there is
        /// one (spectral build), else the RGB albedo. In the RGB build a curve is already stored as its D65 RGB.
        Spectrum reflectanceAt(const SurfaceInteraction& rec, const SampledWavelengths& lambda) const {
#if defined(RAYT_SPECTRAL_MODE)
            if (m_stex) {
                Spectrum s;
                for (int i = 0; i < kSpectralSamples; ++i)
                    s[i] = m_stex->reflectance(rec.uv.x, rec.uv.y, lambda[i]);
                return s;
            }
            if (m_curve) {
                Spectrum s;
                for (int i = 0; i < kSpectralSamples; ++i)
                    s[i] = std::max(m_curve->evaluate(lambda[i], color::Extrapolate::Clamp), Real(0));
                return s;
            }
#endif
            return albedoAt(rec).reflectance(lambda);
        }

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
        Lambertian(const RgbColor& a) : albedo(a) {}

        /// A measured reflectance R(lambda), used as is by the spectral build; the RGB build uses its color under D65.
        Lambertian(std::shared_ptr<const color::SpectralCurve> curve)
            : albedo(RgbColor(color::reflectanceCurveToRgb(*curve))), m_curve(std::move(curve)) {}

        /// A measured reflectance per texel R(u, v, lambda) (Core/SpectralTexture.hpp); the RGB build uses each texel's color under D65.
        explicit Lambertian(std::shared_ptr<const SpectralTexture> tex)
            : albedo(RgbColor(1.0)), m_stex(std::move(tex)) {}

        /// Makes the surface a translucent layer: stopped with probability opacity x texture alpha.
        void setTranslucent(Real opacity) { m_translucent = true; m_opacity = opacity; }

        Vector3 aovAlbedo(const SurfaceInteraction& rec) const override {
            return glm::clamp(albedoAt(rec).rgb(), Vector3(0.0), Vector3(1.0));
        }

        Real opacity(const SurfaceInteraction& rec) const override {
            if (!m_translucent) return Real(1);
            return m_opacity * (m_tex ? m_tex->alpha(rec.uv.x, rec.uv.y) : Real(1));
        }

        /// Albedo = constant color x texture (sampled at the hit point's uv).
        Lambertian(const RgbColor& a, std::shared_ptr<const Texture2D> tex) : albedo(a), m_tex(std::move(tex)) {}

        /**
         * @brief Evaluates the Lambertian BRDF.
         * * Formula: f(wo, wi) = albedo / PI
         * @return The constant spectral reflectance divided by PI.
         */
        Spectrum eval(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda,
            TransportMode mode) const override {

            // Rejects light coming from behind the surface (back-face)
            if (cosNg(rec, wi) <= 0) {  // using geometry normal
                return Spectrum(0.0);
            }

            return reflectanceAt(rec, lambda) * (1.0 / constants::PI);
        }

        /**
         * @brief Evaluates the PDF for cosine-weighted sampling.
         * * Formula: p(wi) = cos(theta) / PI
         * @return The probability density of choosing direction wi.
         */
        Real pdf(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda) const override {
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
            const SampledWavelengths& lambda,
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
            bsdfSample.f = reflectanceAt(rec, lambda) * (Real(1.0) / constants::PI);

            // 6. Set flags indicating a diffuse reflection interaction.
            bsdfSample.flags =
                BxDFFlags::Diffuse |
                BxDFFlags::Reflection;

            return bsdfSample;
        }
    };

} // namespace rayt