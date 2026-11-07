/**
* @ file Materials/DiffuseLight.hpp
*/

#pragma once

#include <algorithm>

#include "Materials/Material.hpp"
#include "Color/SpectralConversion.hpp"
#include "Color/EmissionSpectrum.hpp"

namespace rayt {

    class DiffuseLight : public Material {
    public:
        // color: color and strength of the light (e.g. (10, 10, 10) is a very bright white)
        DiffuseLight(const EmissionSpectrum& color, Real intensity = Real(1), bool twoSided = false)
            : m_emit(color)
            , m_intensity(intensity)
            , m_twoSided(twoSided) {}

        // --------------------------------------------------------
        // Behavior of the BSDF (reflection / transmission)
        // The light source itself does not reflect incident light (it absorbs it).
        // --------------------------------------------------------

        // Value: always black (no reflection)
        Spectrum eval(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda,
            TransportMode mode) const override {
            return Spectrum(0.0);
        }

        // Sampling: there is no reflection direction, so there is nothing to do
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec,
            const Vector3& wo,
            const Point2& u,
            const SampledWavelengths& lambda,
            TransportMode mode) const override {
            return std::nullopt;
        }

        // PDF: 0
        Real pdf(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda) const override {
            return 0.0;
        }

        // --------------------------------------------------------
        // Emission
        // --------------------------------------------------------

        // wo: vector from the surface towards the camera (or the next ray)
        Spectrum emitted(const SurfaceInteraction& rec, const Vector3& wo,
            const SampledWavelengths& lambda) const override {
            // Emits only when seen from the normal's side (the back side is dark)
            if (m_twoSided || rec.frontFace) {
                return m_emit(lambda) * m_intensity;
            }
            return Spectrum(0.0);
        }

        Real emittedLuminance() const override {
            const Vector3& c = m_emit.rgb();
            return (Real(0.2126) * c.x + Real(0.7152) * c.y + Real(0.0722) * c.z) * m_intensity;
        }

        bool emitsTwoSided() const override { return m_twoSided; }

        /// Hue of the emission scaled so that its largest channel is 1 (a light has no reflectance).
        Vector3 aovAlbedo(const SurfaceInteraction& /*rec*/) const override {
            const Vector3& c = m_emit.rgb();
            const Real m = std::max({ c.x, c.y, c.z });
            return m > Real(0) ? c / m : Vector3(1.0);
        }

    private:
        EmissionSpectrum m_emit;
        Real m_intensity = Real(1);
        bool m_twoSided = false;
    };

}