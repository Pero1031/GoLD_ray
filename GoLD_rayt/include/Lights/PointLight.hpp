/**
 * @file Lights/PointLight.hpp
 * @brief Point light implementation (delta-position light).
 *
 * This file defines rayt::PointLight, a delta-position light source that emits
 * from a single point in space. It supports Next Event Estimation (NEE) via
 * Light::sampleLi() by returning the unique direction toward the light and an
 * inverse-square falloff contribution. For delta lights, Light::pdfLi() returns
 * 0 (delta distribution w.r.t. solid angle), while sampleLi() provides a usable
 * sample with pdf = 1 to avoid division-by-zero in estimators and to keep delta
 * handling explicit for MIS.
 *
 * The intensity is currently treated as RGB-scaled strength; it can be extended
 * to wavelength-dependent evaluation using SampledWavelengths in spectral mode.
 */

#pragma once

#include <optional>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Core/Math.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Color/SpectralConversion.hpp"
#include "Color/EmissionSpectrum.hpp"
#include "Lights/Light.hpp"

namespace rayt {

    /**
    * @brief Point light (delta position light).
    *
    * Emits light from a single position in space.
    * - sampleLi: returns a deterministic direction to the light
    * - pdfLi: 0 (delta distribution)
    * - Le(ray): not used (unless you want glowing point in background; typically 0)
    *
    * Spectrum is currently treated as RGB intensity; later it can be wavelength-evaluated.
    */
    class PointLight final : public Light {
    public:
        PointLight() = default;

        PointLight(const Point3& position,
            const EmissionSpectrum& intensity,
            Real scale = 1.0)
            : m_p(position), m_I(intensity), m_scale(scale) {}

        LightType type() const override { return LightType::DeltaPosition; }

        std::optional<LightBounds> bounds() const override {
            LightBounds b;
            b.bounds = AABB(m_p, m_p);
            b.w = Vector3(0, 0, 1);
            const Vector3& c = m_I.rgb();
            b.phi = Real(4) * constants::PI * m_scale * (Real(0.2126) * c.x + Real(0.7152) * c.y + Real(0.0722) * c.z);
            b.cosThetaO = Real(-1);   // emits in all directions
            b.cosThetaE = Real(-1);
            b.twoSided = false;
            return b;
        }

        std::optional<LightLiSample>
        sampleLi(const LightSampleContext& ctx,
                 const Point2& /*u*/,
                 const SampledWavelengths& lambda,
                 bool /*allowIncompletePDF*/ = false) const override
        {

            Vector3 d = m_p - ctx.p;
            Real dist2 = glm::dot(d, d);
            if (dist2 <= Real(0)) return std::nullopt;

            Real invDist = math:: safe_recip(math::safe_sqrt(dist2));
            Vector3 wi = d * invDist;

            LightLiSample ls;
            ls.wi = wi;
            ls.pLight = m_p;
            ls.tMax = math::safe_sqrt(dist2);
            ls.isDelta = true;

            // Delta light: pdf is 1 for the single sampled direction conceptually,
            // but for MIS we typically treat it as delta => pdfLi() returns 0,
            // and sampleLi provides a usable sample with pdf=1.
            // This convention avoids dividing by zero while keeping delta handling explicit.
            ls.pdf = Real(1.0);

            // Radiance arriving from a point light:
            // Often modeled as intensity / r^2 (inverse-square falloff).
            // Here we put the falloff into Li directly (common in simple renderers).
            Real invDist2 = math::safe_recip(dist2);
            ls.Li = (m_I(lambda) * m_scale) * invDist2;

            return ls;
        }

        Real pdfLi(const LightSampleContext& /*ctx*/,
            const Vector3& /*wi*/,
            bool /*allowIncompletePDF*/ = false) const override
        {
            // Delta distribution over directions: pdf w.r.t solid angle is 0 everywhere
            // (except at the single direction, where it is a delta spike).
            // Return 0 as PBRT does for delta lights.
            return 0.0;
        }

    private:
        Point3 m_p = Point3(0.0);
        EmissionSpectrum m_I; // intensity (arbitrary scale; evaluated per wavelength when rendering)
        Real m_scale = Real(1.0);
    };
} // namespace rayt