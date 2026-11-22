/**
 * @file Lights/AreaLight.hpp
 * @brief Diffuse area light implementation bound to a geometric shape.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>

#include "Core/Math.hpp"
#include "Core/Types.hpp"
#include "Geometry/Shape.hpp"
#include "Lights/Light.hpp"
#include "Materials/Material.hpp"
#include "Color/SpectralConversion.hpp"

namespace rayt {

    /**
     * @brief Diffuse area light emitting from a shape surface.
     *
     * This light samples points on the underlying shape and evaluates
     * emitted radiance using the attached material's emitted() function.
     */
    class AreaLight final : public Light {
    public:
        AreaLight() = default;

        AreaLight(std::shared_ptr<Shape> shape,
            std::shared_ptr<Material> material,
            const RgbColor& scale = RgbColor(1.0))
            : m_shape(std::move(shape))
            , m_material(std::move(material))
            , m_scale(scale) {}

        LightType type() const override { return LightType::Area; }

        std::optional<LightBounds> bounds() const override {
            if (!m_shape || !m_material) return std::nullopt;
            LightBounds b;
            b.bounds = m_shape->bounds();
            m_shape->normalBounds(b.w, b.cosThetaO);
            b.twoSided = m_material->emitsTwoSided();
            const Vector3& t = m_scale.rgb();
            const Real tint = Real(0.2126) * t.x + Real(0.7152) * t.y + Real(0.0722) * t.z;
            b.phi = m_material->emittedLuminance() * tint * m_shape->area() * (b.twoSided ? Real(2) : Real(1));
            b.cosThetaE = Real(0);    // diffuse emitter: spreads over the hemisphere (cos(pi/2))
            return b;
        }

        std::optional<LightLiSample>
            sampleLi(const LightSampleContext& ctx, const Point2& u,
                const SampledWavelengths& lambda,
                bool /*allowIncompletePDF*/ = false) const override {
            if (!m_shape || !m_material) return std::nullopt;

            // pdf is already with respect to solid angle at ctx.p; frontFace is set for this view.
            auto shapeSample = m_shape->sampleFrom(ctx.p, u);
            if (!shapeSample || shapeSample->pdf <= Real(0)) return std::nullopt;

            const SurfaceInteraction& pLight = shapeSample->si;

            Vector3 d = pLight.p - ctx.p;
            Real dist2 = glm::dot(d, d);
            if (dist2 <= Real(0)) return std::nullopt;

            Real invDist = math::safe_recip(math::safe_sqrt(dist2));
            Vector3 wi = d * invDist;

            Real pdfW = shapeSample->pdf;

            LightLiSample ls;
            ls.wi = wi;
            ls.pdf = pdfW;
            ls.isDelta = false;
            ls.pLight = pLight.p;
            ls.tMax = math::safe_sqrt(dist2);
            ls.Li = m_material->emitted(pLight, -wi, lambda) * m_scale.reflectance(lambda);

            return ls;
        }

        Real pdfLi(const LightSampleContext& ctx, const Vector3& wi,
            bool /*allowIncompletePDF*/ = false) const override {
            if (!m_shape) return Real(0);

            SurfaceInteraction ref;
            ref.p = ctx.p;
            ref.gn = ctx.gn;
            ref.n = ctx.gn;
            ref.medium = ctx.medium;
            ref.time = ctx.time;

            return m_shape->pdfSurface(ref, wi);
        }

        Spectrum L(const SurfaceInteraction& pLight, const Vector3& w,
            const SampledWavelengths& lambda) const override {
            if (!m_material) return Spectrum(0.0);
            return m_material->emitted(pLight, w, lambda) * m_scale.reflectance(lambda);
        }

        /**
         * @brief A linear polarizer over a flat emitter: a screen (LCD) (docs/POLARIZATION.md).
         *
         * `axis` is the transmission axis (the electric field of the light that leaves along the normal), in the plane of the emitter
         * with unit normal `normal`; `degree` the degree of polarization (0 turns it off). The polarizer is a dichroic sheet of
         * refractive index `filmIndex`: inside it the field is perpendicular to the absorbing axis (cross(normal, axis)) and to the
         * ray, and it leaves through the sheet's surface with the Fresnel transmission amplitudes, so away from the normal the field
         * is not simply the axis projected (filmIndex 1: no refraction). The radiance is left as it is, and the polarization is the
         * same at every wavelength and on both sides.
         */
        void setLinearPolarizer(const Vector3& axis, const Vector3& normal, Real degree = Real(1), Real filmIndex = Real(1.5)) {
            m_polNormal = glm::normalize(normal);
            const Vector3 t = axis - glm::dot(axis, m_polNormal) * m_polNormal;
            m_polDegree = glm::length(t) > Real(1e-9) ? std::clamp(degree, Real(0), Real(1)) : Real(0);
            if (m_polDegree > Real(0)) m_polAbsorb = glm::normalize(glm::cross(m_polNormal, t));
            m_polIndex = std::max(filmIndex, Real(1));
        }
        Real polarizerDegree() const { return m_polDegree; }

        bool polarization(const Vector3& wi, Vector3& xRef, Vector3& s) const override {
            if (!(m_polDegree > Real(0))) return false;
            const Vector3 k = -wi;                                   // the light's direction of travel, away from the emitter
            const Real kn = glm::dot(k, m_polNormal);
            const Vector3 n = kn < Real(0) ? -m_polNormal : m_polNormal;   // the side the light leaves on
            const Vector3 kt = k - glm::dot(k, n) * n;               // tangential part: kept across the surface, divided by the index
            const Real sin2 = glm::dot(kt, kt);
            const Real cosOut = std::sqrt(std::max(Real(0), Real(1) - sin2));
            const Real cosIn = std::sqrt(std::max(Real(0), Real(1) - sin2 / (m_polIndex * m_polIndex)));
            const Vector3 kIn = kt / m_polIndex + cosIn * n;         // the ray inside the sheet
            Vector3 e = glm::cross(kIn, m_polAbsorb);                // the field inside: normal to the absorbing axis and to the ray
            if (sin2 > Real(1e-12)) {
                const Vector3 sAxis = glm::normalize(glm::cross(n, kt));
                const Vector3 pIn = glm::cross(kIn, sAxis), pOut = glm::cross(k, sAxis);
                const Real ts = Real(2) * m_polIndex * cosIn / (m_polIndex * cosIn + cosOut);
                const Real tp = Real(2) * m_polIndex * cosIn / (cosIn + m_polIndex * cosOut);
                e = ts * glm::dot(e, sAxis) * sAxis + tp * glm::dot(e, pIn) * pOut;
            }
            const Real len = glm::length(e);
            if (len < Real(1e-9)) return false;
            xRef = e / len;
            s = Vector3(m_polDegree, Real(0), Real(0));
            return true;
        }

    private:
        std::shared_ptr<Shape> m_shape;
        std::shared_ptr<Material> m_material;
        RgbColor m_scale = RgbColor(1.0);   // tint multiplied onto the emitted spectrum
        Vector3 m_polNormal = Vector3(0, 1, 0), m_polAbsorb = Vector3(1, 0, 0);   // setLinearPolarizer
        Real m_polDegree = Real(0), m_polIndex = Real(1.5);
    };

} // namespace rayt