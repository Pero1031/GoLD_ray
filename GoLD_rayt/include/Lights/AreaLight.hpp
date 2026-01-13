/**
 * @file Lights/AreaLight.hpp
 * @brief Diffuse area light implementation bound to a geometric shape.
 */

#pragma once

#include <memory>
#include <optional>

#include "Core/Math.hpp"
#include "Core/Types.hpp"
#include "Geometry/Shape.hpp"
#include "Lights/Light.hpp"
#include "Materials/Material.hpp"

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
            const Spectrum& scale = Spectrum(1.0))
            : m_shape(std::move(shape))
            , m_material(std::move(material))
            , m_scale(scale) {}

        LightType type() const override { return LightType::Area; }

        std::optional<LightLiSample>
            sampleLi(const LightSampleContext& ctx, const Point2& u,
                const SampledWavelengths& /*lambda*/,
                bool /*allowIncompletePDF*/ = false) const override {
            if (!m_shape || !m_material) return std::nullopt;

            auto shapeSample = m_shape->sampleSurface(u);
            if (!shapeSample || shapeSample->pdf <= Real(0)) return std::nullopt;

            const SurfaceInteraction& pLight = shapeSample->si;

            Vector3 d = pLight.p - ctx.p;
            Real dist2 = glm::dot(d, d);
            if (dist2 <= Real(0)) return std::nullopt;

            Real invDist = math::safe_recip(math::safe_sqrt(dist2));
            Vector3 wi = d * invDist;

            Real cosTheta = std::abs(glm::dot(pLight.gn, -wi));
            if (cosTheta <= Real(0)) return std::nullopt;

            Real pdfW = shapeSample->pdf * dist2 / cosTheta;
            if (pdfW <= Real(0)) return std::nullopt;

            LightLiSample ls;
            ls.wi = wi;
            ls.pdf = pdfW;
            ls.isDelta = false;
            ls.pLight = pLight.p;
            ls.tMax = math::safe_sqrt(dist2);
            ls.Li = m_material->emitted(pLight, -wi) * m_scale;

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
            const SampledWavelengths& /*lambda*/) const override {
            if (!m_material) return Spectrum(0.0);
            return m_material->emitted(pLight, w) * m_scale;
        }

    private:
        std::shared_ptr<Shape> m_shape;
        std::shared_ptr<Material> m_material;
        Spectrum m_scale = Spectrum(1.0);
    };

} // namespace rayt