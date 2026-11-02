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

    private:
        std::shared_ptr<Shape> m_shape;
        std::shared_ptr<Material> m_material;
        RgbColor m_scale = RgbColor(1.0);   // tint multiplied onto the emitted spectrum
    };

} // namespace rayt