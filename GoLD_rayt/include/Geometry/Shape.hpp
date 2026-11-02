/**
 * @file Geometry/Shape.hpp
 * @brief Abstract shape interface for geometric sampling.
 */

#pragma once

#include <cmath>
#include <optional>

#include "Core/Types.hpp"
#include "Core/Interaction.hpp"
#include "Scene/Hittable.hpp"

namespace rayt {

    /**
     * @brief Result of sampling a point on a shape surface.
     *
     * The PDF is expressed with respect to surface area.
     */
    struct ShapeSample {
        SurfaceInteraction si;
        Real pdf = 0;
    };

    /**
     * @brief Base class for geometric shapes that support surface sampling.
     */
    class Shape : public Hittable {
    public:
        ~Shape() override = default;

        /**
         * @brief Total surface area of the shape.
         */
        virtual Real area() const = 0;

        /**
         * @brief Sample a point uniformly on the shape's surface.
         *
         * @param u 2D uniform random sample in [0, 1)^2
         * @return Sampled surface interaction and PDF (area measure)
         */
        virtual std::optional<ShapeSample> sampleSurface(const Point2& u) const = 0;

        /**
         * @brief PDF of sampling a point on the shape visible from a reference point.
         *
         * The PDF is expressed with respect to solid angle at the reference point.
         */
        virtual Real pdfSurface(const SurfaceInteraction& ref, const Vector3& wi) const = 0;

        /**
         * @brief Sample a point on the shape as seen from a reference point.
         *
         * Unlike sampleSurface(), the returned pdf is with respect to SOLID ANGLE at `ref`, and
         * si.frontFace tells whether `ref` is on the side the surface normal points to (a
         * one-sided emitter only emits there). Must be consistent with pdfSurface().
         *
         * Default: uniform area sampling converted to solid angle. Shapes override it with
         * something smarter (Sphere: cone sampling).
         */
        virtual std::optional<ShapeSample> sampleFrom(const Point3& ref, const Point2& u) const {
            auto s = sampleSurface(u);
            if (!s || s->pdf <= Real(0)) return std::nullopt;

            const Vector3 d = s->si.p - ref;
            const Real dist2 = glm::dot(d, d);
            if (dist2 <= Real(0)) return std::nullopt;

            const Vector3 wi = d / std::sqrt(dist2);
            const Real cosL = glm::dot(s->si.gn, -wi);
            if (std::abs(cosL) <= Real(1e-9)) return std::nullopt;

            s->si.frontFace = cosL > Real(0);
            s->pdf = s->pdf * dist2 / std::abs(cosL);
            return s;
        }
    };

} // namespace rayt