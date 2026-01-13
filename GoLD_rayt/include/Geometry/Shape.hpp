/**
 * @file Geometry/Shape.hpp
 * @brief Abstract shape interface for geometric sampling.
 */

#pragma once

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
    };

} // namespace rayt