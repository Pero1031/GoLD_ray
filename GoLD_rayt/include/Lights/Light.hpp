// Lights/Light.hpp
#pragma once

#include <optional>

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Interaction.hpp"

namespace rayt {

    struct LightSample {
        Vector3 wi;     // ref点から光への方向
        Spectrum Li;    // 入射放射輝度
        Real pdf = 0;
        Point3 pLight;  // サンプルした光源上の点（面光源用）
        bool isDelta = false;
    };

    class Light {
    public:
        virtual ~Light() = default;

        // ref（交差点）から見た光をサンプル
        virtual std::optional<LightSample>
            sampleLi(const SurfaceInteraction& ref, const Point2& u) const = 0;

        // その方向をサンプルする確率密度
        virtual Real pdfLi(const SurfaceInteraction& ref, const Vector3& wi) const = 0;

        // 環境光の“背景”（miss時）
        virtual Spectrum Le(const Ray& ray) const { return Spectrum(0.0); }
    };

} // namespace rayt