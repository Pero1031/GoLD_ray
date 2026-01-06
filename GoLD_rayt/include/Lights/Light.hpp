// Lights/Light.hpp
#pragma once

#include <optional>

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Interaction.hpp"
#include "Core/Ray.hpp"

namespace rayt {

    enum class LightType { DeltaPosition, DeltaDirection, Area, Infinite };

    // 将来: 波長サンプルを持つ（PBRTの SampledWavelengths 相当）
    struct SampledWavelengths {
        // 例: N本のlambda[nm] と、そのPDFなど
        // 今は空でもOK（RGBモードでは無視）
    };

    // 「光サンプリングに必要な最小の参照情報」
    struct LightSampleContext {
        Point3 p;
        Normal3 gn;          // geometry normal (推奨)
        const Medium* medium = nullptr;
        Real time = 0;
    };

    struct LightLiSample {
        Vector3 wi;        // ref -> light direction (normalized)
        Spectrum Li;       // radiance arriving from wi (already includes distance/geometry as your convention)
        Real pdf = 0;      // pdf w.r.t. solid angle at ref (IMPORTANT!)
        Point3 pLight;     // sampled point on light (area light)
        Real tMax = 0;     // distance to pLight (for shadow ray)
        bool isDelta = false;  
    };

    // Lightからレイを出す（SampleLe用）
    struct LightLeSample {
        Ray ray;        // emitted ray
        Spectrum Le;    // radiance along ray
        Real pdfPos = 0;
        Real pdfDir = 0;
        bool isDelta = false; // e.g. point light: delta position
    };

    class Light {
    public:
        virtual ~Light() = default;

        virtual LightType type() const = 0;

        // ref（交差点）から見た光をサンプル
        virtual std::optional<LightLiSample>
            sampleLi(const LightSampleContext& ctx, const Point2& u,
                const SampledWavelengths& lambda,
                bool allowIncompletePDF = false) const = 0;

        // その方向をサンプルする確率密度
        virtual Real pdfLi(const LightSampleContext& ctx, const Vector3& wi,
            bool allowIncompletePDF = false) const = 0;

        // InfiniteLights only: miss時の背景放射
        virtual Spectrum Le(const Ray& ray, const SampledWavelengths& lambda) const {
            return Spectrum(0.0);
        }

        // 双方向/ライトトレーシング用（将来）
        virtual std::optional<LightLeSample>
            sampleLe(const Point2& u1, const Point2& u2,
                SampledWavelengths& lambda, Real time) const {
            return std::nullopt;
        }

        virtual void pdfLe(const Ray& ray, Real* pdfPos, Real* pdfDir) const {
            if (pdfPos) *pdfPos = 0;
            if (pdfDir) *pdfDir = 0;
        }

        // AreaLights only: light surfaceからの放射（Material::emittedと同義）
        virtual Spectrum L(const SurfaceInteraction& pLight,
            const Vector3& w,
            const SampledWavelengths& lambda) const {
            return Spectrum(0.0);
        }
    };

} // namespace rayt