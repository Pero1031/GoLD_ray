/*#pragma once

#include "Material.hpp"

namespace rayt {

    class DiffuseLight : public Material {
    public:
        // color: 光の色と強さ (例: (10, 10, 10) なら非常に明るい白)
        DiffuseLight(const Spectrum& color) : m_emit(color) {}

        // 散乱しない (光を反射するのではなく、光源そのものなので false)
        virtual bool scatter(const Ray& r_in, const SurfaceInteraction& rec,
            Spectrum& attenuation, Ray& scattered) const override {
            return false;
        }

        // 発光する！
        virtual Spectrum emitted(const SurfaceInteraction& rec) const override {
            return m_emit;
        }

    private:
        Spectrum m_emit;
    };

}*/

#pragma once

#include "Materials/Material.hpp"

namespace rayt {

    class DiffuseLight : public Material {
    public:
        // color: 光の色と強さ (例: (10, 10, 10) なら非常に明るい白)
        DiffuseLight(const Spectrum& color) : m_emit(color) {}

        // --------------------------------------------------------
        // BSDF (反射・透過) の振る舞い
        // 光源そのものは入射光を反射しない（吸収する）ものとして扱います
        // --------------------------------------------------------

        // 評価値: 常に黒 (反射しない)
        Spectrum eval(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            TransportMode mode) const override {
            return Spectrum(0.0);
        }

        // サンプリング: 反射方向がないので何もしない
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec,
            const Vector3& wo,
            const Point2& u,
            TransportMode mode) const override {
            return std::nullopt;
        }

        // PDF: 0
        Real pdf(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi) const override {
            return 0.0;
        }

        // --------------------------------------------------------
        // 発光 (Emission)
        // --------------------------------------------------------

        // wo: 表面からカメラ(または次のレイ)へ向かうベクトル
        Spectrum emitted(const SurfaceInteraction& rec, const Vector3& wo) const override {
            // 法線と同じ側から見ている場合のみ発光する（裏面は光らない）
            if (glm::dot(rec.n, wo) > 0.0) {
                return m_emit;
            }
            return Spectrum(0.0);
        }

    private:
        Spectrum m_emit;
    };

}