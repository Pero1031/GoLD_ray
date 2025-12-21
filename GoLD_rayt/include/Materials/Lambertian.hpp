#pragma once

#include "Core/Core.hpp"
#include "Materials/Material.hpp"
#include "Core/Utils.hpp"

namespace rayt {

    class Lambertian : public Material {
        Spectrum albedo;
    public:
        // albedo: (R, G, B)
        Lambertian(const Spectrum& a) : albedo(a) {}

        // 1. BSDFの評価 (f = albedo / PI)
        Spectrum eval(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            TransportMode mode) const override {
            // 裏面からの光や、埋もれる光は0にする
            Real cosTheta = glm::dot(rec.n, wi);
            if (cosTheta <= 0) return Spectrum(0.0);

            return albedo * (1.0f / Constants::PI);
        }

        // 2. 確率密度の評価 (pdf = cos(theta) / PI)
        Real pdf(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi) const override {
            Real cosTheta = glm::dot(rec.n, wi);
            if (cosTheta <= 0) return 0.0;

            return cosTheta * (1.0f / Constants::PI);
        }

        // 3. サンプリング
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec,
            const Vector3& wo,
            const Point2& u,
            TransportMode mode) const override {

            BSDFSample bsdfSample;

            // 1. ローカル座標で方向を生成し、ワールド座標へ変換
            //    法線 rec.n を基準にした接空間を作る
            Onb onb(rec.n);
            Vector3 localDir = Utils::randomCosineDirection(u); // ランダムな方向
            bsdfSample.wi = onb.localToWorld(localDir);

            // 2. 幾何学的整合性のチェック（裏面に行っていないか）
            if (glm::dot(rec.n, bsdfSample.wi) <= 0) return std::nullopt;

            // 3. PDF計算: cos(theta) / PI
            //    Integrator側で (f * cos / pdf) を計算するので、正しいPDFを返す
            bsdfSample.pdf = localDir.z / Constants::PI; // localDir.z は cosTheta と同じ

            // 4. BSDF値: albedo / PI
            bsdfSample.f = albedo / Constants::PI;

            // 5. フラグ設定
            bsdfSample.sampledType = BxDFType(BSDF_DIFFUSE | BSDF_REFLECTION);

            return bsdfSample;
        }

        /*virtual bool scatter(const Ray& r_in, const SurfaceInteraction& rec,
            Spectrum& attenuation, Ray& scattered) const override {
            // Lambertian reflection direction: normal + random vector
            Vector3 scatter_direction = rec.n + Utils::randomUnitVector();

            // ランダムベクトルがたまたま法線と真逆でゼロベクトルになるのを防ぐ
            if (glm::length(scatter_direction) < 1e-8)
                scatter_direction = rec.n;

            scattered = Ray(rec.p, glm::normalize(scatter_direction));
            attenuation = m_albedo;
            return true;
        }*/

    /*private:
        Spectrum m_albedo;*/
    };

}