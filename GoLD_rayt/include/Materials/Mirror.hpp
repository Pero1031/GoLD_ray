#pragma once

#include "Core/Core.hpp"
#include "Materials/Material.hpp"

namespace rayt {

    class Mirror : public Material {
    public:
        Spectrum albedo; // 通常は(1,1,1)ですが、色付き鏡のために

        Mirror(const Spectrum& a = Spectrum(1.0)) : albedo(a) {}

        // 鏡面反射は「特定の方向」以外確率0なので、任意の方向に対する評価は常に黒
        Spectrum eval(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            TransportMode mode) const override {
            return Spectrum(0.0);
        }

        // デルタ関数の確率は点評価できないため 0
        Real pdf(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi) const override {
            return 0.0;
        }

        // サンプリングのみ機能する
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec,
            const Vector3& wo,
            const Point2& u,
            TransportMode mode) const override {
            BSDFSample result;

            // 正反射ベクトル R = I - 2(N・I)N
            // ※ woはカメラに向かうベクトルなので、reflectの引数注意
            // 一般的な reflect(v, n) は v - 2*dot(v,n)*n
            // wo (視線) を反転(-wo)して入射として計算するか、公式通りやるか
            // rayt::Reflect の定義によりますが、通常は:
            result.wi = glm::reflect(-wo, rec.n);

            // 鏡面フラグを立てる
            result.sampledType = BxDFType(BSDF_SPECULAR | BSDF_REFLECTION);
            result.pdf = 1.0; // 特異点なのでダミーの1を入れる

            // 積分値（f/pdf に相当する重み）を直接返す
            // 完全鏡面の場合、レンダリング方程式の余弦項が相殺される形で
            // 単純に「反射率(albedo)」だけが残るのが一般的です。
            // ※Integrator側で beta *= f しているので、ここは反射率そのものを返します。

            // ★本来はここにフレネル項 (Fresnel) が掛かります
            result.f = albedo;

            return result;
        }
    };

}