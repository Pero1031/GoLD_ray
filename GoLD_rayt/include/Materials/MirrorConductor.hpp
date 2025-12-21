#pragma once

#include "Core/Core.hpp"
#include "Materials/Material.hpp"
#include "Core/Interaction.hpp"
// #include "Spectrum.hpp" // Spectrumの定義が必要
#include <complex>
#include <algorithm> // std::clampに必要
#include <cmath>

namespace rayt {

    class MirrorConductor : public Material {
        Spectrum eta; // 屈折率の実部 n (RGB)
        Spectrum k;   // 消衰係数 k (複素部) (RGB)

    public:
        // コンストラクタ: 金(Au)や銅(Cu)などのn, kを渡す
        MirrorConductor(const Spectrum& eta, const Spectrum& k) : eta(eta), k(k) {}

        // 鏡面反射なので、散乱確率密度(PDF)と評価値(Eval)は 0 を返す
        Spectrum eval(const SurfaceInteraction&, const Vector3&, const Vector3&, TransportMode) const override {
            return Spectrum(0.0);
        }

        Real pdf(const SurfaceInteraction&, const Vector3&, const Vector3&) const override {
            return 0.0;
        }

        // サンプリング関数
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec,
            const Vector3& wo,
            const Point2& u,
            TransportMode mode) const override {

            BSDFSample bsdfSample;

            // 1. 方向の計算：完全正反射 (Perfect Specular Reflection)
            // wo はカメラ側へのベクトル。入射ベクトルは -wo。
            // 以下の計算は一般的に reflect(inc, n) = inc - 2(n.inc)n
            bsdfSample.wi = glm::reflect(-wo, rec.n);

            // 幾何学的チェック（表面より下に行っていないか）
            Real cosTheta = glm::dot(bsdfSample.wi, rec.n);
            if (cosTheta <= 0) return std::nullopt;

            // 2. フラグ設定
            // Integrator側で「鏡面反射だ」と認識させるために必須
            bsdfSample.sampledType = BxDFType(BSDF_SPECULAR | BSDF_REFLECTION);
            bsdfSample.pdf = 1.0; // 特異点（Delta distribution）のためダミー値

            // 3. BSDF値 (f) の計算
            // Integratorのロジック (beta *= f) に合わせ、フレネル反射率そのものを返す
            // RGB各波長について厳密解を計算
            bsdfSample.f = Spectrum(
                fresnelConductorExact(cosTheta, eta.x, k.x),
                fresnelConductorExact(cosTheta, eta.y, k.y),
                fresnelConductorExact(cosTheta, eta.z, k.z)
            );

            return bsdfSample;
        }

    private:
        // 導体用フレネル計算の厳密解 (Exact solution for Conductors)
        // cosThetaI: 入射角のコサイン (n・wi)
        Real fresnelConductorExact(Real cosThetaI, Real etaVal, Real kVal) const {
            // 数値誤差対策: cosは 0~1 にクランプ
            cosThetaI = std::clamp(cosThetaI, Real(0.0), Real(1.0));

            // 複素屈折率 N = eta + ik
            std::complex<Real> N(etaVal, kVal);

            Real cosThetaI2 = cosThetaI * cosThetaI;
            Real sinThetaI2 = Real(1.0) - cosThetaI2;

            // スネルの法則の複素数版を用いて透過角のsin^2を計算
            std::complex<Real> sinThetaT2 = sinThetaI2 / (N * N);
            std::complex<Real> cosThetaT = std::sqrt(Real(1.0) - sinThetaT2);

            // フレネル係数 (s偏光, p偏光)
            std::complex<Real> r_s = (cosThetaI - N * cosThetaT) / (cosThetaI + N * cosThetaT);
            std::complex<Real> r_p = (N * cosThetaI - cosThetaT) / (N * cosThetaI + cosThetaT);

            // ※注: 一般的な教科書では r_s の分子が (ni cos - nt cos) だったり (cos - N cos) だったり定義揺れがありますが、
            // 最終的に絶対値の二乗を取るため、符号の違いは結果に影響しません。
            // ここでは PBRT v3/v4 の実装に近い形式を採用しています。

            // 非偏光の反射率 = (|r_s|^2 + |r_p|^2) / 2
            return (std::norm(r_s) + std::norm(r_p)) * Real(0.5);
        }
    };

} // namespace rayt