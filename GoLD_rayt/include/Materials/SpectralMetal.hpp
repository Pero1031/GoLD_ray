#pragma once

#include "Materials/Material.hpp"
#include "IO/IORInterpolator.hpp"
//#include "Core/Utils.hpp"
#include "Core/Sampling.hpp"

namespace rayt {

    class SpectralMetal : public Material {
    public:
        // コンストラクタでCSVファイル名を指定
        // roughness: 表面の粗さ (0.0=鏡面, 1.0=ザラザラ)
        SpectralMetal(const std::string& csvPath, Real roughness)
            : m_roughness(std::min(roughness, (Real)1.0))
        {
            // CSVをロード
            IORInterpolator ior;
            if (!ior.loadCSV(csvPath)) {
                std::cerr << "[Error] Failed to load IOR data: " << csvPath << std::endl;
                // エラー時はデフォルトで銀のような値をセット
                m_eta = Spectrum(0.05);
                m_k = Spectrum(3.0);
            }
            else {
                // R, G, B それぞれの代表波長における n, k を取得
                // Red: 650nm, Green: 550nm, Blue: 450nm
                auto complexR = ior.evaluate(650.0);
                auto complexG = ior.evaluate(550.0);
                auto complexB = ior.evaluate(450.0);

                m_eta = Spectrum(complexR.real(), complexG.real(), complexB.real());
                m_k = Spectrum(complexR.imag(), complexG.imag(), complexB.imag());

                std::cout << "[SpectralMetal] Loaded " << csvPath << "\n"
                    << "  R(650nm): n=" << m_eta.r << ", k=" << m_k.r << "\n"
                    << "  G(550nm): n=" << m_eta.g << ", k=" << m_k.g << "\n"
                    << "  B(450nm): n=" << m_eta.b << ", k=" << m_k.b << std::endl;
            }
        }

        // ---------------------------------------------------------------------
        // 1. eval: BSDFの評価
        // ---------------------------------------------------------------------
        // ラフネスがある場合、本当はマイクロファセット分布を評価すべきですが、
        // 簡易実装として「ラフネスがある導体も特異点に近い」とみなして 0 を返すか、
        // あるいは完全鏡面として振る舞わせます。
        // ※今回は「sample」ですべて処理するタイプ（Delta分布扱い）として実装します。
        Spectrum eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, TransportMode mode) const override {
            return Spectrum(0.0);
        }

        // ---------------------------------------------------------------------
        // 2. pdf: 確率密度関数
        // ---------------------------------------------------------------------
        // Delta分布扱いなので 0
        Real pdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi) const override {
            return 0.0;
        }

        // ---------------------------------------------------------------------
        // 3. sample: 次の方向を決定し、重みを計算
        // ---------------------------------------------------------------------
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec,
            const Vector3& wo,
            const Point2& u,
            TransportMode mode) const override {

            BSDFSample bsdfSample;

            // (A) 反射方向の計算
            // 完全鏡面反射ベクトル
            Vector3 reflected = glm::reflect(-wo, rec.n);

            // (B) ラフネスの適用 (旧scatterのロジックを移植)
            // 表面を乱数で揺らす
            if (m_roughness > 0) {
                Vector3 fuzz = rayt::sampling::randomInUnitSphere() * m_roughness;
                reflected = glm::normalize(reflected + fuzz);
            }

            // (C) 幾何チェック (表面の内側に入ってしまったら吸収)
            Real cosTheta = glm::dot(reflected, rec.n);
            if (cosTheta <= 0) return std::nullopt;

            // (D) 結果の格納
            bsdfSample.wi = reflected;

            // フラグ設定:
            // ラフネスがあっても、この簡易実装では「確率的に1方向を選ぶ」ので
            // 数学的には SPECULAR (デルタ分布) として扱ったほうがIntegratorとの相性が良いです。
            // (GLOSSYにすると pdf の計算が必要になるため)
            bsdfSample.sampledType = BxDFType(BSDF_SPECULAR | BSDF_REFLECTION);
            bsdfSample.pdf = 1.0;

            // (E) フレネル反射率の計算 (複素数)
            // 視線ベクトル wo と法線 n の角度を使います (厳密にはハーフベクトルですが、ここでは wo・n で近似)
            Real cosThetaI = std::clamp(glm::dot(wo, rec.n), (Real)0.0, (Real)1.0);

            Spectrum F;
            F.r = fresnelConductorExact(cosThetaI, m_eta.r, m_k.r);
            F.g = fresnelConductorExact(cosThetaI, m_eta.g, m_k.g);
            F.b = fresnelConductorExact(cosThetaI, m_eta.b, m_k.b);

            bsdfSample.f = F;

            return bsdfSample;
        }

    private:
        // 導体(金属)のためのフレネル計算
        // 入射角と複素屈折率 (n + ik) から反射率(Energy)を求める
        // R = |(n - ik - 1)/(n - ik + 1)|^2  (垂直入射の簡易式ではなく、角度依存のフルバージョン)
        Spectrum conductorFresnel(Real cosTheta, const Spectrum& eta, const Spectrum& k) const {

            // 厳密な計算は複雑なため、ここでは以下の近似式(PBRT等で使用)または
            // 各成分ごとの垂直入射近似 + Schlick補正を用いるのが一般的だが、
            // 研究目的なら「角度依存性」も重要。

            // ここでは、実装の容易さと精度のバランスが良い
            // "Fresnel Equations for Conductor" の各RGB成分計算を行う

            Spectrum result;
            result.r = fresnelConductorExact(cosTheta, eta.r, k.r);
            result.g = fresnelConductorExact(cosTheta, eta.g, k.g);
            result.b = fresnelConductorExact(cosTheta, eta.b, k.b);

            return result;
        }

        // 単一波長に対する導体フレネル項の計算
        Real fresnelConductorExact(Real cosTheta, Real eta, Real k) const {
            Real cosTheta2 = cosTheta * cosTheta;
            Real sinTheta2 = 1.0 - cosTheta2;

            Real t0 = eta * eta - k * k - sinTheta2;
            Real a2plusb2 = std::sqrt(t0 * t0 + 4 * eta * eta * k * k);
            Real t1 = a2plusb2 + cosTheta2;
            Real a = std::sqrt(0.5 * (a2plusb2 + t0));
            Real t2 = 2 * cosTheta * a;
            Real Rs = (t1 - t2) / (t1 + t2);

            Real t3 = cosTheta2 * a2plusb2 + sinTheta2 * sinTheta2;
            Real t4 = t2 * sinTheta2;
            Real Rp = Rs * (t3 - t4) / (t3 + t4);

            return 0.5 * (Rp + Rs);
        }

    private:
        Spectrum m_eta; // 屈折率 n
        Spectrum m_k;   // 消衰係数 k
        Real m_roughness;
    };

}