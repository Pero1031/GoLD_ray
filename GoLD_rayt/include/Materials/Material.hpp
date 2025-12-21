#pragma once

#include "Core/Core.hpp"
#include "Core/Interaction.hpp"
#include "Core/Ray.hpp"

#include <optional>

// Include your spectral data handler
// #include "IORInterpolator.hpp" 

namespace rayt {

    // BxDF (BSDF) の性質を表すフラグ
    // 多重重点的サンプリング (MIS) や、鏡面反射の特異点扱いに必要です。
    enum BxDFType {
        BSDF_REFLECTION = 1 << 0,
        BSDF_TRANSMISSION = 1 << 1,
        BSDF_DIFFUSE = 1 << 2,
        BSDF_GLOSSY = 1 << 3,
        BSDF_SPECULAR = 1 << 4,
        BSDF_ALL = 0xFF
    };

    // 放射輝度(Radiance)か、重要度(Importance)か。
    // 屈折を伴う双方向パス（BDPTなど）で非対称なスループットを扱う際に必要。
    enum class TransportMode {
        Radiance,
        Importance
    };

    // サンプリングの結果を格納する構造体
    struct BSDFSample {
        Spectrum f;           // BSDFの値 (throughput)
        Vector3 wi;           // サンプリングされた入射方向 (World Space)
        Real pdf = 0;        // 確率密度
        BxDFType sampledType; // サンプリングされたローブの種類

        bool isSpecular() const { return sampledType & BSDF_SPECULAR; }
    };

    class Material {
    public:
        virtual ~Material() = default;

        /**
         * @brief BSDF (BRDF/BTDF) の値を評価する関数 f(wo, wi)
         * * @param rec 交差表面の情報（法線、接線、UVなど）
         * @param wo 視線方向 (Outdirection direction, World Space)
         * @param wi 光源方向 (Incident direction, World Space)
         * @param mode 輸送モード (通常は Radiance)
         * @return Spectrum BSDFの値 [1/sr]
         * * 注: 厳密なPBRでは、ここで cosine term (n・wi) を含めるかどうかは設計によりますが、
         * 通常は純粋なBSDF値を返し、Integrator側でcosineを掛けます。
         * ただし、デルタ分布（完全鏡面）の場合は0を返します。
         */
        virtual Spectrum eval(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            TransportMode mode = TransportMode::Radiance) const = 0;

        /**
         * @brief BSDFに基づいて新しい入射方向 wi をインポータンスサンプリングする
         * * @param rec 交差表面の情報
         * @param wo 視線方向 (World Space)
         * @param sample RNGから得られたランダムな2次元サンプル値 [0,1)^2
         * @param mode 輸送モード
         * @return std::optional<BSDFSample> サンプリング成功時は結果、全反射などで無効な場合はnullopt
         */
        virtual std::optional<BSDFSample> sample(const SurfaceInteraction& rec,
            const Vector3& wo,
            const Point2& u, // ランダムシード
            TransportMode mode = TransportMode::Radiance) const = 0;

        /**
         * @brief 確率密度関数 (PDF) を評価する
         * MIS (Multiple Importance Sampling) のために必須です。
         * * @param rec 交差表面の情報
         * @param wo 視線方向
         * @param wi 光源方向 (評価したい方向)
         * @return Float pdfの値
         */
        virtual Real pdf(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi) const = 0;

        // -----------------------------------------------------------
        // 以下のヘルパーや、以前の簡易的な関数は必要に応じて残すか削除します
        // -----------------------------------------------------------

        // 自己発光（エリアライト用）
        virtual Spectrum emitted(const SurfaceInteraction& rec,
            const Vector3& wo) const {
            return Spectrum(0.0);
        }

        // ラフネスがあるかどうかなどの判定用
        virtual bool isSpecular() const { return false; }
    };

}