#pragma once

#include "Core/Core.hpp"
#include "renderer/Scene.hpp"   // HittableList, etc.
#include "Renderer/Camera.hpp"
#include "renderer/Film.hpp"
#include "Core/Ray.hpp"
#include "Core/Interaction.hpp"
#include "Materials/Material.hpp"
//#include "Core/Utils.hpp"
#include "Core/Sampling.hpp"

#include <memory>
#include <iostream>

namespace rayt {

    // Abstract base class for all rendering algorithms.
    class Integrator {
    public:
        virtual ~Integrator() = default;

        // 純粋仮想関数
        // ここで Scene や Film の型を使うので、上のincludeが必須です
        virtual void render(const Scene& scene, Film& film) = 0;
    };

    // -------------------------------------------------------------------------

    // Path Tracing Integrator
    class PathIntegrator : public Integrator {
    public:
        // コンストラクタ
        PathIntegrator(std::shared_ptr<Camera> camera, int maxDepth, int spp)
            : m_camera(camera), m_maxDepth(maxDepth), m_spp(spp) {}

        // レンダリングループの実装
        virtual void render(const Scene& scene, Film& film) override {
            int width = film.width();
            int height = film.height();

            std::cout << "[PathIntegrator] Rendering " << width << "x" << height
                << " (" << m_spp << " spp)" << std::endl;

            for (int j = 0; j < height; ++j) {
                // 進捗表示
                std::cout << "\rScanlines remaining: " << (height - j) << " " << std::flush;

                #pragma omp parallel for // 可能なら並列化推奨

                for (int i = 0; i < width; ++i) {
                    Spectrum pixelColor(0.0);

                    for (int s = 0; s < m_spp; ++s) {
                        // アンチエイリアシング用のジッター
                        Real u = (Real(i) + rayt::sampling::Random()) / Real((width));
                        Real v = (Real(j) + rayt::sampling::Random()) / Real((height));

                        Point2 lensSample = sampling::Random2D();

                        Ray r = m_camera->getRay(u, v, lensSample);
                        pixelColor += Li(r, scene);
                    }
                    pixelColor /= Real(m_spp);

                    // NaN除去（デバッグ用）
                    if (HasInvalidValues(pixelColor)) {
                        std::cerr << "NaN detected at " << i << ", " << j << std::endl;
                        pixelColor = Spectrum(0.0);
                    }

                    // 上下反転して保存
                    film.setPixel(i, height - 1 - j, pixelColor);
                }
            }
            std::cout << "\n[PathIntegrator] Done." << std::endl;
        }

        // 放射輝度計算 (Li)
        Spectrum Li(Ray r, const Scene& scene) const {
            Spectrum L(0.0);        // 最終的な放射輝度（Accumulated Radiance）
            Spectrum beta(1.0);     // スループット（Throughput: 経路の重み）
            
            for (int depth = 0; depth < m_maxDepth; ++depth) {
                SurfaceInteraction rec;

                // 1. 交差判定
                if (!scene.hit(r, rec)) {
                    // 背景色 (IBLなどを使う場合はここで計算)
                    Spectrum skyColor(0.0, 0.0, 0.0); 
                    // L += beta * GetSkyColor(r); 
                    L += beta * skyColor;

                    break;
                }

                // 2. 自己発光の加算 (Le)
                // 光源に当たったら、ここまでの減衰(beta)を掛けて足す
                // ※ wo = -r.direction
                L += beta * rec.matPtr->emitted(rec, -r.d);

                // 3. 次の方向をサンプリング (Material::sample)
                // ランダムな乱数を用意 (本来はSamplerクラスから取得すべき)
                Point2 u(rayt::sampling::Random(), rayt::sampling::Random());

                // sample() 呼び出し: wo, uv を渡す
                auto bsdfSample = rec.matPtr->sample(rec, -r.d, u);

                // サンプリング失敗（吸収、全反射角超過など）なら終了
                if (!bsdfSample) {
                    break;
                }

                // 4. スループットの更新 (Beta update)
                // モンテカルロ積分の式: beta_new = beta_old * (f * cos_theta / pdf)

                Spectrum f = bsdfSample->f;
                Real pdf = bsdfSample->pdf;
                Vector3 wi = bsdfSample->wi; // 新しい方向

                // 鏡面反射（デルタ分布）かどうかの判定
                if (bsdfSample->isSpecular()) {
                    // ★ 鏡面反射の場合 (Specular)
                    // PDFは概念上無限大(Dirac Delta)なので、
                    // sample関数内で既に約分された重み (f / pdf) が f に入っているものとして扱う
                    beta *= f;
                }
                else {
                    // ★ 拡散・光沢反射の場合 (Diffuse / Glossy)
                    // 余弦項 (Cosine Term) = (n dot wi)
                    Real cosTheta = std::abs(glm::dot(rec.n, wi));

                    if (pdf > 1e-8f) { // ゼロ除算防止
                        beta *= f * cosTheta / pdf;
                    }
                    else {
                        break;
                    }
                }

                // スループットが0になったら計算打ち切り（ロシアンルーレットもここで入れると良い）
                if (isBlack(beta)) break;

                // 5. レイの更新
                // r = Ray(rec.p + rec.n * constants::RAY_EPSILON, wi);  old
                r = rayt::SpawnRay(rec.p, rec.n, wi);
            }

            return L;
        }

    private:
        std::shared_ptr<Camera> m_camera;
        int m_maxDepth;
        int m_spp;
    };

} // namespace rayt