#pragma once

#include "Core/Core.hpp"
#include "renderer/Scene.hpp"   // HittableList, etc.
#include "Renderer/Camera.hpp"
#include "renderer/Film.hpp"
#include "Core/Ray.hpp"
#include "Core/Interaction.hpp"
#include "Materials/Material.hpp"
#include "Core/Sampling.hpp"
#include "IO/EnvMap.hpp"

#include <memory>
#include <iostream>
#include <algorithm>

namespace rayt {

    // Abstract base class for all rendering algorithms.
    class Integrator {
    public:
        virtual ~Integrator() = default;

        // 純粋仮想関数
        // ここで Scene や Film の型を使うので、上のincludeが必須です
        virtual void render(const Scene& scene, Film& film) = 0;
    };

    // Path Tracing Integrator
    class PathIntegrator : public Integrator {
    public:
        // コンストラクタ
        PathIntegrator(std::shared_ptr<Camera> camera, 
            std::shared_ptr<EnvMap> env, 
            int maxDepth, int spp)
            : m_camera(camera), m_env(env), 
            m_maxDepth(maxDepth), m_spp(spp) {}

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
            Real lastPdf = 0;
            bool lastSpecular = false;
            bool hasLastBsdf = false;
            
            for (int depth = 0; depth < m_maxDepth; ++depth) {
                SurfaceInteraction rec;

                if (!scene.hit(r, rec)) {

                    if (m_env) {
                        Spectrum envL;
                        glm::vec3 rgb = m_env->eval(r.d);
                        envL = Spectrum(rgb.x, rgb.y, rgb.z);

                        if (hasLastBsdf && !lastSpecular) {
                            Real pdfEnv = m_env->pdf(r.d);   // ★ EnvMap に pdf(dir) を用意しておく

                            Real w = 1.0;
                            if (pdfEnv > 0 && lastPdf > 0) {
                                Real a = lastPdf;
                                Real b = pdfEnv;
                                w = (a * a) / (a * a + b * b); // power heuristic
                            }
                            L += beta * envL * w;
                        }
                        else {
                            // カメラレイ直撃 or 鏡面経路は MIS しない
                            L += beta * envL;
                        }
                    }
                    break;
                }


                // 2. 自己発光の加算 (Le)
                // 光源に当たったら、ここまでの減衰(beta)を掛けて足す
                // ※ wo = -r.direction
                L += beta * rec.matPtr->emitted(rec, -r.d);

                // 2.5. Next Event Estimation (Environment Light)
                if (m_env && !rec.matPtr->isSpecular()) {

                    Point2 uLight(sampling::Random(), sampling::Random());

                    Vector3 wi;
                    Real pdfEnv;
                    Vector3 Le = m_env->sample(uLight, wi, pdfEnv);

                    if (pdfEnv > 0 && !isBlack(Le)) {

                        // シャドウレイ
                        Ray shadow = SpawnRay(rec.p, rec.gn, wi);

                        SurfaceInteraction tmp;
                        if (!scene.hit(shadow, tmp)) {

                            // BSDF評価
                            Spectrum f = rec.matPtr->eval(rec, -r.d, wi);
                            if (isBlack(f)) continue; // or continue;

                            // cos項は abs を取る（重要）
                            Real cosTheta = std::abs(glm::dot(rec.n, wi));

                            // BSDF側 pdf
                            Real pdfBsdf = rec.matPtr->pdf(rec, -r.d, wi);

                            // MIS（Power heuristic）
                            Real w = 1.0;
                            if (pdfBsdf > 0) {
                                Real a = pdfEnv;
                                Real b = pdfBsdf;
                                w = (a * a) / (a * a + b * b);
                            }

                            L += beta * f * Spectrum(Le.x, Le.y, Le.z)
                                * cosTheta * (w / pdfEnv);
                        }
                    }
                }


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

                lastPdf = pdf;
                lastSpecular = bsdfSample->isSpecular();
                hasLastBsdf = true;

                // 鏡面反射（デルタ分布）かどうかの判定
                if (bsdfSample->isSpecular()) {
                    // ★ 鏡面反射の場合 (Specular)
                    // PDFは概念上無限大(Dirac Delta)なので、
                    // sample関数内で既に約分された重み (f / pdf) が f に入っているものとして扱う
                    beta *= f;
                }
                else {
                    Real cosTheta = std::abs(glm::dot(rec.n, wi));
                    if (pdf > 1e-8f)
                        beta *= f * cosTheta / pdf;
                    else
                        break;
                }

                // スループットが0になったら計算打ち切り（ロシアンルーレットもここで入れると良い）
                if (isBlack(beta)) break;

                // 5. レイの更新
                // r = Ray(rec.p + rec.n * constants::RAY_EPSILON, wi);  old
                r = rayt::SpawnRay(rec.p, rec.gn, wi);
            }

            return L;
        }

        void renderOnePass(const Scene& scene, Film& film, int sampleIndex) {
            int width = film.width();
            int height = film.height();

            // 並列化 (OpenMP)
            #pragma omp parallel for schedule(dynamic, 1)
            for (int j = 0; j < height; ++j) {
                for (int i = 0; i < width; ++i) {

                    // 1. 1回だけサンプリング
                    Real u = (Real(i) + rayt::sampling::Random()) / Real(width);
                    Real v = (Real(j) + rayt::sampling::Random()) / Real(height);

                    Point2 lensSample = sampling::Random2D();
                    Ray r = m_camera->getRay(u, v, lensSample);

                    Spectrum L = Li(r, scene);

                    // NaNチェック
                    if (HasInvalidValues(L)) L = Spectrum(0.0);

                    // 2. Filmに「加算」する
                    // 注意: Filmクラスの実装によりますが、ここでは「現在の画素値」に足し込む想定です。
                    // もしFilmが「平均値」しか持てない場合は、
                    // new_avg = (old_avg * (N-1) + new_val) / N の計算が必要です。
                    // ここではシンプルに「累積加算」するためのメソッド addSample があると仮定、
                    // なければ setPixel で工夫します。

                    film.addPixel(i, height - 1 - j, L);
                }
            }
        }

    private:
        std::shared_ptr<Camera> m_camera;
        std::shared_ptr<EnvMap> m_env;

        int m_maxDepth;
        int m_spp;

        static bool visible(const Scene& scene, const SurfaceInteraction& ref,
            const Point3& pLight)
        {
            Vector3 toL = pLight - ref.p;
            Real dist = glm::length(toL);
            if (dist <= Real(0)) return false;

            Vector3 wi = toL / dist;

            // 影レイ：ライト手前まで（eps分手前）
            Ray shadow = rayt::SpawnRay(ref.p, ref.gn, wi);
            shadow.tMin = constants::RAY_EPSILON;
            shadow.tMax = dist - constants::RAY_EPSILON;

            SurfaceInteraction tmp;
            return !scene.hit(shadow, tmp);
        }
    };

} // namespace rayt 