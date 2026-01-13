#pragma once

#include <memory>
#include <iostream>
#include <algorithm>

#include "Core/Forward.hpp"
#include "Core/Types.hpp"
#include "Core/Math.hpp"
#include "Core/Ray.hpp"
#include "Core/Sampling.hpp"
#include "Core/Constants.hpp"
#include "Core/Interaction.hpp"
#include "Core/SpectrumUtils.hpp"
#include "Render/Scene.hpp" // HittableList, etc.
#include "Render/Film.hpp"
#include "Render/Camera.hpp"
#include "Lights/Light.hpp"
#include "Lights/EnvLight.hpp"
#include "Lights/PointLight.hpp"
#include "Lights/LightSampler.hpp"
#include "Materials/Material.hpp"
#include "IO/EnvMap.hpp"

namespace rayt {

    // Abstract base class for all rendering algorithms.
    class Integrator {
    public:
        virtual ~Integrator() = default;

        // 純粋仮想関数
        // ここで Scene や Film の型を使うので、上のincludeが必須です
        virtual void render(const Scene& scene, Film& film) = 0;
    };

    /**
     * @brief Path Tracing Integrator with Multiple Importance Sampling.
     *
     * Features:
     * - Supports arbitrary number of lights via LightSampler
     * - Environment/infinite light support
     * - Multiple Importance Sampling (MIS) for variance reduction
     * - Handles both delta and non-delta lights correctly
     */
    class PathIntegrator : public Integrator {
    public:
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
                        Real u = (Real(i) + sampling::Random()) / Real((width));
                        Real v = (Real(j) + sampling::Random()) / Real((height));

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

        /**
         * @brief Compute radiance along a ray using path tracing.
         */
        Spectrum Li(Ray r, const Scene& scene) const {
            Spectrum L(0.0);        // Accumulated Radiance
            Spectrum beta(1.0);     // Throughput

            Real lastPdf = 0;
            bool lastSpecular = false;
            bool hasLastBsdf = false;

            SurfaceInteraction lastRec;
            
            for (int depth = 0; depth < m_maxDepth; ++depth) {
                SurfaceInteraction rec;

                // ========== Ray Miss: Environment Light ==========
                if (!scene.hit(r, rec)) {

                    if (scene.hasEnvLight()) {
                        SampledWavelengths lambda;
                        Spectrum envL = scene.envLight()->Le(r, lambda);

                        if (hasLastBsdf && !lastSpecular) {
                            // MIS with BSDF sampling
                            LightSampleContext ctx{ lastRec.p, lastRec.gn };
                            Real pdfEnv = scene.envLight()->pdfLi(ctx, r.d);

                            Real w = 1.0;
                            if (pdfEnv > 0 && lastPdf > 0) {
                                w = math::powerHeuristic(lastPdf, pdfEnv);
                            }
                            L += beta * envL * w;
                        }
                        else {
                            // Camera ray or specular path: no MIS
                            L += beta * envL;
                        }
                    }
                    break;
                }

                // 現在のヒット点を保存
                lastRec = rec;

                // ========== Emissive Surface Hit ==========
                // ※ wo = -r.direction
                L += beta * rec.matPtr->emitted(rec, -r.d);

                // ========== Next Event Estimation (NEE) ==========
                if (!rec.matPtr->isSpecular()) {
                    L += beta * sampleAllLights(scene, rec, -r.d);
                }

                // ========== BSDF Sampling for Next Path Vertex ==========
                Point2 u(sampling::Random(), sampling::Random());
                auto bsdfSample = rec.matPtr->sample(rec, -r.d, u);

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
                r = rec.SpawnRay(wi);
            }

            return L;
        }

        /**
         * @brief Progressive rendering: render one sample per pixel.
         */
        void renderOnePass(const Scene& scene, Film& film, int sampleIndex) {
            int width = film.width();
            int height = film.height();

            // 並列化 (OpenMP)
            #pragma omp parallel for schedule(dynamic, 1)
            for (int j = 0; j < height; ++j) {
                for (int i = 0; i < width; ++i) {

                    // 1. 1回だけサンプリング
                    Real u = (Real(i) + sampling::Random()) / Real(width);
                    Real v = (Real(j) + sampling::Random()) / Real(height);

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
        int m_maxDepth;
        int m_spp;

        /**
         * @brief Sample all lights (scene lights + environment) with MIS.
         *
         * Strategy:
         * 1. Count total lights (scene + env if present)
         * 2. Uniformly select one light
         * 3. Sample the chosen light
         * 4. Apply MIS weight combining light PDF and BSDF PDF
         * 5. Account for light selection probability
         */
        Spectrum sampleAllLights(const Scene& scene,
            const SurfaceInteraction& rec,
            const Vector3& wo) const {
            Spectrum Ld(0.0);
            SampledWavelengths lambda;
            LightSampleContext ctx{ rec.p, rec.gn };

            // Count available lights
            int numSceneLights = scene.hasLights() ? scene.numLights() : 0;
            size_t numEnvLights = scene.hasEnvLight() ? 1 : 0;
            int totalLights = numSceneLights + numEnvLights;

            if (totalLights == 0) return Ld;

            // ========== Light Selection Strategy ==========
            // Uniformly select one light from all available lights
            Real uSelect = sampling::Random();
            Real lightSelectPdf = Real(1) / Real(totalLights);

            const Light* sampledLight = nullptr;

            if (uSelect < Real(numSceneLights) / Real(totalLights)) {
                // Sample from scene lights using LightSampler
                Real uSampler = uSelect * Real(totalLights) / Real(numSceneLights);
                auto lightSample = scene.lightSampler()->sample(uSampler);

                if (lightSample) {
                    sampledLight = lightSample->light;
                    // Note: lightSample->pdf is the probability of selecting
                    // this light among scene lights, but we need the probability
                    // among ALL lights, so we multiply by the ratio
                    //lightSelectPdf *= lightSample->pdf;
                }
            }
            else if (scene.hasEnvLight()) {
                // Sample environment light
                sampledLight = scene.envLight();
            }

            if (!sampledLight) return Ld;

            // ========== Sample Chosen Light ==========
            Point2 uLight(sampling::Random(), sampling::Random());
            auto ls = sampledLight->sampleLi(ctx, uLight, lambda);

            if (!ls || ls->pdf <= 0 || isBlack(ls->Li)) {
                return Ld;
            }

            Vector3 wi = ls->wi;
            Spectrum Li = ls->Li;
            Real pdfLight = ls->pdf;  // PDF w.r.t. solid angle

            // ========== Visibility Test ==========
            if (!visible(scene, rec, wi, ls->tMax)) {
                return Ld;
            }

            // ========== BSDF Evaluation ==========
            Spectrum f = rec.matPtr->eval(rec, wo, wi);
            if (isBlack(f)) return Ld;

            Real cosTheta = std::abs(glm::dot(rec.n, wi));

            // ========== MIS Weight ==========
            Real w = 1.0;
            if (!ls->isDelta) {
                Real pdfBsdf = rec.matPtr->pdf(rec, wo, wi);
                if (pdfBsdf > 0) {
                    w = math::powerHeuristic(pdfLight, pdfBsdf);
                }
            }

            // ========== Final Contribution ==========
            // Li already includes geometry term from light sampling
            // We need: f * Li * cos / (pdfLight * lightSelectPdf)
            Ld = f * Li * cosTheta * w / (pdfLight * lightSelectPdf);

            return Ld;
        }

        /**
         * @brief Check visibility between surface point and light.
         */
        static bool visible(const Scene& scene,
            const SurfaceInteraction& ref,
            const Vector3& wi,
            Real tMax) {

            Ray shadow = ref.SpawnRay(wi, tMax);

            SurfaceInteraction tmp;
            return !scene.hit(shadow, tmp);
        }

    };

} // namespace rayt 