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

                // #pragma omp parallel for // 可能なら並列化推奨

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

                // 1. 交差判定   
                /*if (!scene.hit(r, rec)) {     旧コード
                    // 背景色 (IBLなどを使う場合はここで計算)
                    //Spectrum skyColor(0.0, 0.0, 0.0); 
                    // L += beta * GetSkyColor(r); 
                    //L += beta * skyColor;

                    //break;

                    Spectrum envL(0.0);

                    if (m_env) {
                        glm::vec3 rgb = m_env->eval(r.d);
                        envL = Spectrum(rgb.x, rgb.y, rgb.z);
                        // 任意：明るさ調整
                        // envL *= 1.0f;
                    }

                    L += beta * envL;
                    break;
                }*

                if (!scene.hit(r, rec)) {

                    Spectrum envL(0.0);
                    if (m_env) {
                        glm::vec3 rgb = m_env->eval(r.d);
                        envL = Spectrum(rgb.x, rgb.y, rgb.z);

                        // MIS weight（BSDFサンプル側）
                        Real pdfEnv = m_env->pdf(r.d);
                        Real pdfBsdf = pdf; // 直前のBSDFサンプルのpdf

                        Real w = 1.0;
                        if (pdfEnv > 0) {
                            Real a = pdfBsdf;
                            Real b = pdfEnv;
                            w = (a * a) / (a * a + b * b);
                        }

                        L += beta * envL * w;
                    }
                    break;
                }*/

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
                /*if (m_env && !rec.matPtr->isSpecular()) {

                    // --- サンプル ---
                    Point2 uLight(sampling::Random(), sampling::Random());

                    Vector3 wi;
                    Real pdfEnv;
                    Vector3 Le = m_env->sample(uLight, wi, pdfEnv);

                    if (pdfEnv > 0 && !isBlack(Le)) {

                        // 幾何的に正しい半球チェック（幾何法線）
                        if (glm::dot(rec.gn, wi) > 0) {

                            // シャドウレイ
                            Ray shadow = SpawnRay(rec.p, rec.gn, wi);

                            SurfaceInteraction tmp;
                            if (!scene.hit(shadow, tmp)) {

                                // BSDF評価
                                Spectrum f = rec.matPtr->eval(rec, -r.d, wi);

                                if (!isBlack(f)) {
                                    Real cosTheta = std::max<Real>(0, glm::dot(rec.n, wi));

                                    // BSDF側pdf（MIS用）
                                    Real pdfBsdf = rec.matPtr->pdf(rec, -r.d, wi);

                                    // MIS（Power Heuristic）
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
                    }
                }*/

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
                /*else {
                    // ★ 拡散・光沢反射の場合 (Diffuse / Glossy)
                     
                    if (glm::dot(rec.gn, wi) <= 0) {
                        break;
                    }

                    // 余弦項 (Cosine Term) = (n dot wi)
                    Real cosTheta = std::max<Real>(0, glm::dot(rec.n, wi));

                    if (pdf > 1e-8f) { // ゼロ除算防止
                        beta *= f * cosTheta / pdf;
                    }
                    else {
                        break;
                    }
                }*/

                // スループットが0になったら計算打ち切り（ロシアンルーレットもここで入れると良い）
                if (isBlack(beta)) break;

                // 5. レイの更新
                // r = Ray(rec.p + rec.n * constants::RAY_EPSILON, wi);  old
                r = rayt::SpawnRay(rec.p, rec.gn, wi);
            }

            return L;
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