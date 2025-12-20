#pragma once

#include "Core.hpp"
#include "Scene.hpp"   // HittableList, etc.
#include "Camera.hpp"
#include "Film.hpp"
#include "Ray.hpp"
#include "Interaction.hpp"
#include "Material.hpp"
#include "Utils.hpp"

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

                for (int i = 0; i < width; ++i) {
                    Spectrum pixelColor(0.0);

                    for (int s = 0; s < m_spp; ++s) {
                        Real u = (Real(i) + Utils::Random()) / Real((width - 1));
                        Real v = (Real(j) + Utils::Random()) / Real((height - 1));

                        Ray r = m_camera->getRay(u, v);
                        pixelColor += Li(r, scene, m_maxDepth);
                    }
                    pixelColor /= Real(m_spp);

                    // 上下反転して保存
                    film.setPixel(i, height - 1 - j, pixelColor);
                }
            }
            std::cout << "\n[PathIntegrator] Done." << std::endl;
        }

        // 放射輝度計算 (Li)
        Spectrum Li(const Ray& r, const Scene& scene, int depth) const {
            SurfaceInteraction rec;

            // 再帰深さ制限
            if (depth <= 0) {
                return Spectrum(0.0);
            }

            // 交差判定
            /*if (scene.hit(r, 0.001, Constants::INFINITY_VAL, rec)) {
                Ray scattered;
                Spectrum attenuation;

                // マテリアル散乱計算
                if (rec.matPtr->scatter(r, rec, attenuation, scattered)) {
                    return attenuation * Li(scattered, scene, depth - 1);
                }
                return Spectrum(0.0);
            }*/

            //---------------------------------------------------------------------------
            // 照明追加ver
            // 1. 交差判定
            if (!scene.hit(r, rec)) {
                // 背景 (Sky)
                // 光源を置くなら、背景は真っ黒(0,0,0)にした方が、照明効果が分かりやすいです
                return Spectrum(0.0, 0.0, 0.0); 
            }

            Ray scattered;
            Spectrum attenuation;

            // ★重要: まず物体からの「発光 (Emitted)」を取得
            Spectrum emitted = rec.matPtr->emitted(rec);

            // 2. マテリアル散乱 (Scatter)
            if (!rec.matPtr->scatter(r, rec, attenuation, scattered)) {
                // 散乱しなかった（光源に当たった場合など）は、発光だけを返す
                return emitted;
            }

            // 3. 再帰 (Recursive)
            // 結果 = 発光 + (減衰率 * 次のレイの光)
            return emitted + attenuation * Li(scattered, scene, depth - 1);

            // 背景色 (Sky)
            /*Vector3 unit_direction = glm::normalize(r.d);
            auto t = 0.5 * (unit_direction.y + 1.0);
            return (1.0 - t) * Spectrum(1.0, 1.0, 1.0) + t * Spectrum(0.5, 0.7, 1.0);*/
        }

    private:
        std::shared_ptr<Camera> m_camera;
        int m_maxDepth;
        int m_spp;
    };

} // namespace rayt