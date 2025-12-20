// GoLD_rayt.cpp 
//
/// @brief 金属をJohnsonの物理データを基にレンダリング
///
/// 

// C++20を仮定している
// 1. Precompiled Header (Must be first)
#include "pch.h"    

// 自作ヘッダー
// 2. Project Headers
#include "IORInterpolator.hpp"
#include "Constants.hpp"
#include "Utils.hpp"
#include "Core.hpp"
#include "Film.hpp"
#include "Camera.hpp"
#include "Hittable.hpp"
#include "HittableList.hpp"
#include "Interaction.hpp"
#include "Material.hpp"
#include "Sphere.hpp"
#include "SpectralMetal.hpp"
#include "Scene.hpp"
#include "Integrator.hpp"
#include "Lambertian.hpp"
#include "DiffuseLight.hpp"
// Note: "Lambertian" or other materials can be added here in the future.

// 画像生成のためのヘッダー
// マクロを書く必要はなし
#include "stb_image_write.h"

using namespace rayt;

// -----------------------------------------------------------------------------
// Scene Configuration
// -----------------------------------------------------------------------------
const int IMAGE_WIDTH = 800;
const int IMAGE_HEIGHT = 450;      // 16:9 Aspect Ratio
const int SAMPLES_PER_PIXEL = 100; // Higher = less noise, slower
const int MAX_DEPTH = 50;          // Max recursion depth for rays

// -----------------------------------------------------------------------------
// The "Integrator" Function
// -----------------------------------------------------------------------------
// Recursively traces rays into the scene to calculate color.
// r: The current ray
// world: The list of all objects in the scene
// depth: Remaining bounce count
/*Spectrum ray_color(const Ray& r, const Hittable& world, int depth) {
    SurfaceInteraction rec;

    // 0. Base case: If ray bounces too many times, no more light is gathered.
    if (depth <= 0) {
        return Spectrum(0.0);
    }

    // 1. Intersection Check
    // tMin is 0.001 to avoid "Shadow Acne" (self-intersection due to float precision)
    if (world.hit(r, 0.001, Constants::INFINITY_VAL, rec)) {
        Ray scattered;
        Spectrum attenuation;

        // 2. Material Scattering
        // If the material scatters the ray (reflects/refracts), trace the new ray.
        if (rec.matPtr->scatter(r, rec, attenuation, scattered)) {
            // Recursive call: Multiply current attenuation by the light from the next bounce.
            return attenuation * ray_color(scattered, world, depth - 1);
        }

        // If the ray is absorbed (e.g., hit a black body), return black.
        return Spectrum(0.0);
    }

    // 3. Background (Sky)
    // If no object is hit, draw a simple gradient sky.
    Vector3 unit_direction = glm::normalize(r.d);
    auto t = 0.5 * (unit_direction.y + 1.0);

    // Linear interpolation: White -> Blue
    // (1.0-t)*White + t*Blue
    return (1.0 - t) * Spectrum(1.0, 1.0, 1.0) + t * Spectrum(0.5, 0.7, 1.0);
}*/

// -----------------------------------------------------------------------------
// Main Entry Point
// -----------------------------------------------------------------------------
int main() {
    // Initialize the Embree device.
    /*RTCDevice device = rtcNewDevice(nullptr);
    if (!device) {
        std::cerr << "Embree init failed\n";
        return -1;
    }
    std::cout << "Embree OK\n";
    rtcReleaseDevice(device);

    // 1. クラスのインスタンスを作成
    IORInterpolator metal_data;

    // 2. CSVファイルを読み込む
    // ※ Johnson.csv はプロジェクト(.vcxproj)と同じフォルダにある前提です
    std::string csv_path = "Johnson.csv";

    std::cout << "Loading: " << csv_path << "..." << std::endl;
    if (!metal_data.loadCSV(csv_path)) {
        std::cerr << "Failed to load CSV file." << std::endl;
        return -1; // エラー終了
    }

    // データが正しく読めたか範囲を表示
    metal_data.printInfo();

    // 3. 任意の波長 (nm) で屈折率を取得してみる
    // 例: 500nm (青緑色) の屈折率を取得
    double target_wl_nm = 550.0;
    std::complex<double> ior = metal_data.evaluate(target_wl_nm);

    std::cout << "\n--- Result at " << target_wl_nm << " nm ---" << std::endl;
    std::cout << "Complex IOR: " << ior << std::endl;
    std::cout << "  n (Refractive Index): " << ior.real() << std::endl;
    std::cout << "  k (Extinction Coeff): " << ior.imag() << std::endl;

    // 4. 垂直入射での反射率 (Reflectance) を計算する
    // 金属の反射率 R = |(n_complex - 1) / (n_complex + 1)|^2
    std::complex<double> num = ior - 1.0;
    std::complex<double> den = ior + 1.0;
    double reflectance = std::norm(num / den); // std::normは絶対値の2乗

    std::cout << "Reflectance: " << reflectance << " (" << reflectance * 100.0 << "%)" << std::endl;

    float roughness = 0.5f;

    // そのまま呼び出すスタイル
    float r2 = Utils::sqr(roughness);
    std::cout << "roughness" << r2 << std::endl;

    const int width = 800;
    const int height = 600;

    // Create the film (sensor)
    rayt::Film film(width, height);

    // ... Rendering Loop ...
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            // Assume we calculated some radiance
            rayt::Spectrum radiance(2.5, 0.5, 0.1); // Bright red light

            // Store it directly
            film.setPixel(x, y, radiance);
        }
    }

    // Save as standard image (Tone mapped & Gamma corrected)
    film.save("result.png");

    // Save as raw data for analysis (Linear & High dynamic range)
    film.save("result_raw.hdr");*/

    //------------------------------------------------------------------------------------------
    // debugの2回目のやつ；integratorを独立させていないものです------------------

    /*std::cout << "[System] Starting Rendering..." << std::endl;

    // --- 1. Scene Setup ---
    HittableList world;

    // Create Gold Material
    // Loads "Johnson.csv" to get physical n, k values for Gold (Au).
    // Roughness 0.2 gives a slightly blurry, realistic metal look.
    auto materialGold = std::make_shared<SpectralMetal>("Johnson.csv", 0.2);

    // Create a generic "Silver-ish" material (using the same class for test, assuming distinct CSV or fallback)
    // If "Silver.csv" doesn't exist, it might fallback or error. 
    // For now, let's make another Gold sphere but smoother.
    auto materialGoldPolished = std::make_shared<SpectralMetal>("Johnson.csv", 0.0);

    // Add Spheres to the world
    // Center Sphere (Rough Gold)
    world.add(std::make_shared<Sphere>(Point3(0, 0, -1), 0.5, materialGold));

    // Ground Sphere (Huge sphere to act as a floor)
    // Note: Since we only have Metal implemented, the floor will be a mirror!
    world.add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, materialGoldPolished));


    // --- 2. Camera Setup ---
    Point3 lookFrom(0, 0, 2);   // Positioned slightly back
    Point3 lookAt(0, 0, -1);    // Looking at the sphere
    Vector3 vUp(0, 1, 0);
    Real distToFocus = glm::length(lookFrom - lookAt);
    Real aperture = 0.1;        // Slight depth of field

    Camera cam(lookFrom, lookAt, vUp, 45.0, Real(IMAGE_WIDTH) / IMAGE_HEIGHT, aperture, distToFocus);


    // --- 3. Film (Image Buffer) Setup ---
    Film film(IMAGE_WIDTH, IMAGE_HEIGHT);


    // --- 4. Rendering Loop ---
    std::cout << "[Render] " << IMAGE_WIDTH << "x" << IMAGE_HEIGHT
        << " (" << SAMPLES_PER_PIXEL << " spp)" << std::endl;

    for (int j = 0; j < IMAGE_HEIGHT; ++j) {
        // Progress Indicator
        std::cout << "\rScanlines remaining: " << (IMAGE_HEIGHT - j) << " " << std::flush;

        for (int i = 0; i < IMAGE_WIDTH; ++i) {
            Spectrum pixelColor(0.0);

            // Anti-Aliasing Loop: Take multiple samples per pixel
            for (int s = 0; s < SAMPLES_PER_PIXEL; ++s) {
                // Randomized UV coordinates inside the pixel
                auto u = (double(i) + Utils::randomFloat()) / (IMAGE_WIDTH - 1);
                auto v = (double(j) + Utils::randomFloat()) / (IMAGE_HEIGHT - 1); // Note: verify your Camera/Texture Y-axis direction

                // Generate Ray
                Ray r = cam.getRay(u, v);

                // Accumulate Radiance
                pixelColor += ray_color(r, world, MAX_DEPTH);
            }

            // Average the color
            pixelColor /= Real(SAMPLES_PER_PIXEL);

            // Store in Film
            // Note: Camera Y might be inverted depending on setup. 
            // If image is upside down, use: film.setPixel(i, IMAGE_HEIGHT - 1 - j, pixelColor);
            film.setPixel(i, IMAGE_HEIGHT - 1 - j, pixelColor);
        }
    }
    std::cout << "\n[Render] Done." << std::endl;


    // --- 5. Save Output ---
    // Save as Tone-mapped PNG
    film.save("render_gold.png");

    // Save as Raw Linear HDR (For analysis)
    film.save("render_gold.hdr");*/

    const int samples_per_pixel = 100; // ← これが samples_per_pixel の定義
    const int max_depth = 50;          // ← これが max_depth の定義
    // -------------------------------------------------------------------------
    // 2. マテリアル作成 (Materials)
    // -------------------------------------------------------------------------
    // Johnson.csv から金 (Au) の物性データ (n, k) を読み込む
    // Roughness: 0.2 (少しざらついたリアルな金属)
    /*auto matRoughGold = std::make_shared<SpectralMetal>("Johnson.csv", 0.2);

    // Roughness: 0.0 (完全な鏡面/研磨された金)
    auto matPolishedGold = std::make_shared<SpectralMetal>("Johnson.csv", 0.0);

    // -------------------------------------------------------------------------
    // 3. 物体配置 (Geometry)
    // -------------------------------------------------------------------------
    // HittableList (物体リスト) をスマートポインタで作成
    auto worldObjects = std::make_shared<HittableList>();

    // 中央の球 (ざらついた金)
    worldObjects->add(std::make_shared<Sphere>(Point3(0, 0, -1), 0.5, matRoughGold));

    // 地面となる巨大な球 (鏡面の金)
    // ※ 現在は金属マテリアルしかないので、床も金になります！
    worldObjects->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matPolishedGold));*/

    // 金 (Roughness 0.2)
    auto matGold = std::make_shared<SpectralMetal>("Johnson.csv", 0.2);

    // 床用 (グレーの拡散反射)
    auto matFloor = std::make_shared<Lambertian>(Spectrum(0.5, 0.5, 0.5));

    // 照明 (非常に明るい白)
    // 値が 1.0 を超えるとHDRとして「輝く」光になります
    auto matLight = std::make_shared<DiffuseLight>(Spectrum(15.0, 15.0, 15.0));


    // -------------------------------------------------------------------------
    // 3. 物体配置
    // -------------------------------------------------------------------------
    auto worldObjects = std::make_shared<HittableList>();

    // メインの金色の球
    worldObjects->add(std::make_shared<Sphere>(Point3(0, 0, -1), 0.5, matGold));

    // 床 (Lambertianに変更)
    worldObjects->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));

    // ★光源を追加★
    // 球の少し右上、手前側に配置
    // 半径を小さくすると鋭いハイライト、大きくすると柔らかいハイライトになります
    worldObjects->add(std::make_shared<Sphere>(Point3(2, 2, 1), 0.5, matLight));

    // -------------------------------------------------------------------------
    // 4. シーン構築 (Scene Construction)
    // -------------------------------------------------------------------------
    // 物体リストを Scene クラスにラップする
    // (将来的に光源ライトなどもここで管理します)
    Scene scene(worldObjects);

    // -------------------------------------------------------------------------
    // 5. カメラ設定 (Camera Setup)
    // -------------------------------------------------------------------------
    Point3 lookFrom(0, 0, 2);   // カメラの位置 (少し手前)
    Point3 lookAt(0, 0, -1);    // 注視点 (球の中心)
    Vector3 vUp(0, 1, 0);       // 上方向

    // 焦点距離と絞り (被写界深度用)
    Real distToFocus = glm::length(lookFrom - lookAt);
    Real aperture = 0.05;       // 0.0ならピンホール(ボケなし)、大きくするとボケる

    auto camera = std::make_shared<Camera>(
        lookFrom, lookAt, vUp,
        45.0,                       // 視野角 (FOV)
        double(IMAGE_WIDTH) / IMAGE_HEIGHT,     // アスペクト比
        aperture,
        distToFocus
    );

    // -------------------------------------------------------------------------
    // 6. フィルム (画像バッファ) 作成
    // -------------------------------------------------------------------------
    Film film(IMAGE_WIDTH, IMAGE_HEIGHT);

    // -------------------------------------------------------------------------
    // 7. 積分器 (Integrator) 作成 & レンダリング実行
    // -------------------------------------------------------------------------
    // パストレーシング法を使用
    auto integrator = std::make_unique<PathIntegrator>(camera, max_depth, samples_per_pixel);

    // ★ ここで全ての計算が走ります ★
    integrator->render(scene, film);

    // -------------------------------------------------------------------------
    // 8. 保存 (Save Output)
    // -------------------------------------------------------------------------
    // PNG (トーンマップ済み・鑑賞用)
    film.save("render_gold_pbrt.png");

    // HDR (生データ・研究解析用)
    film.save("render_gold_pbrt.hdr");
    
    

    return 0;
}