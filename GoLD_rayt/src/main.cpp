// GoLD_rayt.cpp 
//
/// @brief 金属をJohnsonの物理データを基にレンダリング
///
/// 

// C++20を仮定している
// Precompiled Header (Must be first)
#include "pch.h"    

// 自作ヘッダー
// Project Headers
// Core
#include "Core/Constants.hpp"
//#include "Core/Utils.hpp"
#include "Core/Core.hpp"
#include "Core/Interaction.hpp"
#include "Core/Math.hpp"
#include "Core/Sampling.hpp"

// Geometry
#include "Geometry/Hittable.hpp"
#include "Geometry/HittableList.hpp"
#include "Geometry/Sphere.hpp"

// IO
#include "IO/IORInterpolator.hpp"

// Renderer
#include "Renderer/Film.hpp"
#include "Renderer/Camera.hpp"
#include "Renderer/Scene.hpp"
#include "Renderer/Integrator.hpp"

// Materials
#include "Materials/Material.hpp"
//#include "Materials/SpectralMetal.hpp"
#include "Materials/Lambertian.hpp"
#include "Materials/DiffuseLight.hpp"
#include "Materials/Lambertian.hpp"
#include "Materials/MirrorConductor.hpp"
#include "Materials/Mirror.hpp"
// Note: "Lambertian" or other materials can be added here in the future.

#include "IO/ImageLoader.hpp"
#include "IO/EnvMap.hpp"

#include <filesystem>


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

const std::string ENV_HDR_PATH = "assets/env/grace-new.hdr";

// -----------------------------------------------------------------------------
// Main Entry Point
// -----------------------------------------------------------------------------
int main() {


    // -------------------------------------------------------------------------
// EnvMap (HDRI) 読み込み
// -------------------------------------------------------------------------
    const std::string envPath = "assets/env/grace-new.hdr";

    std::cout << "CWD = " << std::filesystem::current_path() << std::endl;

    std::shared_ptr<rayt::EnvMap> env = nullptr;
    try {
        auto envImg = rayt::io::loadHDR(envPath);
        env = std::make_shared<rayt::EnvMap>(std::move(envImg));
        std::cout << "[EnvMap] Loaded: " << envPath << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "[EnvMap] Failed: " << e.what() << "\n";
        std::cerr << "[EnvMap] Fallback to black background.\n";
    }



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

    //const int samples_per_pixel = 100; // ← これが samples_per_pixel の定義
    //const int max_depth = 50;          // ← これが max_depth の定義
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



    /*
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
    film.save("render_gold_pbrt.hdr");*/





    std::cout << "[System] Initializing..." << std::endl;

    // -------------------------------------------------------------------------
    // 1. マテリアルの作成
    // -------------------------------------------------------------------------

    // --- 床用: グレーの拡散反射 ---
    auto matFloor = std::make_shared<Lambertian>(Spectrum(0.5f, 0.5f, 0.5f));

    // --- 光源: 明るい白 (DiffuseLightを使用) ---
    // 値が1.0を超えると発光体として機能します
    auto matLight = std::make_shared<DiffuseLight>(Spectrum(15.0f, 15.0f, 15.0f));

    // --- 主役: 金 (MirrorConductor) ---
    // 波長ごとの屈折率データ (RGB近似値)
    // Red(650nm), Green(550nm), Blue(450nm) 付近の値を設定
    // Au (Gold):
    // n: R=0.16, G=0.42, B=1.45
    // k: R=3.48, G=2.45, B=1.77
    Spectrum n_Au(0.16f, 0.42f, 1.45f);
    Spectrum k_Au(3.48f, 2.45f, 1.77f);

    auto matGold = std::make_shared<MirrorConductor>(n_Au, k_Au);

    // -------------------------------------------------------------------------
    // 2. 物体の配置 (Scene)
    // -------------------------------------------------------------------------
    auto worldObjects = std::make_shared<HittableList>();

    // 床 (巨大な球)
    worldObjects->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));

    // 中央の球 (金)
    worldObjects->add(std::make_shared<Sphere>(Point3(0, 0, -1), 0.5, matGold));

    // 光源 (右上・手前)
    //worldObjects->add(std::make_shared<Sphere>(Point3(1.5, 2.0, 1.0), 0.5, matLight));

    // 補助光源（左側・遠く）
    //worldObjects->add(std::make_shared<Sphere>(Point3(-2.0, 1.0, -2.0), 0.3, matLight));

    Scene scene(worldObjects);

    // -------------------------------------------------------------------------
    // 3. カメラ設定
    // -------------------------------------------------------------------------
    Point3 lookFrom(0, 0.5, 2.5); // 少し高い位置から見下ろす
    Point3 lookAt(0, 0, -1);
    Vector3 vUp(0, 1, 0);

    Real distToFocus = glm::length(lookFrom - lookAt);
    Real aperture = 0.0; // ピンホールカメラ（ボケなし）でテスト

    auto camera = std::make_shared<Camera>(
        lookFrom, lookAt, vUp,
        35.0, // FOV
        double(IMAGE_WIDTH) / IMAGE_HEIGHT,
        aperture,
        distToFocus
    );

    // -------------------------------------------------------------------------
    // 4. レンダリング準備
    // -------------------------------------------------------------------------
    Film film(IMAGE_WIDTH, IMAGE_HEIGHT);

    // 新しいIntegratorを使用
    // max_depth, spp を渡す
    //auto integrator = std::make_unique<PathIntegrator>(camera, MAX_DEPTH, SAMPLES_PER_PIXEL);
    auto integrator = std::make_unique<PathIntegrator>(camera, env, MAX_DEPTH, SAMPLES_PER_PIXEL);

    // -------------------------------------------------------------------------
    // 5. レンダリング実行
    // -------------------------------------------------------------------------
    std::cout << "[Render] Start PBR rendering..." << std::endl;
    integrator->render(scene, film);

    // -------------------------------------------------------------------------
    // 6. 保存
    // -------------------------------------------------------------------------
    std::cout << "[Output] Saving images..." << std::endl;
    film.save("result_gold_pbr.png");
    film.save("result_gold_pbr.hdr");

    std::cout << "[System] Finished." << std::endl;
    

    return 0;
}