// GoLD_rayt.cpp 
//
/// @brief 金属をJohnsonの物理データを基にレンダリング
///
/// 

// C++20を仮定している
// Precompiled Header (Must be first)
#include "pch.h"    

#define GLFW_INCLUDE_NONE
#include <glad/glad.h>
#include <GLFW/glfw3.h>

//#include <imgui.h>
//#include <backends/imgui_impl_glfw.h>
//#include <backends/imgui_impl_opengl3.h>

#include "UI/ImGuiLayer.hpp"
#include "Core/RenderSettings.hpp"

// 自作ヘッダー
// Project Headers
// Core
#include "Core/Constants.hpp"
#include "Core/Core.hpp"
#include "Core/Interaction.hpp"
#include "Core/Math.hpp"
#include "Core/Sampling.hpp"
#include "Core/AABB.hpp"
#include "Core/Assert.hpp"
#include "Core/Image.hpp"

// accelerates
#include "Accelerators/BVH.hpp" 

// Geometry
#include "Geometry/Hittable.hpp"
#include "Geometry/Aggregate.hpp"
#include "Geometry/Sphere.hpp"
#include "Geometry/Frame.hpp"

// Renderer
#include "Render/Film.hpp"
#include "Render/Camera.hpp"
#include "Render/Scene.hpp"
#include "Render/Integrator.hpp"

// Materials
#include "Materials/Material.hpp"
//#include "Materials/SpectralMetal.hpp"
#include "Materials/Lambertian.hpp"
#include "Materials/DiffuseLight.hpp"
#include "Materials/Lambertian.hpp"
//#include "Materials/MirrorConductor.hpp"
//#include "Materials/Mirror.hpp"
#include "Materials/RoughConductor.hpp"
#include "Materials/Dielectric.hpp"
// Note: "Lambertian" or other materials can be added here in the future.

// Microfacet
#include "Microfacet/Distribution.hpp"
#include "Microfacet/GGX.hpp"

// IO
#include "IO/ImageLoader.hpp"
#include "IO/EnvMap.hpp"

#include <filesystem>

#include "DebugTools/FrameDebug.hpp"


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

// env path
const std::string ENV_HDR_PATH = "assets/env/bryanston_park_sunrise_2k.hdr";

// -----------------------------------------------------------------------------
// Main Entry Point
// -----------------------------------------------------------------------------
//int main() {

    // debug frame 
    // rayt::debug::TestFrameRoundTrip();


// -------------------------------------------------------------------------
// EnvMap (HDRI) 読み込み
// -------------------------------------------------------------------------
    //const std::string envPath = "assets/env/bryanston_park_sunrise_2k.hdr";

    /*std::cout << "CWD = " << std::filesystem::current_path() << std::endl;

    std::shared_ptr<rayt::EnvMap> env = nullptr;
    try {
        auto envImg = rayt::io::loadHDR(ENV_HDR_PATH);
        env = std::make_shared<rayt::EnvMap>(std::move(envImg));
        std::cout << "[EnvMap] Loaded: " << ENV_HDR_PATH << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "[EnvMap] Failed: " << e.what() << "\n";
        std::cerr << "[EnvMap] Fallback to black background.\n";
    }

    std::cout << "[System] Initializing..." << std::endl;

    //-------------------------------------------------------------------------------------------
    // ラフのテスト
    // -------------------------------------------------------------------------
    // 1. マテリアルの作成 (Roughness Test)
    // -------------------------------------------------------------------------

    // 床用
    auto matFloor = std::make_shared<Lambertian>(Spectrum(0.5, 0.5, 0.5)); // 少し暗くして反射を目立たせる

    // 金の光学定数 (Au)
    Spectrum n_Au(0.16, 0.42, 1.45);
    Spectrum k_Au(3.48, 2.45, 1.77);

    // ★比較用: 3段階の粗さを作成
    // 0.01: ほぼ鏡 (MirrorConductorと比較用)
    // 0.20: 少しぼやけた金属
    // 0.50: マットな金属（ブラスト仕上げ風）
    auto matGoldSmooth = std::make_shared<RoughConductor>(n_Au, k_Au, 0.01);
    auto matGoldMedium = std::make_shared<RoughConductor>(n_Au, k_Au, 0.20);
    auto matGoldRough = std::make_shared<RoughConductor>(n_Au, k_Au, 0.50);

    //auto matGlass = std::make_shared<Dielectric>(1.5, 0.0); // 粗さ0 = 完全透明
    //auto matFrosted = std::make_shared<Dielectric>(1.5, 0.2); // 粗さ0.2 = すりガラス

    // -------------------------------------------------------------------------
    // 2. 物体の配置 (Scene)
    // -------------------------------------------------------------------------
    //auto worldObjects = std::make_shared<HittableList>();
    // -------------------------------------------------------------------------
    // 2. 物体の配置 (Scene) - 1万個チャレンジ
    // -------------------------------------------------------------------------
    auto worldObjects = std::make_shared<HittableList>();

    // 乱数生成器の準備 (配置用)
    std::mt19937 rng(12345); // シード固定で毎回同じ配置にする
    std::uniform_real_distribution<Real> distPos(-15.0, 15.0); // X, Yの範囲
    std::uniform_real_distribution<Real> distDepth(-30.0, 0.0); // Zの範囲 (手前から奥へ)
    std::uniform_real_distribution<Real> distRadius(0.05, 0.1);  // 半径の範囲
    std::uniform_real_distribution<Real> distMat(0.0, 1.0);     // マテリアル抽選用

    std::cout << "[Scene] Generating 10,000 spheres..." << std::endl;

    // ★ 1万個ループ
    for (int i = 0; i < 10000; ++i) {
        // ランダムな位置とサイズ
        Real x = distPos(rng);
        Real y = distPos(rng) * 0.6; // Y方向は少し平べったく
        Real z = distDepth(rng);
        Real r = distRadius(rng);

        Point3 center(x, y, z);

        // マテリアルをランダムに割り当て
        std::shared_ptr<Material> mat;
        Real roll = distMat(rng);

        if (roll < 0.33) {
            mat = matGoldSmooth; // ピカピカの金
        }
        else if (roll < 0.66) {
            mat = matGoldMedium; // 普通の金
        }
        else {
            mat = matGoldRough;  // マットな金
        }

        worldObjects->add(std::make_shared<Sphere>(center, r, mat));
    }

    // 床も一応置いておく
    // worldObjects->add(std::make_shared<Sphere>(Point3(0, -1000.0 - 5.0, -1), 1000.0, matFloor));

    // ==========================================
    // ★ BVHの構築 (ここが最重要)
    // ==========================================
    std::cout << "[System] Building BVH for " << worldObjects->objects.size() << " objects..." << std::endl;

    // 時間計測しても面白いです
    auto startBuild = std::chrono::high_resolution_clock::now();

    auto bvhRoot = std::make_shared<BVHNode>(worldObjects->objects, 0, worldObjects->objects.size());

    auto endBuild = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> buildTime = endBuild - startBuild;
    std::cout << "[System] BVH Built in " << buildTime.count() << " ms" << std::endl;

    // BVHをシーンに渡す
    Scene scene(bvhRoot);

    // 床
    //worldObjects->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));

    // 球を横に3つ並べる
    // 左: ツルツル
    //worldObjects->add(std::make_shared<Sphere>(Point3(-1.2, 0, -1), 0.5, matGoldSmooth));

    // 中央: 少し粗い
    //worldObjects->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), 0.5, matGoldMedium));

    // 右: かなり粗い
    //worldObjects->add(std::make_shared<Sphere>(Point3(1.2, 0, -1), 0.5, matGoldRough));

    // ガラス
    //worldObjects->add(std::make_shared<Sphere>(Point3(0.7, 0, -1), 0.5, matGlass));

    // すりガラス
    //worldObjects->add(std::make_shared<Sphere>(Point3(-0.7, 0, -1), 0.5, matFrosted));

    // ==========================================
    // ★ BVHの構築 (ここを追加！)
    // ==========================================
    //std::cout << "[System] Building BVH..." << std::endl;

    // HittableListの中にある object (vector) を取り出して、BVHNodeを作る
    // 引数: (オブジェクトのリスト, 開始インデックス, 終了インデックス)
    //auto bvhRoot = std::make_shared<BVHNode>(worldObjects->objects, 0, worldObjects->objects.size());

    // ★ リストではなく、爆速になった bvhRoot をシーンに渡す
    //Scene scene(bvhRoot);

    //Scene scene(worldObjects);

    // -------------------------------------------------------------------------
    // 3. カメラ設定
    // -------------------------------------------------------------------------
    /*Point3 lookFrom(0, 0.5, 2.5); // 少し高い位置から見下ろす
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
    );*/

    // -------------------------------------------------------------------------
    // 3. カメラ設定 (1万個用の広角・俯瞰ショット)
    // -------------------------------------------------------------------------
    // 少し手前・上空から、奥の群集を見下ろす位置
    /*Point3 lookFrom(0, 10.0, 15.0);

    // 群集の中心あたりを見る
    Point3 lookAt(0, 0, -15.0);

    Vector3 vUp(0, 1, 0);

    Real distToFocus = glm::length(lookFrom - lookAt);
    Real aperture = 0.0;

    // FOVを 35 -> 50 に広げて、視界を広くする
    auto camera = std::make_shared<Camera>(
        lookFrom, lookAt, vUp,
        50.0, // FOV: 広角気味にして全体を映す
        double(IMAGE_WIDTH) / IMAGE_HEIGHT,
        aperture,
        distToFocus
    );
    

    //--------------------------------------------------------------------------
    // BVHのデバッグ
    
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
    // film.save("result_gold_pbr.hdr");

    std::cout << "[System] Finished." << std::endl;
    
    
    

    

    return 0;
}*/


#include "Application/Application.hpp"

int main()
{
    App app(1280, 800, "GoLD_rayt Research UI");
    app.run();
    return 0;
}
    