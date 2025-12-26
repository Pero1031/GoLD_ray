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
#include "Core/Core.hpp"
#include "Core/Interaction.hpp"
#include "Core/Math.hpp"
#include "Core/Sampling.hpp"
#include "Core/AABB.hpp"
#include "Core/Assert.hpp"
#include "Core/Image.hpp"

// Geometry
#include "Geometry/Hittable.hpp"
#include "Geometry/HittableList.hpp"
#include "Geometry/Sphere.hpp"
#include "Geometry/Frame.hpp"

// IO
#include "IO/IORInterpolator.hpp"

// Renderer
#include "Renderer/Film.hpp"
#include "Renderer/Camera.hpp"
#include "Renderer/Scene.hpp"
#include "Renderer/Integrator.hpp"
#include "Renderer/BVH.hpp"

// Materials
#include "Materials/Material.hpp"
//#include "Materials/SpectralMetal.hpp"
#include "Materials/Lambertian.hpp"
#include "Materials/DiffuseLight.hpp"
#include "Materials/Lambertian.hpp"
#include "Materials/MirrorConductor.hpp"
#include "Materials/Mirror.hpp"
#include "Materials/RoughConductor.hpp"
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
const std::string ENV_HDR_PATH = "assets/env/grace-new.hdr";

// -----------------------------------------------------------------------------
// Main Entry Point
// -----------------------------------------------------------------------------
int main() {

    // debug frame 
    // rayt::debug::TestFrameRoundTrip();


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

    std::cout << "[System] Initializing..." << std::endl;

    // -------------------------------------------------------------------------
    // 1. マテリアルの作成
    // -------------------------------------------------------------------------

    // --- 床用: グレーの拡散反射 ---
    /*auto matFloor = std::make_shared<Lambertian>(Spectrum(0.5, 0.5, 0.5));

    // --- 光源: 明るい白 (DiffuseLightを使用) ---
    // 値が1.0を超えると発光体として機能します
    auto matLight = std::make_shared<DiffuseLight>(Spectrum(15.0, 15.0, 15.0));

    // --- 主役: 金 (MirrorConductor) ---
    // 波長ごとの屈折率データ (RGB近似値
    // Red(650nm), Green(550nm), Blue(450nm) 付近の値を設定
    // Au (Gold):
    // n: R=0.16, G=0.42, B=1.45
    // k: R=3.48, G=2.45, B=1.77
    Spectrum n_Au(0.16, 0.42, 1.45);
    Spectrum k_Au(3.48, 2.45, 1.77);

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

    Scene scene(worldObjects);*/

    //-------------------------------------------------------------------------------------------
    // ラフのテスト
    // -------------------------------------------------------------------------
    // 1. マテリアルの作成 (Roughness Test)
    // -------------------------------------------------------------------------

    // 床用
    auto matFloor = std::make_shared<Lambertian>(Spectrum(0.2, 0.2, 0.2)); // 少し暗くして反射を目立たせる

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

    // -------------------------------------------------------------------------
    // 2. 物体の配置 (Scene)
    // -------------------------------------------------------------------------
    auto worldObjects = std::make_shared<HittableList>();

    // 床
    worldObjects->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));

    // 球を横に3つ並べる
    // 左: ツルツル
    worldObjects->add(std::make_shared<Sphere>(Point3(-1.2, 0, -1), 0.5, matGoldSmooth));

    // 中央: 少し粗い
    worldObjects->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), 0.5, matGoldMedium));

    // 右: かなり粗い
    worldObjects->add(std::make_shared<Sphere>(Point3(1.2, 0, -1), 0.5, matGoldRough));

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
    // film.save("result_gold_pbr.hdr");

    std::cout << "[System] Finished." << std::endl;
    

    return 0;
}