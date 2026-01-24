/*
* @file src/Application/Application.cpp
* @brief Application implementation
*/

// Precompiled header(PCH)
#include "pch.h"

#include "Application/Application.hpp"

#include <imgui.h>
#include <iostream>
#include <vector>
#include <filesystem>

// Project header
// Core
#include "Core/Core.hpp"
#include "Core/Interaction.hpp"
#include "Core/Sampling.hpp"
#include "Core/Image.hpp"

// Acceleration structures
#include "Accelerators/BVH.hpp" 

// Scene
#include "Scene/Aggregate.hpp"
#include "Scene/Scene.hpp"

// Geometry
#include "Geometry/Sphere.hpp"

// Renderer
#include "Render/Film.hpp"
#include "Render/Camera.hpp"
#include "Render/Integrator.hpp"

// Color
#include "Color/ColorTransform.hpp"

// Lights
#include "Lights/AreaLight.hpp"
#include "Lights/PointLight.hpp"

// Materials
#include "Materials/Material.hpp"
#include "Materials/Lambertian.hpp"
#include "Materials/DiffuseLight.hpp"
#include "Materials/RoughConductor.hpp"
#include "Materials/Dielectric.hpp"

// IO
#include "IO/ImageLoader.hpp"
#include "IO/EnvMap.hpp"

// 複素屈折率関連
#include "IO/ComplexIorLoader.hpp"
#include "Color/ComplexIorTable.hpp"
#include "Color/IorRgbApprox.hpp"

using namespace rayt; // Convenience: avoid qualifying with rayt:: in this .cpp

// -----------------------------------------------------------------------------
// Asset paths
// -----------------------------------------------------------------------------

// Environment map path
const std::string ENV_HDR_PATH = "assets/env/bryanston_park_sunrise_2k.hdr";

// Gold IOR data path (RefractiveIndex.info / Johnson & Christy style CSV)
const std::string AU_IOR_CSV_PATH = "assets/Metal_data/Johnson_Cu.csv";

// -----------------------------------------------------------------------------
// Scene Configuration
// -----------------------------------------------------------------------------
const int IMAGE_WIDTH = 800;
const int IMAGE_HEIGHT = 450;      // 16:9 Aspect Ratio
const int SAMPLES_PER_PIXEL = 256; // Higher = less noise, slower
const int MAX_DEPTH = 10;          // Max recursion depth for rays


App::App(int width, int height, const char* title)
    : m_width(width), m_height(height)
{
    // -----------------------------------------------------
    // GLFW / OpenGL Init
    // -----------------------------------------------------
    // Initialize GLFW (windowing + OpenGL context creation).
    if (!glfwInit()) exit(1);

    // Request an OpenGL 3.3 Core Profile context.
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SCALE_TO_MONITOR, GLFW_TRUE);

    // Create a window + OpenGL context.
    m_window = glfwCreateWindow(width, height, title, nullptr, nullptr);
    if (!m_window) {
        glfwTerminate();
        exit(1);
    }

    // Make the created context current before any OpenGL calls.
    glfwMakeContextCurrent(m_window);

    // V-Sync OFF: maximizes FPS but may cause tearing and higher GPU usage.
    glfwSwapInterval(0);

    // Load OpenGL function pointers via GLAD.
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) exit(1);

    // Project-specific initialization (renderer, resources, UI, etc.)
    init();
}

App::~App() {
    m_uiLayer.shutdown();
    glDeleteTextures(1, &m_renderTexture);
    glfwDestroyWindow(m_window);
    glfwTerminate();
}

void App::init() {
    // --------------------------------------------------
    // 1. UI & OpenGL テクスチャ初期化
    // --------------------------------------------------
    m_uiLayer.init(m_window);

    glGenTextures(1, &m_renderTexture);
    glBindTexture(GL_TEXTURE_2D, m_renderTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, m_width, m_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glBindTexture(GL_TEXTURE_2D, 0);

    // --------------------------------------------------
    // 2. 環境マップ (EnvMap) 読み込み
    // --------------------------------------------------
    std::cout << "CWD = " << std::filesystem::current_path() << std::endl;
    
    std::shared_ptr<EnvMap> envMap = nullptr;  // 環境マップクラスをインスタンス化
    std::shared_ptr<EnvLight> envLight = nullptr;

    try {
        std::cout << "[EnvMap] Loading: " << ENV_HDR_PATH << std::endl;
        auto envImg = io::loadHDR(ENV_HDR_PATH);   // HDRでダウンロード
        envMap = std::make_shared<EnvMap>(std::move(envImg));
        std::cout << "[EnvMap] Success.\n";

        // EnvMapからEnvLightを作成 (スケール1.0)
        envLight = std::make_shared<EnvLight>(envMap, 1.0);
        std::cout << "[EnvLight] Created from environment map.\n";
    }
    catch (const std::exception& e) {
        std::cerr << "[EnvMap] Failed: " << e.what() << "\n";
        std::cerr << "[EnvMap] Fallback to black background.\n";
        envMap = nullptr;
        envLight = nullptr;
    }

    // --------------------------------------------------
    // 3. マテリアル作成 
    // --------------------------------------------------

    // 拡散反射の床
    auto matFloor = std::make_shared<Lambertian>(Spectrum(0.5, 0.5, 0.5));

    // --------------------------------------------------
    // 金の光学定数 (Au) をCSVからロードして RGB近似
    // --------------------------------------------------
    rayt::color::ComplexIorTable auTable;
    std::string auErr;

    // fallback（ロード失敗時の保険）
    Vector3 eta_rgb_fallback(0.16, 0.42, 1.45);
    Vector3 k_rgb_fallback(3.48, 2.45, 1.77);

    Vector3 eta_rgb = eta_rgb_fallback;
    Vector3 k_rgb = k_rgb_fallback;

    if (!rayt::io::loadComplexIor_RefractiveIndexInfoCsv(AU_IOR_CSV_PATH, auTable, &auErr)) {
        std::cerr << "[Au IOR] Failed: " << auErr << "\n";
        std::cerr << "[Au IOR] Fallback to hard-coded eta/k.\n";
    }
    else {
        // 3波長サンプルでRGB近似（暫定。後でCMF積分に置換可能）
        rayt::color::approxIorToRgb_3Wavelengths(auTable, eta_rgb, k_rgb, 450, 550, 650);

        std::cout << "[Au IOR] Loaded: " << AU_IOR_CSV_PATH << "\n";
        std::cout << "[Au IOR] RGB approx eta=("
            << eta_rgb.x << "," << eta_rgb.y << "," << eta_rgb.z << ") "
            << "k=(" << k_rgb.x << "," << k_rgb.y << "," << k_rgb.z << ")\n";
    }

    // Spectrum が今は glm::vec3 相当の想定
    Spectrum n_Au = eta_rgb;
    Spectrum k_Au = k_rgb;

    // 3段階の粗さ
    // 0.01: ほぼ鏡 (MirrorConductorと比較用)
    // 0.20: 少しぼやけた金属
    // 0.50: マットな金属（ブラスト仕上げ風）
    auto matGoldSmooth = std::make_shared<RoughConductor>(n_Au, k_Au, 0.01);
    auto matGoldMedium = std::make_shared<RoughConductor>(n_Au, k_Au, 0.20);
    auto matGoldRough = std::make_shared<RoughConductor>(n_Au, k_Au, 0.50);

    //auto matGlass = std::make_shared<Dielectric>(1.5, 0.0); // 粗さ0 = 完全透明
    //auto matFrosted = std::make_shared<Dielectric>(1.5, 0.2); // 粗さ0.2 = すりガラス

    // --------------------------------------------------
    // 4. シーン構築
    // --------------------------------------------------
    auto worldList = std::make_shared<Aggregate>();

    // リストに対して add する
    worldList->add(std::make_shared<Sphere>(Point3(-1.2, 0, -1), 0.5, matGoldSmooth));
    worldList->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), 0.5, matGoldMedium));
    worldList->add(std::make_shared<Sphere>(Point3(1.2, 0, -1), 0.5, matGoldRough));

    // ガラス
    //worldList->add(std::make_shared<Sphere>(Point3(0.7, 0, -1), 0.5, matGlass));

    // すりガラス
    //worldList->add(std::make_shared<Sphere>(Point3(-0.7, 0, -1), 0.5, matFrosted));


    // 大きな球で床を作製
    worldList->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));

    // --------------------------------------------------
    // 5. Scene作成と光源設定
    // --------------------------------------------------
    m_scene = std::make_shared<Scene>(worldList);

    // 環境光を設定（存在する場合）
    if (envLight) {
        m_scene->setEnvLight(envLight);
        std::cout << "[Scene] Environment light set.\n";
    }

    // オプション: ポイントライトや他の光源を追加
    // 例: シーン内にポイントライトを追加する場合
    
    /*auto pointLight = std::make_shared<PointLight>(
        Point3(-1.2, 3, 0),           // 位置
        Spectrum(0, 0, 100)       // 強度
    );
    m_scene->addLight(pointLight);
    std::cout << "[Scene] Point light added.\n";*/

    // エリアライトの設定
    //m_scene->addLight(areaLight);
    //std::cout << "[Scene] Area light added.\n";

    // IMPORTANT: シーンをファイナライズ（LightSamplerを構築）
    m_scene->finalize("uniform");
    std::cout << "[Scene] Finalized with "
        << m_scene->numLights() << " scene lights";
    if (m_scene->hasEnvLight()) {
        std::cout << " + environment light";
    }
    std::cout << ".\n";

    // --------------------------------------------------
    // 6. カメラ設定 
    // --------------------------------------------------
    Point3 lookFrom(0, 0.5, 2.5); 
    Point3 lookAt(0, 0, -1);
    Vector3 vUp(0, 1, 0);

    Real distToFocus = glm::length(lookFrom - lookAt);
    Real aperture = 0.0; // ピンホール
    Real vfov = 35.0;    
    Real aspect = static_cast<Real>(m_width) / static_cast<Real>(m_height);

    m_camera = std::make_shared<Camera>(
        lookFrom, lookAt, vUp,
        vfov, aspect, aperture, distToFocus
    );

    // --------------------------------------------------
    // 7. レンダラー生成
    // --------------------------------------------------
    m_film = std::make_unique<Film>(m_width, m_height);

    // 新しいPathIntegratorのコンストラクタ: (camera, maxDepth, spp)
    m_integrator = std::make_unique<PathIntegrator>(
        m_camera,
        m_settings.maxDepth,  // 最大深度
        m_settings.targetSpp     // プログレッシブ用のspp
    );

    std::cout << "[Integrator] PathIntegrator created (maxDepth="
        << m_settings.maxDepth << ", spp=" << m_settings.targetSpp << ").\n";
}

// テクスチャ更新 
void App::updateTexture() {
    const auto* pixels = m_film->getData();
    float scale = 1.0f / std::max(1, m_currentSpp);

    static std::vector<unsigned char> displayBuffer;
    if (displayBuffer.size() != m_width * m_height * 4)
        displayBuffer.resize(m_width * m_height * 4);

    #pragma omp parallel for
    for (int i = 0; i < m_width * m_height; ++i) {
        // ... (前のコードと同じ変換処理) ...
        // 省略せずに書く場合は前回のコードを参照
        rayt::Spectrum s = pixels[i] * static_cast<rayt::Real>(scale);
        // Gamma correction & Tone mapping...
        const Real exposureStops = Real(0); // ここはUIや設定から渡す想定
        s = rayt::color::toDisplaySRGB_Reinhard(s, exposureStops);

        displayBuffer[i * 4 + 0] = (unsigned char)(255.0f * rayt::math::saturate((float)s.r));
        displayBuffer[i * 4 + 1] = (unsigned char)(255.0f * rayt::math::saturate((float)s.g));
        displayBuffer[i * 4 + 2] = (unsigned char)(255.0f * rayt::math::saturate((float)s.b));
        displayBuffer[i * 4 + 3] = 255;
    }

    glBindTexture(GL_TEXTURE_2D, m_renderTexture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, m_width, m_height, GL_RGBA, GL_UNSIGNED_BYTE, displayBuffer.data());
}

// メインループ 
void App::run() {
    while (!glfwWindowShouldClose(m_window)) {
        glfwPollEvents();

        if (m_settings.requestReset) {
            m_film->clear();
            m_currentSpp = 0;
            m_settings.requestReset = false;
        }

        if (!m_settings.paused && m_currentSpp < m_settings.targetSpp) {
            m_integrator->renderOnePass(*m_scene, *m_film, m_currentSpp);
            m_currentSpp++;
            updateTexture();
        }

        m_uiLayer.beginFrame();
        m_uiLayer.draw(m_settings);

        {
            ImGui::Begin("Viewport");
            ImVec2 size = ImGui::GetContentRegionAvail();
            ImGui::Image((void*)(intptr_t)m_renderTexture, size, ImVec2(0, 0), ImVec2(1, 1));
            ImGui::End();
        }

        int w, h;
        glfwGetFramebufferSize(m_window, &w, &h);
        glViewport(0, 0, w, h);
        glClear(GL_COLOR_BUFFER_BIT);

        m_uiLayer.endFrame();
        glfwSwapBuffers(m_window);
    }
}