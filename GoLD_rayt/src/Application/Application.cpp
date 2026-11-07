/*
* @file src/Application/Application.cpp
* @brief Application implementation
*/

// Precompiled header(PCH)
#include "pch.h"

#include "Application/Application.hpp"

#include <imgui.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
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
#include "Scene/Primitive.hpp"

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

// Complex refractive index
#include "IO/ComplexIorLoader.hpp"
#include "Color/ComplexIorTable.hpp"
#include "Color/IorRgbApprox.hpp"

using namespace rayt; // Convenience: avoid qualifying with rayt:: in this .cpp

// -----------------------------------------------------------------------------
// Asset paths
// -----------------------------------------------------------------------------

// Environment map path
//const std::string ENV_HDR_PATH = "assets/env/grace-new.hdr";
//const std::string ENV_HDR_PATH = "assets/env/bryanston_park_sunrise_2k.hdr";
const std::string ENV_HDR_PATH = "assets/env/stpeters_probe_latlong.hdr";

// Gold IOR data path (RefractiveIndex.info / Johnson & Christy style CSV)
const std::string AU_IOR_CSV_PATH = "assets/Metal_data/Johnson_Au.csv";

// -----------------------------------------------------------------------------
// Scene Configuration
// -----------------------------------------------------------------------------
const int IMAGE_WIDTH = 800;
const int IMAGE_HEIGHT = 450;      // 16:9 Aspect Ratio
const int SAMPLES_PER_PIXEL = 1024; // Higher = less noise, slower
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
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // req modern
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
    std::cout << "[Render] mode = " << rayt::kRenderModeName << "  (switch: Core/RenderMode.hpp)\n";

    // --------------------------------------------------
    // 1. UI and OpenGL texture
    // --------------------------------------------------
    m_uiLayer.init(m_window);

    glGenTextures(1, &m_renderTexture);
    glBindTexture(GL_TEXTURE_2D, m_renderTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, m_width, m_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glBindTexture(GL_TEXTURE_2D, 0);

    // --------------------------------------------------
    // 2. Materials, environment maps, scene
    // --------------------------------------------------
    std::cout << "CWD = " << std::filesystem::current_path() << std::endl;
    loadMaterials();
    buildScene();

    // --------------------------------------------------
    // 3. Camera, film, integrator
    // --------------------------------------------------
    Point3 lookFrom(0, 0.5, 2.5);
    Point3 lookAt(0, 0, -1);
    Vector3 vUp(0, 1, 0);

    Real distToFocus = glm::length(lookFrom - lookAt);
    Real aperture = 0.0; // pinhole
    Real vfov = 35.0;
    Real aspect = static_cast<Real>(m_width) / static_cast<Real>(m_height);

    m_camera = std::make_shared<Camera>(
        lookFrom, lookAt, vUp,
        vfov, aspect, aperture, distToFocus
    );

    m_film = std::make_unique<Film>(m_width, m_height);

    // PathIntegrator(camera, maxDepth, spp)
    m_integrator = std::make_unique<PathIntegrator>(
        m_camera,
        m_settings.maxDepth,
        m_settings.targetSpp
    );

    std::cout << "[Integrator] PathIntegrator created (maxDepth="
        << m_settings.maxDepth << ", spp=" << m_settings.targetSpp << ").\n";
}

void App::loadMaterials() {
    // Environment maps: every .hdr in assets/env (index 0 = none). The default is the one named by
    // ENV_HDR_PATH when it exists (it is not tracked by Git), otherwise the first one found.
    m_envPaths.assign(1, "");
    m_settings.envNames.assign(1, "None (black)");
    {
        std::error_code ec;
        std::vector<std::filesystem::path> found;
        for (const auto& entry : std::filesystem::directory_iterator("assets/env", ec)) {
            if (entry.is_regular_file() && entry.path().extension() == ".hdr") found.push_back(entry.path());
        }
        std::sort(found.begin(), found.end());
        for (const auto& path : found) {
            m_envPaths.push_back(path.generic_string());
            m_settings.envNames.push_back(path.filename().string());
        }
        m_settings.env = 0;
        for (size_t i = 1; i < m_envPaths.size(); ++i) {
            if (m_envPaths[i] == ENV_HDR_PATH) m_settings.env = static_cast<int>(i);
        }
        if (m_settings.env == 0 && m_envPaths.size() > 1) m_settings.env = 1;
        std::cout << "[EnvMap] " << (m_envPaths.size() - 1) << " map(s) found in assets/env\n";
    }

    // Gold optical constants (Au) from the CSV. The table is handed to the material as-is: spectral mode
    // evaluates it per wavelength, RGB mode reduces it to RGB eta/k by CMF integration.
    m_auTable = std::make_shared<rayt::color::ComplexIorTable>();
    std::string auErr;
    if (!rayt::io::loadComplexIor_RefractiveIndexInfoCsv(AU_IOR_CSV_PATH, *m_auTable, &auErr)) {
        std::cerr << "[Au IOR] Failed: " << auErr << "\n";
        std::cerr << "[Au IOR] Fallback to hard-coded eta/k.\n";

        // Fallback: a tiny table built from the old hard-coded RGB constants (450/550/650 nm).
        m_auTable->assign({
            { 450.0, 1.45, 1.77 },
            { 550.0, 0.42, 2.45 },
            { 650.0, 0.16, 3.48 },
        });
    }
    else {
        std::cout << "[Au IOR] Loaded: " << AU_IOR_CSV_PATH << "\n";
    }

    // Log the RGB reduction of the table (what RGB mode renders with).
    {
        Vector3 eta_rgb(1.0), k_rgb(0.0);
        rayt::color::integrateIorToRgb_CMF(*m_auTable, eta_rgb, k_rgb);
        std::cout << "[Au IOR] RGB (CMF-integrated) eta=("
            << eta_rgb.x << "," << eta_rgb.y << "," << eta_rgb.z << ") "
            << "k=(" << k_rgb.x << "," << k_rgb.y << "," << k_rgb.z << ")\n";
        std::cout << "[Au IOR] F0 (normal incidence) = ("
            << rayt::color::conductorNormalIncidenceReflectance(eta_rgb.x, k_rgb.x) << ","
            << rayt::color::conductorNormalIncidenceReflectance(eta_rgb.y, k_rgb.y) << ","
            << rayt::color::conductorNormalIncidenceReflectance(eta_rgb.z, k_rgb.z) << ")\n";
    }

    // Three gold spheres: almost a mirror (0.01) / the one the roughness slider controls / matte (0.50).
    m_goldSmooth = std::make_shared<RoughConductor>(m_auTable, 0.01);
    m_goldCenter = std::make_shared<RoughConductor>(m_auTable, Real(m_settings.roughness));
    m_goldRough = std::make_shared<RoughConductor>(m_auTable, 0.50);
    m_appliedRoughness = m_settings.roughness;

    m_glass = std::make_shared<Dielectric>(6.0, 0.01);     // roughness ~0 = clear
    m_frosted = std::make_shared<Dielectric>(1.2, 0.2);    // roughness 0.2 = frosted
    m_floor = std::make_shared<Lambertian>(RgbColor(0.5, 0.5, 0.5));
}

std::shared_ptr<EnvLight> App::loadEnv(int index) {
    if (index <= 0 || index >= static_cast<int>(m_envPaths.size())) return nullptr;
    if (const auto it = m_envCache.find(index); it != m_envCache.end()) return it->second;

    std::shared_ptr<EnvLight> light;
    try {
        std::cout << "[EnvMap] Loading: " << m_envPaths[static_cast<size_t>(index)] << std::endl;
        auto envImg = io::loadHDR(m_envPaths[static_cast<size_t>(index)]);
        auto envMap = std::make_shared<EnvMap>(std::move(envImg));
        light = std::make_shared<EnvLight>(envMap, 1.0);
        std::cout << "[EnvMap] Success.\n";
    }
    catch (const std::exception& e) {
        std::cerr << "[EnvMap] Failed: " << e.what() << "\n";
        std::cerr << "[EnvMap] Fallback to black background.\n";
    }
    m_envCache[index] = light;   // a failed load is cached as well, so it is not retried every time
    return light;
}

void App::buildScene() {
    auto worldList = std::make_shared<Aggregate>();

    if (m_settings.scene == 0) {
        worldList->add(std::make_shared<Sphere>(Point3(-1.2, 0, -1), 0.5, m_goldSmooth));
        worldList->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), 0.5, m_goldCenter));
        worldList->add(std::make_shared<Sphere>(Point3(1.2, 0, -1), 0.5, m_goldRough));
    }
    else {
        worldList->add(std::make_shared<Sphere>(Point3(0.7, 0, -1), 0.5, m_glass));
        worldList->add(std::make_shared<Sphere>(Point3(-0.7, 0, -1), 0.5, m_frosted));
    }

    // A big sphere is the floor.
    worldList->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, m_floor));

    m_scene = std::make_shared<Scene>(worldList);

    if (auto envLight = loadEnv(m_settings.env)) {
        m_scene->setEnvLight(envLight);
        std::cout << "[Scene] Environment light set.\n";
    }

    // IMPORTANT: finalize builds the LightSampler (default: spatial importance, PBRT-style).
    m_scene->finalize();
    std::cout << "[Scene] " << (m_settings.scene == 0 ? "gold spheres" : "glass spheres") << ", finalized with "
        << m_scene->numLights() << " scene lights";
    if (m_scene->hasEnvLight()) std::cout << " + environment light";
    std::cout << ".\n";
}

void App::updateTexture() {
    const auto* pixels = m_film->getData();
    float scale = 1.0f / std::max(1, m_currentSpp);

    static std::vector<unsigned char> displayBuffer;
    if (displayBuffer.size() != static_cast<size_t>(m_width) * static_cast<size_t>(m_height) * 4)
        displayBuffer.resize(static_cast<size_t>(m_width) * static_cast<size_t>(m_height) * 4);

    const Real exposureStops = Real(m_settings.exposure);   // from the UI slider
    m_shownExposure = m_settings.exposure;

    #pragma omp parallel for
    for (int i = 0; i < m_width * m_height; ++i) {
        rayt::Vector3 s = pixels[i] * static_cast<rayt::Real>(scale);
        // Exposure -> Reinhard tone mapping -> sRGB encoding
        s = rayt::color::toDisplaySRGB_Reinhard(s, exposureStops);

        displayBuffer[i * 4 + 0] = (unsigned char)(255.0f * rayt::math::saturate((float)s[0]));
        displayBuffer[i * 4 + 1] = (unsigned char)(255.0f * rayt::math::saturate((float)s[1]));
        displayBuffer[i * 4 + 2] = (unsigned char)(255.0f * rayt::math::saturate((float)s[2]));
        displayBuffer[i * 4 + 3] = 255;
    }

    glBindTexture(GL_TEXTURE_2D, m_renderTexture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, m_width, m_height, GL_RGBA, GL_UNSIGNED_BYTE, displayBuffer.data());
}

void App::saveImage() const {
    const int spp = std::max(1, m_currentSpp);
    const std::time_t now = std::time(nullptr);
    std::tm local{};
    localtime_s(&local, &now);
    char stamp[32];
    std::strftime(stamp, sizeof(stamp), "%Y%m%d_%H%M%S", &local);
    const std::string base = std::string("gui_") + stamp + "_" + std::to_string(spp) + "spp";

    const float average = 1.0f / static_cast<float>(spp);
    m_film->save(base + ".png", average * std::pow(2.0f, m_settings.exposure));   // exposure in stops = a factor 2^EV on the radiance
    m_film->save(base + ".hdr", average);                                         // linear, not tone mapped
    std::cout << "[GUI] saved " << base << ".png / .hdr\n";
}

// Main loop
void App::run() {
    while (!glfwWindowShouldClose(m_window)) {
        frame();
    }
}

// One iteration of the main loop: apply the UI requests, render the passes, update the texture, draw the UI.
void App::frame() {
    glfwPollEvents();

    // The scene or the environment map changed: rebuild (the accumulation restarts).
    if (m_settings.requestRebuild) {
        buildScene();
        m_settings.requestRebuild = false;
        m_settings.requestReset = true;
    }

    // The center gold sphere's roughness changed. This happens between two render passes, so no ray is in flight.
    if (m_settings.roughness != m_appliedRoughness) {
        m_goldCenter->setRoughness(Real(m_settings.roughness));
        m_appliedRoughness = m_settings.roughness;
        m_settings.requestReset = true;
    }

    bool textureDirty = false;
    if (m_settings.requestReset) {
        m_film->clear();
        m_currentSpp = 0;
        m_settings.requestReset = false;
        textureDirty = true;
    }

    if (!m_settings.paused) {
        for (int i = 0; i < m_settings.sppPerFrame && m_currentSpp < m_settings.targetSpp; ++i) {
            m_integrator->renderOnePass(*m_scene, *m_film, m_currentSpp);
            ++m_currentSpp;
            textureDirty = true;
        }
    }

    // A new exposure only needs the tone mapping again, not new samples (also while paused or finished).
    if (m_settings.exposure != m_shownExposure) textureDirty = true;
    if (textureDirty) updateTexture();

    if (m_settings.requestSave) {
        saveImage();
        m_settings.requestSave = false;
    }
    m_settings.currentSpp = m_currentSpp;

    m_uiLayer.beginFrame();
    m_uiLayer.draw(m_settings);

    {
        ImGui::Begin("Viewport");
        // Keep the aspect ratio of the image and center it in the window.
        const ImVec2 avail = ImGui::GetContentRegionAvail();
        const float aspect = static_cast<float>(m_width) / static_cast<float>(m_height);
        ImVec2 size = avail;
        if (size.x / aspect > size.y) size.x = size.y * aspect; else size.y = size.x / aspect;
        if (size.x > 1.0f && size.y > 1.0f) {
            ImGui::SetCursorPos(ImVec2(ImGui::GetCursorPosX() + (avail.x - size.x) * 0.5f,
                                       ImGui::GetCursorPosY() + (avail.y - size.y) * 0.5f));
            ImGui::Image((void*)(intptr_t)m_renderTexture, size, ImVec2(0, 0), ImVec2(1, 1));
        }
        ImGui::End();
    }

    int w, h;
    glfwGetFramebufferSize(m_window, &w, &h);
    glViewport(0, 0, w, h);
    glClear(GL_COLOR_BUFFER_BIT);

    m_uiLayer.endFrame();
    glfwSwapBuffers(m_window);
}

bool App::selfTest(bool keepSavedFiles) {
    bool allPassed = true;
    auto check = [&](const std::string& what, bool ok, const std::string& detail = "") {
        std::cout << "[SelfTest] " << what << (detail.empty() ? "" : ": " + detail) << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };
    auto anyPixel = [&]() {
        const Vector3* d = m_film->getData();
        for (int i = 0; i < m_width * m_height; ++i) if (d[i].x > 0 || d[i].y > 0 || d[i].z > 0) return true;
        return false;
    };

    // Every frame() renders sppPerFrame passes of the current scene.
    m_settings.sppPerFrame = 1;
    for (int i = 0; i < 3; ++i) frame();
    check("3 frames of 1 pass accumulate 3 spp", m_currentSpp == 3, std::to_string(m_currentSpp) + " spp");
    check("the film has light in it", anyPixel());

    m_settings.sppPerFrame = 4;
    frame();
    check("SPP / frame = 4 adds 4 passes", m_currentSpp == 7, std::to_string(m_currentSpp) + " spp");

    m_settings.exposure = 1.0f;
    frame();
    check("a new exposure does not restart the accumulation", m_currentSpp == 11 && m_shownExposure == 1.0f,
        std::to_string(m_currentSpp) + " spp, shown exposure " + std::to_string(m_shownExposure));

    m_settings.roughness = 0.6f;
    frame();
    check("a new roughness restarts the accumulation", m_currentSpp == 4 && m_appliedRoughness == 0.6f,
        std::to_string(m_currentSpp) + " spp");

    const auto glassScene = m_scene;
    m_settings.scene = 1;
    m_settings.requestRebuild = true;
    frame();
    check("the glass scene replaces the scene and restarts", m_scene != glassScene && m_currentSpp == 4, std::to_string(m_currentSpp) + " spp");

    const int envIndex = m_settings.env;
    m_settings.env = 0;
    m_settings.requestRebuild = true;
    frame();
    check("environment 'None' removes the environment light", !m_scene->hasEnvLight());
    if (m_envPaths.size() > 1) {
        m_settings.env = std::max(1, envIndex);
        m_settings.requestRebuild = true;
        frame();
        check("an environment map adds the environment light", m_scene->hasEnvLight(), m_settings.envNames[static_cast<size_t>(m_settings.env)]);
    }
    m_settings.scene = 0;
    m_settings.requestRebuild = true;
    frame();

    // Save: gui_<time>_<spp>spp.png and .hdr in the working directory.
    std::error_code ec;
    std::vector<std::filesystem::path> before;
    for (const auto& e : std::filesystem::directory_iterator(".", ec)) before.push_back(e.path());
    m_settings.requestSave = true;
    frame();
    std::vector<std::filesystem::path> created;
    for (const auto& e : std::filesystem::directory_iterator(".", ec))
        if (std::find(before.begin(), before.end(), e.path()) == before.end()) created.push_back(e.path());
    bool png = false, hdr = false;
    for (const auto& f : created) {
        png = png || (f.extension() == ".png" && f.filename().string().rfind("gui_", 0) == 0 && std::filesystem::file_size(f) > 0);
        hdr = hdr || (f.extension() == ".hdr" && f.filename().string().rfind("gui_", 0) == 0 && std::filesystem::file_size(f) > 0);
    }
    check("'Save image' writes a PNG and an HDR", png && hdr);
    for (const auto& f : created) {
        if (keepSavedFiles) std::cout << "[SelfTest] kept " << f.filename().string() << "\n";
        else std::filesystem::remove(f, ec);
    }

    // Pause and the sample limit.
    const int spp = m_currentSpp;
    m_settings.paused = true;
    frame();
    check("pause stops the accumulation", m_currentSpp == spp);
    m_settings.paused = false;
    m_settings.targetSpp = spp;
    frame();
    check("the sample limit stops the accumulation", m_currentSpp == spp);

    std::cout << "[SelfTest] RESULT: " << (allPassed ? "PASS" : "FAIL") << std::endl;
    return allPassed;
}
