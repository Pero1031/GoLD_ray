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

// MMD files (lists for the preview)
#include "Cli/PmxSupport.hpp"
#include "IO/ImageWriter.hpp"

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
    m_preview.reset();   // its OpenGL objects go while the context exists
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
    findMmdFiles();

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

    // Gold optical constants (Au): the table is handed to the materials as it is (spectral mode evaluates it per
    // wavelength, RGB mode reduces it to RGB eta/k).
    m_auTable = rayt::io::loadGoldIorTable(AU_IOR_CSV_PATH);

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

    // MMD preview: load what the settings ask for; the path tracer rests meanwhile
    const bool preview = (m_settings.mode == 1);
    if (preview && m_settings.requestModelLoad) {
        loadPreviewModel();
        m_settings.requestModelLoad = false;
    }
    if (preview) handlePoseRequests();
    if (preview && m_settings.requestCameraReset && m_preview && m_preview->loaded()) {
        m_orbit = m_preview->framingCamera();
        m_settings.requestCameraReset = false;
    }

    bool textureDirty = false;
    if (m_settings.requestReset) {
        m_film->clear();
        m_currentSpp = 0;
        m_settings.requestReset = false;
        textureDirty = true;
    }

    if (!m_settings.paused && !preview) {
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

    drawViewport();

    int w, h;
    glfwGetFramebufferSize(m_window, &w, &h);
    glViewport(0, 0, w, h);
    glClear(GL_COLOR_BUFFER_BIT);

    m_uiLayer.endFrame();
    glfwSwapBuffers(m_window);
}

void App::findMmdFiles() {
    const std::filesystem::path current = (m_settings.pmx >= 0 && m_settings.pmx < static_cast<int>(m_pmxFiles.size())) ? m_pmxFiles[static_cast<size_t>(m_settings.pmx)] : std::filesystem::path();
    m_pmxFiles = cli::FindPmxFiles();
    m_settings.pmxNames.clear();
    size_t largest = 0;
    for (size_t i = 0; i < m_pmxFiles.size(); ++i) {
        m_settings.pmxNames.push_back(cli::PathToUtf8(m_pmxFiles[i].parent_path().filename()) + " / " + cli::PathToUtf8(m_pmxFiles[i].filename()));
        std::error_code ec;
        if (std::filesystem::file_size(m_pmxFiles[i], ec) > std::filesystem::file_size(m_pmxFiles[largest], ec)) largest = i;
    }
    m_settings.pmx = static_cast<int>(largest);   // the largest file, as the ray tracer's PMX scene picks by default
    for (size_t i = 0; i < m_pmxFiles.size(); ++i) if (!current.empty() && m_pmxFiles[i] == current) m_settings.pmx = static_cast<int>(i);   // a rescan keeps the model
    m_poseFiles.assign(1, std::filesystem::path());
    m_settings.poseNames.assign(1, "(rest pose)");
    for (const auto& f : cli::FindMmdFiles("pose", ".vpd")) {
        m_poseFiles.push_back(f);
        m_settings.poseNames.push_back(cli::PathToUtf8(f.parent_path().filename()) + " / " + cli::PathToUtf8(f.filename()));
    }
    std::cout << "[Preview] " << m_pmxFiles.size() << " model(s) and " << (m_poseFiles.size() - 1) << " pose file(s) in MMD_model\n";
}

void App::loadPreviewModel() {
    if (m_pmxFiles.empty()) return;
    if (!m_preview) m_preview = std::make_unique<preview::MmdPreview>();
    const bool first = !m_preview->loaded();
    const auto& pmx = m_pmxFiles[static_cast<size_t>(std::clamp(m_settings.pmx, 0, static_cast<int>(m_pmxFiles.size()) - 1))];
    const auto& pose = m_poseFiles[static_cast<size_t>(std::clamp(m_settings.pose, 0, static_cast<int>(m_poseFiles.size()) - 1))];
    std::string error;
    m_selectedBone = m_dragBone = -1;
    m_settings.selectedBone.clear();
    if (m_preview->load(pmx, pose, &error)) {
        char buf[256];
        std::snprintf(buf, sizeof(buf), "Loaded in %.0f ms (read, pose, convert, upload); height %.2f m.", m_preview->loadMilliseconds(), double(m_preview->modelHeight()));
        m_settings.previewStatus = buf;
        if (first) m_settings.requestCameraReset = true;
    }
    else {
        m_settings.previewStatus = "Could not load: " + error;
        std::cerr << "[Preview] " << m_settings.previewStatus << "\n";
    }
}

void App::handlePoseRequests() {
    RenderSettings& s = m_settings;
    if (!m_preview || !m_preview->loaded()) { s.requestUndo = s.requestResetBone = s.requestResetPose = s.requestSavePose = false; return; }
    if (s.requestUndo) {
        s.poseStatus = m_preview->undo() ? "Undone." : "Nothing to undo.";
        s.requestUndo = false;
    }
    if (s.requestResetBone) {
        if (m_selectedBone >= 0) { m_preview->pushUndo(); m_preview->resetBone(m_selectedBone); s.poseStatus = "Bone reset: " + s.selectedBone; }
        s.requestResetBone = false;
    }
    if (s.requestResetPose) {
        m_preview->pushUndo();
        m_preview->resetPose();
        s.poseStatus = "Every bone at rest.";
        s.requestResetPose = false;
    }
    if (s.requestSavePose) {
        const std::time_t now = std::time(nullptr);
        std::tm local{};
        localtime_s(&local, &now);
        char stamp[32];
        std::strftime(stamp, sizeof(stamp), "%Y%m%d_%H%M%S", &local);
        // MMD_model/pose/gui, next to the models
        std::filesystem::path mmd = "../MMD_model";
        for (std::filesystem::path d = m_pmxFiles.empty() ? std::filesystem::path() : m_pmxFiles.front().parent_path(); !d.empty() && d.has_relative_path(); d = d.parent_path()) {
            if (d.filename() == "MMD_model") { mmd = d; break; }
            if (d == d.parent_path()) break;
        }
        const std::filesystem::path file = mmd / "pose" / "gui" / (std::string("pose_") + stamp + ".vpd");
        std::string error;
        if (m_preview->savePose(file, &error)) {
            // the new file joins the list, selected (without loading it again: it is what is shown)
            findMmdFiles();
            std::error_code ec;
            for (size_t i = 1; i < m_poseFiles.size(); ++i)
                if (std::filesystem::equivalent(m_poseFiles[i], file, ec)) s.pose = static_cast<int>(i);
            s.poseStatus = "Saved MMD_model/pose/gui/" + file.filename().string() + " (ray tracer: --scene pmx --pose-index " + std::to_string(s.pose - 1) + ").";
        }
        else s.poseStatus = "Could not save: " + error;
        s.requestSavePose = false;
    }
}

void App::drawViewport() {
    ImGui::Begin("Viewport");
    const ImVec2 avail = ImGui::GetContentRegionAvail();
    if (m_settings.mode == 1 && m_preview) {
        // The preview fills the window at the framebuffer's resolution; the mouse drives the orbit camera.
        const ImVec2 scale = ImGui::GetIO().DisplayFramebufferScale;
        const int w = std::max(1, static_cast<int>(avail.x * scale.x)), h = std::max(1, static_cast<int>(avail.y * scale.y));
        if (avail.x > 1.0f && avail.y > 1.0f) {
            preview::PreviewOptions o;
            o.lightAzimuth = m_settings.lightAzimuth;
            o.lightElevation = m_settings.lightElevation;
            o.shadows = m_settings.previewShadows;
            o.edges = m_settings.previewEdges;
            o.floor = m_settings.previewFloor;
            o.edgeScale = m_settings.previewEdgeScale;
            m_preview->render(w, h, m_orbit, o);
            ImGui::Image((void*)(intptr_t)m_preview->texture(), avail, ImVec2(0, 1), ImVec2(1, 0));   // OpenGL rows run from the bottom
            const ImVec2 r0 = ImGui::GetItemRectMin(), r1 = ImGui::GetItemRectMax();
            const bool hovered = ImGui::IsItemHovered();
            const ImGuiIO& io = ImGui::GetIO();

            // ---- the bones: dots on the joints, lines to the parents; a click on a dot selects it, a drag turns it ----
            const glm::dmat4 vp = preview::MmdPreview::viewProjection(w, h, m_orbit);
            const auto& bones = m_preview->bones();
            std::vector<ImVec2> screen(bones.size());
            std::vector<char> onScreen(bones.size(), 0);
            for (size_t i = 0; i < bones.size(); ++i) {
                const glm::dvec4 c = vp * glm::dvec4(bones[i].position, 1.0);
                if (c.w <= 1e-6) continue;
                screen[i] = ImVec2(r0.x + float((c.x / c.w * 0.5 + 0.5) * (r1.x - r0.x)), r0.y + float((0.5 - c.y / c.w * 0.5) * (r1.y - r0.y)));
                onScreen[i] = 1;
            }
            const float dotR = std::max(3.0f, ImGui::GetFontSize() * 0.22f);
            int nearest = -1;
            float nearestD2 = (dotR * 2.2f) * (dotR * 2.2f);
            if (m_settings.showBones && hovered) {
                for (size_t i = 0; i < bones.size(); ++i) {
                    if (!bones[i].handle || !onScreen[i]) continue;
                    const float dx = screen[i].x - io.MousePos.x, dy = screen[i].y - io.MousePos.y;
                    if (dx * dx + dy * dy < nearestD2) { nearestD2 = dx * dx + dy * dy; nearest = static_cast<int>(i); }
                }
            }
            if (hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                m_dragBone = nearest;
                if (nearest >= 0) {
                    m_selectedBone = nearest;
                    m_settings.selectedBone = bones[static_cast<size_t>(nearest)].name;
                    m_preview->pushUndo();
                }
            }
            if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) m_dragBone = -1;
            if (m_dragBone >= 0 && (io.MouseDelta.x != 0.0f || io.MouseDelta.y != 0.0f) && onScreen[static_cast<size_t>(m_dragBone)]) {
                const Vector3 fwd = m_orbit.forward();
                const Vector3 right = glm::normalize(glm::cross(fwd, Vector3(0, 1, 0)));
                const Vector3 up = glm::cross(right, fwd);
                const Vector3 toViewer = -fwd;
                glm::dquat q(1, 0, 0, 0);
                if (io.KeyShift) {
                    // tilt: a trackball about the joint; the side facing the viewer follows the mouse
                    const Vector3 d = right * Real(io.MouseDelta.x) - up * Real(io.MouseDelta.y);
                    const Vector3 axis = glm::cross(toViewer, d);
                    const Real len = glm::length(axis);
                    if (len > Real(1e-9)) q = glm::angleAxis(glm::radians(Real(0.5)) * glm::length(d), axis / len);
                }
                else {
                    // in the picture plane: the angle the mouse sweeps round the joint's dot (screen y runs down)
                    const ImVec2 c = screen[static_cast<size_t>(m_dragBone)];
                    const ImVec2 a(io.MousePos.x - io.MouseDelta.x - c.x, io.MousePos.y - io.MouseDelta.y - c.y), b(io.MousePos.x - c.x, io.MousePos.y - c.y);
                    if (a.x * a.x + a.y * a.y > 16.0f && b.x * b.x + b.y * b.y > 16.0f) {
                        const Real angle = std::atan2(Real(a.x) * Real(-b.y) - Real(-a.y) * Real(b.x), Real(a.x) * Real(b.x) + Real(a.y) * Real(b.y));
                        q = glm::angleAxis(angle, toViewer);
                    }
                }
                m_preview->rotateBone(m_dragBone, q);
                char buf[160];
                std::snprintf(buf, sizeof(buf), "Turning %s: the pose updates in %.1f ms.", bones[static_cast<size_t>(m_dragBone)].name.c_str(), m_preview->poseMilliseconds());
                m_settings.poseStatus = buf;
            }
            if (m_settings.showBones) {
                ImDrawList* dl = ImGui::GetWindowDrawList();
                dl->PushClipRect(r0, r1, true);
                for (size_t i = 0; i < bones.size(); ++i) {
                    if (!bones[i].handle || !onScreen[i]) continue;
                    const int p = bones[i].parent;
                    if (p >= 0 && onScreen[static_cast<size_t>(p)]) dl->AddLine(screen[static_cast<size_t>(p)], screen[i], IM_COL32(40, 60, 120, 160), 1.5f);
                }
                for (size_t i = 0; i < bones.size(); ++i) {
                    if (!bones[i].handle || !onScreen[i]) continue;
                    const bool sel = static_cast<int>(i) == m_selectedBone, hot = static_cast<int>(i) == nearest;
                    const ImU32 fill = sel ? IM_COL32(255, 150, 30, 255) : hot ? IM_COL32(120, 200, 255, 255) : IM_COL32(60, 110, 230, 220);
                    dl->AddCircleFilled(screen[i], sel ? dotR * 1.4f : dotR, fill);
                    dl->AddCircle(screen[i], sel ? dotR * 1.4f : dotR, IM_COL32(255, 255, 255, 200), 0, 1.0f);
                }
                if (nearest >= 0 && m_dragBone < 0) ImGui::SetTooltip("%s", bones[static_cast<size_t>(nearest)].name.c_str());
                dl->PopClipRect();
            }
            if (hovered && io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_Z)) m_settings.requestUndo = true;

            if (hovered && m_dragBone < 0) {
                if (ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
                    m_orbit.yaw -= Real(io.MouseDelta.x) * Real(0.35);
                    m_orbit.pitch = std::clamp(m_orbit.pitch + Real(io.MouseDelta.y) * Real(0.35), Real(-89), Real(89));
                }
                if (ImGui::IsMouseDragging(ImGuiMouseButton_Right) || ImGui::IsMouseDragging(ImGuiMouseButton_Middle)) {
                    const Vector3 fwd = m_orbit.forward();
                    const Vector3 right = glm::normalize(glm::cross(fwd, Vector3(0, 1, 0)));
                    const Vector3 up = glm::cross(right, fwd);
                    const Real k = m_orbit.distance * Real(2) * std::tan(glm::radians(m_orbit.vfov * Real(0.5))) / Real(std::max(1.0f, avail.y));
                    m_orbit.target += (-right * Real(io.MouseDelta.x) + up * Real(io.MouseDelta.y)) * k;
                }
                if (io.MouseWheel != 0.0f) m_orbit.distance = std::max(Real(0.05), m_orbit.distance * std::pow(Real(0.88), Real(io.MouseWheel)));
            }
        }
    }
    else {
        // Keep the aspect ratio of the image and center it in the window.
        const float aspect = static_cast<float>(m_width) / static_cast<float>(m_height);
        ImVec2 size = avail;
        if (size.x / aspect > size.y) size.x = size.y * aspect; else size.y = size.x / aspect;
        if (size.x > 1.0f && size.y > 1.0f) {
            ImGui::SetCursorPos(ImVec2(ImGui::GetCursorPosX() + (avail.x - size.x) * 0.5f,
                                       ImGui::GetCursorPosY() + (avail.y - size.y) * 0.5f));
            ImGui::Image((void*)(intptr_t)m_renderTexture, size, ImVec2(0, 0), ImVec2(1, 1));
        }
    }
    ImGui::End();
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

    // The MMD preview: the first switch loads the default model; the picture differs from the plain background.
    if (!m_pmxFiles.empty()) {
        m_settings.mode = 1;
        m_settings.requestModelLoad = true;
        const auto t0 = std::chrono::steady_clock::now();
        frame();
        check("the MMD preview loads a model", m_preview && m_preview->loaded(), m_settings.previewStatus);
        if (m_preview && m_preview->loaded()) {
            for (int i = 0; i < 3; ++i) frame();
            const int frames = 30;
            const auto t1 = std::chrono::steady_clock::now();
            for (int i = 0; i < frames; ++i) frame();
            const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t1).count() / frames;
            // read the preview's texture: some texels must differ clearly from the background
            GLint w = 0, h = 0;
            glBindTexture(GL_TEXTURE_2D, m_preview->texture());
            glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
            glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);
            std::vector<unsigned char> px(static_cast<size_t>(w) * static_cast<size_t>(h) * 4);
            glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, px.data());
            glBindTexture(GL_TEXTURE_2D, 0);
            const preview::PreviewOptions defaults;
            size_t differ = 0;
            for (size_t i = 0; i + 3 < px.size(); i += 4) {
                const int dr = std::abs(int(px[i]) - int(defaults.background.x * 255)), dg = std::abs(int(px[i + 1]) - int(defaults.background.y * 255));
                if (dr + dg > 40) ++differ;
            }
            check("the preview draws the model", differ > px.size() / 4 / 50, std::to_string(w) + " x " + std::to_string(h) + ", "
                + std::to_string(100.0 * double(differ) / double(std::max<size_t>(1, px.size() / 4))) + " % of the texels differ from the background");
            std::cout << "[SelfTest] preview: " << ms << " ms per frame (incl. ImGui), first frame with loading "
                << std::chrono::duration<double, std::milli>(t1 - t0).count() << " ms\n";
            if (keepSavedFiles) {   // the preview as a PNG (rows flipped: OpenGL starts at the bottom)
                std::vector<unsigned char> rgb(static_cast<size_t>(w) * static_cast<size_t>(h) * 3);
                for (int y = 0; y < h; ++y)
                    for (int x = 0; x < w; ++x)
                        for (int c = 0; c < 3; ++c)
                            rgb[(static_cast<size_t>(y) * w + x) * 3 + c] = px[(static_cast<size_t>(h - 1 - y) * w + x) * 4 + c];
                rayt::io::WriteImage8("gui_preview_selftest.png", w, h, 3, rgb);
                std::cout << "[SelfTest] kept gui_preview_selftest.png\n";
            }

            // Posing: turn the left arm 30 degrees about the axis toward the viewer; its child joint must land where that
            // rotation about the shoulder puts it (the mirrored z of MMD's space and the parent's frame both enter here).
            const auto& bones = m_preview->bones();
            int arm = -1, child = -1;
            for (size_t i = 0; i < bones.size(); ++i) if (bones[i].name == "\xE5\xB7\xA6\xE8\x85\x95") arm = static_cast<int>(i);   // 左腕
            for (size_t i = 0; i < bones.size() && arm >= 0; ++i)
                if (bones[i].parent == arm && glm::length(bones[i].position - bones[static_cast<size_t>(arm)].position) > Real(0.05)) { child = static_cast<int>(i); break; }
            if (arm >= 0 && child >= 0) {
                const Vector3 pivot = bones[static_cast<size_t>(arm)].position, before = bones[static_cast<size_t>(child)].position;
                const glm::dquat q = glm::angleAxis(glm::radians(Real(30)), -m_orbit.forward());
                m_preview->pushUndo();
                m_preview->rotateBone(arm, q);
                const Vector3 expected = pivot + q * (before - pivot);
                const Real err = glm::length(m_preview->bones()[static_cast<size_t>(child)].position - expected);
                check("turning the left arm moves its child joint as the rotation says", err < Real(1e-6),
                    bones[static_cast<size_t>(child)].name + " off by " + std::to_string(err * 1000) + " mm, pose update " + std::to_string(m_preview->poseMilliseconds()) + " ms");
                frame();
                if (keepSavedFiles) {
                    GLint pw = 0, ph = 0;
                    glBindTexture(GL_TEXTURE_2D, m_preview->texture());
                    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &pw);
                    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &ph);
                    std::vector<unsigned char> pp(static_cast<size_t>(pw) * static_cast<size_t>(ph) * 4), rgb(static_cast<size_t>(pw) * static_cast<size_t>(ph) * 3);
                    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, pp.data());
                    glBindTexture(GL_TEXTURE_2D, 0);
                    for (int y = 0; y < ph; ++y)
                        for (int x = 0; x < pw; ++x)
                            for (int c = 0; c < 3; ++c)
                                rgb[(static_cast<size_t>(y) * pw + x) * 3 + c] = pp[(static_cast<size_t>(ph - 1 - y) * pw + x) * 4 + c];
                    rayt::io::WriteImage8("gui_pose_selftest.png", pw, ph, 3, rgb);
                    std::cout << "[SelfTest] kept gui_pose_selftest.png\n";
                }

                // save -> load gives the same pose; undo gives the rest back
                const Vector3 posed = m_preview->bones()[static_cast<size_t>(child)].position;
                const std::filesystem::path vpd = "gui_pose_selftest.vpd";
                std::string err2;
                const bool saved = m_preview->savePose(vpd, &err2);
                const auto& pmx = m_pmxFiles[static_cast<size_t>(m_settings.pmx)];
                const bool reloaded = saved && m_preview->load(pmx, vpd, &err2);
                const Real roundTrip = reloaded ? glm::length(m_preview->bones()[static_cast<size_t>(child)].position - posed) : Real(1);
                check("a saved pose loads back the same", reloaded && roundTrip < Real(1e-5), std::to_string(roundTrip * 1000) + " mm " + err2);
                if (!keepSavedFiles) { std::error_code ec; std::filesystem::remove(vpd, ec); }
                else std::cout << "[SelfTest] kept gui_pose_selftest.vpd\n";
                m_preview->load(pmx, std::filesystem::path(), &err2);
                m_preview->pushUndo();
                m_preview->rotateBone(arm, q);
                m_preview->undo();
                check("undo puts the arm back", glm::length(m_preview->bones()[static_cast<size_t>(child)].position - before) < Real(1e-9));
            }
            else check("the model has a left arm bone with a child joint", false);
        }
        m_settings.mode = 0;
    }

    std::cout << "[SelfTest] RESULT: " << (allPassed ? "PASS" : "FAIL") << std::endl;
    return allPassed;
}
