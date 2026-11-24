/**
 * @file ImGuiLayer.cpp
 * @brief Implementation of the ImGuiLayer class.
 *
 * This file handles the initialization, rendering loop, and resource management
 * for the Dear ImGui interface. It includes support for High-DPI scaling,
 * Docking architecture, and OpenGL/GLFW backend integration.
 */

#include "pch.h"

#include "ImGuiLayer.hpp"
#include "Core/RenderSettings.hpp"
#include "Core/RenderMode.hpp"

#include <algorithm>
#include <string>

#include <GLFW/glfw3.h>
#include <filesystem>
#include <imgui.h>
#include <imgui_internal.h>   // DockBuilder: the default arrangement of the windows
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>

 /**
  * @brief Initializes the ImGui context and backends.
  *
  * Sets up the docking configuration, applies the dark style, initializes
  * GLFW and OpenGL3 backends, and calculates the initial High-DPI scale factor.
  *
  * @param window Pointer to the GLFW window handle.
  */
void ImGuiLayer::init(GLFWwindow* window) {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;   // Enable Docking
    // io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable; // Enable Multi-Viewport / Platform Windows
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableSetMousePos;

    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    // The monitor's scale (Windows' "scale and layout": 1 at 100 %, 1.5 at 150 % ...)
    float xs = 1.0f, ys = 1.0f;
    glfwGetWindowContentScale(window, &xs, &ys);
    m_dpiScale = std::max(1.0f, (xs + ys) * 0.5f);

    // A real font, rasterized at the size it is shown (the built-in one is a small bitmap font that blurs when scaled):
    // Meiryo (Latin, kana, kanji) with Microsoft YaHei merged in for the simplified Chinese of model and material names.
    m_fontScale = m_dpiScale;
    const float px = 17.0f * m_fontScale;
    bool haveFont = false;
    if (std::filesystem::exists("C:/Windows/Fonts/meiryo.ttc")) {
        haveFont = io.Fonts->AddFontFromFileTTF("C:/Windows/Fonts/meiryo.ttc", px, nullptr, io.Fonts->GetGlyphRangesJapanese()) != nullptr;
        if (haveFont && std::filesystem::exists("C:/Windows/Fonts/msyh.ttc")) {
            ImFontConfig merge;
            merge.MergeMode = true;
            io.Fonts->AddFontFromFileTTF("C:/Windows/Fonts/msyh.ttc", px, &merge, io.Fonts->GetGlyphRangesChineseSimplifiedCommon());
        }
    }
    if (!haveFont) {
        ImFontConfig cfg;
        cfg.SizePixels = 13.0f * m_fontScale;
        io.Fonts->AddFontDefault(&cfg);
    }
    io.Fonts->Build();

    m_baseStyle = new ImGuiStyle(ImGui::GetStyle());
    applyUiScale();
}

void ImGuiLayer::applyUiScale() {
    ImGuiStyle& style = ImGui::GetStyle();
    style = *m_baseStyle;
    style.ScaleAllSizes(m_dpiScale * m_userScale);
    ImGui::GetIO().FontGlobalScale = m_dpiScale * m_userScale / m_fontScale;
}

/**
 * @brief Callback for handling window content scale changes (DPI changes).
 *
 * This is usually called by the GLFW callback when the window is moved
 * to a monitor with a different DPI.
 *
 * @param xscale The new x-axis scale factor.
 * @param yscale The new y-axis scale factor.
 */
void ImGuiLayer::onContentScaleChanged(float xscale, float yscale) {
    m_dpiScale = std::max(1.0f, (xscale + yscale) * 0.5f);
    applyUiScale();
}

/**
 * @brief Starts a new ImGui frame.
 *
 * Must be called before issuing any ImGui drawing commands.
 */
void ImGuiLayer::beginFrame() {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

/**
 * @brief Renders the UI layout.
 *
 * This function sets up the main DockSpace that covers the entire viewport,
 * allowing other windows to dock into it. It also renders the specific
 * application windows (e.g., Render Settings).
 *
 * @param s Reference to the RenderSettings object to manipulate via UI.
 */
void ImGuiLayer::draw(RenderSettings& s) {
    // 1. Setup DockSpace (Boilerplate code)
    // Create a "Background Window" that covers the entire screen
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(viewport->WorkPos);
    ImGui::SetNextWindowSize(viewport->WorkSize);
    ImGui::SetNextWindowViewport(viewport->ID);

    // Style settings (no rounding, no borders)
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);

    // Zero padding to ensure the window fits perfectly
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));

    // Configure flags for the background window
    ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoDocking;
    window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
    window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;
    // Note: use NoBackground if you want transparency, but since we have a DockSpace, default is fine.

    // Begin "Background"
    ImGui::Begin("MainDockSpace", nullptr, window_flags);

    // Restore style (subsequent windows will use standard style)
    ImGui::PopStyleVar(3);

    // Actual DockSpace ID
    ImGuiID dockspace_id = ImGui::GetID("MyDockSpace");
    // The arrangement: the settings on the left, the viewport on the right (at start, and on "Reset layout"; a layout saved
    // in imgui.ini from another screen could leave windows larger than the window or outside it)
    if (m_resetLayout || s.requestLayoutReset) {
        ImGui::DockBuilderRemoveNode(dockspace_id);
        ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
        ImGui::DockBuilderSetNodeSize(dockspace_id, viewport->WorkSize);
        ImGuiID left = 0, right = 0;
        ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 0.30f, &left, &right);
        ImGui::DockBuilderDockWindow("Render Settings", left);
        ImGui::DockBuilderDockWindow("Viewport", right);
        ImGui::DockBuilderFinish(dockspace_id);
        m_resetLayout = false;
        s.requestLayoutReset = false;
    }
    ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), ImGuiDockNodeFlags_None);

    // 2. Draw Content Windows

    // Render Settings Window
    ImGui::Begin("Render Settings"); 

    ImGui::Text("Application Average: %.3f ms/frame (%.1f FPS)",
        1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
    ImGui::Text("Render mode: %s  (compile-time switch: Core/RenderMode.hpp)", rayt::kRenderModeName);
    ImGui::SliderFloat("UI scale", &s.uiScale, 0.6f, 2.5f, "x %.2f");
    if (ImGui::IsItemDeactivatedAfterEdit() && s.uiScale != m_userScale) { m_userScale = s.uiScale; applyUiScale(); }
    ImGui::SameLine();
    if (ImGui::Button("Reset layout")) s.requestLayoutReset = true;
    ImGui::Separator();

    if (ImGui::Combo("View", &s.mode, "Path tracing (spheres)\0MMD preview (OpenGL, toon)\0") && s.mode == 1 && !s.pmxNames.empty())
        s.requestModelLoad = s.requestModelLoad || s.previewStatus.empty();   // the first switch loads the model
    ImGui::Separator();

    if (s.mode == 1) {
        // --- MMD preview ---
        if (s.pmxNames.empty()) ImGui::TextWrapped("No .pmx file was found in MMD_model/.");
        else {
            auto combo = [&](const char* label, std::vector<std::string>& names, int& index) {
                bool changed = false;
                if (ImGui::BeginCombo(label, names[static_cast<size_t>(index)].c_str())) {
                    for (int i = 0; i < static_cast<int>(names.size()); ++i) {
                        const bool selected = (i == index);
                        if (ImGui::Selectable(names[static_cast<size_t>(i)].c_str(), selected) && !selected) { index = i; changed = true; }
                        if (selected) ImGui::SetItemDefaultFocus();
                    }
                    ImGui::EndCombo();
                }
                return changed;
            };
            if (combo("Model", s.pmxNames, s.pmx)) { s.requestModelLoad = true; s.requestCameraReset = true; }
            if (combo("Pose", s.poseNames, s.pose)) s.requestModelLoad = true;
        }
        ImGui::SliderFloat("Light azimuth", &s.lightAzimuth, -180.0f, 180.0f, "%.0f deg");
        ImGui::SliderFloat("Light elevation", &s.lightElevation, 0.0f, 90.0f, "%.0f deg");
        ImGui::Checkbox("Shadows", &s.previewShadows);
        ImGui::SameLine();
        ImGui::Checkbox("Outlines", &s.previewEdges);
        ImGui::SameLine();
        ImGui::Checkbox("Floor", &s.previewFloor);
        ImGui::SliderFloat("Outline width", &s.previewEdgeScale, 0.0f, 4.0f, "x %.1f");
        if (ImGui::Button("Frame the model")) s.requestCameraReset = true;
        ImGui::Separator();

        // --- posing ---
        ImGui::Checkbox("Show bones", &s.showBones);
        ImGui::Text("Bone: %s", s.selectedBone.empty() ? "(click a dot in the viewport)" : s.selectedBone.c_str());
        if (ImGui::Button("Undo (Ctrl+Z)")) s.requestUndo = true;
        ImGui::SameLine();
        if (ImGui::Button("Reset bone")) s.requestResetBone = true;
        ImGui::SameLine();
        if (ImGui::Button("Reset pose")) s.requestResetPose = true;
        if (ImGui::Button("Save pose (.vpd)")) s.requestSavePose = true;
        if (!s.poseStatus.empty()) ImGui::TextWrapped("%s", s.poseStatus.c_str());
        ImGui::Separator();
        ImGui::TextWrapped("Viewport: drag a bone's dot to turn it in the picture plane (Shift: tilt it toward or away from you); "
            "drag elsewhere to turn the camera, right or middle drag moves it, the wheel zooms.");
        if (!s.previewStatus.empty()) ImGui::TextWrapped("%s", s.previewStatus.c_str());
        ImGui::End(); // End Render Settings
        ImGui::End(); // End MainDockSpace
        return;
    }

    // --- Scene ---
    if (ImGui::Combo("Scene", &s.scene, "Gold spheres\0Glass spheres\0")) s.requestRebuild = true;
    if (ImGui::BeginCombo("Environment", s.envNames[static_cast<size_t>(s.env)].c_str())) {
        for (int i = 0; i < static_cast<int>(s.envNames.size()); ++i) {
            const bool selected = (i == s.env);
            if (ImGui::Selectable(s.envNames[static_cast<size_t>(i)].c_str(), selected) && !selected) {
                s.env = i;
                s.requestRebuild = true;
            }
            if (selected) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
    if (s.scene == 0) {
        // The material is changed between two render passes by the application loop (it also restarts the accumulation).
        ImGui::SliderFloat("Center gold roughness", &s.roughness, 0.01f, 1.0f, "%.2f");
    }
    ImGui::Separator();

    // --- Rendering ---
    ImGui::Checkbox("Paused", &s.paused);
    ImGui::SliderInt("SPP / frame", &s.sppPerFrame, 1, 64);
    ImGui::InputInt("Max SPP", &s.targetSpp, 64, 1024);
    s.targetSpp = std::max(1, s.targetSpp);
    ImGui::ProgressBar(std::min(1.0f, static_cast<float>(s.currentSpp) / static_cast<float>(s.targetSpp)), ImVec2(-1.0f, 0.0f),
        (std::to_string(s.currentSpp) + " / " + std::to_string(s.targetSpp) + " spp").c_str());

    if (ImGui::Button("Reset accumulation")) {
        s.requestReset = true;
    }
    ImGui::Separator();

    // --- Display and output ---
    ImGui::SliderFloat("Exposure (EV)", &s.exposure, -6.0f, 6.0f, "%.1f");
    if (ImGui::Button("Save image")) {
        s.requestSave = true;   // writes gui_<time>_<spp>spp.png (with the exposure above) and .hdr (linear, raw)
    }

    ImGui::End(); // End Render Settings


    // (Future: Add Viewport windows here)


    // 3. End DockSpace

    ImGui::End(); // End MainDockSpace
}

/**
 * @brief Finalizes the frame and renders the draw data to the OpenGL context.
 */
void ImGuiLayer::endFrame() {
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

/**
 * @brief Shuts down ImGui and releases backend resources.
 */
void ImGuiLayer::shutdown() {
    delete m_baseStyle;
    m_baseStyle = nullptr;
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}