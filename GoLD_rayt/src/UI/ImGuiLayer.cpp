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

#include <GLFW/glfw3.h>
#include <imgui.h>
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

    // Determine scale from DPI (fallback if unavailable)
    float xs = 1.0f, ys = 1.0f;
    glfwGetWindowContentScale(window, &xs, &ys);
    float scale = (xs + ys) * 0.5f;

    // Fix for environments that return 1.0 even on high-DPI screens:
    // Ensure minimum readability.
    if (scale < 1.25f) scale = 2.0f;  // Preference: could be 1.5f, etc.

    applyUiScale(scale); 
}

/**
 * @brief Applies a global UI scale factor to fonts and style elements.
 *
 * Calculates the ratio between the new scale and the previously applied scale
 * to prevent cumulative scaling errors when resizing or changing monitors.
 *
 * @param newScale The target scale factor (e.g., 2.0 for 200% scaling).
 */
void ImGuiLayer::applyUiScale(float newScale) {
    if (newScale <= 0.0f) return;

    // Calculate ratio to avoid accumulation issues
    float ratio = newScale / m_uiScale;
    m_uiScale = newScale;

    ImGuiIO& io = ImGui::GetIO();
    io.FontGlobalScale *= ratio;
    ImGui::GetStyle().ScaleAllSizes(ratio);

    // Build default font if none exists
    if (io.Fonts->Fonts.empty()) {
        io.Fonts->AddFontDefault();
        io.Fonts->Build();
    }
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
    float scale = (xscale + yscale) * 0.5f;
    if (scale < 1.25f) scale = 2.0f; // Fallback logic
    applyUiScale(scale);
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
    ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), ImGuiDockNodeFlags_None);

    // 2. Draw Content Windows

    // Render Settings Window
    ImGui::Begin("Render Settings"); 

    ImGui::Text("Application Average: %.3f ms/frame (%.1f FPS)",
        1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
    ImGui::Separator();

    ImGui::Checkbox("Paused", &s.paused);
    ImGui::SliderInt("SPP / frame", &s.sppPerFrame, 1, 64);
    // ... 他の項目 ...

    if (ImGui::Button("Reset accumulation")) {
        s.requestReset = true;
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
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}