/*
* @file include/Application/Application.hpp
* @brief Defines the main application entry class (window loop, rendering, UI).
*/

#pragma once

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Render/Film.hpp"
#include "Scene/Scene.hpp"
#include "Render/Integrator.hpp" 
#include "../src/UI/ImGuiLayer.hpp"
#include "Core/RenderSettings.hpp"

namespace rayt {
    class EnvLight;
    class RoughConductor;
    namespace color { class ComplexIorTable; }
}

/**
 * @class App
 * @brief The main application class responsible for managing the window, rendering loop, and UI integration.
 *
 * This class handles the initialization of GLFW and OpenGL, manages the main event loop,
 * and coordinates data transfer between the CPU-based ray tracer (Integrator/Film)
 * and the GPU-based display (OpenGL Texture).
 */
class App {
public:
    /**
     * @brief Constructs the Application.
     *
     * @param width The initial width of the window in pixels.
     * @param height The initial height of the window in pixels.
     * @param title The title text displayed on the window border.
     */
    App(int width, int height, const char* title);

    /**
     * @brief Destructor.
     *
     * Cleans up GLFW resources, OpenGL textures, and UI contexts.
     */
    ~App();

    /**
     * @brief Starts the main application loop.
     *
     * This method enters an infinite loop that processes input, updates the scene,
     * performs rendering (integration), and draws the UI until the window is closed.
     */
    void run();

    /**
     * @brief Drives the application loop through scripted settings changes and checks the results.
     *
     * Used by `GoLD_rayt.exe --selftest`: scene / environment switches, exposure, roughness, SPP per
     * frame, pause, sample limit and saving, without any user input. The window stays open meanwhile.
     * @param keepSavedFiles keep the image written by the "Save image" check instead of deleting it
     * @return true if every check passed
     */
    bool selfTest(bool keepSavedFiles = false);

private:
    /// One iteration of the main loop (events, UI requests, render passes, texture, UI drawing).
    void frame();

    /**
     * @brief Initializes the application resources.
     *
     * Sets up the OpenGL context, ImGui context, creates the initial scene,
     * and initializes the renderer components (Camera, Integrator, Film).
     */
    void init();

    /**
     * @brief Updates the application state per frame.
     *
     * Handles user input (keyboard/mouse) and updates camera parameters or settings.
     */
    void update();

    /**
     * @brief Renders the current frame to the screen.
     *
     * This involves binding the OpenGL texture, drawing a fullscreen quad (or ImGui window)
     * with the ray-traced image, and rendering the ImGui overlay.
     */
    void render();

    /**
     * @brief Updates the OpenGL texture with the latest pixel data.
     *
     * This method uploads the accumulated pixel data from the CPU-side Film buffer
     * to the GPU-side OpenGL texture (`m_renderTexture`).
     * It should be called whenever the film has been updated by the integrator.
     */
    void updateTexture();

    /**
     * @brief Builds the scene selected in the settings (gold or glass spheres, floor, environment light).
     *
     * Called at startup and whenever the scene or the environment map changes; the materials are
     * created once (loadMaterials) and shared by every scene.
     */
    void buildScene();

    /// Creates the materials (gold from the Johnson Au table, glass, floor) and the list of environment maps.
    void loadMaterials();

    /// Environment light of settings index `index` (0 = none), loaded on first use and cached.
    std::shared_ptr<rayt::EnvLight> loadEnv(int index);

    /// Writes the accumulated image: PNG (with the exposure of the UI) and a linear HDR.
    void saveImage() const;

private:
    // Window management
    GLFWwindow* m_window = nullptr;   ///< Pointer to the GLFW window handle.
    int m_width, m_height;            

    // Rendering components
    std::shared_ptr<rayt::Scene> m_scene;
    std::shared_ptr<rayt::Camera> m_camera;
    std::unique_ptr<rayt::PathIntegrator> m_integrator;
    std::unique_ptr<rayt::Film> m_film;

    // Materials shared by the scenes
    std::shared_ptr<rayt::color::ComplexIorTable> m_auTable;
    std::shared_ptr<rayt::RoughConductor> m_goldSmooth, m_goldCenter, m_goldRough;
    std::shared_ptr<rayt::Material> m_glass, m_frosted, m_floor;

    // Environment maps found in assets/env (index 0 = none)
    std::vector<std::string> m_envPaths;
    std::map<int, std::shared_ptr<rayt::EnvLight>> m_envCache;

    // Display resources
    GLuint m_renderTexture = 0;
    int m_currentSpp = 0;
    float m_shownExposure = 0.0f;       ///< exposure the texture was last converted with
    float m_appliedRoughness = 0.2f;    ///< roughness the center gold sphere currently has

    // UI & Settings
    ImGuiLayer m_uiLayer;
    RenderSettings m_settings;
};