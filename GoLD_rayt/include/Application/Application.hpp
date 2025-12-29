#pragma once

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <memory>

#include "Renderer/Film.hpp"
#include "Renderer/Scene.hpp"
#include "Renderer/Integrator.hpp" 
#include "../src/UI/ImGuiLayer.hpp"
#include "Core/RenderSettings.hpp"

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

private:
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

private:
    // Window management
    GLFWwindow* m_window = nullptr;   ///< Pointer to the GLFW window handle.
    int m_width, m_height;            

    // Rendering components
    std::shared_ptr<rayt::Scene> m_scene;
    std::shared_ptr<rayt::Camera> m_camera;
    std::unique_ptr<rayt::PathIntegrator> m_integrator;
    std::unique_ptr<rayt::Film> m_film;

    // Display resources
    GLuint m_renderTexture = 0;
    int m_currentSpp = 0;

    // UI & Settings
    ImGuiLayer m_uiLayer;
    RenderSettings m_settings;
};