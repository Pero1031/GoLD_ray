/**
 * @file ImGuiLayer.hpp
 * @brief Manages the user interface layer using Dear ImGui.
 */

#pragma once

struct GLFWwindow;
struct RenderSettings;

/**
 * @brief Handles the lifecycle and rendering of the ImGui-based user interface.
 * 
 * This class encapsulates the initialization, frame management, and rendering
 * of the GUI. It interacts with the GLFW window and modifies the renderer settings.
 */
class ImGuiLayer {
public:
    /**
     * @brief Initializes the ImGui context and backends.
     * Sets up the ImGui style, connects to the GLFW window, initializes the
     * OpenGL backend, and determines the initial DPI scaling.
     * @param window Pointer to the GLFW window handle.
     */
    void init(GLFWwindow* window);

     /**
      * @brief Starts a new ImGui frame.
      * Must be called at the beginning of the render loop, before issuing any
      * ImGui commands. Handles input polling and time updates.
      */
    void beginFrame();

    /**
     * @brief Renders the UI components to control renderer settings.
     * Draws widgets (sliders, checkboxes, etc.) that allow the user to
     * modify the render settings in real-time.
     * @param[in,out] settings Reference to the render settings structure to be modified by the UI.
     */
    void draw(RenderSettings& settings);

    /**
     * @brief Finalizes the frame and executes the draw data.
     * Calls ImGui::Render() and dispatches the draw commands to the GPU.
     */
    void endFrame();

    /**
     * @brief Cleans up ImGui resources.
     * Shuts down the ImGui backends and destroys the context.
     */
    void shutdown();

    /**
     * @brief Callback to handle changes in the window content scale (DPI).
     * This should be called from the GLFW window content scale callback.
     * It recalculates and applies the UI scaling factor.
     * @param xscale The new x-axis content scale.
     * @param yscale The new y-axis content scale.
     */
    void onContentScaleChanged(float xscale, float yscale);

private:
    /**
     * @brief Internal helper to apply the global font and style scaling.
     * @param newScale The target scaling factor (e.g., 1.0 for 100%, 1.5 for 150%).
     */
    void applyUiScale(float newScale);

    /// Current UI scaling factor.
    float m_uiScale = 1.0f;
};