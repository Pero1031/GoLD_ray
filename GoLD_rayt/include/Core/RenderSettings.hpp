/*
* @file include/Core/RenderSettings.hpp
*
* Settings shared between the GUI (ImGuiLayer) and the application loop (App).
*/

#pragma once

#include <string>
#include <vector>

// For the view modes (add more as needed)
enum class ViewMode {
    Beauty = 0,
    Normal,
    Depth,
    Albedo
};

struct RenderSettings {
    // Render control
    bool paused = false;
    bool requestReset = false;
    bool requestSave = false;        // "Save image" button: write the accumulated image once
    bool requestRebuild = false;     // scene or environment map changed: rebuild the scene and restart

    // Quality settings
    int sppPerFrame = 1;             // render passes (1 spp each) per displayed frame
    int targetSpp = 1024;            // stop after this many samples per pixel
    int maxDepth = 8;

    // --- Scene ---
    int scene = 0;                   // 0 = gold spheres, 1 = glass spheres
    int env = 0;                     // index into envNames (0 = no environment map)
    std::vector<std::string> envNames{ "None (black)" };   // filled at startup from assets/env/*.hdr

    // --- Material tweak ---
    float roughness = 0.2f;          // [0,1] roughness of the CENTER gold sphere (the outer ones stay at 0.01 / 0.5)

    // Camera and display
    float exposure = 0.0f;           // display exposure in stops (EV); the PNG written by "Save image" uses it too
    ViewMode viewMode = ViewMode::Beauty;

    // Resolution (initial value)
    int width = 1280;
    int height = 800;

    // Read-only information for the UI
    int currentSpp = 0;
};
