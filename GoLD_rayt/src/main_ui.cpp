// GoLD_rayt.cpp 
//
/// @brief 金属をJohnsonの物理データを基にレンダリング
///
/// 

// C++20を仮定している
// Precompiled Header (Must be first)
#include "pch.h"    

#define GLFW_INCLUDE_NONE
#include <glad/glad.h>
#include <GLFW/glfw3.h>

//#include <imgui.h>
//#include <backends/imgui_impl_glfw.h>
//#include <backends/imgui_impl_opengl3.h>

using namespace rayt;

#include "Application/Application.hpp"

int main(int argc, char** argv)
{
    App app(1280, 800, "GoLD_rayt Research UI");

    // --selftest: scripted check of the UI-driven features (scene, environment, roughness, exposure, save ...);
    // exits with code 2 on failure instead of waiting for the window to be closed.
    // --selftest-keep also keeps the saved image (gui_<time>_<spp>spp.png / .hdr) so it can be looked at.
    if (argc > 1 && (std::string(argv[1]) == "--selftest" || std::string(argv[1]) == "--selftest-keep")) {
        return app.selfTest(std::string(argv[1]) == "--selftest-keep") ? 0 : 2;
    }

    app.run();
    return 0;
}