#pragma once
// Rendering a scene: from the scene table entry to the saved image.

#include <memory>

#include "Cli/Options.hpp"

namespace rayt {
    class Camera;
}

namespace rayt::cli {

struct SceneSetup;

struct ImageSize {
    int width;
    int height;
};

// Picture size: --width, and the height from the format of the scene (16:9, or 3:4 for a standing figure).
ImageSize ComputeImageSize(const CliOptions& opt);

// The camera a render of this scene builds (its lens, the stop and the focus), without the notes it prints; null if it cannot be built.
// For the checks of the kernel renderer (Cli/KernelScenePack.hpp).
std::shared_ptr<Camera> BuildSceneCamera(const CliOptions& opt, const SceneSetup& setup, int imageWidth, int imageHeight);

// Builds the scene named by opt.scene and renders it. Returns the exit code (0 = fine, 2 = a check failed).
int RenderScene(CliOptions& opt);

} // namespace rayt::cli
