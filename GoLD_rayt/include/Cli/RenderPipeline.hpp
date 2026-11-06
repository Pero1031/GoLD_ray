#pragma once
// Rendering a scene: from the scene table entry to the saved image.

#include "Cli/Options.hpp"

namespace rayt::cli {

struct ImageSize {
    int width;
    int height;
};

// Picture size: --width, and the height from the format of the scene (16:9, or 3:4 for a standing figure).
ImageSize ComputeImageSize(const CliOptions& opt);

// Builds the scene named by opt.scene and renders it. Returns the exit code (0 = fine, 2 = a check failed).
int RenderScene(CliOptions& opt);

} // namespace rayt::cli
