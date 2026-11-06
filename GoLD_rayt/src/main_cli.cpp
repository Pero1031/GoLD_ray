// GoLD_rayt command line renderer.
//
// main() parses the options, runs a command that does not render (self-checks, measurement tools), or renders the scene.
//   include/Cli, src/Cli        the CLI: Options, Commands, Checks/, Tools/, Scenes/ (one builder per scene), RenderPipeline
//   include/*, src/*            the renderer (header-only engine)

#include "pch.h"

#include "Cli/Commands.hpp"
#include "Cli/Options.hpp"
#include "Cli/RenderPipeline.hpp"
#include "Core/Precision.hpp"

using namespace rayt;
using namespace rayt::cli;

int main(int argc, char** argv) {
    CliOptions opt;
    if (const std::optional<int> exitCode = ParseArgs(argc, argv, opt)) return *exitCode;
    const ImageSize size = ComputeImageSize(opt);
    std::cout << "[Render] mode = " << kRenderModeName << "  (switch: Core/RenderMode.hpp)\n";
    std::cout << "[Options] scene=" << opt.scene << " " << size.width << "x" << size.height
        << " spp=" << opt.spp << " maxDepth=" << opt.maxDepth << " out=" << opt.out << "\n";
    if (opt.accel == "embree" && opt.scene != "mesh" && opt.scene != "bench") {
        // Embree is wired up for single-material triangle meshes only (the mesh and bench scenes); the PMX scene
        // (many materials, alpha cut-outs) and the primitive scenes always use their own BVH / lists.
        std::cout << "[Accel] note: --accel embree has no effect in the '" << opt.scene << "' scene\n";
    }

    if (const std::optional<int> exitCode = RunNonRenderCommand(opt)) return *exitCode;
    return RenderScene(opt);
}
