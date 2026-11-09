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
#include "Core/Sampling.hpp"

#include <random>

using namespace rayt;
using namespace rayt::cli;

int main(int argc, char** argv) {
    CliOptions opt;
    if (const std::optional<int> exitCode = ParseArgs(argc, argv, opt)) return *exitCode;
    const ImageSize size = ComputeImageSize(opt);
    // Every run has a seed, so that a picture can be made again from what its log says (--seed N).
    if (opt.seed < 0) opt.seed = static_cast<long long>(std::random_device{}() & 0x7fffffffu);
    sampling::SetRunSeed(static_cast<std::uint64_t>(opt.seed));
    sampling::SeedStream(0);   // this thread builds the scene
    std::cout << "[Render] mode = " << kRenderModeName << "  (switch: Core/RenderMode.hpp)\n";
    std::cout << "[Options] scene=" << opt.scene << " " << size.width << "x" << size.height
        << " spp=" << opt.spp << " maxDepth=" << opt.maxDepth << " seed=" << opt.seed << " out=" << opt.out << "\n";
    if (opt.accel == "embree" && opt.scene != "mesh" && opt.scene != "bench") {
        // Embree is wired up for single-material triangle meshes only (the mesh and bench scenes); the PMX scene
        // (many materials, alpha cut-outs) and the primitive scenes always use their own BVH / lists.
        std::cout << "[Accel] note: --accel embree has no effect in the '" << opt.scene << "' scene\n";
    }

    if (const std::optional<int> exitCode = RunNonRenderCommand(opt)) return *exitCode;
    return RenderScene(opt);
}
