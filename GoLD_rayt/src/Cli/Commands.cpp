// The commands of the CLI that do not render a picture.
#include "pch.h"
#include "Cli/Commands.hpp"
#include "Cli/Checks.hpp"
#include "Cli/PmxSupport.hpp"
#include "Cli/Tools.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

std::optional<int> RunNonRenderCommand(const CliOptions& opt) {
    if (opt.scene == "lights") return RunLightTests(opt);
    if (opt.scene == "aov") return RunAovTests();
    if (opt.scene == "sss") return RunSssTests(opt);
    if (opt.scene == "skinlayers") return RunSkinLayers(opt);
    if (opt.scene == "volume") return RunVolumeTests(opt);
    if (opt.scene == "caustic") return RunCausticTests(opt);
    if (opt.scene == "manifoldtest") return RunManifoldTests(opt, false);
    if (opt.scene == "manifoldcup") return RunManifoldTests(opt, true);
    if (opt.scene == "metals") return RunMetalTests(opt);
    if (opt.scene == "bench") return RunBvhBenchmark(opt);
    if (opt.scene == "pmx" && opt.poseList) {
        const auto files = FindMmdFiles("pose", ".vpd");
        for (size_t i = 0; i < files.size(); ++i)
            std::cout << "[Pose] " << i << ": " << PathToUtf8(files[i].parent_path().filename()) << "/" << PathToUtf8(files[i].filename()) << "\n";
        if (files.empty()) std::cout << "[Pose] no .vpd files found in MMD_model/pose\n";
        return 0;
    }
    if (opt.scene == "pmx" && opt.pmxList) {
        const auto files = FindPmxFiles();
        for (size_t i = 0; i < files.size(); ++i)
            std::cout << "[PMX] " << i << ": " << PathToUtf8(files[i].filename()) << "  (" << std::filesystem::file_size(files[i]) / 1024 << " KB)\n";
        if (files.empty()) std::cout << "[PMX] no .pmx files found in MMD_model/\n";
        return 0;
    }
    return std::nullopt;
}

} // namespace rayt::cli
