// The commands of the CLI that do not render a picture.
#include "pch.h"
#include "Cli/Commands.hpp"
#include "Cli/Checks.hpp"
#include "Cli/PmxSupport.hpp"
#include "Cli/Tools.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

namespace {

struct Command {
    CommandDef def;
    std::function<int(const CliOptions&)> run;
};

const std::vector<Command>& Commands() {
    static const std::vector<Command> commands = {
        { { "lights", "light source checks against closed forms and numerical integrals" }, [](const CliOptions& o) { return RunLightTests(o); } },
        { { "aov", "denoiser feature (albedo / normal) checks" }, [](const CliOptions&) { return RunAovTests(); } },
        { { "dispersion", "Sellmeier presets of glasses (n_d, Abbe number) and a diamond in a white furnace" }, [](const CliOptions&) { return RunDispersionTests(); } },
        { { "thinfilm", "thin-film interference checks (Airy's formula, sample() weights, white furnace)" }, [](const CliOptions&) { return RunThinFilmTests(); } },
        { { "sss", "subsurface scattering checks" }, [](const CliOptions& o) { return RunSssTests(o); } },
        { { "sssrw", "the BSSRDF against a volumetric random walk on a sphere and thin discs (--sssrw-albedo, --sssrw-g, --sssrw-ell, --sssrw-A)" }, [](const CliOptions& o) { return RunSssRandomWalk(o); } },
        { { "skinlayers", "Monte Carlo radial profile of a layered skin (--melanin, --blood, --skin-out ...)" }, [](const CliOptions& o) { return RunSkinLayers(o); } },
        { { "volume", "volumetric random walk checks" }, [](const CliOptions& o) { return RunVolumeTests(o); } },
        { { "layered", "layered coat (clear coat over diffuse / gold): closed form, eval vs sample, energy, reciprocity" }, [](const CliOptions& o) { return RunLayeredTests(o); } },
        { { "aniso", "anisotropic metal: tangent of the shapes, anisotropic GGX, rotation of the roughness axes" }, [](const CliOptions& o) { return RunAnisotropyTests(o); } },
        { { "caustic", "caustic photon map check against plain path tracing" }, [](const CliOptions& o) { return RunCausticTests(o); } },
        { { "lvc", "bidirectional path tracing with a cache of light paths in a fog (prototype) against path tracing" }, [](const CliOptions& o) { return RunLvcTests(o); } },
        { { "manifoldtest", "prism caustic: specular manifold NEE against the photon map" }, [](const CliOptions& o) { return RunManifoldTests(o, false); } },
        { { "manifoldcup", "tumbler caustic: specular manifold NEE against the photon map" }, [](const CliOptions& o) { return RunManifoldTests(o, true); } },
        { { "metals", "metal data checks (file layouts, reflectance input, edge tint)" }, [](const CliOptions& o) { return RunMetalTests(o); } },
        { { "bench", "BVH benchmark: recursive reference vs SAH BVH" }, [](const CliOptions& o) { return RunBvhBenchmark(o); } },
    };
    return commands;
}

} // namespace

const CommandDef* FindCommand(const std::string& name) {
    for (const Command& c : Commands())
        if (name == c.def.name) return &c.def;
    return nullptr;
}

std::vector<const CommandDef*> AllCommands() {
    std::vector<const CommandDef*> all;
    for (const Command& c : Commands()) all.push_back(&c.def);
    return all;
}

std::optional<int> RunNonRenderCommand(const CliOptions& opt) {
    for (const Command& c : Commands())
        if (opt.scene == c.def.name) return c.run(opt);
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
