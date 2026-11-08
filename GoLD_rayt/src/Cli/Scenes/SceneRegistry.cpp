// The scene table: name, description, picture format, default path depth, kind of environment and the builder of every scene.
#include "pch.h"
#include "Cli/Scenes.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

namespace {

constexpr double kDefaultHeightPerWidth = double(IMAGE_HEIGHT) / IMAGE_WIDTH;   // 16:9

const SceneDef kScenes[] = {
    { "gold", "three gold spheres of different roughness", kDefaultHeightPerWidth, 0,
      [](const CliOptions&) { return EnvKind::Hdr; }, BuildGoldSpheres },
    { "goldrough", "rough gold spheres (--no-ms: without the multiple-scattering compensation)", kDefaultHeightPerWidth, 0,
      [](const CliOptions&) { return EnvKind::Hdr; }, BuildGoldRough },
    { "glass", "clear, dispersive (--abbe) and frosted glass spheres", kDefaultHeightPerWidth, 0,
      [](const CliOptions&) { return EnvKind::Hdr; }, BuildGlassSpheres },
    { "studio", "gold, glass and a red sphere under a tungsten softbox and an LED", kDefaultHeightPerWidth, 0,
      [](const CliOptions&) { return EnvKind::None; }, BuildStudio },
    { "mesh", "triangle meshes (a gold torus, flat and smooth spheres); runs the mesh checks first", kDefaultHeightPerWidth, 0,
      [](const CliOptions&) { return EnvKind::None; }, BuildMesh },
    { "spd", "colored spheres and gold under one light spectrum (--light, --wb)", kDefaultHeightPerWidth, 0,
      [](const CliOptions&) { return EnvKind::None; }, BuildLightSpectra },
    { "goldbowl", "a gold bowl with a gold ball: where RGB and spectral rendering differ (--light)", kDefaultHeightPerWidth, 0,
      [](const CliOptions&) { return EnvKind::None; }, BuildGoldBowl },
    { "metalprobe", "a smooth metal sphere in a uniform shell, for the color per angle of incidence (--metal, --light)", kDefaultHeightPerWidth, 0,
      [](const CliOptions&) { return EnvKind::None; }, BuildMetalProbe },
    { "pmx", "an MMD character from MMD_model/ (--pmx-list, --pose-list, --expr, --sss)", 4.0 / 3.0, 0,
      [](const CliOptions&) { return EnvKind::Hdr; }, BuildPmx },
    { "skindemo", "a skin sphere and ring for subsurface scattering (--sss)", kDefaultHeightPerWidth, 0,
      [](const CliOptions&) { return EnvKind::None; }, BuildSkinDemo },
    { "blob", "a glass blob filled with a scattering medium", 0.75, 256,
      [](const CliOptions&) { return EnvKind::Hdr; }, BuildBlob },
    { "water", "a tumbler of water (--water-lamp, --water-sky, --caustics)", kDefaultHeightPerWidth, 64,
      [](const CliOptions& o) { return o.waterSky ? EnvKind::HdrWithPeak : EnvKind::None; }, BuildWater },
    { "sea", "a Gerstner-wave sea over sand (--sea-steep, --sea-depth)", kDefaultHeightPerWidth, 48,
      [](const CliOptions& o) { return o.seaHdr ? EnvKind::HdrWithPeak : EnvKind::ProceduralSky; }, BuildSea },
    { "pool", "a swimming pool with a caustic network (--caustics)", kDefaultHeightPerWidth, 48,
      [](const CliOptions&) { return EnvKind::Hdr; }, BuildPool },
    { "manifold", "the caustic of a point light through a prism (--caustic-method)", kDefaultHeightPerWidth, 32,
      [](const CliOptions&) { return EnvKind::None; }, BuildManifold },
    { "fogbeam", "a slit of light through a dispersive prism in fog: volumetric caustics by light tracing (--fog-sigma, --abbe)", kDefaultHeightPerWidth, 32,
      [](const CliOptions&) { return EnvKind::None; }, BuildFogBeam },
    { "crystal", "faceted glass on a dark stage with dispersion and caustics (--crystal-inside)", kDefaultHeightPerWidth, 48,
      [](const CliOptions&) { return EnvKind::None; }, BuildCrystal },
    { "furnace", "white furnace: objects in a uniform white environment must keep their energy", kDefaultHeightPerWidth, 0,
      [](const CliOptions&) { return EnvKind::Uniform; }, BuildFurnace },
    { "file", "a scene written as data in a JSON file: GoLD_rayt scene.json (materials, objects, lights, camera)", kDefaultHeightPerWidth, 0,
      [](const CliOptions& o) {
          if (!o.sceneDescription) return EnvKind::None;
          switch (o.sceneDescription->environment.type) {
          case EnvironmentSpec::Type::Hdr: return EnvKind::Hdr;
          case EnvironmentSpec::Type::Uniform: return EnvKind::Uniform;
          default: return EnvKind::None;
          }
      }, BuildFromDescription },
};

} // namespace

const SceneDef* FindScene(const std::string& name) {
    for (const SceneDef& s : kScenes)
        if (name == s.name) return &s;
    return nullptr;
}

std::vector<const SceneDef*> AllScenes() {
    std::vector<const SceneDef*> all;
    for (const SceneDef& s : kScenes) all.push_back(&s);
    return all;
}

} // namespace rayt::cli
