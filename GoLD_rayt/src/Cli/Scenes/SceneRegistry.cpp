// The scene table: name, picture format, default path depth, kind of environment and the builder of every scene.
#include "pch.h"
#include "Cli/Scenes.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

namespace {

constexpr double kDefaultHeightPerWidth = double(IMAGE_HEIGHT) / IMAGE_WIDTH;   // 16:9

const SceneDef kScenes[] = {
    { "gold", kDefaultHeightPerWidth, 0, [](const CliOptions&) { return EnvKind::Hdr; },
      BuildGoldSpheres },
    { "goldrough", kDefaultHeightPerWidth, 0, [](const CliOptions&) { return EnvKind::Hdr; },
      BuildGoldRough },
    { "glass", kDefaultHeightPerWidth, 0, [](const CliOptions&) { return EnvKind::Hdr; },
      BuildGlassSpheres },
    { "studio", kDefaultHeightPerWidth, 0, [](const CliOptions&) { return EnvKind::None; },
      BuildStudio },
    { "mesh", kDefaultHeightPerWidth, 0, [](const CliOptions&) { return EnvKind::None; },
      BuildMesh },
    { "spd", kDefaultHeightPerWidth, 0, [](const CliOptions&) { return EnvKind::None; },
      BuildLightSpectra },
    { "goldbowl", kDefaultHeightPerWidth, 0, [](const CliOptions&) { return EnvKind::None; },
      BuildGoldBowl },
    { "metalprobe", kDefaultHeightPerWidth, 0, [](const CliOptions&) { return EnvKind::None; },
      BuildMetalProbe },
    { "pmx", 4.0 / 3.0, 0, [](const CliOptions&) { return EnvKind::Hdr; },
      BuildPmx },
    { "skindemo", kDefaultHeightPerWidth, 0, [](const CliOptions&) { return EnvKind::None; },
      BuildSkinDemo },
    { "blob", 0.75, 256, [](const CliOptions&) { return EnvKind::Hdr; },
      BuildBlob },
    { "water", kDefaultHeightPerWidth, 64, [](const CliOptions& o) { return o.waterSky ? EnvKind::HdrWithPeak : EnvKind::None; },
      BuildWater },
    { "sea", kDefaultHeightPerWidth, 48, [](const CliOptions& o) { return o.seaHdr ? EnvKind::HdrWithPeak : EnvKind::ProceduralSky; },
      BuildSea },
    { "pool", kDefaultHeightPerWidth, 48, [](const CliOptions&) { return EnvKind::Hdr; },
      BuildPool },
    { "manifold", kDefaultHeightPerWidth, 32, [](const CliOptions&) { return EnvKind::None; },
      BuildManifold },
    { "crystal", kDefaultHeightPerWidth, 48, [](const CliOptions&) { return EnvKind::None; },
      BuildCrystal },
    { "furnace", kDefaultHeightPerWidth, 0, [](const CliOptions&) { return EnvKind::Uniform; },
      BuildFurnace },
};

} // namespace

const SceneDef* FindScene(const std::string& name) {
    for (const SceneDef& s : kScenes)
        if (name == s.name) return &s;
    return nullptr;
}

} // namespace rayt::cli
