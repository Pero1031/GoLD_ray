#pragma once
// The scene builders of the CLI (one function per --scene that renders an image; see Scene.hpp for the contract).

#include "Cli/Scene.hpp"

namespace rayt::cli {

std::optional<int> BuildGoldSpheres(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildGoldRough(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildLightSpectra(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildGoldBowl(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildMetalProbe(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildGlassSpheres(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildFurnace(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildBlob(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildManifold(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildFogBeam(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildFogSpot(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildCausticGallery(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildCrystal(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildWater(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildSea(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildPool(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildStudio(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildMesh(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildSplats(SceneContext& ctx, SceneSetup& out);   // a synthetic cloud of Gaussians around a gold sphere; runs the splat checks first
std::optional<int> BuildSkinDemo(SceneContext& ctx, SceneSetup& out);
std::optional<int> BuildPmx(SceneContext& ctx, SceneSetup& out);

std::optional<int> BuildFromDescription(SceneContext& ctx, SceneSetup& out);   // the scene "file": a JSON scene written as data

} // namespace rayt::cli
