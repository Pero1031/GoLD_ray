#pragma once
// MMD (PMX) helpers of the CLI: finding model / pose files, printing, converting a model to a mesh and materials.

#include <filesystem>
#include <memory>
#include <string>
#include <vector>
#include "Core/Texture.hpp"
#include "Geometry/TriangleMesh.hpp"
#include "IO/PmxLoader.hpp"
#include "Materials/Material.hpp"
#include "Materials/Subsurface.hpp"

namespace rayt::cli {

// -----------------------------------------------------------------------------
// Triangle mesh checks (run at the start of --scene mesh)
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// MMD (PMX) models  (--scene pmx)
// -----------------------------------------------------------------------------
// The files live in MMD_model/ (git-ignored, third-party assets). Names are CJK, which a Windows
// command line cannot carry reliably, so a model is chosen by its index in --pmx-list.
std::string PathToUtf8(const std::filesystem::path& p);

std::filesystem::path PathFromUtf8(const std::string& s);

std::vector<std::filesystem::path> FindMmdFiles(const std::string& subdir, const std::string& extension);

std::vector<std::filesystem::path> FindPmxFiles();

Real SrgbToLinear(Real c);

void PrintPmxInfo(const io::PmxModel& m, const std::filesystem::path& file);

/// PMX material -> renderer material.
///
/// MMD materials are an NPR recipe (diffuse color, a Phong specular term, an additive sphere map for
/// highlights, a toon ramp for the shadow side). A path tracer produces shadows and highlights by
/// itself, so only the intent is kept:
///  - PMX specular > 0 : a glossy dielectric coat (weight = specular strength, roughness from the
///    shininess via Blinn-Phong exponent -> GGX alpha = sqrt(2 / (n + 2))).
///  - this model has specular 0 everywhere and fakes gloss with sphere-mapped overlay layers, so
///    its gloss is restored from the material NAMES (hair, eyes, ornaments, skin). That table is
///    specific to this model; --pmx-matte turns every preset off.
///  - toon ramps, sphere maps and edge lines are not used (shading comes from the lighting).
std::shared_ptr<Material> MakePmxMaterial(const io::PmxMaterial& pm, const RgbColor& base, const std::shared_ptr<Texture2D>& tex, bool matte, const std::shared_ptr<const SubsurfaceProfile>& skinProfile = nullptr);

/// PMX -> TriangleMesh in the renderer's convention (right-handed, +Y up, metres, feet on y = 0, x/z centered).
/// MMD is left-handed with 1 unit ~ 8 cm; mirroring z fixes the handedness, yaw (degrees) turns the model about +Y.
std::shared_ptr<TriangleMesh> PmxToMesh(const io::PmxModel& m, Real yawDegrees, Real* outHeight, Vector3* outSize = nullptr);
} // namespace rayt::cli
