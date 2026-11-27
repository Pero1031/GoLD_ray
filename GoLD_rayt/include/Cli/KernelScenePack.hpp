#pragma once
// The kernel renderer's scene (Kernel/KernelScene.hpp) built from a scene file, and the --device cpu-kernel | gpu render
// (docs/GPU.md). The kernel renderer models a subset of the engine (quads, spheres and triangle meshes; diffuse surfaces with a
// constant albedo, isotropic rough conductors, glass, clear coats over a diffuse or a metal base, emitting quads and spheres; an
// environment map; a pinhole, a thin lens or a real lens of spherical surfaces; the RGB build); anything else is refused with a list of
// what is missing, instead of being rendered differently.

#include <functional>
#include <map>
#include <string>
#include <vector>

#include "Cli/Options.hpp"
#include "Cli/SceneDescription.hpp"
#include "Color/EmissionSpectrum.hpp"
#include "Kernel/KernelScene.hpp"

namespace rayt {
    class Camera;
    class Film;
    class Material;
    class EnvLight;
}

namespace rayt::cli {

struct SceneSetup;
struct Environment;

/// The arrays of a kernel scene, owned on the CPU, and the header (counts, camera, environment, render settings).
struct KernelSceneHost {
    std::vector<kernel::KFloat3> positions, normals, faceNormals;
    std::vector<uint32_t> indices, triMaterial;
    std::vector<int32_t> triLight;
    std::vector<uint8_t> triSmooth;
    std::vector<kernel::KMaterial> materials;
    std::vector<kernel::KSphere> spheres;
    std::vector<kernel::KLight> lights;
    std::vector<float> lightCdf;
    std::vector<float> envPixels, envMarginalFunc, envMarginalCdf, envCondFunc, envCondCdf, envCondInt;
    std::vector<kernel::LensSurface<float>> lensSurfaces;
    std::vector<float> lensPupil, lensBladeNormals;
    kernel::KScene header{};   ///< the fields that are not pointers

    /// The scene with its pointers on these CPU arrays (valid while this object lives and is not changed).
    kernel::KScene view() const;
};

/// Builds the kernel scene of a scene file seen through `camera` (which must be a plain pinhole). `engineMaterials` are the engine's
/// materials of the scene file by name (SceneSetup::namedMaterials): materials whose parameters the engine derives (a conductor's
/// RGB index from its table) are read from them, so the kernel never derives them a second way; so is the environment (`envLight`,
/// the engine's environment light: its map, its sampling distribution and its scale). False: `why` lists every part of
/// the scene or of the options that the kernel renderer does not model.
bool PackKernelScene(const SceneDescription& d, const std::map<std::string, const Material*>& engineMaterials, const EnvLight* envLight, const Camera& camera,
    const CliOptions& opt, int width, int height, KernelSceneHost& out, std::string& why);

/// The engine's environment light (its map, its sampling distribution, its scale) into the kernel scene; false if there is none.
bool PackEnvironment(const EnvLight* envLight, KernelSceneHost& out);

/// A scene file as a kernel scene of width x height pixels, built as the CLI builds it (the engine's materials and camera included) but
/// without rendering; `edit` (optional) changes the scene's description first; `engine` (optional) is called with the engine's build of
/// the scene and its camera once the kernel scene is packed (for checks that compare the two). For the checks (Checks/GpuChecks.cpp).
bool PackSceneFile(const std::string& path, const CliOptions& opt, int width, int height, const std::function<void(SceneDescription&)>& edit,
    KernelSceneHost& out, std::string& why, const std::function<void(const SceneSetup&, const Environment&, const Camera&)>& engine = {});

/// The kernel render of a scene file (--device cpu-kernel: the CPU, gpu: the GPU module) into `film`. Returns 0, or an exit code
/// after printing why it cannot render.
int RenderOnKernel(const CliOptions& opt, const std::map<std::string, const Material*>& engineMaterials, const EnvLight* envLight, const Camera& camera,
    Film& film);

/// Renders a kernel scene on the CPU (Embree in float, OpenMP over rows): `mean` receives width x height x 3 floats, the mean of
/// each pixel (row 0 = the bottom). `albedo`, `normal` (optional, both or neither): the denoiser's features, the mean over the
/// samples of each pixel as the engine averages them (a sample that recorded none adds zero). False with the reason in `why`.
bool RenderKernelCpu(const kernel::KScene& s, std::vector<float>& mean, std::string& why, bool progress,
    std::vector<float>* albedo = nullptr, std::vector<float>* normal = nullptr);

/// The same on the GPU (the GPU module, OptiX). `unavailable` (optional) is set when there is no module or no usable GPU, as
/// opposed to a failed render; `device` (optional) receives a description of the GPU.
bool RenderKernelGpu(const kernel::KScene& s, std::vector<float>& mean, std::string& why, bool progress, bool* unavailable = nullptr,
    std::string* device = nullptr, std::vector<float>* albedo = nullptr, std::vector<float>* normal = nullptr);

/// Follows the given camera rays through the kernel renderer on the CPU (ray i draws from the stream of pixel i, sample 0) and
/// returns the denoiser's features of each, and the material of the first surface it hits (kernel::kNoHit: none). For the checks.
bool KernelRayAovs(const kernel::KScene& s, const std::vector<kernel::KRay>& rays, std::vector<kernel::KAov>& aovs, std::vector<uint32_t>& firstMaterial,
    std::string& why);

/// Renders a kernel scene on the CPU, one thread, and writes down every ray it casts with Embree's answer: closest-hit queries
/// (shadow 0) as they are, shadow queries (shadow 1) as prim 0 (blocked) or kNoHit. For comparing the GPU's tracer with Embree's.
bool RecordPathRays(const kernel::KScene& s, std::vector<kernel::KRay>& rays, std::vector<kernel::KHit>& hits, std::vector<uint8_t>& shadow, std::string& why);

/// The spectrum of an emission material (Scenes/FileScene.cpp).
bool MakeEmissionSpectrum(const SpectrumSpec& s, EmissionSpectrum& out);

} // namespace rayt::cli
