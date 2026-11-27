#pragma once
// The kernel renderer's scene (Kernel/KernelScene.hpp) built from a scene file, and the --device cpu-kernel | gpu render
// (docs/GPU.md). The kernel renderer models a subset of the engine (quads and triangle meshes; diffuse surfaces with a constant
// albedo, isotropic rough conductors, emitting quads; a uniform environment; a pinhole camera; the RGB build); anything else is
// refused with a list of what is missing, instead of being rendered differently.

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
}

namespace rayt::cli {

/// The arrays of a kernel scene, owned on the CPU, and the header (counts, camera, environment, render settings).
struct KernelSceneHost {
    std::vector<kernel::KFloat3> positions, normals, faceNormals;
    std::vector<uint32_t> indices, triMaterial;
    std::vector<int32_t> triLight;
    std::vector<uint8_t> triSmooth;
    std::vector<kernel::KMaterial> materials;
    std::vector<kernel::KQuadLight> lights;
    std::vector<float> lightCdf;
    kernel::KScene header{};   ///< the fields that are not pointers

    /// The scene with its pointers on these CPU arrays (valid while this object lives and is not changed).
    kernel::KScene view() const;
};

/// Builds the kernel scene of a scene file seen through `camera` (which must be a plain pinhole). `engineMaterials` are the engine's
/// materials of the scene file by name (SceneSetup::namedMaterials): materials whose parameters the engine derives (a conductor's
/// RGB index from its table) are read from them, so the kernel never derives them a second way. False: `why` lists every part of
/// the scene or of the options that the kernel renderer does not model.
bool PackKernelScene(const SceneDescription& d, const std::map<std::string, const Material*>& engineMaterials, const Camera& camera, const CliOptions& opt,
    int width, int height, KernelSceneHost& out, std::string& why);

/// A scene file as a kernel scene of width x height pixels, built as the CLI builds it (the engine's materials included) but without
/// rendering; `edit` (optional) changes the scene's description first. For the checks (Checks/GpuChecks.cpp).
bool PackSceneFile(const std::string& path, const CliOptions& opt, int width, int height, const std::function<void(SceneDescription&)>& edit,
    KernelSceneHost& out, std::string& why);

/// The kernel render of a scene file (--device cpu-kernel: the CPU, gpu: the GPU module) into `film`. Returns 0, or an exit code
/// after printing why it cannot render.
int RenderOnKernel(const CliOptions& opt, const std::map<std::string, const Material*>& engineMaterials, const Camera& camera, Film& film);

/// Renders a kernel scene on the CPU (Embree in float, OpenMP over rows): `mean` receives width x height x 3 floats, the mean of
/// each pixel (row 0 = the bottom). False with the reason in `why`.
bool RenderKernelCpu(const kernel::KScene& s, std::vector<float>& mean, std::string& why, bool progress);

/// The same on the GPU (the GPU module, OptiX). `unavailable` (optional) is set when there is no module or no usable GPU, as
/// opposed to a failed render; `device` (optional) receives a description of the GPU.
bool RenderKernelGpu(const kernel::KScene& s, std::vector<float>& mean, std::string& why, bool progress, bool* unavailable = nullptr,
    std::string* device = nullptr);

/// The spectrum of an emission material (Scenes/FileScene.cpp).
bool MakeEmissionSpectrum(const SpectrumSpec& s, EmissionSpectrum& out);

} // namespace rayt::cli
