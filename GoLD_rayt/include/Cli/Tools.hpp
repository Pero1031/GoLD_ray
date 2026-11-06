#pragma once
// Measurement and debugging tools of the CLI (not self-checks, not scenes).

#include <string>
#include <utility>
#include <vector>
#include "Materials/Material.hpp"
#include "Render/Camera.hpp"
#include "Render/Film.hpp"
#include "Scene/Scene.hpp"
#include "Cli/Options.hpp"

namespace rayt::cli {

// -----------------------------------------------------------------------------
// Layered skin profile (--scene skinlayers)
// -----------------------------------------------------------------------------
int RunSkinLayers(const CliOptions& opt);

// -----------------------------------------------------------------------------
// BVH benchmark (--scene bench): the reference recursive BVH vs the SAH flat BVH
// -----------------------------------------------------------------------------
int RunBvhBenchmark(const CliOptions& opt);

// metalprobe report: mean film color of the metal sphere per angle of incidence (one specular
// bounce off a uniform emitter, so the value is F(theta) x light). Bins are +- 1.5 degrees wide.
void ReportMetalProbe(const Scene& scene, const Camera& camera, const Film& film, const Material* metal, const std::string& label);

// Debug aid: follow a few light paths of one pixel and print which object each bounce hits.
//   [in]/[out] = the ray hits the surface from outside / from inside, (R)/(T) = reflected / transmitted.
void TracePixelDebug(const Scene& scene, const Camera& camera, int px, int py, int w, int h, int maxDepth, int numPaths, const std::vector<std::pair<std::string, const Material*>>& names);
} // namespace rayt::cli
