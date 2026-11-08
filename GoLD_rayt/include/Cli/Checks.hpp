#pragma once
// Self-checks of the CLI (--scene furnace | lights | aov | sss | volume | caustic | manifoldtest | manifoldcup | metals | mesh, --pmx-check).
// Each prints OK / FAIL per case; the Run* functions return the exit code (2 = a check failed).

#include <vector>
#include "IO/PmxLoader.hpp"
#include "Render/Camera.hpp"
#include "Render/Film.hpp"
#include "Scene/Scene.hpp"
#include "Cli/Options.hpp"

namespace rayt::cli {

// One object of the furnace test with the mean color it must reach.
//   expected: 1.0 for lossless materials. Models that lose energy by construction use the value
//   of an independent numerical integration (conductors) or a measured regression value (frosted
//   glass, which has no closed form here).
struct FurnaceObject {
    std::string label;
    const Material* material;   // nullptr = the background (environment)
    Real expected;
    Real tolerance;             // allowed |mean - expected| per color channel
};

// Mean film color per object seen through the pixel centers (furnace test report).
// Returns false when any object is outside its tolerance (the caller turns that into an exit code).
bool CheckFurnace(const Scene& scene, const Camera& camera, const Film& film, const std::vector<FurnaceObject>& objects);

int RunLightTests(const CliOptions& opt);

// -----------------------------------------------------------------------------
// AOV tests (--scene aov): the albedo / normal features recorded for the denoiser
// -----------------------------------------------------------------------------
int RunAovTests();

// -----------------------------------------------------------------------------
// Thin-film interference tests (--scene thinfilm): Airy's formula against closed forms, sample() weights, furnace
// -----------------------------------------------------------------------------
int RunThinFilmTests();

// -----------------------------------------------------------------------------
// Dispersion tests (--scene dispersion): the Sellmeier presets (n_d, Abbe number) and a diamond ball in a white furnace
// -----------------------------------------------------------------------------
int RunDispersionTests();

// -----------------------------------------------------------------------------
// Subsurface scattering tests (--scene sss)
// -----------------------------------------------------------------------------
int RunSssTests(const CliOptions& opt);

int RunVolumeTests(const CliOptions& opt);

// Layered coat checks (--scene layered): closed form of a smooth coat over a diffuse base, eval() against sample(), energy, reciprocity.
int RunLayeredTests(const CliOptions& opt);

// Anisotropy checks (--scene aniso): the tangent of the shapes, the anisotropic GGX, sample() against eval(), the rotation of the axes.
int RunAnisotropyTests(const CliOptions& opt);

// Ocean checks (--scene ocean): the wave spectrum and the set of waves drawn from it, the dispersion relation, the mesh.
int RunOceanTests(const CliOptions& opt);

// -----------------------------------------------------------------------------
// Caustic photon map test (--scene caustic)
// -----------------------------------------------------------------------------
// A glass sphere (radius 5 cm) on a diffuse floor under a uniform white sky. The floor radiance at a few points near the sphere is
// estimated twice: by plain path tracing (caustics found by BSDF sampling: unbiased, a reference because a uniform sky is easy
// to hit) and by path tracing + caustic photons (L S+ D from the photon map, the same paths dropped from the BSDF sampling).
// The two must agree. Photons emitted = --spp x 8000 (about 2 M by default).
int RunCausticTests(const CliOptions& opt);

// --scene lvc: VolumeBPT (bidirectional path tracing with a cache of light paths, in a gray fog) against plain path tracing with next event
// estimation in the medium. The same small image is rendered by both; the strategies alone (next event estimation, light tracing) and all
// together must agree with the reference in the total and block by block.
int RunLvcTests(const CliOptions& opt);

// --scene manifoldtest: the caustic of a point light through the prism, computed by the specular manifold NEE and by the caustic photon
// map (the existing, independent estimator), at a few lit floor points. The difference to the plain render (direct light only) is the
// caustic part of each; the two must agree (the photon map blurs, so points in smooth parts of the pattern are used).
int RunManifoldTests(const CliOptions& opt, bool cup);

int RunMetalTests(const CliOptions& opt);

bool RunMeshChecks(const CliOptions& opt);

/// Skeleton summary and self-tests of the posing code (--pmx-check).
bool RunPmxChecks(const io::PmxModel& model);
} // namespace rayt::cli
