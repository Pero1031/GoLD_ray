#pragma once
// How the CLI describes a demo scene.
//
//   SceneContext  what every scene builder gets (options, the shared materials, the environment)
//   SceneSetup    what a scene builder hands back (objects, lights, camera, caustic target, hooks after rendering)
//   SceneDef      one line of the scene table: name, picture format, default depth, kind of environment, builder
//
// A scene builder only places objects and decides how the scene is seen. The common steps (registering lights, caustics,
// the integrator, rendering, saving) are in RenderPipeline.hpp. The pipeline does not know the scenes by name: what a scene
// needs from it (photon counts, the fog's lighting) is in its SceneSetup.

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "Core/Sampling.hpp"
#include "Color/ComplexIorTable.hpp"
#include "IO/EnvMap.hpp"
#include "Lights/Light.hpp"
#include "Materials/Dielectric.hpp"
#include "Materials/Lambertian.hpp"
#include "Materials/RoughConductor.hpp"
#include "Render/Camera.hpp"
#include "Render/Film.hpp"
#include "Render/Integrator.hpp"
#include "Render/SpecularNee.hpp"
#include "Scene/Aggregate.hpp"
#include "Scene/Scene.hpp"
#include "Cli/Options.hpp"
#include "Cli/SceneDescription.hpp"

namespace rayt::cli {

// Materials and data that several scenes use (built once, see MakeStandardMaterials).
struct StdMaterials {
    std::shared_ptr<Lambertian> floor;
    std::shared_ptr<color::ComplexIorTable> auTable;                              // Johnson & Christy gold
    std::shared_ptr<RoughConductor> goldSmooth, goldMedium, goldRough;            // roughness 0.01 / 0.20 / 0.50
    std::shared_ptr<Dielectric> glass, glassDisp, frosted;                        // clear / dispersive (--abbe) / frosted (0.2)
    std::shared_ptr<Lambertian> white;                                            // albedo 1 (furnace reference)
};

// The environment map and its light (either may be null: black background).
struct Environment {
    std::shared_ptr<EnvMap> map;
    std::shared_ptr<EnvLight> light;
    Vector3 peakDir = Vector3(0, 1, 0);   // direction of the brightest spot (water / sea scenes follow it)
};

enum class EnvKind {
    None,          // black background
    Hdr,           // the --env map
    HdrWithPeak,   // the --env map, and the direction of its brightest spot is searched (water / sea scenes)
    Uniform,       // uniform white (furnace)
    ProceduralSky  // blue gradient and a low sun (sea scene)
};

struct SceneContext {
    const CliOptions& opt;
    const StdMaterials& mats;
    Environment& env;
};

// Bounding sphere of the specular objects: where the caustic photons are aimed.
struct CausticTarget {
    Point3 center = Point3(0.0);
    Real radius = 0;
};

// How a scene uses the caustic photon map and the specular manifold NEE. The scene states its defaults here; the pipeline
// applies them together with the options (--caustics, --caustic-radius, --caustic-method, --caustic-trace).
struct CausticSetup {
    bool photonMap = false;                      // the scene has a caustic photon map (--caustics N)
    long long defaultPhotons = -1;               // photons when --caustics is not given (-1: no photon map unless it is asked for)
    long long defaultPhotonsWithManifold = -1;   // the same while the manifold NEE solves the T T chains (-1: defaultPhotons)
    bool photonMethodNeedsPhotons = false;       // --caustic-method photons with --caustics 0 still traces defaultPhotons (the map is then the only source of the caustic)
    Real defaultRadius = Real(0.0015);           // gather radius (m) when --caustic-radius is not given
    std::string defaultMethod = "photons";       // the caustic method when --caustic-method is not given
    std::string manifoldUnavailable;             // why the manifold NEE cannot be used in this set-up (empty: it can); the photon map is used instead
    bool lightTraceOnly = false;                 // the caustics exist only as light-traced photons: --caustic-trace is always on
    bool throughMedia = false;                   // the photons go on through the fog and are connected to the camera at every collision
};

// Everything the pipeline hands to the hooks after rendering.
struct RenderResult {
    const Scene& scene;
    const Camera& camera;
    const Film& film;
    PathIntegrator& integrator;
};

struct SceneSetup {
    std::shared_ptr<Aggregate> world = std::make_shared<Aggregate>();
    std::vector<std::shared_ptr<Light>> lights;                  // lights to register for next-event estimation
    CameraSpec camera;
    std::optional<Vector3> background;                           // color seen by camera rays that hit nothing (instead of the environment)
    std::string lightSampler;                                    // how next-event estimation picks a light, when the scene needs a specific one (empty = --light-sampler)
    CausticTarget causticTarget;
    CausticSetup caustics;                                       // the scene's defaults for the photon map and the manifold NEE
    bool mediaNee = false;                                       // the fog is lit by next event estimation at every collision, and its boundary does not block shadow rays (as --media-nee)
    const PointLight* bptLight = nullptr;                        // fogspot: the lamp and the fog, for --lvc (VolumeBPT)
    const Medium* bptMedium = nullptr;
    Real bptAlbedo = Real(0.5);                                  // albedo of the diffuse surfaces (all the same)
    std::vector<SpecularNee::Body> manifoldBodies;               // flat-faced glass bodies the specular manifold NEE can solve
    std::vector<std::shared_ptr<void>> keepAlive;                // owners of things the scene only points to
    std::function<void(const RenderResult&)> afterRender;        // right after the render, before saving
    std::function<bool(const RenderResult&)> afterSave;          // after saving; false = a check failed (exit code 2)
    bool checksPassed = true;                                    // self-checks that run while the scene is built (mesh)
};

// A scene builder returns std::nullopt when the scene was built, or an exit code when the program should stop (an error, or an
// information request such as --pmx-info).
using SceneBuilder = std::optional<int> (*)(SceneContext& ctx, SceneSetup& out);

struct SceneDef {
    const char* name;
    const char* description;                   // one line, shown by --list-scenes
    double heightPerWidth;                     // picture format: height / width
    int defaultMaxDepth;                       // used when --depth is not given (0 = the global default)
    EnvKind (*envKind)(const CliOptions&);
    SceneBuilder build;
};

// The scene table (null for a name that is not a scene, e.g. a self-check).
const SceneDef* FindScene(const std::string& name);

// Every scene, in table order.
std::vector<const SceneDef*> AllScenes();

// The materials shared by the scenes (prints the gold table it loaded).
StdMaterials MakeStandardMaterials(const CliOptions& opt);

// Loads or builds the environment the scene wants.
Environment MakeEnvironment(const CliOptions& opt, EnvKind kind);

} // namespace rayt::cli
