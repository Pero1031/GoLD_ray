// GoLD_rayt.cpp 
//
/// @brief 金属をJohnsonの物理データを基にレンダリング
///
/// 

// C++20を仮定している
// Precompiled Header (Must be first)
#include "pch.h"    

#define GLFW_INCLUDE_NONE
#include <glad/glad.h>
#include <GLFW/glfw3.h>

//#include <imgui.h>
//#include <backends/imgui_impl_glfw.h>
//#include <backends/imgui_impl_opengl3.h>

#include "UI/ImGuiLayer.hpp"
#include "Core/RenderSettings.hpp"

// 自作ヘッダー
// Project Headers
// Core
#include "Core/Core.hpp"
#include "Core/Interaction.hpp"
#include "Core/Sampling.hpp"
#include "Core/Image.hpp"

// accelerates
#include "Accelerators/BVH.hpp" 

// Geometry
#include "Scene/Hittable.hpp"
#include "Scene/Aggregate.hpp"
#include "Geometry/Sphere.hpp"
#include "Geometry/Quad.hpp"
#include "Geometry/TriangleMesh.hpp"
#include "Accelerators/EmbreeMesh.hpp"
#include "IO/ObjLoader.hpp"
#include "IO/PmxLoader.hpp"
#include "IO/PmxPose.hpp"
#include "IO/PmxMorph.hpp"
#include <glm/gtc/matrix_transform.hpp>
#include "Geometry/Frame.hpp"
#include "Scene/Primitive.hpp"
#include "Scene/Scene.hpp"

// Renderer
#include "Render/Film.hpp"
#include "Render/Camera.hpp"
#include "Scene/Scene.hpp"
#include "Render/Integrator.hpp"
#include "Render/Denoiser.hpp"

// Materials
#include "Materials/Material.hpp"
#include "Materials/Lambertian.hpp"
#include "Materials/CoatedDiffuse.hpp"
#include "Materials/DiffuseLight.hpp"
#include "Materials/RoughConductor.hpp"
#include "Materials/Dielectric.hpp"

// Lights
#include "Lights/AreaLight.hpp"
#include "Lights/PointLight.hpp"
#include "Color/EmissionSpectrum.hpp"
#include "Color/ConductorFromReflectance.hpp"
#include "Core/Fresnel.hpp"

// IO
#include "IO/ImageLoader.hpp"
#include "IO/EnvMap.hpp"
#include "IO/ComplexIorLoader.hpp"
#include "IO/SpectralCurveLoader.hpp"
#include "Color/IorRgbApprox.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <functional>
#include <map>
#include <stdexcept>
#include <vector>

#include "DebugTools/FrameDebug.hpp"


// 画像生成のためのヘッダー
// マクロを書く必要はなし
#include "stb_image_write.h"

using namespace rayt;

// -----------------------------------------------------------------------------
// Scene Configuration
// -----------------------------------------------------------------------------
const int IMAGE_WIDTH = 800;
const int IMAGE_HEIGHT = 450;      // 16:9 Aspect Ratio
const int SAMPLES_PER_PIXEL = 256; // Higher = less noise, slower
const int MAX_DEPTH = 16;          // Max recursion depth for rays

// env path
const std::string ENV_HDR_PATH = "assets/env/bryanston_park_sunrise_2k.hdr";

const std::string AU_IOR_CSV_PATH = "assets/Metal_data/Johnson_Au.csv";

// -----------------------------------------------------------------------------
// Command line
//   main_cli [--scene gold|glass|furnace] [--width N] [--spp N] [--depth N]
//            [--env file.hdr] [--abbe V] [--out file.png] [--trace x,y]
//
//   --trace x,y : debug aid. Prints a few light paths of the pixel (x,y) (image coordinates,
//                 origin at the top-left) with the object hit at every bounce, then exits.
//
//   gold    : three gold spheres (default, as before)
//   glass   : clear / dispersive / frosted glass spheres on a floor
//   furnace : objects in a uniform white environment. The mean color of every object must match
//             its expected value (1.0 for lossless materials; a value below/above means energy is
//             lost/created, i.e. a bug). Prints OK/FAIL per object and exits with code 2 on FAIL.
//             Debug builds also assert that all ray directions have unit length.
//   spd     : colored diffuse spheres + gold under ONE sphere light whose spectrum is chosen with
//             --light (d65 | tungsten = illuminant A | candle = 1850 K blackbody | led = synthetic
//             SPD loaded from assets/Illuminants/*.csv); --wb white-balances the film to that light.
//             No environment light (black background). Spectral mode shows metamerism best.
//   pmx     : an MMD model from MMD_model/ (--pmx-list, --pmx-index N, --pmx-info, --pmx-yaw deg). Diffuse textures (uv, repeat, bilinear)
//             and alpha cut-outs; "+" overlay materials are skipped unless --pmx-all. HDR environment lighting.
//   goldrough: gold spheres of roughness 0.4 / 0.7 / 1.0 under the HDR environment; --no-ms = without the
//             multiple-scattering compensation (to compare).
//   studio  : gold / dispersive glass / red diffuse sphere lit by a tungsten rectangular softbox
//             and a small cool-white LED sphere light (SPDs from EmissionSpectrum / CSV). No environment.
//   mesh    : triangle meshes: runs mesh checks first (geometry, OBJ, BVH vs brute force, furnace,
//             emissive triangles), then renders a tilted gold torus (assets/Meshes/torus.obj) and a
//             flat- vs smooth-shaded low-poly sphere pair under the studio lights.
//   bench   : BVH benchmark (recursive reference vs SAH flat BVH) on icospheres of 20k / 82k / 328k triangles;
//             rays = --spp x 2000. Use a Release build for meaningful speeds.
//   metals  : metal spectral data tests (no image): file layouts / units, reflectance input,
//             edge-tint accuracy against the Johnson tables, rendered gold. Samples = --spp x 1024.
//   lights  : light source tests (no image). A Lambertian ground is lit by one light and the
//             radiance of one ray is compared with a closed-form / independent numerical value
//             (point light, sphere area light with and without NEE, colored emission, environment
//             maps). Samples per test = --spp x 4096. Exits with code 2 on FAIL.
//   aov     : denoiser feature tests (no image): the albedo / normal recorded for a camera path on
//             a diffuse surface, a coated sphere, a light, through smooth glass and for escaping rays.
// Without options the behavior is the same as before (gold scene, constants above).
// -----------------------------------------------------------------------------
struct CliOptions {
    std::string scene = "gold";
    int width = IMAGE_WIDTH;           // height follows the 16:9 aspect ratio of the constants
    int spp = SAMPLES_PER_PIXEL;
    int maxDepth = MAX_DEPTH;
    std::string env = ENV_HDR_PATH;
    Real abbe = 20.0;                  // Abbe number of the dispersive glass (smaller = stronger)
    std::string out;                   // default: result_<scene>.png ("gold" keeps result_gold_pbr.png)
    int pmxIndex = -1;                     // --pmx-index N (see --pmx-list); -1 = the largest .pmx
    bool pmxList = false;                  // --pmx-list : print the .pmx files found in MMD_model/ and exit
    Real pmxZoom = 1;                      // --pmx-zoom Z : camera distance / Z
    Real pmxFocus = 0.5;                   // --pmx-focus f : look at this fraction of the model height (0.9 = face)
    int poseIndex = -1;                    // --pose-index N (see --pose-list): apply a .vpd pose from MMD_model/pose
    bool poseList = false;                 // --pose-list : print the .vpd files found and exit
    bool morphList = false;                // --morph-list : print the model's morphs and exit
    std::string morphSpec;                 // --morph "12=1.0,15=0.5" : morph weights by index
    std::string exprFile;                  // --expr file.txt : "name = weight" lines (UTF-8) from MMD_model/expression
    bool pmxCheck = false;                 // --pmx-check : skeleton summary + pose self-tests, then exit
    std::string studioLook = "default";    // --studio-look default|warm|cool|dramatic|soft : lighting of --pmx-studio
    bool aov = false;                      // --aov : also write <out>_albedo.png and <out>_normal.png (denoiser features)
    bool denoise = false;                  // --denoise : also write <out>_denoised.* (Intel OIDN, uses albedo + normal)
    bool cleanAux = false;                 // --clean-aux : tell OIDN the albedo / normal are noise-free (skips its prefilter)
    bool pmxStudio = false;                // --pmx-studio : studio lighting (key/fill/rim soft boxes, dim backdrop), framing by pose width
    bool pmxLight = false;                 // --pmx-light : add a white softbox in front-left (highlights show the gloss)
    bool pmxMatte = false;                 // --pmx-matte : every material plain diffuse (no gloss presets)
    bool pmxAll = false;                   // --pmx-all : also draw the "+" overlay materials
    bool pmxInfoOnly = false;              // --pmx-info : load, print the model summary, exit (no render)
    Real pmxYaw = 0;                       // --pmx-yaw deg : turn the model about +Y
    bool noMs = false;                     // --no-ms : single-scattering GGX only (goldrough scene)
    std::string accel = "own";             // own | embree (triangle mesh intersection)
    std::string lightSampler = "spatial";  // uniform | power | spatial (light selection for NEE)
    std::string light = "tungsten";    // spd scene: d65 | tungsten | candle | led
    bool whiteBalance = false;         // spd scene: white-balance the film to the light
    int traceX = -1, traceY = -1;      // --trace x,y : print a few light paths of one pixel and exit
};

static bool ParseArgs(int argc, char** argv, CliOptions& o) {
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            auto value = [&](const char* name) -> std::string {
                if (i + 1 >= argc) throw std::runtime_error(std::string("missing value for ") + name);
                return argv[++i];
            };
            if (a == "--scene")      o.scene = value("--scene");
            else if (a == "--width") o.width = std::stoi(value("--width"));
            else if (a == "--spp")   o.spp = std::stoi(value("--spp"));
            else if (a == "--depth") o.maxDepth = std::stoi(value("--depth"));
            else if (a == "--env")   o.env = value("--env");
            else if (a == "--abbe")  o.abbe = std::stod(value("--abbe"));
            else if (a == "--out")   o.out = value("--out");
            else if (a == "--light") o.light = value("--light");
            else if (a == "--light-sampler") o.lightSampler = value("--light-sampler");
            else if (a == "--accel") o.accel = value("--accel");
            else if (a == "--no-ms") o.noMs = true;
            else if (a == "--pmx-index") o.pmxIndex = std::stoi(value("--pmx-index"));
            else if (a == "--pmx-list") o.pmxList = true;
            else if (a == "--pmx-info") o.pmxInfoOnly = true;
            else if (a == "--pmx-all") o.pmxAll = true;
            else if (a == "--pmx-matte") o.pmxMatte = true;
            else if (a == "--pmx-light") o.pmxLight = true;
            else if (a == "--pmx-studio") o.pmxStudio = true;
            else if (a == "--studio-look") o.studioLook = value("--studio-look");
            else if (a == "--aov") o.aov = true;
            else if (a == "--denoise") o.denoise = true;
            else if (a == "--clean-aux") o.cleanAux = true;
            else if (a == "--pose-index") o.poseIndex = std::stoi(value("--pose-index"));
            else if (a == "--pose-list") o.poseList = true;
            else if (a == "--pmx-check") o.pmxCheck = true;
            else if (a == "--morph-list") o.morphList = true;
            else if (a == "--morph") o.morphSpec = value("--morph");
            else if (a == "--expr") o.exprFile = value("--expr");
            else if (a == "--pmx-zoom") o.pmxZoom = std::stod(value("--pmx-zoom"));
            else if (a == "--pmx-focus") o.pmxFocus = std::stod(value("--pmx-focus"));
            else if (a == "--pmx-yaw") o.pmxYaw = std::stod(value("--pmx-yaw"));
            else if (a == "--wb")    o.whiteBalance = true;
            else if (a == "--trace") {
                const std::string xy = value("--trace");
                const size_t comma = xy.find(',');
                if (comma == std::string::npos) throw std::runtime_error("--trace needs x,y");
                o.traceX = std::stoi(xy.substr(0, comma));
                o.traceY = std::stoi(xy.substr(comma + 1));
            }
            else throw std::runtime_error("unknown option: " + a);
        }
        if (o.scene != "gold" && o.scene != "glass" && o.scene != "furnace" && o.scene != "lights" && o.scene != "aov" && o.scene != "spd" && o.scene != "metals" && o.scene != "studio" && o.scene != "mesh" && o.scene != "bench" && o.scene != "goldrough" && o.scene != "pmx")
            throw std::runtime_error("--scene must be gold, glass, furnace, lights, aov, spd, metals, studio, mesh, bench, goldrough or pmx");
        if (o.accel != "own" && o.accel != "embree")
            throw std::runtime_error("--accel must be own or embree");
        if (o.studioLook != "default" && o.studioLook != "warm" && o.studioLook != "cool" && o.studioLook != "dramatic" && o.studioLook != "soft")
            throw std::runtime_error("--studio-look must be default, warm, cool, dramatic or soft");
        if (o.lightSampler != "uniform" && o.lightSampler != "power" && o.lightSampler != "spatial" && o.lightSampler != "all")
            throw std::runtime_error("--light-sampler must be uniform, power, spatial or all");
        if (o.light != "d65" && o.light != "tungsten" && o.light != "candle" && o.light != "led")
            throw std::runtime_error("--light must be d65, tungsten, candle or led");
        if (o.width < 16 || o.spp < 1 || o.maxDepth < 1)
            throw std::runtime_error("--width must be >= 16, --spp and --depth >= 1");
    }
    catch (const std::exception& e) {
        std::cerr << "[Args] " << e.what() << "\n"
            << "usage: main_cli [--scene gold|glass|furnace|lights|aov|spd|metals|studio|mesh|bench|goldrough|pmx] [--light d65|tungsten|candle|led] [--wb] [--light-sampler uniform|power|spatial|all] [--accel own|embree] [--no-ms] [--pmx-list] [--pmx-index N] [--pmx-info] [--pmx-all] [--pmx-matte] [--pmx-light] [--pmx-studio] [--studio-look default|warm|cool|dramatic|soft] [--aov] [--denoise] [--clean-aux] [--pose-list] [--pose-index N] [--pmx-check] [--morph-list] [--morph i=w,..] [--expr file] [--pmx-yaw deg] [--pmx-zoom Z] [--pmx-focus f] [--width N] [--spp N] [--depth N]"
            << " [--env file.hdr] [--abbe V] [--out file.png] [--trace x,y]\n";
        return false;
    }
    if (o.out.empty()) {
        o.out = (o.scene == "gold") ? "result_gold_pbr.png" : "result_" + o.scene + ".png";
        if (o.scene == "spd") o.out = "result_spd_" + o.light + (o.whiteBalance ? "_wb" : "") + ".png";
    }
    return true;
}

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
static bool CheckFurnace(const Scene& scene, const Camera& camera, const Film& film,
    const std::vector<FurnaceObject>& objects)
{
    struct Acc { Vector3 sum = Vector3(0.0); long long n = 0; };
    std::map<const Material*, Acc> acc;   // key nullptr = background (environment)

    const int w = film.width();
    const int h = film.height();
    for (int j = 0; j < h; ++j) {
        for (int i = 0; i < w; ++i) {
            // Same pixel convention as PathIntegrator::render (film row = h - 1 - j).
            const Ray r = camera.getRay((Real(i) + 0.5) / Real(w), (Real(j) + 0.5) / Real(h), Point2(0.5, 0.5));
            SurfaceInteraction rec;
            const Material* m = scene.hit(r, rec) ? rec.matPtr : nullptr;
            Acc& a = acc[m];
            a.sum += film.getData()[static_cast<size_t>(h - 1 - j) * static_cast<size_t>(w) + static_cast<size_t>(i)];
            a.n++;
        }
    }

    std::cout << "\n[Furnace] mean linear RGB per object vs the expected value (+- tolerance, per channel)\n";
    bool allPassed = true;
    for (const FurnaceObject& obj : objects) {
        const auto it = acc.find(obj.material);
        if (it == acc.end() || it->second.n == 0) {
            std::cout << "[Furnace]   " << obj.label << ": no pixels found -> FAIL\n";
            allPassed = false;
            continue;
        }
        const Vector3 mean = it->second.sum / Real(it->second.n);
        const bool ok = std::abs(mean.x - obj.expected) <= obj.tolerance
            && std::abs(mean.y - obj.expected) <= obj.tolerance
            && std::abs(mean.z - obj.expected) <= obj.tolerance;
        allPassed = allPassed && ok;
        std::cout << "[Furnace]   " << obj.label << ": (" << mean.x << ", " << mean.y << ", " << mean.z
            << ")  expected " << obj.expected << " +- " << obj.tolerance
            << "  [" << it->second.n << " px]  " << (ok ? "OK" : "FAIL") << "\n";
    }
    std::cout << "[Furnace] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed;
}

// Debug aid: follow a few light paths of one pixel and print which object each bounce hits.
//   [in]/[out] = the ray hits the surface from outside / from inside, (R)/(T) = reflected / transmitted.
static void TracePixelDebug(const Scene& scene, const Camera& camera, int px, int py, int w, int h,
    int maxDepth, int numPaths,
    const std::vector<std::pair<std::string, const Material*>>& names)
{
    auto label = [&](const Material* m) {
        for (const auto& [n, p] : names) if (p == m) return n;
        return std::string("?");
    };
    const SampledWavelengths lambda;   // RGB-mode default (wavelength-independent materials)
    std::cout << "[Trace] pixel (" << px << ", " << py << ") of " << w << "x" << h << ", " << numPaths << " paths\n";
    for (int k = 0; k < numPaths; ++k) {
        // Same convention as PathIntegrator: image row 0 is the top, camera t = 1 at the top.
        Ray r = camera.getRay((Real(px) + 0.5) / Real(w), Real(1) - (Real(py) + 0.5) / Real(h), Point2(0.5, 0.5));
        std::cout << "  path " << k << ":";
        for (int depth = 0; depth < maxDepth; ++depth) {
            SurfaceInteraction rec;
            if (!scene.hit(r, rec)) {
                std::cout << " -> ENV [ray o=(" << r.o.x << "," << r.o.y << "," << r.o.z << ") d=("
                    << r.d.x << "," << r.d.y << "," << r.d.z << ") tMin=" << r.tMin << " tMax=" << r.tMax << "]";
                break;
            }
            std::cout << " -> " << label(rec.matPtr) << (rec.frontFace ? "[in]" : "[out]")
                << "@(" << rec.p.x << "," << rec.p.y << "," << rec.p.z << ")";
            const auto s = rec.matPtr->sample(rec, -r.d, Point2(sampling::Random(), sampling::Random()), lambda);
            if (!s) { std::cout << " (no sample: " << (rec.frontFace ? "?" : "TIR") << ")"; break; }
            std::cout << (s->isTransmission() ? "(T)" : "(R)");
            r = rec.SpawnRay(s->wi);
        }
        std::cout << "\n";
    }
}

// -----------------------------------------------------------------------------
// Light tests (--scene lights)
// -----------------------------------------------------------------------------

// Mean linear RGB (and its standard error) of Li() along one fixed ray.
static void EstimateRadiance(const Scene& scene, const Ray& ray, int n, int maxDepth,
    Vector3& mean, Vector3& se, int& badSamples)
{
    PathIntegrator integrator(nullptr, maxDepth, 1);
    double sx = 0, sy = 0, sz = 0, qx = 0, qy = 0, qz = 0;
    int bad = 0;

    // OpenMP 2.0 (MSVC): scalar reductions and a signed int loop counter only.
#pragma omp parallel for reduction(+:sx,sy,sz,qx,qy,qz,bad)
    for (int i = 0; i < n; ++i) {
#if defined(RAYT_SPECTRAL_MODE)
        SampledWavelengths lambda = SampledWavelengths::sampleVisible(sampling::Random());
#else
        SampledWavelengths lambda;
#endif
        const Spectrum L = integrator.Li(ray, scene, lambda);
        const Vector3 c = SpectrumToRGB(L, lambda);
        if (math::hasNonFinite(c)) { ++bad; continue; }
        sx += c.x; sy += c.y; sz += c.z;
        qx += c.x * c.x; qy += c.y * c.y; qz += c.z * c.z;
    }

    const double inv = 1.0 / double(n);
    mean = Vector3(sx, sy, sz) * inv;
    const Vector3 var(std::max(0.0, qx * inv - mean.x * mean.x),
        std::max(0.0, qy * inv - mean.y * mean.y),
        std::max(0.0, qz * inv - mean.z * mean.z));
    se = Vector3(std::sqrt(var.x * inv), std::sqrt(var.y * inv), std::sqrt(var.z * inv));
    badSamples = bad;
}

// Passes when every channel is within 4 standard errors (+0.5 % of the value) of the expectation.
static bool ReportLightCase(const std::string& name, const Vector3& expected,
    const Vector3& mean, const Vector3& se, int bad)
{
    bool ok = (bad == 0);
    for (int c = 0; c < 3; ++c) {
        ok = ok && std::abs(mean[c] - expected[c]) <= 4.0 * se[c] + 0.005 * std::abs(expected[c]) + 1e-6;
    }
    std::cout << "[Lights]   " << name << "\n"
        << "             expected (" << expected.x << ", " << expected.y << ", " << expected.z << ")"
        << "   got (" << mean.x << ", " << mean.y << ", " << mean.z << ") +- ("
        << se.x << ", " << se.y << ", " << se.z << ")"
        << "   ratio (" << mean.x / expected.x << ", " << mean.y / expected.y << ", " << mean.z / expected.z << ")"
        << (bad ? "  NaN samples: " + std::to_string(bad) : std::string())
        << "   " << (ok ? "OK" : "FAIL") << "\n";
    return ok;
}

// Lambertian ground (a huge sphere, locally flat) with whatever `addLights` puts above it.
static std::shared_ptr<Scene> MakeGroundScene(Real albedo,
    const std::function<void(Aggregate&, Scene&)>& addLights, const char* lightSampler = "spatial")
{
    auto world = std::make_shared<Aggregate>();
    world->add(std::make_shared<Sphere>(Point3(0, -1000, 0), 1000.0,
        std::make_shared<Lambertian>(RgbColor(albedo))));
    auto scene = std::make_shared<Scene>(world);
    addLights(*world, *scene);
    scene->finalize(lightSampler);
    return scene;
}

// Spherical emitter. registerAsLight = false leaves it out of the light list (no NEE: pure BSDF sampling).
static void AddSphereLight(Aggregate& world, Scene& scene, const Point3& center, Real radius,
    const EmissionSpectrum& color, Real intensity, bool registerAsLight)
{
    auto mat = std::make_shared<DiffuseLight>(color, intensity);
    auto shape = std::make_shared<Sphere>(center, radius, mat);
    if (registerAsLight) {
        auto light = std::make_shared<AreaLight>(shape, mat);
        world.add(std::make_shared<Primitive>(shape, mat, light));
        scene.addLight(light);
    }
    else {
        world.add(std::make_shared<Primitive>(shape, mat));
    }
}

// Rectangular emitter (corner + two edges; emits toward cross(edge1, edge2)).
static std::shared_ptr<Light> AddQuadLight(Aggregate& world, Scene* scene, const Point3& corner,
    const Vector3& e1, const Vector3& e2, const EmissionSpectrum& color, Real intensity,
    bool registerAsLight, bool twoSided = false)
{
    auto mat = std::make_shared<DiffuseLight>(color, intensity, twoSided);
    auto shape = std::make_shared<Quad>(corner, e1, e2, mat);
    std::shared_ptr<Light> light;
    if (registerAsLight) {
        light = std::make_shared<AreaLight>(shape, mat);
        if (scene) scene->addLight(light);
    }
    world.add(std::make_shared<Primitive>(shape, mat, light));
    return light;
}

// Independent reference for a quad light: radiance of the ground point at the origin (normal +y)
//   Lo = albedo / pi * L * Integral cos(receiver) cos(light) / r^2 dA     (midpoint rule on the quad)
static Real QuadGroundRadiance(const Point3& corner, const Vector3& e1, const Vector3& e2,
    bool twoSided, Real L, Real albedo)
{
    const Vector3 nn = glm::cross(e1, e2);
    const Real area = glm::length(nn);
    const Vector3 n = nn / area;
    const int N = 600;
    Real E = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            const Point3 q = corner + e1 * ((Real(i) + 0.5) / N) + e2 * ((Real(j) + 0.5) / N);
            const Vector3 d = q;                         // from the origin to the light point
            const Real r2 = glm::dot(d, d), r = std::sqrt(r2);
            const Real cosR = d.y / r;                   // receiver normal +y
            const Real cosL = glm::dot(n, -d) / r;       // > 0: the emitting side faces the receiver
            if (cosR <= 0 || (cosL <= 0 && !twoSided)) continue;
            E += cosR * std::abs(cosL) / r2;
        }
    }
    E *= L * area / Real(N * N);
    return albedo / constants::PI * E;
}

// Environment map: dim bluish sky that brightens upward + a small bright sun (elevation 50 deg).
static std::shared_ptr<EnvMap> MakeSunSkyEnv(Vector3& sunDirOut) {
    const int W = 256, H = 128;
    const Real elev = glm::radians(50.0), azim = glm::radians(30.0);
    const Vector3 sun(std::cos(elev) * std::cos(azim), std::sin(elev), std::cos(elev) * std::sin(azim));
    sunDirOut = sun;

    std::vector<Vector3> px(static_cast<size_t>(W) * H);
    for (int y = 0; y < H; ++y) {
        const Real theta = constants::PI * (Real(y) + 0.5) / Real(H);     // row 0 = zenith
        for (int x = 0; x < W; ++x) {
            const Real phi = 2.0 * constants::PI * (Real(x) + 0.5) / Real(W) - constants::PI;
            const Vector3 d(std::sin(theta) * std::cos(phi), std::cos(theta), std::sin(theta) * std::sin(phi));
            Vector3 c = Vector3(0.8, 0.9, 1.0) * (0.2 + 0.3 * std::max(d.y, 0.0));
            if (std::acos(std::clamp(glm::dot(d, sun), -1.0, 1.0)) < 0.05) c += Vector3(500.0, 450.0, 400.0);
            px[static_cast<size_t>(y) * W + x] = c;
        }
    }
    return std::make_shared<EnvMap>(Image(W, H, std::move(px)));
}

// Independent reference: radiance leaving a white-ish Lambertian ground point (normal +y) under
// an environment map = albedo / pi * Integral_{upper hemisphere} Le(w) cos(theta) dw, evaluated
// on a dense (theta, phi) grid with EnvMap::eval (no importance sampling involved).
static Vector3 GroundRadianceFromEnv(const EnvMap& env, Real albedo) {
    const int nTheta = 1024, nPhi = 2048;
    Vector3 E(0.0);
    for (int i = 0; i < nTheta; ++i) {
        const Real theta = 0.5 * constants::PI * (Real(i) + 0.5) / Real(nTheta);
        const Real dTheta = 0.5 * constants::PI / Real(nTheta);
        const Real dPhi = 2.0 * constants::PI / Real(nPhi);
        for (int j = 0; j < nPhi; ++j) {
            const Real phi = 2.0 * constants::PI * (Real(j) + 0.5) / Real(nPhi);
            const Vector3 d(std::sin(theta) * std::cos(phi), std::cos(theta), std::sin(theta) * std::sin(phi));
            E += env.eval(d) * (std::cos(theta) * std::sin(theta) * dTheta * dPhi);
        }
    }
    return E * (albedo / constants::PI);
}

// Consistency of EnvMap::sample()/pdf(): the pdf integrates to 1, sample() and pdf() agree, and
// Le/pdf estimates the same integral as the grid.
static bool CheckEnvSampling(const EnvMap& env, const std::string& label, int n) {
    // Integral of pdf over the sphere on a deterministic (theta, phi) grid (a small bright sun
    // would make a Monte Carlo estimate of this far too noisy).
    double pdfIntegral = 0.0;
    Vector3 gridLe(0.0);   // Integral Le dw by the same grid
    {
        const int nTheta = 1024, nPhi = 2048;
        const double dTheta = constants::PI / nTheta, dPhi = 2.0 * constants::PI / nPhi;
        for (int i = 0; i < nTheta; ++i) {
            const Real theta = (Real(i) + 0.5) * dTheta;
            for (int j = 0; j < nPhi; ++j) {
                const Real phi = (Real(j) + 0.5) * dPhi;
                const Vector3 d(std::sin(theta) * std::cos(phi), std::cos(theta), std::sin(theta) * std::sin(phi));
                pdfIntegral += env.pdf(d) * std::sin(theta) * dTheta * dPhi;
                gridLe += env.eval(d) * (std::sin(theta) * dTheta * dPhi);
            }
        }
    }

    Vector3 isLe(0.0);
    int mismatches = 0;   // sample() and pdf() disagree (> 0.1 %)
    int isCount = 0;

    for (int i = 0; i < n; ++i) {
        Vector3 wi(0.0);
        Real pdf = 0;
        const Vector3 Le = env.sample(Point2f(float(sampling::Random()), float(sampling::Random())), wi, pdf);
        if (pdf > 0) {
            isLe += Le / pdf;
            if (std::abs(env.pdf(wi) - pdf) > 1e-3 * pdf) ++mismatches;   // only at texel borders (float round-off)
            ++isCount;
        }
    }
    isLe /= double(n);
    const double mismatchRate = double(mismatches) / std::max(isCount, 1);

    const bool okPdf = std::abs(pdfIntegral - 1.0) < 0.01;
    const bool okMatch = mismatchRate < 1e-3;
    const bool okLe = std::abs(isLe.y - gridLe.y) <= 0.03 * std::abs(gridLe.y);
    std::cout << "[Lights]   env sampling (" << label << "): Integral pdf dw = " << pdfIntegral
        << (okPdf ? " OK" : " FAIL") << ";  sample()/pdf() mismatch rate = " << mismatchRate
        << (okMatch ? " OK" : " FAIL") << ";  Integral Le dw (G): grid " << gridLe.y
        << " vs importance " << isLe.y << (okLe ? " OK" : " FAIL") << "\n";
    return okPdf && okMatch && okLe;
}

// -----------------------------------------------------------------------------
// AOV tests (--scene aov): the albedo / normal features recorded for the denoiser
// -----------------------------------------------------------------------------
static int RunAovTests() {
    bool allPassed = true;
    const Vector3 groundAlbedo(0.3, 0.5, 0.7), coatedAlbedo(0.8, 0.2, 0.1), lightColor(1.0, 0.4, 0.1);
    const int maxDepth = 8;

    auto buildScene = [&](bool withSky) {
        auto world = std::make_shared<Aggregate>();
        world->add(std::make_shared<Sphere>(Point3(0, -1000, 0), 1000.0, std::make_shared<Lambertian>(RgbColor(groundAlbedo))));
        world->add(std::make_shared<Sphere>(Point3(-1.5, 0.5, 0), 0.5, std::make_shared<Dielectric>(1.5, 0.0)));
        world->add(std::make_shared<Sphere>(Point3(1.5, 0.5, 0), 0.5,
            std::make_shared<CoatedDiffuse>(RgbColor(coatedAlbedo), Real(0.4))));
        auto scene = std::make_shared<Scene>(world);
        AddSphereLight(*world, *scene, Point3(0, 3, -3), 0.5, EmissionSpectrum(RgbColor(lightColor)), 10.0, true);
        if (withSky) {
            auto sky = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(0.5))));
            scene->setEnvLight(std::make_shared<EnvLight>(sky, 1.0));
        }
        scene->finalize("spatial");
        return scene;
    };
    const auto sceneSky = buildScene(true);
    const auto sceneNoSky = buildScene(false);

    struct Result { int n = 0, valid = 0; Vector3 albedo = Vector3(0.0), normal = Vector3(0.0); int nearGround = 0; };
    auto trace = [&](const Scene& scene, const Ray& ray, int n) {
        Result res;
        PathIntegrator integrator(nullptr, maxDepth, 1);
        for (int i = 0; i < n; ++i) {
#if defined(RAYT_SPECTRAL_MODE)
            SampledWavelengths lambda = SampledWavelengths::sampleVisible(sampling::Random());
#else
            SampledWavelengths lambda;
#endif
            PathAOV aov;
            integrator.Li(ray, scene, lambda, &aov);
            ++res.n;
            if (!aov.valid) continue;
            ++res.valid;
            res.albedo += aov.albedo;
            res.normal += aov.normal;
            if (glm::length(aov.albedo - groundAlbedo) < 1e-6) ++res.nearGround;
        }
        if (res.valid > 0) { res.albedo /= Real(res.valid); res.normal /= Real(res.valid); }
        return res;
    };
    auto report = [&](const std::string& name, bool ok) {
        std::cout << "[Aov]   " << name << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };
    auto close = [](const Vector3& a, const Vector3& b, Real tol) { return glm::length(a - b) <= tol; };

    std::cout << "[Aov] albedo / normal features of a camera path\n";

    {   // diffuse ground (a huge sphere, so its normal is +y only to within 0.002)
        const Result r = trace(*sceneSky, Ray(Point3(0.5, 2, 2), Vector3(0, -1, 0)), 200);
        report("diffuse ground: albedo = material color, normal = +y",
            r.valid == r.n && close(r.albedo, groundAlbedo, 1e-6) && close(r.normal, Vector3(0, 1, 0), 5e-3));
    }
    {   // coated diffuse sphere, hit head-on (normal points back at the ray origin)
        const Result r = trace(*sceneSky, Ray(Point3(1.5, 0.5, 3), Vector3(0, 0, -1)), 200);
        report("coated diffuse sphere: albedo = base color, normal = +z",
            r.valid == r.n && close(r.albedo, coatedAlbedo, 1e-6) && close(r.normal, Vector3(0, 0, 1), 1e-6));
    }
    {   // emitter: the hue of the emission, scaled to a maximum of 1
        const Result r = trace(*sceneSky, Ray(Point3(0, 3, 0), Vector3(0, 0, -1)), 200);
        report("sphere light: albedo = emission hue (max 1), normal = +z",
            r.valid == r.n && close(r.albedo, lightColor, 1e-6) && close(r.normal, Vector3(0, 0, 1), 1e-6));
    }
    {   // smooth glass sphere over the ground: the ray is followed through it, so nearly every path ends
        // with the ground color (the rest are Fresnel reflections that escape to the sky). If the glass
        // itself were recorded, the albedo would be white (the Material default) and this fraction ~0.
        const Result r = trace(*sceneSky, Ray(Point3(-1.5, 3, 0), Vector3(0, -1, 0)), 2000);
        const double groundFraction = double(r.nearGround) / double(r.n);
        std::cout << "[Aov]     (through glass: " << 100.0 * groundFraction << " % of the paths record the ground color)\n";
        // (about 1 path in 2000 bounces around inside the sphere until the depth limit and records nothing)
        report("smooth glass is skipped: albedo is the ground behind it", r.valid >= 0.99 * r.n && groundFraction > 0.85);
    }
    {   // a ray that leaves the scene
        const Result sky = trace(*sceneSky, Ray(Point3(0, 5, 5), Vector3(0, 1, 0)), 2000);
        report("escaping ray with an environment light: albedo = environment color, normal = 0",
            sky.valid == sky.n && close(sky.albedo, Vector3(0.5), 0.03) && close(sky.normal, Vector3(0.0), 1e-9));
        const Result none = trace(*sceneNoSky, Ray(Point3(0, 5, 5), Vector3(0, 1, 0)), 50);
        report("escaping ray without environment: no feature recorded", none.valid == 0);
    }

    std::cout << "[Aov] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

static int RunLightTests(const CliOptions& opt) {
    const int n = opt.spp * 4096;   // samples per test (default --spp 256 -> about 1M)
    const int depth = std::max(opt.maxDepth, 4);
    std::cout << "[Lights] " << n << " samples per test, albedo 0.5 ground, probe ray hits the origin\n";
    bool allPassed = true;
    const Real albedo = 0.5;

    const Ray probe(Point3(-3, 3, 0), glm::normalize(Vector3(1, -1, 0)));
    auto run = [&](const std::string& name, const std::shared_ptr<Scene>& scene, const Vector3& expected) {
        Vector3 mean, se;
        int bad = 0;
        EstimateRadiance(*scene, probe, n, depth, mean, se, bad);
        allPassed = ReportLightCase(name, expected, mean, se, bad) && allPassed;
    };

    // Overhead light at distance d = 2 from the probe point.
    // Point light:   Lo = albedo / pi * I / d^2
    run("point light I=10, d=2",
        MakeGroundScene(albedo, [](Aggregate&, Scene& s) {
            s.addLight(std::make_shared<PointLight>(Point3(0, 2, 0), RgbColor(1.0), 10.0)); }),
        Vector3(albedo / constants::PI * 10.0 / 4.0));

    // Sphere light (radius R, radiance L, center distance d, angle theta to the normal):
    //   E = pi L (R/d)^2 cos(theta)   ->   Lo = albedo L (R/d)^2 cos(theta)
    const Real overhead = albedo * 10.0 * (0.5 * 0.5) / 4.0;                       // theta = 0, d = 2
    const Real oblique = albedo * 10.0 * (0.25 / 4.5) * (1.5 / std::sqrt(4.5));    // center (1.5, 1.5, 0)
    for (int variant = 0; variant < 2; ++variant) {
        const bool nee = (variant == 0);
        const std::string tag = nee ? "NEE + MIS" : "BSDF sampling only";
        run("sphere light overhead, R=0.5 L=10  [" + tag + "]",
            MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                AddSphereLight(w, s, Point3(0, 2, 0), 0.5, RgbColor(1.0), 10.0, nee); }),
            Vector3(overhead));
        run("sphere light oblique,  R=0.5 L=10  [" + tag + "]",
            MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                AddSphereLight(w, s, Point3(1.5, 1.5, 0), 0.5, RgbColor(1.0), 10.0, nee); }),
            Vector3(oblique));
        run("sphere light overhead, big R=1.2 L=10  [" + tag + "]",
            MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                AddSphereLight(w, s, Point3(0, 2.5, 0), 1.2, RgbColor(1.0), 10.0, nee); }),
            Vector3(albedo * 10.0 * (1.2 * 1.2) / (2.5 * 2.5)));
        run("colored sphere light (1, .4, .1)*10  [" + tag + "]",
            MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                AddSphereLight(w, s, Point3(0, 2, 0), 0.5, RgbColor(1.0, 0.4, 0.1), 10.0, nee); }),
            Vector3(overhead) * Vector3(1.0, 0.4, 0.1));
    }

    // Rectangular (quad) lights: reference = numerical integration over the quad.
    //   (a) facing the ground, (b) tilted toward it, (c) facing away (must give 0), (d) two-sided, facing away.
    {
        struct QuadCase { const char* name; Point3 corner; Vector3 e1, e2; bool twoSided; };
        const QuadCase cases[] = {
            { "quad light overhead 1x1 at height 2, facing down", Point3(-0.5, 2, -0.5), Vector3(1, 0, 0), Vector3(0, 0, 1), false },
            { "quad light tilted toward the receiver",             Point3(1.0, 1.0, -0.7), Vector3(0, 0, 1.4), Vector3(-0.5, 0.7, 0), false },
            { "quad light CLOSE strip 0.6x1.8 at height 0.3",      Point3(0.1, 0.3, -0.9), Vector3(0.6, 0, 0), Vector3(0, 0, 1.8), false },
            { "quad light facing AWAY (one-sided, expect 0)",      Point3(-0.5, 2, -0.5), Vector3(0, 0, 1), Vector3(1, 0, 0), false },
            { "quad light facing away but two-sided",              Point3(-0.5, 2, -0.5), Vector3(0, 0, 1), Vector3(1, 0, 0), true },
        };
        for (const QuadCase& c : cases) {
            const Real expected = QuadGroundRadiance(c.corner, c.e1, c.e2, c.twoSided, 10.0, albedo);
            for (int variant = 0; variant < 2; ++variant) {
                const bool nee = (variant == 0);
                run(std::string(c.name) + (nee ? "  [NEE + MIS]" : "  [BSDF sampling only]"),
                    MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                        AddQuadLight(w, &s, c.corner, c.e1, c.e2, EmissionSpectrum(RgbColor(1.0)), 10.0, nee, c.twoSided); }),
                    Vector3(expected));
            }
        }
    }

    // Several lights of very different brightness at once. Every selection strategy must give the
    // same (sum of the single-light closed forms) value; the standard error shows how efficient it is.
    {
        const Real dA = 6.0;
        const Real sphereFar = albedo * 400.0 * (0.5 * 0.5) / (dA * dA);                         // bright far sphere
        const Real pointNear = albedo / constants::PI * 2.0 * (0.9 / 1.5) / 2.25;               // dim near point (1.2, 0.9, 0)
        const Point3 qc(0.4, 1.2, -0.4); const Vector3 qe1(1, 0, 0), qe2(0, 0, 0.8);
        const Real quadNear = QuadGroundRadiance(qc, qe1, qe2, false, 3.0, albedo);              // dim near quad
        const Real envLevel = albedo * 0.1;                                                      // uniform sky 0.1
        const Real expected = sphereFar + pointNear + quadNear + envLevel;
        auto dimEnv = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(0.1))));
        std::cout << "[Lights]   (multi-light scene: far sphere " << sphereFar << " + point " << pointNear
            << " + quad " << quadNear << " + sky " << envLevel << ")\n";
        for (const bool withSky : { true, false }) {
            for (const char* sampler : { "uniform", "power", "spatial", "all" }) {
                run(std::string(withSky ? "4 lights (with sky)" : "3 lights (no sky)") + ", light sampler = " + sampler,
                    MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                        AddSphereLight(w, s, Point3(0, dA, 0), 0.5, EmissionSpectrum(RgbColor(1.0)), 400.0, true);
                        s.addLight(std::make_shared<PointLight>(Point3(1.2, 0.9, 0), RgbColor(1.0), 2.0));
                        AddQuadLight(w, &s, qc, qe1, qe2, EmissionSpectrum(RgbColor(1.0)), 3.0, true);
                        if (withSky) s.setEnvLight(std::make_shared<EnvLight>(dimEnv, 1.0)); }, sampler),
                    Vector3(withSky ? expected : expected - envLevel));
            }
        }
    }

    // Light source spectra (EmissionSpectrum). Unit-luminance lights: Y of the result equals the
    // D65 point-light value above, only the color differs.
    {
        const Real pointY = albedo / constants::PI * 10.0 / 4.0;
        const EmissionSpectrum tungsten = EmissionSpectrum::IlluminantA();
        auto pointScene = [&](const EmissionSpectrum& e) {
            return MakeGroundScene(albedo, [&](Aggregate&, Scene& s) {
                s.addLight(std::make_shared<PointLight>(Point3(0, 2, 0), e, 10.0)); });
        };

        // (a) the SPD itself: CIE illuminant A has chromaticity xy = (0.4476, 0.4074).
        color::CIE_XYZ sum;
        for (int l = 380; l <= 780; ++l) {
            const Real S = color::evaluateBlackbody(Real(l), 2856.0);
            const color::CIE_XYZ c = color::evaluateCIE_XYZ_1931_MultiLobe(Real(l));
            sum.X += c.X * S; sum.Y += c.Y * S; sum.Z += c.Z * S;
        }
        const Real x = sum.X / (sum.X + sum.Y + sum.Z), y = sum.Y / (sum.X + sum.Y + sum.Z);
        const bool okA = std::abs(x - 0.4476) < 0.005 && std::abs(y - 0.4074) < 0.005;
        std::cout << "[Lights]   blackbody 2856 K chromaticity: xy = (" << x << ", " << y
            << ") vs CIE A (0.4476, 0.4074)   " << (okA ? "OK" : "FAIL") << "\n";
        allPassed = okA && allPassed;

        // (b) point light with the SPD: expected = sensor RGB of the SPD (Y = 1) x D65-light value.
        run("point light, illuminant A (tungsten)", pointScene(tungsten), tungsten.rgb() * pointY);
        run("point light, 1850 K candle blackbody", pointScene(EmissionSpectrum::Blackbody(1850.0)),
            EmissionSpectrum::Blackbody(1850.0).rgb() * pointY);

        // (c) CSV round trip: write a 10 nm SPD table (arbitrary scale), load it, same light.
        const std::string csvPath = (std::filesystem::temp_directory_path() / "rayt_test_spd.csv").string();
        {
            std::ofstream f(csvPath);
            f << "wavelength_nm,power\n";
            for (int l = 300; l <= 900; l += 10) f << l << "," << 7.0 * color::evaluateBlackbody(Real(l), 2856.0) << "\n";
        }
        color::SpectralCurve curve;
        std::string err;
        if (io::loadSpectralCurveCsv(csvPath, curve, &err)) {
            const EmissionSpectrum fromCsv = EmissionSpectrum::FromCurve(curve);
            const Vector3 d = glm::abs(fromCsv.rgb() / tungsten.rgb() - Vector3(1.0));
            const bool okC = std::max({ d.x, d.y, d.z }) < 0.005;
            std::cout << "[Lights]   SPD from CSV (" << curve.size() << " rows) vs generated A: max RGB difference "
                << std::max({ d.x, d.y, d.z }) * 100.0 << " %   " << (okC ? "OK" : "FAIL") << "\n";
            allPassed = okC && allPassed;
            run("point light, SPD loaded from CSV", pointScene(fromCsv), fromCsv.rgb() * pointY);
        }
        else {
            std::cout << "[Lights]   SPD CSV load failed: " << err << "   FAIL\n";
            allPassed = false;
        }
        std::error_code ec;
        std::filesystem::remove(csvPath, ec);

        // (d) the same SPD through a DiffuseLight sphere (NEE + MIS): expected = rgb x overhead value.
        run("sphere light, illuminant A  [NEE + MIS]",
            MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                AddSphereLight(w, s, Point3(0, 2, 0), 0.5, tungsten, 10.0, true); }),
            tungsten.rgb() * overhead);

        // (e) white balance on the tungsten light: the result must be neutral with the same luminance.
        color::setWhiteBalance(tungsten);
        run("point light, illuminant A with white balance (neutral)", pointScene(tungsten), Vector3(pointY));
        color::resetWhiteBalance();
    }

    // Environment maps: reference = dense grid integration of EnvMap::eval.
    Vector3 sunDir(0.0);
    auto sunEnv = MakeSunSkyEnv(sunDir);
    allPassed = CheckEnvSampling(*sunEnv, "synthetic sun + sky", n / 4) && allPassed;
    run("environment: synthetic sun + sky (importance sampled)",
        MakeGroundScene(albedo, [&](Aggregate&, Scene& s) {
            s.setEnvLight(std::make_shared<EnvLight>(sunEnv, 1.0)); }),
        GroundRadianceFromEnv(*sunEnv, albedo));
    run("environment: synthetic sun + sky, scale 2",
        MakeGroundScene(albedo, [&](Aggregate&, Scene& s) {
            s.setEnvLight(std::make_shared<EnvLight>(sunEnv, 2.0)); }),
        GroundRadianceFromEnv(*sunEnv, albedo) * 2.0);

    try {
        auto hdr = std::make_shared<EnvMap>(io::loadHDR(opt.env));
        allPassed = CheckEnvSampling(*hdr, opt.env, n / 4) && allPassed;
        run("environment: " + opt.env,
            MakeGroundScene(albedo, [&](Aggregate&, Scene& s) {
                s.setEnvLight(std::make_shared<EnvLight>(hdr, 1.0)); }),
            GroundRadianceFromEnv(*hdr, albedo));
    }
    catch (const std::exception& e) {
        std::cout << "[Lights]   (skipped HDR test: " << e.what() << ")\n";
    }

    std::cout << "[Lights] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

// -----------------------------------------------------------------------------
// Metal data tests (--scene metals)
// -----------------------------------------------------------------------------

// Largest |difference| of eta and k between two tables on a wavelength grid (relative to max(1, value)).
static Real MaxNkDifference(const color::ComplexIorTable& a, const color::ComplexIorTable& b) {
    Real worst = 0;
    for (int l = 400; l <= 700; l += 25) {
        const auto [ea, ka] = a.evaluateNK(Real(l));
        const auto [eb, kb] = b.evaluateNK(Real(l));
        worst = std::max({ worst, std::abs(ea - eb) / std::max(Real(1), std::abs(ea)),
                                  std::abs(ka - kb) / std::max(Real(1), std::abs(ka)) });
    }
    return worst;
}

// Normal-incidence reflectance curve R(l) of a table, sampled every 5 nm over 300..900 nm.
static color::SpectralCurve ReflectanceCurveOf(const color::ComplexIorTable& t) {
    color::SpectralCurve c;
    for (int l = 300; l <= 900; l += 5) {
        const auto [e, k] = t.evaluateNK(Real(l));
        c.add(Real(l), color::conductorNormalIncidenceReflectance(e, k));
    }
    c.finalize();
    return c;
}

static int RunMetalTests(const CliOptions& opt) {
    bool allPassed = true;
    auto report = [&](const std::string& what, Real value, Real limit) {
        const bool ok = value <= limit;
        std::cout << "[Metals]   " << what << ": " << value << " (limit " << limit << ")   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };

    color::ComplexIorTable ref;
    std::string err;
    if (!io::loadComplexIor_RefractiveIndexInfoCsv(AU_IOR_CSV_PATH, ref, &err)) {
        std::cerr << "[Metals] cannot load " << AU_IOR_CSV_PATH << ": " << err << "\n";
        return 1;
    }
    std::cout << "[Metals] reference: " << AU_IOR_CSV_PATH << " (" << ref.size() << " rows)\n";

    // ---- 1. loaders: the same data written in other layouts / units must give the same table ----
    const std::filesystem::path dir = std::filesystem::temp_directory_path();
    auto write = [&](const char* name, const std::function<void(std::ofstream&)>& body) {
        const std::string path = (dir / name).string();
        std::ofstream f(path);
        f.precision(12);
        body(f);
        return path;
    };
    const auto& rows = ref.etaCurve().samples();
    const auto& kRows = ref.kCurve().samples();

    const std::string pathNm = write("rayt_au_nm.csv", [&](std::ofstream& f) {
        f << "# wavelength_nm, n, k\n";
        for (size_t i = 0; i < rows.size(); ++i) f << rows[i].lambda_nm << "," << rows[i].value << "," << kRows[i].value << "\n";
    });
    const std::string pathUm = write("rayt_au_um.tsv", [&](std::ofstream& f) {
        f << "wl\tn\tk\n";
        for (size_t i = 0; i < rows.size(); ++i) f << rows[i].lambda_nm / 1000.0 << "\t" << rows[i].value << "\t" << kRows[i].value << "\n";
    });
    const std::string pathEv = write("rayt_au_ev.txt", [&](std::ofstream& f) {
        f << "eV n k\n";
        for (size_t i = rows.size(); i-- > 0;) f << 1239.841984 / rows[i].lambda_nm << " " << rows[i].value << " " << kRows[i].value << "\n";
    });
    const std::string pathN = write("rayt_au_n.csv", [&](std::ofstream& f) {
        for (const auto& s : rows) f << s.lambda_nm << "," << s.value << "\n";
    });
    const std::string pathK = write("rayt_au_k.csv", [&](std::ofstream& f) {
        for (const auto& s : kRows) f << s.lambda_nm << "," << s.value << "\n";
    });

    auto checkLoad = [&](const char* label, bool loaded, const color::ComplexIorTable& t, const std::string& e) {
        if (!loaded) { std::cout << "[Metals]   " << label << ": load failed (" << e << ")   FAIL\n"; allPassed = false; return; }
        report(std::string("loader ") + label + ", max |d(eta,k)| vs reference", MaxNkDifference(ref, t), 1e-6);
    };
    {
        color::ComplexIorTable t; std::string e;
        bool ok = io::loadComplexIor_ThreeColumn(pathNm, t, io::WavelengthUnit::Nm, &e);   checkLoad("3 columns, nm", ok, t, e);
        ok = io::loadComplexIor_ThreeColumn(pathUm, t, io::WavelengthUnit::Auto, &e);     checkLoad("3 columns, um (tab, auto)", ok, t, e);
        ok = io::loadComplexIor_ThreeColumn(pathEv, t, io::WavelengthUnit::Ev, &e);       checkLoad("3 columns, eV (space)", ok, t, e);
        ok = io::loadComplexIor_SeparateFiles(pathN, pathK, t, io::WavelengthUnit::Nm, &e); checkLoad("separate n / k files", ok, t, e);
        ok = io::loadComplexIor(pathUm, t, io::WavelengthUnit::Auto, &e);                 checkLoad("auto-detect (3 columns)", ok, t, e);
        ok = io::loadComplexIor(AU_IOR_CSV_PATH, t, io::WavelengthUnit::Auto, &e);        checkLoad("auto-detect (RefractiveIndex.info)", ok, t, e);
    }

    // ---- 2. reflectance input: F0 must be reproduced exactly ----
    const color::SpectralCurve R = ReflectanceCurveOf(ref);
    const color::ComplexIorTable fromR = color::complexIorTableFromReflectance(R);
    {
        Real worst = 0;
        for (int l = 360; l <= 830; l += 7) {   // off the 5 nm grid on purpose (interpolation error included)
            const auto [e, k] = fromR.evaluateNK(Real(l));
            worst = std::max(worst, std::abs(color::conductorNormalIncidenceReflectance(e, k) - R.evaluate(Real(l))));
        }
        report("reflectance -> (eta,k): max |F0 - R|, 360..830 nm", worst, 2e-3);

        // percent-valued CSV gives the same table
        const std::string pathR = write("rayt_au_R_percent.csv", [&](std::ofstream& f) {
            f << "wavelength_nm,reflectance_percent\n";
            for (const auto& s : R.samples()) f << s.lambda_nm << "," << s.value * 100.0 << "\n";
        });
        color::SpectralCurve Rloaded;
        if (io::loadReflectanceCsv(pathR, Rloaded, &err)) {
            report("reflectance CSV in percent: max |R - R_loaded|", [&] {
                Real w = 0; for (int l = 400; l <= 700; l += 10) w = std::max(w, std::abs(R.evaluate(Real(l)) - Rloaded.evaluate(Real(l)))); return w; }(), 1e-6);
        }
        else { std::cout << "[Metals]   reflectance CSV load failed: " << err << "   FAIL\n"; allPassed = false; }
        std::error_code ec;
        for (const std::string& p : { pathNm, pathUm, pathEv, pathN, pathK, pathR }) std::filesystem::remove(p, ec);
    }

    // ---- 3. angular error of the reflectance-only model against the true (eta, k) ----
    //   (this is the approximation the edge tint stands for; thresholds are the calibrated values + margin)
    for (const char* metal : { "Au", "Ag", "Cu" }) {
        color::ComplexIorTable truth;
        const std::string path = std::string("assets/Metal_data/Johnson_") + metal + ".csv";
        if (!io::loadComplexIor(path, truth, io::WavelengthUnit::Auto, &err)) {
            std::cout << "[Metals]   cannot load " << path << ": " << err << "   FAIL\n";
            allPassed = false;
            continue;
        }
        const color::ComplexIorTable approx = color::complexIorTableFromReflectance(ReflectanceCurveOf(truth));
        std::cout << "[Metals]   " << metal << " max |F(theta) - F_true(theta)| over 400..700 nm:";
        struct { int deg; Real limit; } checks[] = { {0, 0.003}, {45, 0.03}, {60, 0.06}, {70, 0.10}, {80, 0.20} };
        bool ok = true;
        for (const auto& c : checks) {
            Real worst = 0;
            const Real cosT = std::cos(glm::radians(Real(c.deg)));
            for (int l = 400; l <= 700; l += 10) {
                const auto [e1, k1] = truth.evaluateNK(Real(l));
                const auto [e2, k2] = approx.evaluateNK(Real(l));
                worst = std::max(worst, std::abs(fresnel::fresnelConductor(cosT, e1, k1) - fresnel::fresnelConductor(cosT, e2, k2)));
            }
            std::cout << "  " << c.deg << "deg " << worst;
            ok = ok && worst <= c.limit;
        }
        std::cout << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    }

    // ---- 4. through the renderer: true gold vs gold rebuilt from its reflectance ----
    //   A rough-conductor sphere in a uniform white environment; one ray at the center
    //   (near-normal incidence) and one near the rim (about 70 deg).
    {
        auto envMap = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        auto makeScene = [&](const color::ComplexIorTable& t) {
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Sphere>(Point3(0, 0, 0), 1.0,
                std::make_shared<RoughConductor>(std::make_shared<color::ComplexIorTable>(t), 0.2)));
            auto scene = std::make_shared<Scene>(world);
            scene->setEnvLight(std::make_shared<EnvLight>(envMap, 1.0));
            scene->finalize();
            return scene;
        };
        const auto sceneTrue = makeScene(ref);
        const auto sceneFromR = makeScene(color::complexIorTableFromReflectance(ReflectanceCurveOf(ref)));
        const int n = opt.spp * 1024;
        for (const Real x : { Real(0), Real(0.94) }) {
            const Ray ray(Point3(x, 0, 5), Vector3(0, 0, -1));
            Vector3 m1, s1, m2, s2;
            int bad = 0;
            EstimateRadiance(*sceneTrue, ray, n, 8, m1, s1, bad);
            EstimateRadiance(*sceneFromR, ray, n, 8, m2, s2, bad);
            const Vector3 d = glm::abs(m1 - m2);
            const Real limit = (x == 0) ? Real(0.01) : Real(0.15);   // rim: known edge-tint approximation, regression guard
            const Real worst = std::max({ d.x, d.y, d.z });
            std::cout << "[Metals]   rendered gold, " << (x == 0 ? "center (~0 deg)" : "rim (~70 deg)")
                << ": true (" << m1.x << ", " << m1.y << ", " << m1.z << ")  from R (" << m2.x << ", " << m2.y << ", " << m2.z
                << ")  max diff " << worst << " (limit " << limit << ")   " << (worst <= limit + 4 * std::max({ s1.x, s2.x }) ? "OK" : "FAIL") << "\n";
            allPassed = allPassed && (worst <= limit + 4 * std::max({ s1.x, s2.x }));
        }
    }

    // ---- 5. multiple-scattering compensation (Kulla-Conty) ----
    {
        // 5a. near-perfect conductor in a uniform white environment: must return all the light.
        auto envMap = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        auto perfect = std::make_shared<color::ComplexIorTable>();
        perfect->assign({ { 360.0, 1.0, 100.0 }, { 830.0, 1.0, 100.0 } });
        const int n = opt.spp * 512;
        for (const Real roughness : { Real(0.2), Real(0.5), Real(0.8), Real(1.0) }) {
            std::string line = "[Metals]   perfect conductor, roughness " + std::to_string(roughness).substr(0, 3) + ": albedo (center | rim)";
            bool ok = true;
            for (const bool ms : { false, true }) {
                auto mat = std::make_shared<RoughConductor>(perfect, roughness);
                mat->setMultipleScattering(ms);
                auto world = std::make_shared<Aggregate>();
                world->add(std::make_shared<Sphere>(Point3(0, 0, 0), 1.0, mat));
                auto scene = std::make_shared<Scene>(world);
                scene->setEnvLight(std::make_shared<EnvLight>(envMap, 1.0));
                scene->finalize();
                Vector3 mc, sc, mr, sr;
                int bad = 0;
                EstimateRadiance(*scene, Ray(Point3(0, 0, 5), Vector3(0, 0, -1)), n, 16, mc, sc, bad);
                EstimateRadiance(*scene, Ray(Point3(0.94, 0, 5), Vector3(0, 0, -1)), n, 16, mr, sr, bad);
                line += std::string(ms ? "   with MS: " : "   single: ") + std::to_string(mc.x).substr(0, 5) + " | " + std::to_string(mr.x).substr(0, 5);
                if (ms) ok = std::abs(mc.x - 1.0) <= 4 * sc.x + 0.02 && std::abs(mr.x - 1.0) <= 4 * sr.x + 0.03 && bad == 0;
            }
            std::cout << line << "   " << (ok ? "OK" : "FAIL") << "\n";
            allPassed = allPassed && ok;
        }

        // 5b. gold, roughness 0.6: reciprocity, pdf normalization, and sample() vs uniform integration.
        auto gold = std::make_shared<RoughConductor>(std::make_shared<color::ComplexIorTable>(ref), 0.6);
        SurfaceInteraction rec;
        rec.p = Point3(0.0); rec.n = Vector3(0, 0, 1); rec.gn = Vector3(0, 0, 1);
#if defined(RAYT_SPECTRAL_MODE)
        const SampledWavelengths lambda0 = SampledWavelengths::sampleVisible(Real(0.37));   // fixed wavelengths: a meaningful, repeatable albedo
#else
        const SampledWavelengths lambda0;
#endif
        auto randomDir = [&]() {
            const Real z = sampling::Random(0.05, 1.0), phi = 2 * constants::PI * sampling::Random();
            const Real r = std::sqrt(1 - z * z);
            return Vector3(r * std::cos(phi), r * std::sin(phi), z);
        };
        Real worstRecip = 0;
        for (int i = 0; i < 2000; ++i) {
            const Vector3 a = randomDir(), b = randomDir();
            const Spectrum f1 = gold->eval(rec, a, b, lambda0, TransportMode::Radiance);
            const Spectrum f2 = gold->eval(rec, b, a, lambda0, TransportMode::Radiance);
            const Vector3 v1 = SpectrumToRGB(f1, lambda0), v2 = SpectrumToRGB(f2, lambda0);
            worstRecip = std::max(worstRecip, std::abs(v1.x - v2.x) / std::max(Real(1e-9), std::max(v1.x, v2.x)));
        }
        report("gold r=0.6: reciprocity f(wo,wi) = f(wi,wo), max relative difference", worstRecip, 1e-6);

        for (const Real cosO : { Real(0.9), Real(0.5), Real(0.15) }) {
            const Vector3 wo(std::sqrt(1 - cosO * cosO), 0, cosO);
            double pdfInt = 0, albedoU = 0, albedoS = 0, sqU = 0, sqS = 0;
            const int m = 200000;
            for (int i = 0; i < m; ++i) {
                // uniform hemisphere: integral of pdf, and integral of f cos
                const Real z = sampling::Random(), phi = 2 * constants::PI * sampling::Random(), r = std::sqrt(std::max(Real(0), 1 - z * z));
                const Vector3 wi(r * std::cos(phi), r * std::sin(phi), z);
                pdfInt += gold->pdf(rec, wo, wi, lambda0) * 2 * constants::PI;
                const double u = SpectrumToRGB(gold->eval(rec, wo, wi, lambda0, TransportMode::Radiance), lambda0).x * wi.z * 2 * constants::PI;
                albedoU += u; sqU += u * u;
                // importance sampled through the material's own sample()
                const auto s = gold->sample(rec, wo, Point2(sampling::Random(), sampling::Random()), lambda0, TransportMode::Radiance);
                const double w = s ? SpectrumToRGB(s->f, lambda0).x * std::abs(glm::dot(rec.n, s->wi)) / s->pdf : 0.0;
                albedoS += w; sqS += w * w;
            }
            pdfInt /= m; albedoU /= m; albedoS /= m;
            const double seU = std::sqrt(std::max(0.0, sqU / m - albedoU * albedoU) / m), seS = std::sqrt(std::max(0.0, sqS / m - albedoS * albedoS) / m);
            // The pdf integrates to <= 1: GGX reflections that leave the hemisphere are lost (sample() fails for them).
            const bool ok = pdfInt <= 1.02 && pdfInt > 0.8 && std::abs(albedoU - albedoS) <= 4 * std::hypot(seU, seS) + 1e-3;
            std::cout << "[Metals]   gold r=0.6, cos(theta_o)=" << cosO << ": integral of pdf = " << pdfInt << " (<= 1), albedo (R) by uniform integration "
                << albedoU << " vs by sample() " << albedoS << "   " << (ok ? "OK" : "FAIL") << "\n";
            allPassed = allPassed && ok;
        }
    }

    // ---- 6. CoatedDiffuse (diffuse under a GGX coat, metallic blend) ----
    {
        SurfaceInteraction rec;
        rec.p = Point3(0.0); rec.n = Vector3(0, 0, 1); rec.gn = Vector3(0, 0, 1);
#if defined(RAYT_SPECTRAL_MODE)
        const SampledWavelengths lam = SampledWavelengths::sampleVisible(Real(0.37));
#else
        const SampledWavelengths lam;
#endif
        auto randomDir = [&]() {
            const Real z = sampling::Random(0.05, 1.0), phi = 2 * constants::PI * sampling::Random();
            const Real r = std::sqrt(1 - z * z);
            return Vector3(r * std::cos(phi), r * std::sin(phi), z);
        };
        struct Case { const char* name; Real albedo, roughness, weight, metallic; };
        for (const Case& c : { Case{ "plastic r=0.3", 0.8, 0.3, 1.0, 0.0 }, Case{ "plastic r=0.7", 0.8, 0.7, 1.0, 0.0 },
                               Case{ "half coat, 30% metal r=0.4", 0.8, 0.4, 0.5, 0.3 }, Case{ "metal r=0.4", 0.8, 0.4, 1.0, 1.0 } }) {
            auto mat = std::make_shared<CoatedDiffuse>(RgbColor(c.albedo), c.roughness, c.weight, c.metallic);
            Real worstRecip = 0;
            for (int i = 0; i < 2000; ++i) {
                const Vector3 a = randomDir(), b = randomDir();
                const Vector3 v1 = SpectrumToRGB(mat->eval(rec, a, b, lam, TransportMode::Radiance), lam);
                const Vector3 v2 = SpectrumToRGB(mat->eval(rec, b, a, lam, TransportMode::Radiance), lam);
                worstRecip = std::max(worstRecip, std::abs(v1.x - v2.x) / std::max(Real(1e-9), std::max(v1.x, v2.x)));
            }
            bool ok = worstRecip < 1e-6;
            std::ostringstream line;
            line << "[Metals]   CoatedDiffuse " << c.name << ": reciprocity " << worstRecip;
            for (const Real cosO : { Real(0.8), Real(0.3) }) {
                const Vector3 wo(std::sqrt(1 - cosO * cosO), 0, cosO);
                double albedoU = 0, albedoS = 0, sqU = 0, sqS = 0, pdfInt = 0;
                const int m = 200000;
                for (int i = 0; i < m; ++i) {
                    const Real z = sampling::Random(), phi = 2 * constants::PI * sampling::Random(), r = std::sqrt(std::max(Real(0), 1 - z * z));
                    const Vector3 wi(r * std::cos(phi), r * std::sin(phi), z);
                    pdfInt += mat->pdf(rec, wo, wi, lam) * 2 * constants::PI;
                    const double u = SpectrumToRGB(mat->eval(rec, wo, wi, lam, TransportMode::Radiance), lam).x * wi.z * 2 * constants::PI;
                    albedoU += u; sqU += u * u;
                    const auto s = mat->sample(rec, wo, Point2(sampling::Random(), sampling::Random()), lam, TransportMode::Radiance);
                    const double w = s ? SpectrumToRGB(s->f, lam).x * std::abs(glm::dot(rec.n, s->wi)) / s->pdf : 0.0;
                    albedoS += w; sqS += w * w;
                }
                pdfInt /= m; albedoU /= m; albedoS /= m;
                const double seU = std::sqrt(std::max(0.0, sqU / m - albedoU * albedoU) / m), seS = std::sqrt(std::max(0.0, sqS / m - albedoS * albedoS) / m);
                const bool okHere = pdfInt <= 1.02 && std::abs(albedoU - albedoS) <= 4 * std::hypot(seU, seS) + 1e-3;
                ok = ok && okHere;
                line << ";  cos " << cosO << ": albedo " << albedoU << " vs sample() " << albedoS << ", pdf integral " << pdfInt;
            }
            std::cout << line.str() << "   " << (ok ? "OK" : "FAIL") << "\n";
            allPassed = allPassed && ok;
        }

        // A white (albedo 1) plastic sphere in a uniform white environment can only return <= 1 of the light.
        auto envMap = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        for (const Real rough : { Real(0.2), Real(0.5) }) {
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Sphere>(Point3(0, 0, 0), 1.0, std::make_shared<CoatedDiffuse>(RgbColor(1.0), rough, 1.0, 0.0)));
            auto scene = std::make_shared<Scene>(world);
            scene->setEnvLight(std::make_shared<EnvLight>(envMap, 1.0));
            scene->finalize();
            Vector3 mc, sc;
            int bad = 0;
            EstimateRadiance(*scene, Ray(Point3(0, 0, 5), Vector3(0, 0, -1)), opt.spp * 512, 16, mc, sc, bad);
            const bool ok = mc.x <= 1.0 + 4 * sc.x + 0.005 && mc.x > 0.9 && bad == 0;
            std::cout << "[Metals]   white plastic r=" << rough << " in a white furnace: " << mc.x << " (<= 1)   " << (ok ? "OK" : "FAIL") << "\n";
            allPassed = allPassed && ok;
        }
    }

    std::cout << "[Metals] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

// -----------------------------------------------------------------------------
// Triangle mesh checks (run at the start of --scene mesh)
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// MMD (PMX) models  (--scene pmx)
// -----------------------------------------------------------------------------
// The files live in MMD_model/ (git-ignored, third-party assets). Names are CJK, which a Windows
// command line cannot carry reliably, so a model is chosen by its index in --pmx-list.

static std::string PathToUtf8(const std::filesystem::path& p) {
    const std::u8string s = p.u8string();
    return std::string(s.begin(), s.end());
}

static std::filesystem::path PathFromUtf8(const std::string& s) {
    return std::filesystem::path(std::u8string(s.begin(), s.end()));
}

static std::vector<std::filesystem::path> FindMmdFiles(const std::string& subdir, const std::string& extension) {
    namespace fs = std::filesystem;
    std::vector<fs::path> found;
    for (const std::string base : { "../MMD_model", "MMD_model" }) {
        const fs::path dir = subdir.empty() ? fs::path(base) : fs::path(base) / subdir;
        std::error_code ec;
        if (!fs::exists(dir, ec)) continue;
        for (const auto& e : fs::recursive_directory_iterator(dir, ec)) {
            if (e.is_regular_file() && e.path().extension() == extension) found.push_back(e.path());
        }
        if (!found.empty()) break;
    }
    std::sort(found.begin(), found.end());
    return found;
}

static std::vector<std::filesystem::path> FindPmxFiles() { return FindMmdFiles("", ".pmx"); }

static Real SrgbToLinear(Real c) { return std::pow(std::clamp(c, Real(0), Real(1)), Real(2.2)); }

static void PrintPmxInfo(const io::PmxModel& m, const std::filesystem::path& file) {
    namespace fs = std::filesystem;
    std::cout << "[PMX] " << PathToUtf8(file.filename()) << "   model name: " << m.name << "  (PMX " << m.version
        << ", " << m.vertices.size() << " vertices, " << m.indices.size() / 3 << " triangles, " << m.textures.size()
        << " textures, " << m.materials.size() << " materials, " << m.additionalUv << " extra UV sets)\n";

    int skin[5] = {};
    for (const auto& v : m.vertices) ++skin[std::clamp(v.weightType, 0, 4)];
    std::cout << "[PMX]   skinning: BDEF1 " << skin[0] << ", BDEF2 " << skin[1] << ", BDEF4 " << skin[2] << ", SDEF " << skin[3] << ", QDEF " << skin[4]
        << "   (bones/morphs are in the file but not read: this is a static render)\n";

    int existing = 0;
    for (const auto& t : m.textures) {
        std::string rel = t;
        std::replace(rel.begin(), rel.end(), '\\', '/');
        std::error_code ec;
        if (fs::exists(file.parent_path() / PathFromUtf8(rel), ec)) ++existing;
    }
    std::cout << "[PMX]   texture files found next to the model: " << existing << " of " << m.textures.size() << "\n";

    for (size_t i = 0; i < m.materials.size(); ++i) {
        const auto& mat = m.materials[i];
        std::cout << "[PMX]   material " << i << ": " << mat.name << "  tris " << mat.triangleCount()
            << "  diffuse (" << mat.diffuse.x << ", " << mat.diffuse.y << ", " << mat.diffuse.z << ") alpha " << mat.diffuse.w
            << (mat.doubleSided() ? "  double-sided" : "")
            << (mat.texture >= 0 ? "  tex " + m.textures[static_cast<size_t>(mat.texture)] : std::string("  no texture"))
            << "  spec (" << mat.specular.x << "," << mat.specular.y << "," << mat.specular.z << ") shin " << mat.shininess
            << (mat.sphereTexture >= 0 ? "  sphere" + std::to_string(mat.sphereMode) + " " + m.textures[static_cast<size_t>(mat.sphereTexture)] : std::string())
            << (mat.sharedToon ? "  toon#" + std::to_string(mat.toon) : std::string("  own toon"))
            << "\n";
    }
}

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
static std::shared_ptr<Material> MakePmxMaterial(const io::PmxMaterial& pm, const RgbColor& base,
    const std::shared_ptr<Texture2D>& tex, bool matte)
{
    if (matte) return std::make_shared<Lambertian>(base, tex);

    const Real specLum = Real(0.2126) * pm.specular.x + Real(0.7152) * pm.specular.y + Real(0.0722) * pm.specular.z;
    if (specLum > Real(0.01)) {
        const Real alpha = std::sqrt(Real(2) / (std::max<Real>(pm.shininess, Real(1)) + Real(2)));
        return std::make_shared<CoatedDiffuse>(base, std::sqrt(alpha), std::min(Real(1), Real(2) * specLum), Real(0), tex);
    }

    const auto has = [&](const char* s) { return pm.name.find(s) != std::string::npos; };
    if (has("髮") || has("髪")) return std::make_shared<CoatedDiffuse>(base, Real(0.30), Real(1.0), Real(0), tex);   // hair: soft shine
    if (has("目")) return std::make_shared<CoatedDiffuse>(base, Real(0.08), Real(1.0), Real(0), tex);                 // eyes: wet cornea
    if (pm.name == "饰") return std::make_shared<CoatedDiffuse>(base, Real(0.30), Real(1.0), Real(0.85), tex);        // gold ornaments
    if (has("体") || has("肌") || has("颜") || has("首")) return std::make_shared<CoatedDiffuse>(base, Real(0.50), Real(0.4), Real(0), tex);   // skin: faint sheen
    return std::make_shared<Lambertian>(base, tex);
}

/// PMX -> TriangleMesh in the renderer's convention (right-handed, +Y up, metres, feet on y = 0, x/z centered).
/// MMD is left-handed with 1 unit ~ 8 cm; mirroring z fixes the handedness, yaw (degrees) turns the model about +Y.
static std::shared_ptr<TriangleMesh> PmxToMesh(const io::PmxModel& m, Real yawDegrees, Real* outHeight, Vector3* outSize = nullptr) {
    constexpr Real kScale = 0.08;
    auto mesh = std::make_shared<TriangleMesh>();
    const Real yaw = glm::radians(yawDegrees), cy = std::cos(yaw), sy = std::sin(yaw);
    auto rotate = [&](const Vector3& v) { return Vector3(cy * v.x + sy * v.z, v.y, -sy * v.x + cy * v.z); };

    Vector3 lo(1e30), hi(-1e30);
    for (const auto& v : m.vertices) {
        const Vector3 p = rotate(Vector3(v.position.x, v.position.y, -v.position.z) * kScale);
        mesh->p.push_back(Point3(p));
        mesh->n.push_back(glm::normalize(rotate(Vector3(v.normal.x, v.normal.y, -v.normal.z))));
        mesh->uv.push_back(v.uv);
        lo = glm::min(lo, p);
        hi = glm::max(hi, p);
    }
    const Vector3 shift(-(lo.x + hi.x) / 2, -lo.y, -(lo.z + hi.z) / 2);
    for (auto& p : mesh->p) p += shift;
    mesh->indices = m.indices;
    if (outHeight) *outHeight = hi.y - lo.y;
    if (outSize) *outSize = hi - lo;
    return mesh;
}

/// Skeleton summary and self-tests of the posing code (--pmx-check).
static bool RunPmxChecks(const io::PmxModel& model) {
    bool allPassed = true;
    auto check = [&](const std::string& what, bool ok, const std::string& detail) {
        std::cout << "[Pose]   " << what << ": " << detail << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };

    const size_t nb = model.bones.size();
    int ik = 0, grantRot = 0, grantTrans = 0, afterPhysics = 0, roots = 0;
    for (const auto& b : model.bones) {
        ik += b.isIk(); grantRot += b.inheritsRotation(); grantTrans += b.inheritsTranslation(); afterPhysics += b.afterPhysics();
        roots += (b.parent < 0);
    }
    int unskinned = 0;
    for (const auto& v : model.vertices) if (v.bones[0] < 0) ++unskinned;
    std::cout << "[Pose] skeleton: " << nb << " bones (" << roots << " roots, " << ik << " IK, " << grantRot << " rotation-grant, " << grantTrans
        << " translation-grant, " << afterPhysics << " deformed after physics), " << unskinned << " vertices without bones\n";
    {
        std::map<int, int> usage;   // bones that actually move vertices
        for (const auto& v : model.vertices) for (int k = 0; k < 4; ++k) if (v.bones[k] >= 0 && v.weights[k] > 0) ++usage[v.bones[k]];
        int usedGrant = 0, usedIk = 0;
        for (const auto& [bone, n] : usage) { usedGrant += model.bones[static_cast<size_t>(bone)].inheritsRotation(); usedIk += model.bones[static_cast<size_t>(bone)].isIk(); }
        std::cout << "[Pose] skinning uses " << usage.size() << " bones, " << usedGrant << " of them grant-rotation bones, " << usedIk << " IK bones\n";
    }

    // 1. identity pose must reproduce the rest pose
    {
        io::VpdPose none;
        io::VpdBone dummy; dummy.name = "__none__";
        none.bones.push_back(dummy);
        const io::PmxModel same = io::applyPose(model, none);
        Real worst = 0;
        for (size_t i = 0; i < model.vertices.size(); ++i) {
            worst = std::max(worst, glm::length(Vector3(same.vertices[i].position) - Vector3(model.vertices[i].position)));
            worst = std::max(worst, glm::length(Vector3(same.vertices[i].normal) - Vector3(model.vertices[i].normal)));
        }
        check("rest pose reproduced (max vertex / normal difference; normals are renormalized)", worst < 1e-5, std::to_string(worst));
    }

    // 2. rotate the root-most "center" bone by 90 degrees about +Y: everything under it must turn rigidly about its rest position
    {
        int center = -1;
        for (size_t i = 0; i < nb; ++i) if (model.bones[i].name == "センター") { center = static_cast<int>(i); break; }
        if (center < 0) { check("center bone test", false, "no bone named センター"); }
        else {
            io::VpdPose p;
            io::VpdBone b; b.name = "センター"; b.qy = std::sin(constants::PI / 4); b.qw = std::cos(constants::PI / 4);
            p.bones.push_back(b);
            const io::PmxModel moved = io::applyPose(model, p);
            const Vector3 c(model.bones[static_cast<size_t>(center)].position);
            const glm::dmat3 R = glm::dmat3(glm::mat4_cast(glm::angleAxis(constants::PI / 2, glm::dvec3(0, 1, 0))));
            // descendants of the center bone
            std::vector<char> under(nb, 0);
            for (size_t i = 0; i < nb; ++i) { int cur = static_cast<int>(i); for (size_t guard = 0; guard <= nb && cur >= 0; ++guard) { if (cur == center) { under[i] = 1; break; } cur = model.bones[static_cast<size_t>(cur)].parent; } }
            int total = 0, matching = 0;
            for (size_t i = 0; i < model.vertices.size(); ++i) {
                const auto& v = model.vertices[i];
                bool all = v.bones[0] >= 0;
                for (int k = 0; k < 4 && all; ++k) if (v.bones[k] >= 0 && v.weights[k] > 0 && !under[static_cast<size_t>(v.bones[k])]) all = false;
                if (!all) continue;
                ++total;
                const Vector3 expected = c + R * (Vector3(v.position) - c);
                if (glm::length(Vector3(moved.vertices[i].position) - expected) < 1e-6) ++matching;
            }
            check("rotating センター by 90 deg about +Y turns " + std::to_string(total) + " vertices rigidly about its rest position", total > 0 && matching == total,
                std::to_string(matching) + " of " + std::to_string(total) + " match");
        }
    }

    // 3. grant bones: with an identity own rotation, a bone inheriting rotation with weight w must end up with q^w of its grant parent
    {
        int tested = 0, passed = 0;
        for (size_t i = 0; i < nb && tested < 6; ++i) {
            const auto& b = model.bones[i];
            if (!b.inheritsRotation() || b.grantParent < 0) continue;
            const auto& gp = model.bones[static_cast<size_t>(b.grantParent)];
            if (gp.inheritsRotation()) continue;   // keep the test simple: a plain FK source
            io::VpdPose p;
            io::VpdBone vb; vb.name = gp.name; vb.qx = 0.3; vb.qy = 0.2; vb.qz = 0.1; vb.qw = std::sqrt(1 - 0.09 - 0.04 - 0.01);
            p.bones.push_back(vb);
            io::PoseReport rep;
            io::applyPose(model, p, &rep);
            const glm::dquat src = rep.effectiveRotation[static_cast<size_t>(b.grantParent)];
            const glm::dquat expected = io::detail::quatPow(src, b.grantWeight);
            const glm::dquat got = rep.effectiveRotation[i];
            ++tested;
            if (std::abs(std::abs(glm::dot(expected, got)) - 1.0) < 1e-9) ++passed;
            std::cout << "[Pose]     grant " << model.bones[i].name << " <- " << gp.name << " x" << b.grantWeight << ": "
                << (std::abs(std::abs(glm::dot(expected, got)) - 1.0) < 1e-9 ? "ok" : "MISMATCH") << "\n";
        }
        check("grant (append) rotation of " + std::to_string(tested) + " sample bones", tested == 0 || passed == tested, std::to_string(passed) + " of " + std::to_string(tested));
    }
    std::cout << "[Pose] checks: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed;
}

// Triangle mesh -> Hittable with the chosen intersection backend (own SAH BVH or Embree).
static std::shared_ptr<Hittable> MakeMeshHittable(const std::shared_ptr<const TriangleMesh>& mesh,
    const std::shared_ptr<Material>& material, const std::string& accel)
{
    if (accel == "embree") return std::make_shared<EmbreeMesh>(mesh, material);
    return BuildMeshBVH(mesh, material);
}

// -----------------------------------------------------------------------------
// BVH benchmark (--scene bench): the reference recursive BVH vs the SAH flat BVH
// -----------------------------------------------------------------------------
static int RunBvhBenchmark(const CliOptions& opt) {
    using clock = std::chrono::steady_clock;
    auto ms = [](clock::time_point a, clock::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    const int nRays = opt.spp * 2000;
    bool ok = true;
    std::cout << "[Bench] " << nRays << " random rays per mesh, single thread (use a Release build)\n";

    auto white = std::make_shared<Lambertian>(RgbColor(1.0));
    for (const int subdiv : { 5, 6, 7 }) {
        const auto mesh = std::make_shared<const TriangleMesh>(TriangleMesh::Icosphere(subdiv, 1.0, Point3(0.0), true));
        std::vector<std::shared_ptr<Hittable>> prims;
        for (size_t i = 0; i < mesh->numTriangles(); ++i) prims.push_back(std::make_shared<Triangle>(mesh, i, white));

        // Rays from a sphere of radius 3 toward random points of a ball of radius 1.3 (about half hit).
        std::vector<Ray> rays(static_cast<size_t>(nRays));
        for (Ray& r : rays) {
            const Vector3 from = glm::normalize(Vector3(sampling::Random(-1, 1), sampling::Random(-1, 1), sampling::Random(-1, 1))) * 3.0;
            Vector3 to(0.0);
            do { to = Vector3(sampling::Random(-1.3, 1.3), sampling::Random(-1.3, 1.3), sampling::Random(-1.3, 1.3)); } while (glm::dot(to, to) > 1.69);
            r = Ray(from, glm::normalize(to - from));
        }

        auto v1 = prims;
        const auto t0 = clock::now();
        BVHNode legacy(v1, 0, v1.size());
        const auto t1 = clock::now();
        BVH flat(prims);
        const auto t2 = clock::now();

        SurfaceInteraction rec;
        int hitsLegacy = 0, hitsFlat = 0, mismatches = 0;
        std::vector<Real> tLegacy(rays.size()), tFlat(rays.size());
        const auto c0 = clock::now();
        for (size_t i = 0; i < rays.size(); ++i) { const bool h = legacy.hit(rays[i], rec); tLegacy[i] = h ? rec.t : -1; hitsLegacy += h; }
        const auto c1 = clock::now();
        for (size_t i = 0; i < rays.size(); ++i) { const bool h = flat.hit(rays[i], rec); tFlat[i] = h ? rec.t : -1; hitsFlat += h; }
        const auto c2 = clock::now();
        int occ = 0;
        for (size_t i = 0; i < rays.size(); ++i) occ += flat.occluded(rays[i]) ? 1 : 0;
        const auto c3 = clock::now();

        for (size_t i = 0; i < rays.size(); ++i) {
            if (std::abs(tLegacy[i] - tFlat[i]) > 1e-9) ++mismatches;
        }
        ok = ok && mismatches == 0 && occ == hitsFlat;

        // Embree (float traversal; the hit record is built in Real precision by the Triangle).
        const auto e0 = clock::now();
        EmbreeMesh embree(mesh, white);
        const auto e1 = clock::now();
        int hitsEmbree = 0, embreeMismatch = 0;
        std::vector<Real> tEmbree(rays.size());
        const auto e2 = clock::now();
        for (size_t i = 0; i < rays.size(); ++i) { const bool h = embree.hit(rays[i], rec); tEmbree[i] = h ? rec.t : -1; hitsEmbree += h; }
        const auto e3 = clock::now();
        int occEmbree = 0;
        for (size_t i = 0; i < rays.size(); ++i) occEmbree += embree.occluded(rays[i]) ? 1 : 0;
        const auto e4 = clock::now();
        for (size_t i = 0; i < rays.size(); ++i) {
            if (std::abs(tEmbree[i] - tFlat[i]) > 1e-6) ++embreeMismatch;
        }
        const bool embreeOk = embreeMismatch <= std::max(2, nRays / 10000);   // a ray grazing an edge may differ in float vs double
        ok = ok && embreeOk;
        std::cout << "[Bench]     Embree: build " << ms(e0, e1) << " ms (incl. Triangle objects, all threads);  closest hit "
            << double(nRays) / 1e6 / (ms(e2, e3) / 1000.0) << " Mray/s (x" << ms(c1, c2) / ms(e2, e3) << " vs flat BVH);  any-hit "
            << double(nRays) / 1e6 / (ms(e3, e4) / 1000.0) << " Mray/s;  " << hitsEmbree << " hits, " << embreeMismatch
            << " differ from the double BVH (|dt|>1e-6)   " << (embreeOk ? "OK" : "FAIL") << "\n";
        (void)occEmbree;

        if (subdiv == 6) {
            // Coherent rays (camera-like scanline order from one point) on the same BVH: shows how much
            // of the random-ray cost is cache misses.
            const int side = static_cast<int>(std::sqrt(double(nRays)));
            std::vector<Ray> cam;
            cam.reserve(static_cast<size_t>(side) * static_cast<size_t>(side));
            for (int y = 0; y < side; ++y)
                for (int x = 0; x < side; ++x) {
                    const Vector3 target(-1.2 + 2.4 * (x + 0.5) / side, -1.2 + 2.4 * (y + 0.5) / side, 0.0);
                    cam.emplace_back(Point3(0, 0, 3), glm::normalize(target - Vector3(0, 0, 3)));
                }
            int h = 0;
            const auto s0 = clock::now();
            for (const Ray& r : cam) h += flat.hit(r, rec) ? 1 : 0;
            const auto s1 = clock::now();
            std::cout << "[Bench]     coherent camera rays (" << cam.size() << "): " << double(cam.size()) / 1e6 / (ms(s0, s1) / 1000.0)
                << " Mray/s (" << h << " hits)\n";

            // Leaf size sensitivity (closest hit, same rays).
            for (const int leaf : { 1, 2, 4, 8 }) {
                BVH b(prims, leaf);
                int h = 0;
                const auto s0 = clock::now();
                for (size_t i = 0; i < rays.size(); ++i) h += b.hit(rays[i], rec) ? 1 : 0;
                const auto s1 = clock::now();
                std::cout << "[Bench]     leaf size " << leaf << ": " << double(nRays) / 1e6 / (ms(s0, s1) / 1000.0) << " Mray/s, "
                    << b.numNodes() << " nodes (" << h << " hits)\n";
            }
        }

        const double mr = double(nRays) / 1e6;
        std::cout << "[Bench] " << mesh->numTriangles() << " triangles:  build  recursive " << ms(t0, t1) << " ms,  SAH flat " << ms(t1, t2)
            << " ms  (" << flat.numNodes() << " nodes, depth " << flat.maxDepth() << ")\n"
            << "[Bench]     closest hit  recursive " << mr / (ms(c0, c1) / 1000.0) << " Mray/s,  flat " << mr / (ms(c1, c2) / 1000.0)
            << " Mray/s  (x" << ms(c0, c1) / ms(c1, c2) << ");  any-hit flat " << mr / (ms(c2, c3) / 1000.0) << " Mray/s\n"
            << "[Bench]     results: " << hitsLegacy << " / " << hitsFlat << " hits, " << mismatches << " distance mismatches, any-hit = "
            << occ << "   " << ((mismatches == 0 && occ == hitsFlat) ? "OK" : "FAIL") << "\n";
    }
    std::cout << "[Bench] RESULT: " << (ok ? "PASS" : "FAIL") << "\n";
    return ok ? 0 : 2;
}

static bool RunMeshChecks(const CliOptions& opt) {
    bool allPassed = true;
    auto check = [&](const std::string& what, bool pass, const std::string& detail) {
        std::cout << "[Mesh]   " << what << ": " << detail << "   " << (pass ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && pass;
    };
    auto white = std::make_shared<Lambertian>(RgbColor(1.0));

    // a. geometry of a procedural icosphere (an inscribed polyhedron: a bit smaller than the sphere)
    auto sphere = std::make_shared<TriangleMesh>(TriangleMesh::Icosphere(3, 1.0, Point3(0.0), true));
    {
        const Real a = sphere->surfaceArea() / (4 * constants::PI);
        const Real v = sphere->signedVolume() / (4.0 / 3.0 * constants::PI);
        check("icosphere (" + std::to_string(sphere->numTriangles()) + " triangles) area / sphere, volume / sphere",
            a > 0.97 && a <= 1.0 && v > 0.95 && v <= 1.0, std::to_string(a) + ", " + std::to_string(v) + "  (volume > 0 = outward winding)");
    }

    // b. one ray: hit distance and the interpolated (smooth) normal
    const std::shared_ptr<const TriangleMesh> sphereC = sphere;
    auto bvh = BuildMeshBVH(sphereC, white);
    {
        SurfaceInteraction rec;
        const bool hit = bvh->hit(Ray(Point3(0, 0, 5), Vector3(0, 0, -1)), rec);
        const bool ok = hit && rec.t > 3.98 && rec.t <= 4.0 + 1e-9 && glm::dot(rec.n, Vector3(0, 0, 1)) > 0.999 && glm::dot(rec.gn, Vector3(0, 0, 1)) > 0.99;
        check("ray at the pole: t, shading normal . z, geometric normal . z", ok,
            hit ? std::to_string(rec.t) + ", " + std::to_string(glm::dot(rec.n, Vector3(0, 0, 1))) + ", " + std::to_string(glm::dot(rec.gn, Vector3(0, 0, 1))) : std::string("no hit"));
    }

    // c. OBJ round trip
    {
        const std::string path = (std::filesystem::temp_directory_path() / "rayt_test_mesh.obj").string();
        {
            std::ofstream f(path);
            f.precision(12);
            for (const auto& p : sphere->p) f << "v " << p.x << " " << p.y << " " << p.z << "\n";
            for (const auto& n : sphere->n) f << "vn " << n.x << " " << n.y << " " << n.z << "\n";
            for (size_t i = 0; i < sphere->numTriangles(); ++i)
                f << "f " << sphere->indices[3 * i] + 1 << "//" << sphere->indices[3 * i] + 1 << " "
                  << sphere->indices[3 * i + 1] + 1 << "//" << sphere->indices[3 * i + 1] + 1 << " "
                  << sphere->indices[3 * i + 2] + 1 << "//" << sphere->indices[3 * i + 2] + 1 << "\n";
        }
        TriangleMesh back;
        std::string err;
        const bool loaded = io::loadObj(path, back, &err);
        const bool ok = loaded && back.numTriangles() == sphere->numTriangles() && back.hasNormals()
            && std::abs(back.signedVolume() - sphere->signedVolume()) < 1e-9;
        check("OBJ write + load (triangles, normals, volume)", ok,
            loaded ? std::to_string(back.numTriangles()) + " triangles" : err);
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }

    // d. BVH vs brute force on the torus file (2304 triangles): same hit / miss / distance for random rays
    {
        TriangleMesh torus;
        std::string err;
        if (!io::loadObj("assets/Meshes/torus.obj", torus, &err)) {
            check("torus.obj", false, err);
        }
        else {
            auto mesh = std::make_shared<const TriangleMesh>(std::move(torus));
            std::vector<std::shared_ptr<Triangle>> tris;
            auto torusBvh = BuildMeshBVH(mesh, white, &tris);
            Aggregate brute;
            for (const auto& t : tris) brute.add(t);
            int mismatches = 0, hits = 0;
            for (int i = 0; i < 20000; ++i) {
                const Vector3 target(sampling::Random(-0.7, 0.7), sampling::Random(-0.25, 0.25), sampling::Random(-0.7, 0.7));
                const Vector3 origin = glm::normalize(Vector3(sampling::Random(-1, 1), sampling::Random(-1, 1), sampling::Random(-1, 1))) * 3.0;
                const Ray ray(origin, glm::normalize(target - origin));
                SurfaceInteraction a, b;
                const bool ha = torusBvh->hit(ray, a), hb = brute.hit(ray, b);
                hits += hb ? 1 : 0;
                if (ha != hb || (ha && std::abs(a.t - b.t) > 1e-9)) ++mismatches;
                // any-hit with a random shortened range must agree with the closest-hit distance
                Ray shortRay = ray;
                shortRay.tMax = sampling::Random(0.5, 5.0);
                if (torusBvh->occluded(shortRay) != brute.occluded(shortRay) || brute.occluded(shortRay) != (hb && b.t <= shortRay.tMax)) ++mismatches;
            }
            check("torus.obj " + std::to_string(mesh->numTriangles()) + " triangles, volume " + std::to_string(mesh->signedVolume())
                + ": BVH vs brute force on 20000 rays (" + std::to_string(hits) + " hits), mismatches",
                mismatches == 0 && hits > 1000 && mesh->signedVolume() > 0, std::to_string(mismatches));

            // Embree (float traversal, Real hit record) vs brute force on the same kind of rays.
            {
                EmbreeMesh embree(mesh, white);
                int diff = 0, normalDiff = 0, n = 0;
                for (int i = 0; i < 20000; ++i) {
                    const Vector3 target(sampling::Random(-0.7, 0.7), sampling::Random(-0.25, 0.25), sampling::Random(-0.7, 0.7));
                    const Vector3 origin = glm::normalize(Vector3(sampling::Random(-1, 1), sampling::Random(-1, 1), sampling::Random(-1, 1))) * 3.0;
                    const Ray ray(origin, glm::normalize(target - origin));
                    SurfaceInteraction a, b;
                    const bool ha = embree.hit(ray, a), hb = brute.hit(ray, b);
                    ++n;
                    if (ha != hb || (ha && std::abs(a.t - b.t) > 1e-6)) { ++diff; continue; }
                    if (ha && (glm::dot(a.n, b.n) < 0.9999 || a.matPtr != b.matPtr)) ++normalDiff;
                    Ray shortRay = ray;
                    shortRay.tMax = sampling::Random(0.5, 5.0);
                    if (embree.occluded(shortRay) != brute.occluded(shortRay)) ++diff;
                }
                check("Embree vs brute force on " + std::to_string(n) + " rays (hit/miss, distance, any-hit; normals + material)",
                    diff <= 5 && normalDiff == 0, std::to_string(diff) + " ray(s) differ at float/double edge cases, " + std::to_string(normalDiff) + " normal/material differences");
            }
        }
    }

    // e. white (albedo 1) diffuse mesh sphere in a uniform white environment must stay exactly white
    {
        auto envMap = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        for (const bool smooth : { true, false }) {
            auto mesh = std::make_shared<const TriangleMesh>(TriangleMesh::Icosphere(2, 1.0, Point3(0.0), smooth));
            auto world = std::make_shared<Aggregate>();
            world->add(BuildMeshBVH(mesh, white));
            auto scene = std::make_shared<Scene>(world);
            scene->setEnvLight(std::make_shared<EnvLight>(envMap, 1.0));
            scene->finalize();
            double sum = 0;
            const int n = opt.spp * 256;
            Vector3 m, se;
            int bad = 0;
            EstimateRadiance(*scene, Ray(Point3(0.3, 0.2, 5), Vector3(0, 0, -1)), n, 16, m, se, bad);
            const bool ok = std::abs(m.x - 1.0) <= 4 * se.x + 0.01 && bad == 0;
            check(std::string("furnace, white diffuse mesh sphere, ") + (smooth ? "smooth normals" : "flat normals"), ok,
                "mean " + std::to_string(m.x) + " +- " + std::to_string(se.x));
            (void)sum;
        }
    }

    // f. two emissive triangles = the Quad light of --scene lights (numerical reference)
    {
        const Point3 corner(-0.5, 2, -0.5);
        const Vector3 e1(1, 0, 0), e2(0, 0, 1);
        auto mesh = std::make_shared<TriangleMesh>();
        mesh->p = { corner, corner + e1, corner + e1 + e2, corner + e2 };
        mesh->indices = { 0, 1, 2, 0, 2, 3 };
        const std::shared_ptr<const TriangleMesh> meshC = mesh;
        const Real albedo = 0.5;
        auto scene = MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
            auto mat = std::make_shared<DiffuseLight>(EmissionSpectrum(RgbColor(1.0)), 10.0);
            for (size_t i = 0; i < meshC->numTriangles(); ++i) {
                auto tri = std::make_shared<Triangle>(meshC, i, mat);
                auto light = std::make_shared<AreaLight>(tri, mat);
                w.add(std::make_shared<Primitive>(tri, mat, light));
                s.addLight(light);
            }
        });
        Vector3 m, se;
        int bad = 0;
        EstimateRadiance(*scene, Ray(Point3(-3, 3, 0), glm::normalize(Vector3(1, -1, 0))), opt.spp * 1024, 8, m, se, bad);
        const Real expected = QuadGroundRadiance(corner, e1, e2, false, 10.0, albedo);
        check("emissive triangle pair (area lights, NEE+MIS) vs quad reference " + std::to_string(expected),
            std::abs(m.x - expected) <= 4 * se.x + 0.005 * expected && bad == 0,
            std::to_string(m.x) + " +- " + std::to_string(se.x));
    }

    std::cout << "[Mesh] checks: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed;
}

// -----------------------------------------------------------------------------
// Main Entry Point
// -----------------------------------------------------------------------------
int main(int argc, char** argv) {

    CliOptions opt;
    if (!ParseArgs(argc, argv, opt)) return 1;
    const int imageWidth = opt.width;
    // The pmx scene is a standing figure: 3:4 portrait instead of 16:9.
    const double heightPerWidth = (opt.scene == "pmx") ? 4.0 / 3.0 : double(IMAGE_HEIGHT) / IMAGE_WIDTH;
    const int imageHeight = std::max(1, static_cast<int>(std::lround(double(opt.width) * heightPerWidth)));
    std::cout << "[Render] mode = " << kRenderModeName << "  (switch: Core/RenderMode.hpp)\n";
    std::cout << "[Options] scene=" << opt.scene << " " << imageWidth << "x" << imageHeight
        << " spp=" << opt.spp << " maxDepth=" << opt.maxDepth << " out=" << opt.out << "\n";

    if (opt.scene == "lights") return RunLightTests(opt);
    if (opt.scene == "aov") return RunAovTests();
    if (opt.scene == "metals") return RunMetalTests(opt);
    if (opt.scene == "bench") return RunBvhBenchmark(opt);
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
    bool meshChecksPassed = true;
    if (opt.scene == "mesh") meshChecksPassed = RunMeshChecks(opt);

    // --------------------------------------------------
    // 2. 環境マップ (EnvMap) 読み込み
    // --------------------------------------------------
    std::cout << "CWD = " << std::filesystem::current_path() << std::endl;

    std::shared_ptr<EnvMap> envMap = nullptr;  // 環境マップクラスをインスタンス化
    std::shared_ptr<EnvLight> envLight = nullptr;

    if (opt.scene == "furnace") {
        // Uniform white environment: radiance (1,1,1) in every direction.
        envMap = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        envLight = std::make_shared<EnvLight>(envMap, 1.0);
        std::cout << "[EnvMap] Uniform white environment (furnace test).\n";
    }
    else if (opt.scene == "spd" || opt.scene == "studio" || opt.scene == "mesh") {
        std::cout << "[EnvMap] none (spd scene: black background, one sphere light).\n";
    }
    else try {
        std::cout << "[EnvMap] Loading: " << opt.env << std::endl;
        auto envImg = io::loadHDR(opt.env);   // HDRでダウンロード
        envMap = std::make_shared<EnvMap>(std::move(envImg));
        std::cout << "[EnvMap] Success.\n";

        // EnvMapからEnvLightを作成 (スケール1.0)
        envLight = std::make_shared<EnvLight>(envMap, 1.0);
        std::cout << "[EnvLight] Created from environment map.\n";
    }
    catch (const std::exception& e) {
        std::cerr << "[EnvMap] Failed: " << e.what() << "\n";
        std::cerr << "[EnvMap] Fallback to black background.\n";
        envMap = nullptr;
        envLight = nullptr;
    }

    //-------------------------------------------------------------------------------------------
    // ラフのテスト
    // -------------------------------------------------------------------------
    // 1. マテリアルの作成 (Roughness Test)
    // -------------------------------------------------------------------------

    // 床用
    auto matFloor = std::make_shared<Lambertian>(RgbColor(0.5, 0.5, 0.5)); // 少し暗くして反射を目立たせる

    // --------------------------------------------------
    // 金の光学定数 (Au) をCSVからロード
    // The table is handed to the material as-is: spectral mode evaluates it per wavelength,
    // RGB mode reduces it to RGB eta/k by CMF integration.
    // --------------------------------------------------
    auto auTable = std::make_shared<rayt::color::ComplexIorTable>();
    std::string auErr;

    if (!rayt::io::loadComplexIor_RefractiveIndexInfoCsv(AU_IOR_CSV_PATH, *auTable, &auErr)) {
        std::cerr << "[Au IOR] Failed: " << auErr << "\n";
        std::cerr << "[Au IOR] Fallback to hard-coded eta/k.\n";

        // Fallback: a tiny table built from the old hard-coded RGB constants (450/550/650 nm).
        auTable->assign({
            { 450.0, 1.45, 1.77 },
            { 550.0, 0.42, 2.45 },
            { 650.0, 0.16, 3.48 },
        });
    }
    else {
        std::cout << "[Au IOR] Loaded: " << AU_IOR_CSV_PATH << "\n";
    }

    // Log the RGB reduction of the table (what RGB mode renders with).
    {
        Vector3 eta_rgb(1.0), k_rgb(0.0);
        rayt::color::integrateIorToRgb_CMF(*auTable, eta_rgb, k_rgb);
        std::cout << "[Au IOR] RGB (CMF-integrated) eta=("
            << eta_rgb.x << "," << eta_rgb.y << "," << eta_rgb.z << ") "
            << "k=(" << k_rgb.x << "," << k_rgb.y << "," << k_rgb.z << ")\n";
        std::cout << "[Au IOR] F0 (normal incidence) = ("
            << rayt::color::conductorNormalIncidenceReflectance(eta_rgb.x, k_rgb.x) << ","
            << rayt::color::conductorNormalIncidenceReflectance(eta_rgb.y, k_rgb.y) << ","
            << rayt::color::conductorNormalIncidenceReflectance(eta_rgb.z, k_rgb.z) << ")\n";
    }

    // ★比較用: 3段階の粗さを作成
    // 0.01: ほぼ鏡 (MirrorConductorと比較用)
    // 0.20: 少しぼやけた金属
    // 0.50: マットな金属（ブラスト仕上げ風）
    auto matGoldSmooth = std::make_shared<RoughConductor>(auTable, 0.01);
    auto matGoldMedium = std::make_shared<RoughConductor>(auTable, 0.20);
    auto matGoldRough = std::make_shared<RoughConductor>(auTable, 0.50);

    // Glass and reference materials for the "glass" and "furnace" scenes.
    auto matGlass = std::make_shared<Dielectric>(1.5, 0.0);                      // 粗さ0 = 完全透明 (no dispersion)
    auto matGlassDisp = std::make_shared<Dielectric>(1.5, 0.0, 0.0, opt.abbe);   // clear glass with dispersion (spectral mode)
    auto matFrosted = std::make_shared<Dielectric>(1.5, 0.2);                    // 粗さ0.2 = すりガラス
    auto matWhite = std::make_shared<Lambertian>(RgbColor(1.0));                 // albedo 1 (furnace reference)

    // -------------------------------------------------------------------------
    // 2. 物体の配置 (Scene)
    // -------------------------------------------------------------------------
    auto worldList = std::make_shared<Aggregate>();

    EmissionSpectrum spdEmission;                                 // spd scene: spectrum of its sphere light
    Point3 pmxCamFrom(0.0), pmxCamTarget(0.0);                    // pmx scene: camera framing the model
    bool pmxCamera = false;
    std::vector<std::shared_ptr<Light>> studioLights;             // studio scene lights, registered below
    std::shared_ptr<Light> spdLight;                              // spd scene: registered for NEE below
    std::vector<FurnaceObject> furnaceObjects;                    // what the furnace test checks
    std::vector<std::shared_ptr<Material>> furnaceMaterials;      // owns materials only used there

    if (opt.scene == "gold") {
        // リストに対して add する
        worldList->add(std::make_shared<Sphere>(Point3(-1.2, 0, -1), 0.5, matGoldSmooth));
        worldList->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), 0.5, matGoldMedium));
        worldList->add(std::make_shared<Sphere>(Point3(1.2, 0, -1), 0.5, matGoldRough));

        // 大きな球で床を作製
        worldList->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));
    }
    else if (opt.scene == "pmx") {
        const auto files = FindPmxFiles();
        if (files.empty()) { std::cerr << "[PMX] no .pmx files found in MMD_model/\n"; return 1; }
        size_t pick = 0;
        if (opt.pmxIndex >= 0) pick = std::min(static_cast<size_t>(opt.pmxIndex), files.size() - 1);
        else for (size_t i = 1; i < files.size(); ++i) if (std::filesystem::file_size(files[i]) > std::filesystem::file_size(files[pick])) pick = i;

        io::PmxModel model;
        std::string pmxErr;
        if (!io::loadPmx(files[pick], model, &pmxErr)) { std::cerr << "[PMX] load failed: " << pmxErr << "\n"; return 1; }
        PrintPmxInfo(model, files[pick]);
        if (opt.pmxCheck) return RunPmxChecks(model) ? 0 : 2;
        if (opt.morphList) {
            static const char* kTypes[] = { "group", "vertex", "bone", "uv", "uv1", "uv2", "uv3", "uv4", "material", "flip", "impulse" };
            static const char* kPanels[] = { "-", "brow", "eye", "mouth", "other" };
            for (size_t i = 0; i < model.morphs.size(); ++i) {
                const auto& m = model.morphs[i];
                std::cout << "[Morph] " << i << ": " << m.name << "  [" << kPanels[std::clamp(m.panel, 0, 4)] << ", " << kTypes[std::clamp(m.type, 0, 10)]
                    << ", " << m.offsetCount() << " offsets]\n";
            }
            return 0;
        }

        // ---- pose file (bones, and expression entries) ----
        io::VpdPose vpd;
        bool havePose = false;
        std::string poseLabel;
        if (opt.poseIndex >= 0) {
            const auto poses = FindMmdFiles("pose", ".vpd");
            if (poses.empty()) { std::cerr << "[Pose] no .vpd files found in MMD_model/pose\n"; return 1; }
            const auto& poseFile = poses[std::min(static_cast<size_t>(opt.poseIndex), poses.size() - 1)];
            std::string poseErr;
            if (!io::loadVpd(poseFile, vpd, &poseErr)) { std::cerr << "[Pose] load failed: " << poseErr << "\n"; return 1; }
            havePose = true;
            poseLabel = PathToUtf8(poseFile.parent_path().filename()) + "/" + PathToUtf8(poseFile.filename());
        }

        // ---- morph weights: --morph i=w,..  /  --expr name  /  Morph blocks of the pose file ----
        std::map<int, Real> morphWeights;
        auto morphByName = [&](const std::string& name) {
            for (size_t i = 0; i < model.morphs.size(); ++i) if (model.morphs[i].name == name) return static_cast<int>(i);
            return -1;
        };
        auto setMorph = [&](const std::string& key, Real w, const char* from) {
            int idx = -1;
            const bool numeric = !key.empty() && std::all_of(key.begin(), key.end(), [](unsigned char c) { return std::isdigit(c) != 0; });
            if (numeric) idx = std::stoi(key);
            else idx = morphByName(key);
            if (idx < 0 || idx >= static_cast<int>(model.morphs.size())) { std::cout << "[Morph] " << from << ": no morph " << key << " in this model\n"; return; }
            morphWeights[idx] = w;
        };
        if (!opt.morphSpec.empty()) {
            std::stringstream ss(opt.morphSpec);
            std::string item;
            while (std::getline(ss, item, ',')) {
                const size_t eq = item.find('=');
                if (eq != std::string::npos) setMorph(item.substr(0, eq), std::stod(item.substr(eq + 1)), "--morph");
            }
        }
        if (!opt.exprFile.empty()) {
            const auto files = FindMmdFiles("expression", ".txt");
            const auto it = std::find_if(files.begin(), files.end(), [&](const std::filesystem::path& f) { return PathToUtf8(f.stem()) == opt.exprFile; });
            std::vector<std::pair<std::string, Real>> items;
            std::string exprErr;
            if (it == files.end() || !io::loadExpression(*it, items, &exprErr)) { std::cerr << "[Morph] expression " << opt.exprFile << " not found in MMD_model/expression\n"; return 1; }
            for (const auto& [name, w] : items) setMorph(name, w, "--expr");
        }
        if (havePose) for (const auto& [name, w] : vpd.morphs) setMorph(name, w, "pose file");

        std::vector<io::VpdBone> morphBones;
        if (!morphWeights.empty()) {
            model = io::applyMorphs(model, morphWeights, &morphBones);
            std::cout << "[Morph] applied:";
            for (const auto& [idx, w] : morphWeights) std::cout << " " << model.morphs[static_cast<size_t>(idx)].name << "=" << w;
            std::cout << "\n";
        }
        if (havePose || !morphBones.empty()) {
            io::mergePoseBones(vpd, morphBones);
            io::PoseReport rep;
            model = io::applyPose(model, vpd, &rep);
            std::cout << "[Pose] " << (havePose ? poseLabel : std::string("(bone morphs only)")) << "  made for " << vpd.parentModel
                << ": " << vpd.bones.size() << " bones, " << rep.applied << " applied";
            if (!rep.missing.empty()) {
                std::cout << ", " << rep.missing.size() << " not in this model:";
                for (const auto& n : rep.missing) std::cout << " " << n;
            }
            std::cout << "\n";
        }
        if (opt.pmxInfoOnly) return 0;

        Real height = 1;
        Vector3 modelSize(1.0);
        const std::shared_ptr<TriangleMesh> editableMesh = PmxToMesh(model, opt.pmxYaw, &height, &modelSize);
        {
            // A translucent layer (blush ...) is a copy of the surface it covers, lying exactly on it. A ray that
            // passes through the layer must still find the surface below, so the layer is lifted 0.5 mm along
            // its normals (far more than the 1e-5 ray offset, invisible in the picture).
            std::vector<char> lifted(editableMesh->p.size(), 0);
            for (const auto& pm : model.materials) {
                const bool overlayName = !pm.name.empty() && pm.name.back() == '+';
                if (pm.diffuse.w <= Real(0.001) || !(pm.alphaFromMorph && (overlayName || pm.diffuse.w < Real(0.999)))) continue;
                for (int i = pm.firstIndex; i < pm.firstIndex + pm.indexCount; ++i) {
                    const size_t v = static_cast<size_t>(editableMesh->indices[static_cast<size_t>(i)]);
                    if (!lifted[v]) { editableMesh->p[v] += editableMesh->n[v] * Real(5e-4); lifted[v] = 1; }
                }
            }
        }
        const std::shared_ptr<const TriangleMesh> mesh = editableMesh;
        std::vector<std::shared_ptr<Hittable>> prims;
        std::map<int, std::shared_ptr<Texture2D>> textureCache;   // PMX texture index -> loaded texture (nullptr = failed)
        int skippedOverlay = 0, skippedHidden = 0, alphaMaterials = 0, missingTextures = 0, translucentLayers = 0;
        for (const auto& pm : model.materials) {
            if (pm.diffuse.w <= Real(0.001)) { ++skippedHidden; continue; }   // fully transparent material
            // MMD models stack extra layers on top of the real surface and switch them with morphs; names
            // ending in '+' mark them (sphere-map highlights, blush ...). A layer a morph has switched on
            // (diffuse alpha changed by a morph) is drawn as a translucent layer; the others would cover the
            // real surface as opaque diffuse, so they are left out unless --pmx-all.
            const bool overlayName = !pm.name.empty() && pm.name.back() == '+';
            bool translucent = pm.alphaFromMorph && (overlayName || pm.diffuse.w < Real(0.999));
            if (overlayName && !pm.alphaFromMorph && !opt.pmxAll) { ++skippedOverlay; continue; }

            std::shared_ptr<Texture2D> tex;
            if (pm.texture >= 0) {
                auto it = textureCache.find(pm.texture);
                if (it == textureCache.end()) {
                    std::string rel = model.textures[static_cast<size_t>(pm.texture)];
                    std::replace(rel.begin(), rel.end(), '\\', '/');
                    it = textureCache.emplace(pm.texture, Texture2D::load(files[pick].parent_path() / PathFromUtf8(rel))).first;
                    if (!it->second) { ++missingTextures; std::cout << "[PMX] could not load texture " << rel << "\n"; }
                }
                tex = it->second;
            }

            const RgbColor baseColor(SrgbToLinear(pm.diffuse.x), SrgbToLinear(pm.diffuse.y), SrgbToLinear(pm.diffuse.z));
            std::shared_ptr<Material> mat;
            if (translucent) {
                auto layer = std::make_shared<Lambertian>(baseColor, tex);
                layer->setTranslucent(std::clamp(pm.diffuse.w, Real(0), Real(1)));
                mat = layer;
                ++translucentLayers;
            }
            else {
                mat = MakePmxMaterial(pm, baseColor, tex, opt.pmxMatte);
            }
            const bool cutout = tex && tex->hasAlpha() && !translucent;   // a translucent layer uses the texture alpha as opacity instead
            if (cutout) ++alphaMaterials;
            for (int t = pm.firstIndex / 3; t < (pm.firstIndex + pm.indexCount) / 3; ++t) {
                const Vector3 e1 = mesh->p[static_cast<size_t>(mesh->indices[3 * t + 1])] - mesh->p[static_cast<size_t>(mesh->indices[3 * t])];
                const Vector3 e2 = mesh->p[static_cast<size_t>(mesh->indices[3 * t + 2])] - mesh->p[static_cast<size_t>(mesh->indices[3 * t])];
                if (glm::length(glm::cross(e1, e2)) < 1e-12) continue;   // degenerate triangle
                auto tri = std::make_shared<Triangle>(mesh, static_cast<size_t>(t), mat);
                if (cutout) tri->setAlphaTexture(tex, Real(0.5));
                if (translucent) tri->setCastsShadow(false);
                prims.push_back(tri);
            }
        }
        size_t texBytes = 0;
        for (const auto& [idx, t] : textureCache) if (t) texBytes += static_cast<size_t>(t->width()) * static_cast<size_t>(t->height()) * 4;
        std::cout << "[PMX] " << prims.size() << " triangles in the scene, height " << height << " m;  " << textureCache.size() - missingTextures
            << " textures loaded (" << texBytes / (1024 * 1024) << " MB), " << alphaMaterials << " materials with alpha cut-outs, skipped "
            << skippedOverlay << " overlay ('+') and " << skippedHidden << " fully transparent materials, " << translucentLayers << " translucent layers\n";
        worldList->add(std::make_shared<BVH>(std::move(prims)));
        worldList->add(std::make_shared<Sphere>(Point3(0, -100, 0), 100.0, matFloor));
        // Framing: the picture is 3:4 with a 35 degree vertical field of view; a wide pose (a raised leg) is
        // fitted by its width as well as its height.
        const Real zoom = std::max(opt.pmxZoom, Real(0.1));
        const Real extent = std::max(height, modelSize.x / Real(0.75));
        pmxCamTarget = Point3(0, height * opt.pmxFocus, 0);
        pmxCamFrom = Point3(0, height * (opt.pmxFocus + 0.05 / zoom), extent * 1.22 / (2 * std::tan(glm::radians(17.5))) / zoom);
        pmxCamera = true;
        if (opt.pmxStudio) {
            // Studio: a dim backdrop, key / fill soft boxes and two rim lights; the HDR only adds a faint ambient.
            // Looks (--studio-look): color temperature (K; 0 = neutral white) and radiance of each light,
            // the ambient HDR scale and the backdrop reflectance. Positions and sizes are shared.
            struct Lamp { Real kelvin, radiance; };
            struct StudioLook { Lamp key, fill, rimLeft, rimRight; Real env; RgbColor wall; };
            static const std::map<std::string, StudioLook> looks = {
                { "default",  { { 5200, 13 },   { 0, 3.2 },     { 0, 28 },      { 7500, 22 },    0.10, RgbColor(0.16, 0.17, 0.21) } },
                { "warm",     { { 3200, 14 },   { 4500, 2.0 },  { 2500, 30 },   { 3800, 22 },    0.08, RgbColor(0.20, 0.16, 0.13) } },
                { "cool",     { { 8500, 8 },    { 10000, 3.0 }, { 12000, 28 },  { 9000, 22 },    0.12, RgbColor(0.13, 0.16, 0.22) } },
                { "dramatic", { { 4300, 20 },   { 0, 0.4 },     { 6000, 35 },   { 6500, 28 },    0.02, RgbColor(0.07, 0.07, 0.09) } },
                { "soft",     { { 5600, 7 },    { 5600, 6.0 },  { 5600, 10 },   { 5600, 8 },     0.35, RgbColor(0.45, 0.46, 0.50) } },
            };
            const StudioLook& look = looks.at(opt.studioLook);
            auto spectrum = [](const Lamp& l) {
                return l.kelvin > 0 ? EmissionSpectrum::Blackbody(l.kelvin) : EmissionSpectrum(RgbColor(1.0));
            };
            if (envLight) envLight->setScale(look.env);
            const Point3 target(0, height * 0.55, 0);
            auto softbox = [&](const Point3& center, Real w, Real h, const EmissionSpectrum& e, Real radiance) {
                const Vector3 n = glm::normalize(target - center);
                const Vector3 e1 = glm::normalize(glm::cross(Vector3(0, 1, 0), n)) * w;
                const Vector3 e2 = glm::normalize(glm::cross(n, e1)) * h;
                studioLights.push_back(AddQuadLight(*worldList, nullptr, center - e1 * 0.5 - e2 * 0.5, e1, e2, e, radiance, true));
            };
            const Real H = height;
            softbox(Point3(-1.7, H * 1.05, 2.3), 1.3, 1.3, spectrum(look.key), look.key.radiance);               // key: front-left, above
            softbox(Point3(2.3, H * 0.60, 2.4), 1.8, 1.2, spectrum(look.fill), look.fill.radiance);              // fill: front-right, soft
            softbox(Point3(-1.5, H * 1.10, -1.7), 0.5, 1.8, spectrum(look.rimLeft), look.rimLeft.radiance);      // rim: back-left strip
            softbox(Point3(1.8, H * 1.00, -1.6), 0.5, 1.8, spectrum(look.rimRight), look.rimRight.radiance);     // rim: back-right strip
            const Real wallZ = -(modelSize.z * 0.5 + 1.8);
            worldList->add(std::make_shared<Quad>(Point3(-14, 0, wallZ), Vector3(28, 0, 0), Vector3(0, 12, 0),
                std::make_shared<Lambertian>(look.wall)));
            // Dark glossy studio floor (normal +y: cross(z, x)) lying 1 mm above the plain floor: soft reflections of the lights and the figure.
            worldList->add(std::make_shared<Quad>(Point3(-14, 0.001, wallZ), Vector3(0, 0, 16), Vector3(28, 0, 0),
                std::make_shared<CoatedDiffuse>(RgbColor(0.09, 0.09, 0.11), Real(0.20), Real(1.0), Real(0))));
        }
        if (opt.pmxLight) {
            // Softbox 1.0 x 0.7 m, 2.2 m from the focus point, front-left and a little above, facing it.
            const Point3 target = pmxCamTarget;
            const Point3 center = target + Vector3(-1.3, 0.9, 1.6);
            const Vector3 n = glm::normalize(target - center);
            const Vector3 e1 = glm::normalize(glm::cross(Vector3(0, 1, 0), n)) * 1.0;
            const Vector3 e2 = glm::normalize(glm::cross(n, e1)) * 0.7;
            studioLights.push_back(AddQuadLight(*worldList, nullptr, center - e1 * 0.5 - e2 * 0.5, e1, e2,
                EmissionSpectrum(RgbColor(1.0)), 8.0, true));
        }
    }
    else if (opt.scene == "goldrough") {
        // Rough gold spheres (roughness 0.4 / 0.7 / 1.0): multiple scattering matters here (--no-ms turns it off).
        for (const auto& [x, roughness] : { std::pair<Real, Real>{ -1.2, 0.4 }, std::pair<Real, Real>{ 0.0, 0.7 }, std::pair<Real, Real>{ 1.2, 1.0 } }) {
            auto m = std::make_shared<RoughConductor>(auTable, roughness);
            m->setMultipleScattering(!opt.noMs);
            worldList->add(std::make_shared<Sphere>(Point3(x, 0, -1), 0.5, m));
        }
        worldList->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));
    }
    else if (opt.scene == "glass") {
        // left: clear glass, center: dispersive glass, right: frosted glass
        worldList->add(std::make_shared<Sphere>(Point3(-1.2, 0, -1), 0.5, matGlass));
        worldList->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), 0.5, matGlassDisp));
        worldList->add(std::make_shared<Sphere>(Point3(1.2, 0, -1), 0.5, matFrosted));
        worldList->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));
    }
    else if (opt.scene == "studio") {
        const Real r = 0.45;
        worldList->add(std::make_shared<Sphere>(Point3(-1.15, 0, -1), r, matGoldMedium));
        worldList->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), r, matGlassDisp));
        worldList->add(std::make_shared<Sphere>(Point3(1.15, 0, -1), r, std::make_shared<Lambertian>(RgbColor(0.75, 0.2, 0.18))));
        worldList->add(std::make_shared<Sphere>(Point3(0, -100.45, -1), 100.0, std::make_shared<Lambertian>(RgbColor(0.6, 0.6, 0.6))));

        // Tungsten softbox above and behind the spheres (emits downward: cross(x, z) = -y).
        studioLights.push_back(AddQuadLight(*worldList, nullptr, Point3(-1.8, 2.6, -2.0),
            Vector3(2.4, 0, 0), Vector3(0, 0, 1.8), EmissionSpectrum::IlluminantA(), 7.0, true));

        // Small cool-white LED (synthetic SPD from CSV) on the right.
        color::SpectralCurve ledCurve;
        std::string ledErr;
        if (!io::loadSpectralCurveCsv("assets/Illuminants/led_synthetic_coolwhite.csv", ledCurve, &ledErr)) {
            std::cerr << "[SPD] " << ledErr << "\n";
            return 1;
        }
        auto ledMat = std::make_shared<DiffuseLight>(EmissionSpectrum::FromCurve(ledCurve), 30.0);
        auto ledShape = std::make_shared<Sphere>(Point3(2.4, 1.4, 0.2), 0.3, ledMat);
        auto ledLight = std::make_shared<AreaLight>(ledShape, ledMat);
        worldList->add(std::make_shared<Primitive>(ledShape, ledMat, ledLight));
        studioLights.push_back(ledLight);
    }
    else if (opt.scene == "mesh") {
        // Gold torus from an OBJ file, tilted toward the camera.
        TriangleMesh torus;
        std::string objErr;
        if (!io::loadObj("assets/Meshes/torus.obj", torus, &objErr)) {
            std::cerr << "[Mesh] " << objErr << "\n";
            return 1;
        }
        torus.transform(glm::dmat3(glm::rotate(glm::dmat4(1.0), glm::radians(Real(65)), Vector3(1, 0, 0))), Vector3(0.0, 0.22, -1.2));
        worldList->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(torus)), matGoldMedium, opt.accel));

        // The same 80-triangle sphere twice: flat shading (left) vs interpolated vertex normals (right).
        auto red = std::make_shared<Lambertian>(RgbColor(0.75, 0.2, 0.18));
        worldList->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(TriangleMesh::Icosphere(1, 0.4, Point3(-1.45, 0, -0.9), false)), red, opt.accel));
        worldList->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(TriangleMesh::Icosphere(1, 0.4, Point3(1.45, 0, -0.9), true)), red, opt.accel));
        worldList->add(std::make_shared<Sphere>(Point3(0, -100.45, -1), 100.0, std::make_shared<Lambertian>(RgbColor(0.6, 0.6, 0.6))));

        studioLights.push_back(AddQuadLight(*worldList, nullptr, Point3(-1.8, 2.6, -2.0),
            Vector3(2.4, 0, 0), Vector3(0, 0, 1.8), EmissionSpectrum::IlluminantA(), 7.0, true));
        color::SpectralCurve ledCurve;
        std::string ledErr;
        if (!io::loadSpectralCurveCsv("assets/Illuminants/led_synthetic_coolwhite.csv", ledCurve, &ledErr)) {
            std::cerr << "[SPD] " << ledErr << "\n";
            return 1;
        }
        auto ledMat = std::make_shared<DiffuseLight>(EmissionSpectrum::FromCurve(ledCurve), 30.0);
        auto ledShape = std::make_shared<Sphere>(Point3(2.4, 1.4, 0.2), 0.3, ledMat);
        auto ledLight = std::make_shared<AreaLight>(ledShape, ledMat);
        worldList->add(std::make_shared<Primitive>(ledShape, ledMat, ledLight));
        studioLights.push_back(ledLight);
    }
    else if (opt.scene == "spd") {
        // Five diffuse spheres with distinct reflectance spectra + gold, lit by one off-screen sphere light.
        const Real r = 0.3;
        auto lambertian = [](Real R, Real G, Real B) { return std::make_shared<Lambertian>(RgbColor(R, G, B)); };
        worldList->add(std::make_shared<Sphere>(Point3(-1.4, 0, -1), r, lambertian(0.80, 0.15, 0.12)));   // red
        worldList->add(std::make_shared<Sphere>(Point3(-0.7, 0, -1), r, lambertian(0.80, 0.80, 0.80)));   // white
        worldList->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), r, matGoldMedium));                    // gold
        worldList->add(std::make_shared<Sphere>(Point3(0.7, 0, -1), r, lambertian(0.12, 0.30, 0.80)));    // blue
        worldList->add(std::make_shared<Sphere>(Point3(1.4, 0, -1), r, lambertian(0.15, 0.60, 0.20)));    // green
        worldList->add(std::make_shared<Sphere>(Point3(0, -100.3, -1), 100.0, lambertian(0.7, 0.7, 0.7)));

        if (opt.light == "d65")           spdEmission = EmissionSpectrum(RgbColor(1.0));
        else if (opt.light == "tungsten") spdEmission = EmissionSpectrum::IlluminantA();
        else if (opt.light == "candle")   spdEmission = EmissionSpectrum::Blackbody(1850.0);
        else {
            color::SpectralCurve curve;
            std::string err;
            if (!io::loadSpectralCurveCsv("assets/Illuminants/led_synthetic_coolwhite.csv", curve, &err)) {
                std::cerr << "[SPD] " << err << "\n";
                return 1;
            }
            spdEmission = EmissionSpectrum::FromCurve(curve);
        }
        std::cout << "[SPD] light = " << opt.light << ", sensor RGB of the light = ("
            << spdEmission.rgb().x << ", " << spdEmission.rgb().y << ", " << spdEmission.rgb().z << ")"
            << (opt.whiteBalance ? ", film white-balanced to it" : "") << "\n";
        if (opt.whiteBalance) color::setWhiteBalance(spdEmission);

        // Radiance 20 at distance ~3.5 gives a floor radiance of about 0.6.
        auto lightMat = std::make_shared<DiffuseLight>(spdEmission, 20.0);
        auto lightShape = std::make_shared<Sphere>(Point3(-2.2, 3.0, 1.5), 0.8, lightMat);
        spdLight = std::make_shared<AreaLight>(lightShape, lightMat);
        worldList->add(std::make_shared<Primitive>(lightShape, lightMat, spdLight));
    }
    else {
        // furnace: no floor, so the only light path is environment <-> object.
        // The two metals use a near-perfect conductor (k = 100, F ~ 0.9996 at every angle), so
        // only the microfacet model's own single-scattering energy loss is left.
        auto perfectMetal = std::make_shared<rayt::color::ComplexIorTable>();
        perfectMetal->assign({ { 360.0, 1.0, 100.0 }, { 830.0, 1.0, 100.0 } });
        auto matMetalMid = std::make_shared<RoughConductor>(perfectMetal, 0.2);
        auto matMetalRough = std::make_shared<RoughConductor>(perfectMetal, 0.5);
        furnaceMaterials = { matMetalMid, matMetalRough };   // keep them alive until the report

        const Real r = 0.28;
        worldList->add(std::make_shared<Sphere>(Point3(-1.55, 0, -1), r, matGlass));
        worldList->add(std::make_shared<Sphere>(Point3(-0.93, 0, -1), r, matGlassDisp));
        worldList->add(std::make_shared<Sphere>(Point3(-0.31, 0, -1), r, matFrosted));
        worldList->add(std::make_shared<Sphere>(Point3(0.31, 0, -1), r, matWhite));
        worldList->add(std::make_shared<Sphere>(Point3(0.93, 0, -1), r, matMetalMid));
        worldList->add(std::make_shared<Sphere>(Point3(1.55, 0, -1), r, matMetalRough));
        // expected / tolerance:  1.0 for lossless materials;  the near-perfect conductors use an
        // independent numerical integration of single-scattering GGX (0.994 / 0.882);  frosted
        // glass uses the measured value 0.986 (no closed form).
        furnaceObjects = {
            { "background      ", nullptr,              1.000, 0.010 },
            { "clear glass     ", matGlass.get(),       1.000, 0.030 },
            { "dispersive glass", matGlassDisp.get(),   1.000, 0.030 },
            { "white diffuse   ", matWhite.get(),       1.000, 0.030 },
            { "frosted glass   ", matFrosted.get(),     0.986, 0.030 },
            { "metal rough 0.2 ", matMetalMid.get(),    1.000, 0.020 },
            { "metal rough 0.5 ", matMetalRough.get(),  1.000, 0.020 },
        };
    }

    // Area Light (メッシュに紐付いた光源) を作成してシーンに配置
    /*auto lightMat = std::make_shared<DiffuseLight>(Spectrum(50.0, 0.0, 0.0));
    auto lightShape = std::make_shared<Sphere>(Point3(5.0, 10.0, 6.0), 4.0, lightMat);
    auto areaLight = std::make_shared<AreaLight>(lightShape, lightMat, Spectrum(1.0));
    auto lightPrim = std::make_shared<Primitive>(lightShape, lightMat, areaLight);
    worldList->add(lightPrim);*/

    // --------------------------------------------------
    // 5. Scene作成と光源設定
    // --------------------------------------------------
    auto scene = std::make_shared<Scene>(worldList);

    // 環境光を設定（存在する場合）
    if (envLight) {
        scene->setEnvLight(envLight);
        std::cout << "[Scene] Environment light set.\n";
    }

    // オプション: ポイントライトや他の光源を追加
    // 例: シーン内にポイントライトを追加する場合

    /*auto pointLight = std::make_shared<PointLight>(
        Point3(-1.2, 3, 0),           // 位置
        Spectrum(0, 0, 100)       // 強度
    );
    scene->addLight(pointLight);
    std::cout << "[Scene] Point light added.\n";*/

    // エリアライトの設定
    //scene->addLight(areaLight);
    //std::cout << "[Scene] Area light added.\n";

    if (spdLight) scene->addLight(spdLight);
    for (const auto& l : studioLights) scene->addLight(l);

    // IMPORTANT: シーンをファイナライズ（LightSamplerを構築）
    scene->finalize(opt.lightSampler);
    std::cout << "[Scene] Finalized (light sampler: " << opt.lightSampler << ") with "
        << scene->numLights() << " scene lights";
    if (scene->hasEnvLight()) {
        std::cout << " + environment light";
    }
    std::cout << ".\n";

    // -------------------------------------------------------------------------
    // 3. カメラ設定
    // -------------------------------------------------------------------------
    Point3 lookFrom(0, 0.5, 2.5); // 少し高い位置から見下ろす
    Point3 lookAt(0, 0, -1);
    Vector3 vUp(0, 1, 0);
    if (pmxCamera) { lookFrom = pmxCamFrom; lookAt = pmxCamTarget; }

    Real distToFocus = glm::length(lookFrom - lookAt);
    Real aperture = 0.0; // ピンホールカメラ（ボケなし）でテスト

    auto camera = std::make_shared<Camera>(
        lookFrom, lookAt, vUp,
        35.0, // FOV
        double(imageWidth) / imageHeight,
        aperture,
        distToFocus
    );

    if (opt.traceX >= 0) {
        TracePixelDebug(*scene, *camera, opt.traceX, opt.traceY, imageWidth, imageHeight, 12, 24,
            { { "clear", matGlass.get() }, { "disp", matGlassDisp.get() }, { "frosted", matFrosted.get() },
              { "white", matWhite.get() }, { "floor", matFloor.get() },
              { "goldS", matGoldSmooth.get() }, { "goldM", matGoldMedium.get() }, { "goldR", matGoldRough.get() } });
        return 0;
    }

    // --------------------------------------------------
    // 7. レンダラー生成
    // --------------------------------------------------
    Film film(imageWidth, imageHeight);
    if (opt.aov || opt.denoise) film.enableAOVs();   // the denoiser uses them even if they are not written out

    // 新しいPathIntegratorのコンストラクタ: (camera, maxDepth, spp)
    auto integrator = std::make_unique<PathIntegrator>(
        camera,
        opt.maxDepth,  // 最大深度
        opt.spp        // プログレッシブ用のspp
    );

    std::cout << "[Integrator] PathIntegrator created (maxDepth="
        << opt.maxDepth << ", spp=" << opt.spp << ").\n";

    // -------------------------------------------------------------------------
    // 6. レンダリング実行 & 保存
    // -------------------------------------------------------------------------
    std::cout << "[Render] Start PBR rendering..." << std::endl;
    integrator->render(*scene, film);

    std::cout << "[Output] Saving images..." << std::endl;
    film.save(opt.out);
    // film.save("result_gold_pbr.hdr");
    const std::string outStem = std::filesystem::path(opt.out).replace_extension("").string();
    const std::string outExt = std::filesystem::path(opt.out).extension().string();
    if (opt.aov) {
        const std::string aovExt = (outExt == ".hdr") ? ".hdr" : ".png";
        film.saveAlbedo(outStem + "_albedo" + aovExt);
        film.saveNormal(outStem + "_normal" + aovExt);
    }
    if (opt.denoise) {
        Film denoised(imageWidth, imageHeight);
        if (DenoiseFilm(film, denoised, opt.cleanAux)) denoised.save(outStem + "_denoised" + outExt);
    }

    bool furnacePassed = true;
    if (opt.scene == "furnace") {
        furnacePassed = CheckFurnace(*scene, *camera, film, furnaceObjects);
    }

    std::cout << "[System] Finished." << std::endl;

    // Furnace test: exit code 2 when an object is outside its expected value, so a script can check it.
    return (furnacePassed && meshChecksPassed) ? 0 : 2;
}
