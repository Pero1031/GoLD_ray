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
#include "Geometry/Frame.hpp"
#include "Scene/Primitive.hpp"
#include "Scene/Scene.hpp"

// Renderer
#include "Render/Film.hpp"
#include "Render/Camera.hpp"
#include "Scene/Scene.hpp"
#include "Render/Integrator.hpp"

// Materials
#include "Materials/Material.hpp"
#include "Materials/Lambertian.hpp"
#include "Materials/DiffuseLight.hpp"
#include "Materials/RoughConductor.hpp"
#include "Materials/Dielectric.hpp"

// Lights
#include "Lights/AreaLight.hpp"
#include "Lights/PointLight.hpp"
#include "Color/EmissionSpectrum.hpp"

// IO
#include "IO/ImageLoader.hpp"
#include "IO/EnvMap.hpp"
#include "IO/ComplexIorLoader.hpp"
#include "IO/SpectralCurveLoader.hpp"
#include "Color/IorRgbApprox.hpp"

#include <filesystem>
#include <fstream>
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
//   lights  : light source tests (no image). A Lambertian ground is lit by one light and the
//             radiance of one ray is compared with a closed-form / independent numerical value
//             (point light, sphere area light with and without NEE, colored emission, environment
//             maps). Samples per test = --spp x 4096. Exits with code 2 on FAIL.
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
            else if (a == "--trace") {
                const std::string xy = value("--trace");
                const size_t comma = xy.find(',');
                if (comma == std::string::npos) throw std::runtime_error("--trace needs x,y");
                o.traceX = std::stoi(xy.substr(0, comma));
                o.traceY = std::stoi(xy.substr(comma + 1));
            }
            else throw std::runtime_error("unknown option: " + a);
        }
        if (o.scene != "gold" && o.scene != "glass" && o.scene != "furnace" && o.scene != "lights")
            throw std::runtime_error("--scene must be gold, glass, furnace or lights");
        if (o.width < 16 || o.spp < 1 || o.maxDepth < 1)
            throw std::runtime_error("--width must be >= 16, --spp and --depth >= 1");
    }
    catch (const std::exception& e) {
        std::cerr << "[Args] " << e.what() << "\n"
            << "usage: main_cli [--scene gold|glass|furnace|lights] [--width N] [--spp N] [--depth N]"
            << " [--env file.hdr] [--abbe V] [--out file.png] [--trace x,y]\n";
        return false;
    }
    if (o.out.empty()) o.out = (o.scene == "gold") ? "result_gold_pbr.png" : "result_" + o.scene + ".png";
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
    const std::function<void(Aggregate&, Scene&)>& addLights)
{
    auto world = std::make_shared<Aggregate>();
    world->add(std::make_shared<Sphere>(Point3(0, -1000, 0), 1000.0,
        std::make_shared<Lambertian>(RgbColor(albedo))));
    auto scene = std::make_shared<Scene>(world);
    addLights(*world, *scene);
    scene->finalize("uniform");
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
// Main Entry Point
// -----------------------------------------------------------------------------
int main(int argc, char** argv) {

    CliOptions opt;
    if (!ParseArgs(argc, argv, opt)) return 1;
    const int imageWidth = opt.width;
    const int imageHeight = std::max(1, static_cast<int>(std::lround(double(opt.width) * IMAGE_HEIGHT / IMAGE_WIDTH)));
    std::cout << "[Render] mode = " << kRenderModeName << "  (switch: Core/RenderMode.hpp)\n";
    std::cout << "[Options] scene=" << opt.scene << " " << imageWidth << "x" << imageHeight
        << " spp=" << opt.spp << " maxDepth=" << opt.maxDepth << " out=" << opt.out << "\n";

    if (opt.scene == "lights") return RunLightTests(opt);

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
        auTable->samples() = {
            { 450.0, 1.45, 1.77 },
            { 550.0, 0.42, 2.45 },
            { 650.0, 0.16, 3.48 },
        };
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
    else if (opt.scene == "glass") {
        // left: clear glass, center: dispersive glass, right: frosted glass
        worldList->add(std::make_shared<Sphere>(Point3(-1.2, 0, -1), 0.5, matGlass));
        worldList->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), 0.5, matGlassDisp));
        worldList->add(std::make_shared<Sphere>(Point3(1.2, 0, -1), 0.5, matFrosted));
        worldList->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));
    }
    else {
        // furnace: no floor, so the only light path is environment <-> object.
        // The two metals use a near-perfect conductor (k = 100, F ~ 0.9996 at every angle), so
        // only the microfacet model's own single-scattering energy loss is left.
        auto perfectMetal = std::make_shared<rayt::color::ComplexIorTable>();
        perfectMetal->samples() = { { 360.0, 1.0, 100.0 }, { 830.0, 1.0, 100.0 } };
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
            { "metal rough 0.2 ", matMetalMid.get(),    0.994, 0.030 },
            { "metal rough 0.5 ", matMetalRough.get(),  0.882, 0.030 },
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

    // IMPORTANT: シーンをファイナライズ（LightSamplerを構築）
    scene->finalize("uniform");
    std::cout << "[Scene] Finalized with "
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

    bool furnacePassed = true;
    if (opt.scene == "furnace") {
        furnacePassed = CheckFurnace(*scene, *camera, film, furnaceObjects);
    }

    std::cout << "[System] Finished." << std::endl;

    // Furnace test: exit code 2 when an object is outside its expected value, so a script can check it.
    return furnacePassed ? 0 : 2;
}
