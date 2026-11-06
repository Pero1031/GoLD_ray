// Scenes about metals: gold spheres of several roughnesses, the color of the light, a gold bowl, the metal probe.
#include "pch.h"
#include "Cli/Scenes.hpp"
#include "Cli/Checks.hpp"
#include "Cli/LightPresets.hpp"
#include "Cli/PmxSupport.hpp"
#include "Cli/SceneParts.hpp"
#include "Cli/Tools.hpp"
#include "Geometry/ProceduralMesh.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

// --scene gold
std::optional<int> BuildGoldSpheres(SceneContext& ctx, SceneSetup& out) {
    auto& worldList = out.world;
    const auto& matFloor = ctx.mats.floor;
    const auto& matGoldSmooth = ctx.mats.goldSmooth;
    const auto& matGoldMedium = ctx.mats.goldMedium;
    const auto& matGoldRough = ctx.mats.goldRough;

    // リストに対して add する
    worldList->add(std::make_shared<Sphere>(Point3(-1.2, 0, -1), 0.5, matGoldSmooth));
    worldList->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), 0.5, matGoldMedium));
    worldList->add(std::make_shared<Sphere>(Point3(1.2, 0, -1), 0.5, matGoldRough));

    // 大きな球で床を作製
    worldList->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));
    return std::nullopt;
}

// --scene goldrough
std::optional<int> BuildGoldRough(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    const auto& matFloor = ctx.mats.floor;
    const auto& auTable = ctx.mats.auTable;

    // Rough gold spheres (roughness 0.4 / 0.7 / 1.0): multiple scattering matters here (--no-ms turns it off).
    for (const auto& [x, roughness] : { std::pair<Real, Real>{ -1.2, 0.4 }, std::pair<Real, Real>{ 0.0, 0.7 }, std::pair<Real, Real>{ 1.2, 1.0 } }) {
        auto m = std::make_shared<RoughConductor>(auTable, roughness);
        m->setMultipleScattering(!opt.noMs);
        worldList->add(std::make_shared<Sphere>(Point3(x, 0, -1), 0.5, m));
    }
    worldList->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));
    return std::nullopt;
}

// --scene spd
std::optional<int> BuildLightSpectra(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    const auto& matGoldMedium = ctx.mats.goldMedium;
    EmissionSpectrum spdEmission;                                 // spectrum of the sphere light
    std::shared_ptr<Light> spdLight;                              // registered for next-event estimation at the end

    // Five diffuse spheres with distinct reflectance spectra + gold, lit by one off-screen sphere light.
    const Real r = 0.3;
    auto lambertian = [](Real R, Real G, Real B) { return std::make_shared<Lambertian>(RgbColor(R, G, B)); };
    worldList->add(std::make_shared<Sphere>(Point3(-1.4, 0, -1), r, lambertian(0.80, 0.15, 0.12)));   // red
    worldList->add(std::make_shared<Sphere>(Point3(-0.7, 0, -1), r, lambertian(0.80, 0.80, 0.80)));   // white
    worldList->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), r, matGoldMedium));                    // gold
    worldList->add(std::make_shared<Sphere>(Point3(0.7, 0, -1), r, lambertian(0.12, 0.30, 0.80)));    // blue
    worldList->add(std::make_shared<Sphere>(Point3(1.4, 0, -1), r, lambertian(0.15, 0.60, 0.20)));    // green
    worldList->add(std::make_shared<Sphere>(Point3(0, -100.3, -1), 100.0, lambertian(0.7, 0.7, 0.7)));

    if (!MakeLightEmission(opt.light, spdEmission)) return 1;
    std::cout << "[SPD] light = " << opt.light << ", sensor RGB of the light = ("
        << spdEmission.rgb().x << ", " << spdEmission.rgb().y << ", " << spdEmission.rgb().z << ")"
        << (opt.whiteBalance ? ", film white-balanced to it" : "") << "\n";
    if (opt.whiteBalance) color::setWhiteBalance(spdEmission);

    // Radiance 20 at distance ~3.5 gives a floor radiance of about 0.6.
    auto lightMat = std::make_shared<DiffuseLight>(spdEmission, 20.0);
    auto lightShape = std::make_shared<Sphere>(Point3(-2.2, 3.0, 1.5), 0.8, lightMat);
    spdLight = std::make_shared<AreaLight>(lightShape, lightMat);
    worldList->add(std::make_shared<Primitive>(lightShape, lightMat, spdLight));
    if (spdLight) out.lights.push_back(spdLight);
    return std::nullopt;
}

// --scene goldbowl
std::optional<int> BuildGoldBowl(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    const auto& auTable = ctx.mats.auTable;
    EmissionSpectrum spdEmission;                                 // spectrum of the sphere light
    std::shared_ptr<Light> spdLight;                              // registered for next-event estimation at the end

    // Gold where RGB and spectral rendering disagree most: a gold bowl (light bounces between its walls: the error grows with every
    // reflection) and a gold ball inside it, under one sphere light whose spectrum is chosen with --light (d65 | tungsten | candle | led).
    // Metres. A white ball and a gray floor show what the same light does to a neutral surface. Render it in both builds and compare.
    const Real R = Real(0.08), wall = Real(0.003), tEnd = glm::radians(Real(78));
    std::vector<Point2> profile;
    const int steps = 40;
    for (int i = 0; i <= steps; ++i) {   // outer surface from the bottom up to the rim ...
        const Real t = tEnd * Real(i) / Real(steps);
        profile.push_back(Point2(R * std::sin(t), R - R * std::cos(t)));
    }
    for (int i = steps; i >= 0; --i) {   // ... and the inner surface back down to the axis
        const Real t = tEnd * Real(i) / Real(steps);
        profile.push_back(Point2((R - wall) * std::sin(t), R - (R - wall) * std::cos(t)));
    }
    auto bowlMat = std::make_shared<RoughConductor>(auTable, 0.12);
    worldList->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(MakeLatheMesh(profile, 96)), bowlMat, opt.accel));
    worldList->add(std::make_shared<Sphere>(Point3(0.0, wall + 0.0285, 0.0), Real(0.028), bowlMat));          // gold ball in the bowl
    worldList->add(std::make_shared<Sphere>(Point3(0.13, 0.025, 0.03), Real(0.025), std::make_shared<Lambertian>(RgbColor(0.8))));   // white reference
    worldList->add(std::make_shared<Sphere>(Point3(0.0, -100.0, 0.0), 100.0, std::make_shared<Lambertian>(RgbColor(0.5))));          // gray floor

    if (!MakeLightEmission(opt.light, spdEmission)) return 1;
    std::cout << "[GoldBowl] light = " << opt.light << ", sensor RGB of the light = (" << spdEmission.rgb().x << ", " << spdEmission.rgb().y << ", "
        << spdEmission.rgb().z << ")" << (opt.whiteBalance ? ", film white-balanced to it" : "") << "\n";
    if (opt.whiteBalance) color::setWhiteBalance(spdEmission);
    // a sphere light of radius 3 cm, 29 cm from the bowl: about 3 units of irradiance
    auto lightMat = std::make_shared<DiffuseLight>(spdEmission, 90.0);
    auto lightShape = std::make_shared<Sphere>(Point3(-0.18, 0.20, 0.10), Real(0.03), lightMat);
    spdLight = std::make_shared<AreaLight>(lightShape, lightMat);
    worldList->add(std::make_shared<Primitive>(lightShape, lightMat, spdLight));
    // a dim, uniform glow of the same spectrum all around (not a registered light: found by BSDF sampling): without it the inside of the
    // bowl reflects black and the reflections between the gold surfaces cannot be seen
    auto shellMat = std::make_shared<DiffuseLight>(spdEmission, 0.25, true);
    worldList->add(std::make_shared<Primitive>(std::make_shared<Sphere>(Point3(0, 0, 0), 20.0, shellMat), shellMat));
    if (spdLight) out.lights.push_back(spdLight);
    {
        out.camera.from = Point3(0.02, 0.27, 0.33);
        out.camera.at = Point3(0.035, 0.035, 0.0);
        out.camera.fovDeg = 30.0;
    }
    return std::nullopt;
}

// --scene metalprobe
std::optional<int> BuildMetalProbe(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    EmissionSpectrum spdEmission;
    std::shared_ptr<Material> probeMetal;                         // the metal sphere's material

    // A smooth metal sphere (radius 1) seen from far away (nearly orthographic) inside a large
    // two-sided emissive shell: pixel ring radius rho = sin(theta of incidence), one specular bounce.
    // The shell is not a registered light (no NEE needed: the metal is almost a mirror).
    if (!MakeLightEmission(opt.light, spdEmission)) return 1;
    auto table = std::make_shared<rayt::color::ComplexIorTable>();
    std::string err;
    const std::string path = "assets/Metal_data/Johnson_" + opt.metal + ".csv";
    if (!rayt::io::loadComplexIor_RefractiveIndexInfoCsv(path, *table, &err)) {
        std::cerr << "[MetalProbe] cannot load " << path << ": " << err << "\n";
        return 1;
    }
    probeMetal = std::make_shared<RoughConductor>(table, opt.probeRoughness);
    std::cout << "[MetalProbe] metal = " << opt.metal << ", light = " << opt.light
        << ", sensor RGB of the light = (" << spdEmission.rgb().x << ", " << spdEmission.rgb().y
        << ", " << spdEmission.rgb().z << ")\n";

    worldList->add(std::make_shared<Sphere>(Point3(0, 0, 0), 1.0, probeMetal));
    auto shellMat = std::make_shared<DiffuseLight>(spdEmission, 1.0, true);
    auto shell = std::make_shared<Sphere>(Point3(0, 0, 0), 300.0, shellMat);
    worldList->add(std::make_shared<Primitive>(shell, shellMat));

    out.lightSampler = "uniform";   // there is no light to sample
    {
        // Distance 100, unit sphere filling 96 % of the image height: the view direction differs from
        // the optical axis by < 0.6 degrees, so the incidence angle is read off the hit normal.
        out.camera.from = Point3(0, 0, 100);
        out.camera.at = Point3(0, 0, 0);
        out.camera.fovDeg = Real(2) * std::atan(Real(1) / Real(0.96) / Real(100)) * Real(180) / constants::PI;
    }
    out.afterSave = [probeMetal, label = opt.metal + " under " + opt.light + " r=" + std::to_string(opt.probeRoughness).substr(0, 4) +
                                         std::string(" [") + kRenderModeName + "]"](const RenderResult& r) {
        ReportMetalProbe(r.scene, r.camera, r.film, probeMetal.get(), label);
        return true;
    };
    return std::nullopt;
}

} // namespace rayt::cli
