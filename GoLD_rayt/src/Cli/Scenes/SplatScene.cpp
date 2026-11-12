// --scene splats: a synthetic captured room made of Gaussian splats (no data needed) around a gold sphere on a floor; runs the splat checks first.
#include "pch.h"
#include "Cli/Scenes.hpp"
#include "Cli/Checks.hpp"
#include "Cli/SceneParts.hpp"
#include "Geometry/GaussianSplats.hpp"
#include "Geometry/Sphere.hpp"

namespace rayt::cli {

std::optional<int> BuildSplats(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    out.checksPassed = RunSplatChecks(opt);

    // A dome of splats (radius 3 m) with a color that changes with the direction and a few bright "windows": the stand-in for a
    // captured room. Degree 0 harmonics: the color is f_dc alone.
    constexpr Real C0 = 0.28209479177387814;
    std::vector<SplatSource> src;
    std::uint32_t state = 12345;
    auto rnd = [&]() { state = state * 1664525u + 1013904223u; return Real(state >> 8) / Real(1 << 24); };
    for (int i = 0; i < 40000; ++i) {
        const Real u = rnd(), v = rnd();
        const Real phi = 2 * constants::PI * u;
        const Real cosT = 0.05 + 0.95 * v;                          // the upper part of the sphere (cos theta from 0.05 to 1)
        const Real sinT = std::sqrt(std::max(Real(0), 1 - cosT * cosT));
        const Vector3 dir(sinT * std::cos(phi), cosT, sinT * std::sin(phi));
        SplatSource s;
        s.mean = Point3(dir * Real(3.0 + 0.1 * (rnd() - 0.5)));
        s.scale = Vector3(0.06 + 0.06 * rnd(), 0.06 + 0.06 * rnd(), 0.02);   // flat, lying on the dome
        // orient the flat axis (z of the particle) along the dome's normal: a quaternion that turns z to dir
        const Vector3 z(0, 0, 1);
        const Vector3 axis = glm::cross(z, dir);
        const Real w = 1 + glm::dot(z, dir);
        s.rot = glm::normalize(glm::dvec4(w, axis.x, axis.y, axis.z));
        s.opacity = 0.85 + 0.1 * rnd();
        Vector3 rgb(0.35 + 0.25 * cosT, 0.30 + 0.15 * (0.5 + 0.5 * std::sin(phi)), 0.45 - 0.1 * cosT);   // a dim bluish room, warmer at the top
        if (std::abs(phi - 0.9) < 0.25 && cosT > 0.25 && cosT < 0.55) rgb = Vector3(0.98, 0.95, 0.85);   // a window
        if (std::abs(phi - 3.6) < 0.18 && cosT > 0.35 && cosT < 0.6) rgb = Vector3(0.95, 0.80, 0.55);    // a lamp
        s.sh = { float((rgb.x - 0.5) / C0), float((rgb.y - 0.5) / C0), float((rgb.z - 0.5) / C0) };
        src.push_back(std::move(s));
    }
    SplatParams params;
    params.linearize = true;
    params.intensity = 3.0;
    out.splats = std::make_shared<const GaussianSplats>(src, 0, params);

    // the physical objects: a gold sphere that mirrors the dome, a rough glass sphere, and a floor
    auto& world = out.world;
    world->add(std::make_shared<Sphere>(Point3(-0.55, 0.5, -0.6), 0.5, ctx.mats.goldMedium));
    world->add(std::make_shared<Sphere>(Point3(0.75, 0.4, -0.2), 0.4, std::make_shared<Dielectric>(1.5, 0.15)));
    world->add(std::make_shared<Sphere>(Point3(0, -100, 0), 100.0, std::make_shared<Lambertian>(RgbColor(0.5, 0.5, 0.5))));
    // a small lamp, so the floor also has a sharp shadow and the splats are not the only light
    out.lights.push_back(AddQuadLight(*world, nullptr, Point3(-0.6, 2.4, 0.4), Vector3(0.8, 0, 0), Vector3(0, 0, 0.8), EmissionSpectrum::IlluminantA(), 6.0, true));
    out.camera = CameraSpec{ Point3(0.1, 0.9, 2.2), Point3(0.0, 0.45, -0.4), Vector3(0, 1, 0), 42.0 };   // well inside the dome (radius 3)
    return std::nullopt;
}

} // namespace rayt::cli
