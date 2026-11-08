// The BSSRDF (Burley profile + probe rays) against a volumetric random walk on a sphere and on thin discs (--scene sssrw).
//
// Both objects are index matched: the volume is bounded by a null interface (ior 1, no reflection, no refraction) and the BSSRDF material has no coat
// (weight 0, so no Fresnel factors at the entry or the exit). What is left of the difference is what the profile model cannot do: the curvature of a
// surface, the thickness of a body (light that leaves through the other side), the shape of the profile itself.
//
// The BSSRDF parameters come from the medium: the scatter distance l is the reduced mean free path 1 / (sigma_t (1 - g alpha)) (the engine's convention),
// the diffuse color A is the total diffuse reflectance of the semi-infinite medium (docs/SSS.md, 9.1: --sssrw-A, or Monte Carlo of the engine's layered-skin
// tool, which uses the similarity relation). A point light with a narrow cone lights the object; both versions are path traced.
#include "pch.h"
#include "Cli/Tools.hpp"
#include "Cli/EngineHeaders.hpp"
#include "Cli/SceneParts.hpp"
#include "Geometry/ProceduralMesh.hpp"

namespace rayt::cli {

namespace {

struct Pair { double volume = 0, bssrdf = 0; };

// One rendering: the object (a hittable) lit by a spot, seen by a pinhole camera; returns the film.
Film Render(const std::shared_ptr<Hittable>& object, const Point3& lampPos, const Vector3& lampAxis, Real lampHalfAngle, Real lampIntensity,
    const Point3& from, const Point3& at, const Vector3& up, Real fov, int W, int H, int spp, int depth, int volumeCap)
{
    auto world = std::make_shared<Aggregate>();
    world->add(object);
    auto lamp = std::make_shared<PointLight>(lampPos, RgbColor(1.0, 1.0, 1.0), lampIntensity);
    lamp->setCone(lampAxis, lampHalfAngle);
    auto scene = std::make_shared<Scene>(world);
    scene->addLight(lamp);
    scene->finalize("uniform");
    auto cam = std::make_shared<Camera>(from, at, up, fov, double(W) / H, 0.0, glm::length(from - at));
    PathIntegrator pt(cam, depth, spp);
    pt.setMediaNee(true);
    pt.setMaxVolumeEvents(volumeCap);
    Film film(W, H);
    pt.render(*scene, film);
    return film;
}

} // namespace

int RunSssRandomWalk(const CliOptions& opt) {
    const double ell = opt.sssrwEll;                               // the reduced mean free path (m)
    const double alpha = opt.sssrwAlbedo, g = opt.sssrwG;
    const double sigmaT = 1.0 / (ell * (1.0 - g * alpha));         // extinction (1/m)
    double A = opt.sssrwA;
    if (A <= 0) {                                                   // the engine's own Monte Carlo of a semi-infinite medium (similarity relation)
        const double mua = sigmaT * (1 - alpha), musp = sigmaT * alpha * (1 - g);
        const auto r = skin::simulateRadialProfile({ skin::Layer{ 0, mua, musp } }, 64, 1e-6, 0.05, 400000, 3);
        A = r.reflected;
    }
    const int spp = std::max(16, opt.spp), depth = 64, cap = opt.volumeCap;
    std::cout << "[SssRw] medium: albedo " << alpha << ", g " << g << ", reduced mean free path " << ell * 1000 << " mm (sigma_t " << sigmaT
        << " /m); BSSRDF: A " << A << ", scatter distance " << ell * 1000 << " mm; " << spp << " spp\n";

    auto medium = std::make_shared<HomogeneousMedium>(Real(sigmaT), Real(alpha), PhaseFunction(Real(g)));
    auto boundary = std::make_shared<Dielectric>(Real(1.0), Real(0.0));   // a null interface: the medium is index matched
    boundary->setMedium(medium);
    auto profile = std::make_shared<SubsurfaceProfile>(Vector3(ell), Real(1.4));
    auto skin = std::make_shared<CoatedDiffuse>(RgbColor(A), Real(0.5), Real(0.0), Real(0.0));   // no coat: w = 0
    skin->setSubsurface(profile);
    std::shared_ptr<Material> boundaryM = boundary, skinM = skin;

    // ------------------------------------------------------------------------------------------------
    // 1. A sphere (radius 5 l by default), lit from the side at 60 degrees: the radiance along the horizontal diameter
    // ------------------------------------------------------------------------------------------------
    {
        const double R = opt.sssrwSphere * ell;
        const int N = 96;
        const Point3 lamp(0.5 * std::sin(glm::radians(60.0)), 0.0, 0.5 * std::cos(glm::radians(60.0)));
        const Real I = Real(3.14159265 * 0.25);   // an irradiance of pi: a white Lambertian surface facing the lamp has the radiance 1
        const Real fov = Real(2.0 * std::atan(1.25 * R / 0.5) * 180.0 / 3.14159265);
        const Real half = Real(std::atan(1.2 * R / 0.5));
        Film fV = Render(std::make_shared<Sphere>(Point3(0.0), Real(R), boundaryM), lamp, -Vector3(lamp), half, I,
            Point3(0, 0, 0.5), Point3(0.0), Vector3(0, 1, 0), fov, N, N, spp, depth, cap);
        Film fB = Render(std::make_shared<Sphere>(Point3(0.0), Real(R), skinM), lamp, -Vector3(lamp), half, I,
            Point3(0, 0, 0.5), Point3(0.0), Vector3(0, 1, 0), fov, N, N, spp, depth, cap);
        const double pixelPerR = (N * 0.5) / 1.25;          // pixels per sphere radius
        double sumV = 0, sumB = 0;
        for (int i = 0; i < N * N; ++i) { sumV += fV.getData()[i].y; sumB += fB.getData()[i].y; }
        std::cout << "[SssRw] sphere, radius " << opt.sssrwSphere << " l = " << R * 1000 << " mm, lit from the side (60 degrees from the viewing direction)\n";
        std::cout << "[SssRw]   total radiance (sum over the picture), volume " << sumV << ", BSSRDF " << sumB << ", ratio BSSRDF / volume " << sumB / sumV << "\n";
        std::cout << "[SssRw]   radiance along the horizontal diameter (x in sphere radii; the lit side is on the right)\n";
        std::cout << "[SssRw]      x/R    volume   BSSRDF   ratio\n";
        for (int k = -8; k <= 8; ++k) {
            const double x = k / 8.0 * 0.96;
            const int px = static_cast<int>(N * 0.5 + x * pixelPerR);
            double v = 0, b = 0; int n = 0;
            for (int dy = -2; dy <= 2; ++dy) for (int dx = -1; dx <= 1; ++dx) {
                const int xx = px + dx, yy = N / 2 + dy;
                if (xx < 0 || xx >= N) continue;
                v += fV.getData()[yy * N + xx].y; b += fB.getData()[yy * N + xx].y; ++n;
            }
            v /= n; b /= n;
            std::printf("[SssRw]    %+5.2f   %6.3f   %6.3f   %5.2f\n", x, v, b, v > 1e-4 ? b / v : 0.0);
        }
    }

    // ------------------------------------------------------------------------------------------------
    // 2. Thin discs of thickness t (radius 15 l), a spot on the axis: reflection (lit and seen from above), transmission (lit from below)
    // ------------------------------------------------------------------------------------------------
    std::cout << "[SssRw] discs of thickness t (radius 15 l); the spot covers the center; the lamp " << opt.sssrwLampDeg << " degrees from the axis; the mean radiance of the central region of the picture\n";
    std::cout << "[SssRw]   t / l   reflection: volume  BSSRDF  ratio    transmission: volume  BSSRDF  ratio\n";
    for (double tl : { 0.25, 0.5, 1.0, 2.0, 4.0, 8.0 }) {
        const double t = tl * ell, Rd = 15.0 * ell;
        TriangleMesh mesh = MakeLatheMesh({ Point2(0, -t / 2), Point2(Rd, -t / 2), Point2(Rd, t / 2), Point2(0, t / 2) }, 96);
        auto meshP = std::make_shared<const TriangleMesh>(std::move(mesh));
        const int N = 64;
        const double seeHalf = 4.0 * ell;                                      // the picture covers +-4 l around the axis
        const Real fov = Real(2.0 * std::atan(seeHalf / 0.4) * 180.0 / 3.14159265);
        const Real half = Real(std::atan(8.0 * ell / 0.4));                    // the spot covers +-8 l
        const Real I = Real(3.14159265 * 0.16);                                // an irradiance of pi at 0.4 m
        Pair refl, trans;
        for (int mode = 0; mode < 2; ++mode) {
            // mode 0: above (reflection), mode 1: below (transmission); tilted by --sssrw-lamp-deg from the axis, away from the line of the camera
            const double th = glm::radians(opt.sssrwLampDeg);
            const Point3 lamp(0.4 * std::sin(th), (mode == 0 ? 0.4 : -0.4) * std::cos(th), 0.0);
            Film fV = Render(MakeMeshHittable(meshP, boundaryM, opt.accel), lamp, -Vector3(lamp), half, I,
                Point3(0, 0.4, 1e-4), Point3(0.0), Vector3(0, 0, -1), fov, N, N, spp, depth, cap);
            Film fB = Render(MakeMeshHittable(meshP, skinM, opt.accel), lamp, -Vector3(lamp), half, I,
                Point3(0, 0.4, 1e-4), Point3(0.0), Vector3(0, 0, -1), fov, N, N, spp, depth, cap);
            double v = 0, b = 0; int n = 0;
            for (int y = N / 2 - 8; y < N / 2 + 8; ++y) for (int x = N / 2 - 8; x < N / 2 + 8; ++x) { v += fV.getData()[y * N + x].y; b += fB.getData()[y * N + x].y; ++n; }
            (mode == 0 ? refl : trans) = Pair{ v / n, b / n };
        }
        std::printf("[SssRw]   %5.2f          %7.4f  %7.4f  %5.2f              %7.4f  %7.4f  %5.2f\n", tl, refl.volume, refl.bssrdf,
            refl.volume > 1e-6 ? refl.bssrdf / refl.volume : 0.0, trans.volume, trans.bssrdf, trans.volume > 1e-6 ? trans.bssrdf / trans.volume : 0.0);
    }
    return 0;
}

} // namespace rayt::cli
