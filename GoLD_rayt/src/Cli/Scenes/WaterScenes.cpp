// Water scenes: a tumbler of water, a sea, a swimming pool.
#include "pch.h"
#include "Cli/Scenes.hpp"
#include "Cli/Checks.hpp"
#include "Cli/LightPresets.hpp"
#include "Cli/PmxSupport.hpp"
#include "Cli/SceneParts.hpp"
#include "Cli/Tools.hpp"
#include "Geometry/ProceduralMesh.hpp"
#include "Geometry/Quad.hpp"
#include "Geometry/Sphere.hpp"
#include "Lights/AreaLight.hpp"
#include "Lights/PointLight.hpp"
#include "Materials/DiffuseLight.hpp"
#include "Scene/Primitive.hpp"

namespace rayt::cli {

// --scene water
std::optional<int> BuildWater(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    Vector3& envPeakDir = ctx.env.peakDir;
    auto& studioLights = out.lights;
    Point3& causticCenter = out.causticTarget.center;
    Real& causticTargetRadius = out.causticTarget.radius;
    out.caustics.photonMap = true;   // only with --caustics N: there is no default number of photons

    // Everything is in metres: a tumbler (10 cm tall, 6-7.6 cm wide) filled to 8 cm, standing 0.5 mm above a table, in front
    // of a wall of vertical stripes (their mirror-reversed, squeezed image is the classic sign of a water-filled cylinder).
    // The water is a separate closed mesh 0.2 mm inside the glass, so every interface is air | dielectric and no nested
    // dielectric bookkeeping is needed; the thin air gap is invisible at this scale.
    const Real y0 = Real(0.0005);                       // glass base above the table
    const Real gH = 0.100, gBase = 0.008, fill = 0.080, gap = 0.0002;
    const Real rOutBot = 0.030, rOutTop = 0.038, rInBot = 0.0275, rInTop = 0.035;
    auto rInAt = [&](Real y) { return rInBot + (rInTop - rInBot) * (y - gBase) / (gH - gBase); };
    const int rings = 192;

    // Glass: outer wall, flat rim, inner wall, floor of the cavity (counter-clockwise, solid on the left).
    TriangleMesh glassMesh = MakeLatheMesh({ Point2(0, 0), Point2(rOutBot, 0), Point2(rOutTop, gH), Point2(rInTop, gH),
        Point2(rInBot, gBase), Point2(0, gBase) }, rings);
    // Water: bottom, side, flat free surface.
    const Real yb = gBase + gap;
    TriangleMesh waterMesh = MakeLatheMesh({ Point2(0, yb), Point2(rInAt(yb) - gap, yb), Point2(rInAt(fill) - gap, fill), Point2(0, fill) }, rings);
    glassMesh.transform(glm::dmat3(1.0), Vector3(0, y0, 0));
    waterMesh.transform(glm::dmat3(1.0), Vector3(0, y0, 0));
    std::cout << "[Water] glass " << glassMesh.numTriangles() << " triangles, volume " << glassMesh.signedVolume() * 1e6
        << " cm^3; water " << waterMesh.numTriangles() << " triangles, volume " << waterMesh.signedVolume() * 1e6
        << " cm^3 (both must be > 0: closed, outward)\n";

    auto glass = std::make_shared<Dielectric>(Real(1.5), Real(0.0));
    auto water = std::make_shared<Dielectric>(Real(1.333), Real(0.0));
    const Real waterAbsorb = opt.absorbScale >= 0 ? opt.absorbScale : Real(40);
    if (waterAbsorb > 0) {
        color::SpectralCurve sigmaA;
        if (!LoadWaterAbsorption(sigmaA)) return 1;
        auto absorbing = std::make_shared<AbsorbingMedium>(sigmaA, waterAbsorb);
        water->setMedium(absorbing);
        const Real d = Real(2) * rInAt(0.04);   // light path across the glass at half height
        std::cout << "[Water] absorption x" << waterAbsorb << " (1 = physical). Transmittance over the " << d * 100 << " cm path across the glass: "
            << "450 nm " << std::exp(-absorbing->sigmaA(450) * d) << ", 550 nm " << std::exp(-absorbing->sigmaA(550) * d)
            << ", 650 nm " << std::exp(-absorbing->sigmaA(650) * d) << "   (RGB-mode coefficients, 1/m: "
            << absorbing->sigmaRgb().x << ", " << absorbing->sigmaRgb().y << ", " << absorbing->sigmaRgb().z << ")\n";
    }
    else std::cout << "[Water] no absorption (--absorb-scale 0): clear refractive water\n";
    worldList->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(glassMesh)), glass, opt.accel));
    worldList->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(waterMesh)), water, opt.accel));

    // Table and the striped wall behind (white / dark stripes, 12 mm wide).
    worldList->add(std::make_shared<Quad>(Point3(-2, 0, -2), Vector3(0, 0, 4), Vector3(4, 0, 0),
        std::make_shared<Lambertian>(RgbColor(0.55))));                                         // normal +y
    if (!opt.waterSky && !opt.waterLamp) {   // under the open sky there is no wall: it would shade the table from half of the environment
        auto stripeWhite = std::make_shared<Lambertian>(RgbColor(0.80));
        auto stripeDark = std::make_shared<Lambertian>(RgbColor(0.04));
        const Real sw = 0.012;
        for (int k = -30; k < 30; ++k)
            worldList->add(std::make_shared<Quad>(Point3(Real(k) * sw, 0, -0.14), Vector3(sw, 0, 0), Vector3(0, 0.35, 0),
                (k % 2 == 0) ? stripeWhite : stripeDark));
    }

    const Point3 target(0.0, 0.05, 0.0);
    auto softbox = [&](const Point3& center, Real w, Real h, const EmissionSpectrum& e, Real radiance) {
        const Vector3 n = glm::normalize(target - center);
        const Vector3 e1 = glm::normalize(glm::cross(Vector3(0, 1, 0), n)) * w;
        const Vector3 e2 = glm::normalize(glm::cross(n, e1)) * h;
        studioLights.push_back(AddQuadLight(*worldList, nullptr, center - e1 * 0.5 - e2 * 0.5, e1, e2, e, radiance, true));
    };
    if (opt.waterLamp) {
        const Point3 lamp(-0.20, 0.22, 0.16);
        studioLights.push_back(std::make_shared<PointLight>(lamp, RgbColor(1.0, 0.97, 0.92), Real(0.8)));
        envPeakDir = glm::normalize(Vector3(lamp - target));   // the camera azimuth follows the light, like --water-sky
        if (opt.causticPhotons <= 0) std::cout << "[Water] note: a point light casts no caustic without --caustics (it cannot be hit by BSDF sampling)\n";
    }
    else if (!opt.waterSky) {
        softbox(Point3(-0.30, 0.28, 0.22), 0.22, 0.22, EmissionSpectrum(RgbColor(1.0, 0.97, 0.92)), 9.0);   // key: above-left-front
        softbox(Point3(0.35, 0.12, 0.30), 0.25, 0.20, EmissionSpectrum(RgbColor(1.0)), 1.2);               // fill: front-right, dim
        softbox(Point3(0.30, 0.09, 0.00), 0.04, 0.18, EmissionSpectrum(RgbColor(1.0)), 20.0);              // strip on the right: rim light on the glass edge
    }
    causticCenter = Point3(0.0, y0 + gH * 0.5, 0.0);                                  // the tumbler is the only specular object
    causticTargetRadius = std::sqrt(rOutTop * rOutTop + gH * gH * 0.25) * Real(1.02);
    {
        out.camera.from = Point3(0.0, 0.075, 0.40);
        out.camera.at = Point3(0.0, 0.048, 0.0);
        out.camera.fovDeg = 22.0;
        if (opt.waterSky || opt.waterLamp) {   // looking down at the table beside / behind the tumbler, where the caustic lands
            const Real elev = glm::radians(opt.camElev >= 0 ? opt.camElev : Real(38));
            const Real yaw = (opt.camYaw < 1e8) ? glm::radians(opt.camYaw)
                : std::atan2(envPeakDir.x, envPeakDir.z) + glm::radians(Real(35));   // on the sun's side, turned 35 degrees: the caustic falls behind
            out.camera.at = Point3(0.0, 0.03, 0.0);
            out.camera.from = out.camera.at + Vector3(std::cos(elev) * std::sin(yaw), std::sin(elev), std::cos(elev) * std::cos(yaw)) * Real(0.38);
            out.camera.fovDeg = 30.0;
            std::cout << "[Water] camera: elevation " << glm::degrees(elev) << " deg, azimuth " << glm::degrees(yaw) << " deg\n";
        }
    }
    return std::nullopt;
}

// --scene sea
std::optional<int> BuildSea(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    Vector3& envPeakDir = ctx.env.peakDir;
    Real seaViewAz = 0;                                           // azimuth of the view direction (about +y, 0 = +z)

    // Metres. The water is ONE open sheet: a ray that refracts into it carries the absorbing medium (Ray::medium) until it
    // crosses a water surface again, so the seabed under it needs no closed boundary. The camera looks along the map's sun.
    const Real sunAz = std::atan2(envPeakDir.x, envPeakDir.z);
    seaViewAz = (opt.camYaw < 1e8) ? glm::radians(opt.camYaw) : sunAz - glm::radians(Real(15));   // sun slightly right of center
    SeaParams P;
    P.steep = opt.seaSteep;
    P.seed = opt.seaSeed;
    P.windDeg = glm::degrees(std::atan2(-std::cos(seaViewAz), -std::sin(seaViewAz)));   // waves run toward the camera
    TriangleMesh seaMesh = MakeSeaMesh(P);
    std::cout << "[Sea] surface " << seaMesh.numTriangles() << " triangles, camera azimuth " << glm::degrees(seaViewAz) << " deg, wave steepness "
        << opt.seaSteep << ", depth " << opt.seaDepth << " m\n";

    auto water = std::make_shared<Dielectric>(Real(1.333), Real(0.0));
    const Real seaAbsorb = opt.absorbScale >= 0 ? opt.absorbScale : Real(1);
    if (seaAbsorb > 0) {
        color::SpectralCurve sigmaA;
        if (!LoadWaterAbsorption(sigmaA)) return 1;
        auto absorbing = std::make_shared<AbsorbingMedium>(sigmaA, seaAbsorb);
        water->setMedium(absorbing);
        const Real d = Real(2) * opt.seaDepth;   // down and up again
        std::cout << "[Sea] absorption x" << seaAbsorb << ": transmittance over " << d << " m  450 nm " << std::exp(-absorbing->sigmaA(450) * d)
            << ", 550 nm " << std::exp(-absorbing->sigmaA(550) * d) << ", 650 nm " << std::exp(-absorbing->sigmaA(650) * d) << "\n";
    }
    worldList->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(seaMesh)), water, opt.accel));

    // Seabed: sand, with dark rocks (half-buried spheres) that show the refraction and the absorption with depth.
    worldList->add(std::make_shared<Quad>(Point3(-3000, -opt.seaDepth, -3000), Vector3(0, 0, 6000), Vector3(6000, 0, 0),
        std::make_shared<Lambertian>(RgbColor(0.65, 0.58, 0.42))));   // normal +y
    auto rockA = std::make_shared<Lambertian>(RgbColor(0.16, 0.15, 0.14));
    auto rockB = std::make_shared<Lambertian>(RgbColor(0.30, 0.27, 0.22));
    std::mt19937 rockRng(opt.seaSeed + 77);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (int i = 0; i < 120; ++i) {
        const Real dist = Real(2.5) + Real(70) * Real(U(rockRng)) * Real(U(rockRng));
        const Real ang = Real(2) * constants::PI * Real(U(rockRng));
        const Real rad = Real(0.15) + Real(0.5) * Real(U(rockRng)) * Real(U(rockRng));
        worldList->add(std::make_shared<Sphere>(Point3(dist * std::cos(ang), -opt.seaDepth + Real(0.25) * rad, dist * std::sin(ang)), rad,
            (i % 3 == 0) ? rockB : rockA));
    }
    {
        const Real pitch = glm::radians(Real(-2));
        out.camera.from = Point3(0.0, 2.0, 0.0);
        out.camera.at = out.camera.from + Vector3(std::sin(seaViewAz) * std::cos(pitch), std::sin(pitch), std::cos(seaViewAz) * std::cos(pitch)) * Real(10);
        out.camera.fovDeg = 40.0;
    }
    return std::nullopt;
}

// --scene pool
std::optional<int> BuildPool(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    const auto& envLight = ctx.env.light;
    auto& studioLights = out.lights;
    Point3& causticCenter = out.causticTarget.center;
    Real& causticTargetRadius = out.causticTarget.radius;
    out.caustics.photonMap = true;
    out.caustics.defaultPhotons = 8000000;
    out.caustics.defaultRadius = Real(0.015);   // a pool is metres wide: the network on its floor is coarse

    // Metres. Basin 3.2 x 5.2 m, water level y = 0, deck 12 cm above it, the floor at -depth. The camera stands on the deck at +z
    // looking toward -z. The water is one open rippled sheet carrying the absorbing medium; the sun is an area light that
    // takes its caustics from the photon map, the HDR sky (flagged off below) keeps BSDF sampling.
    const Real W = 3.2, D = 5.2, H = opt.poolDepth, deckH = 0.12, L = 300;
    TriangleMesh surface = MakePoolSurface(W, D, opt.poolSteep, opt.seaSeed, Real(0.015));
    std::cout << "[Pool] surface " << surface.numTriangles() << " triangles, depth " << H << " m, ripple steepness " << opt.poolSteep << "\n";
    auto water = std::make_shared<Dielectric>(Real(1.333), Real(0.0));
    const Real poolAbsorb = opt.absorbScale >= 0 ? opt.absorbScale : Real(1);
    if (poolAbsorb > 0) {
        color::SpectralCurve sigmaA;
        if (!LoadWaterAbsorption(sigmaA)) return 1;
        auto absorbing = std::make_shared<AbsorbingMedium>(sigmaA, poolAbsorb);
        water->setMedium(absorbing);
        const Real d = Real(2) * H;
        std::cout << "[Pool] absorption x" << poolAbsorb << ": transmittance over " << d << " m  450 nm " << std::exp(-absorbing->sigmaA(450) * d)
            << ", 550 nm " << std::exp(-absorbing->sigmaA(550) * d) << ", 650 nm " << std::exp(-absorbing->sigmaA(650) * d) << "\n";
    }
    worldList->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(surface)), water, opt.accel));

    // Basin: checkered floor (0.4 m tiles), plain walls, a deck all around.
    auto tileA = std::make_shared<Lambertian>(RgbColor(0.80, 0.86, 0.90));
    auto tileB = std::make_shared<Lambertian>(RgbColor(0.55, 0.68, 0.78));
    const Real tile = 0.4;
    const int nTx = static_cast<int>(std::lround(W / tile)), nTz = static_cast<int>(std::lround(D / tile));
    for (int i = 0; i < nTx; ++i)
        for (int j = 0; j < nTz; ++j)
            worldList->add(std::make_shared<Quad>(Point3(-W / 2 + Real(i) * (W / nTx), -H, -D / 2 + Real(j) * (D / nTz)),
                Vector3(0, 0, D / nTz), Vector3(W / nTx, 0, 0), ((i + j) % 2 == 0) ? tileA : tileB));   // normal +y
    auto wallMat = std::make_shared<Lambertian>(RgbColor(0.70, 0.82, 0.88));
    const Real wallH = H + deckH;
    worldList->add(std::make_shared<Quad>(Point3(-W / 2, -H, -D / 2), Vector3(0, 0, D), Vector3(0, wallH, 0), wallMat));
    worldList->add(std::make_shared<Quad>(Point3(W / 2, -H, -D / 2), Vector3(0, 0, D), Vector3(0, wallH, 0), wallMat));
    worldList->add(std::make_shared<Quad>(Point3(-W / 2, -H, -D / 2), Vector3(W, 0, 0), Vector3(0, wallH, 0), wallMat));
    worldList->add(std::make_shared<Quad>(Point3(-W / 2, -H, D / 2), Vector3(W, 0, 0), Vector3(0, wallH, 0), wallMat));
    auto deck = std::make_shared<Lambertian>(RgbColor(0.42, 0.40, 0.37));
    worldList->add(std::make_shared<Quad>(Point3(-L, deckH, -L), Vector3(0, 0, L - D / 2), Vector3(2 * L, 0, 0), deck));   // beyond -z
    worldList->add(std::make_shared<Quad>(Point3(-L, deckH, D / 2), Vector3(0, 0, L - D / 2), Vector3(2 * L, 0, 0), deck));   // beyond +z
    worldList->add(std::make_shared<Quad>(Point3(-L, deckH, -D / 2), Vector3(0, 0, D), Vector3(L - W / 2, 0, 0), deck));     // beyond -x
    worldList->add(std::make_shared<Quad>(Point3(W / 2, deckH, -D / 2), Vector3(0, 0, D), Vector3(L - W / 2, 0, 0), deck));  // beyond +x

    // Sun: a sphere light 300 m away, 1.2 degree radius (larger and dimmer than the real one: see the sea scene).
    const Real sunAz = glm::radians(opt.poolSunAz), sunEl = glm::radians(opt.poolSunElev);
    const Vector3 sunDir(std::cos(sunEl) * std::sin(sunAz), std::sin(sunEl), std::cos(sunEl) * std::cos(sunAz));
    const Real sunDist = 300, sunRadius = sunDist * std::tan(glm::radians(Real(1.2)));
    auto sunMat = std::make_shared<DiffuseLight>(EmissionSpectrum(RgbColor(1.0, 0.95, 0.85)), Real(4000));   // irradiance about 4.3 on a horizontal surface
    auto sunShape = std::make_shared<Sphere>(Point3(sunDir * sunDist), sunRadius, sunMat);
    auto sunLight = std::make_shared<AreaLight>(sunShape, sunMat);
    worldList->add(std::make_shared<Primitive>(sunShape, sunMat, sunLight));
    studioLights.push_back(sunLight);

    causticCenter = Point3(0.0, 0.0, 0.0);   // the water sheet is the only specular object
    causticTargetRadius = std::sqrt(W * W + D * D) * Real(0.5) * Real(1.02);
    {
        out.camera.from = Point3(0.5, 1.35, 3.1);
        out.camera.at = Point3(0.0, -0.1, -2.0);   // the horizon of the map stays in the upper part of the frame
        out.camera.fovDeg = 52.0;
    }
    // the sun takes its caustics from the photon map, the sky of the map does not (it is estimated better by BSDF sampling)
    if (envLight) {
        envLight->setPhotonCaustics(false);
        envLight->setScale(Real(0.5));   // the sun should dominate (sun : sky about 4 : 1), as in a real pool: that is what makes the caustic network
    }
    return std::nullopt;
}

} // namespace rayt::cli
