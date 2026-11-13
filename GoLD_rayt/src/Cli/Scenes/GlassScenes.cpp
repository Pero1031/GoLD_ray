// Scenes about glass and media: glass spheres, the white furnace, a glass blob, a prism, faceted crystal.
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
#include "IO/ObjLoader.hpp"
#include "Lights/PointLight.hpp"
#include "Materials/CoatedDiffuse.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Integrator.hpp"

namespace rayt::cli {

// --scene glass
std::optional<int> BuildGlassSpheres(SceneContext& ctx, SceneSetup& out) {
    auto& worldList = out.world;
    const auto& matFloor = ctx.mats.floor;
    const auto& matGlass = ctx.mats.glass;
    const auto& matGlassDisp = ctx.mats.glassDisp;
    const auto& matFrosted = ctx.mats.frosted;

    // left: clear glass, center: dispersive glass, right: frosted glass
    worldList->add(std::make_shared<Sphere>(Point3(-1.2, 0, -1), 0.5, matGlass));
    worldList->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), 0.5, matGlassDisp));
    worldList->add(std::make_shared<Sphere>(Point3(1.2, 0, -1), 0.5, matFrosted));
    worldList->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));
    return std::nullopt;
}

// --scene furnace
std::optional<int> BuildFurnace(SceneContext& ctx, SceneSetup& out) {
    auto& worldList = out.world;
    const auto& matGlass = ctx.mats.glass;
    const auto& matGlassDisp = ctx.mats.glassDisp;
    const auto& matFrosted = ctx.mats.frosted;
    const auto& matWhite = ctx.mats.white;
    std::vector<FurnaceObject> furnaceObjects;                    // what the furnace test checks
    std::vector<std::shared_ptr<Material>> furnaceMaterials;      // owns materials only used here

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
    // expected / systematic term (the check adds 4 standard errors from the pixels): 1.0 for lossless materials; the frosted glass is a
    // single-scattering rough dielectric (roughness 0.2) and keeps 0.970 (measured at 256 spp in both builds, 2026-11-13; no closed form).
    // Measured at 256 spp and left as documented losses inside the terms: white diffuse 0.9985, metal 0.2 0.9972, clear glass 0.9995.
    furnaceObjects = {
        { "background      ", nullptr,              1.000, 0.002 },
        { "clear glass     ", matGlass.get(),       1.000, 0.005 },
        { "dispersive glass", matGlassDisp.get(),   1.000, 0.005 },
        { "white diffuse   ", matWhite.get(),       1.000, 0.005 },
        { "frosted glass   ", matFrosted.get(),     0.970, 0.005 },
        { "metal rough 0.2 ", matMetalMid.get(),    1.000, 0.005 },
        { "metal rough 0.5 ", matMetalRough.get(),  1.000, 0.005 },
    };
    out.afterSave = [furnaceObjects, furnaceMaterials](const RenderResult& r) {
        return CheckFurnace(r.scene, r.camera, r.film, furnaceObjects);
    };
    return std::nullopt;
}

// --scene blob
std::optional<int> BuildBlob(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    std::shared_ptr<Material> blobMaterial;                       // the glass + medium material

    // A closed mesh of glass (outward normals) filled with a homogeneous scattering medium: the light enters through the
    // glass, random-walks inside and leaves through the glass again. Lit by the environment light only.
    TriangleMesh blob;
    if (!opt.mesh.empty()) {
        std::string objErr;
        if (!io::loadObj(opt.mesh, blob, &objErr)) { std::cerr << "[Blob] " << objErr << "\n"; return 1; }
    }
    else blob = TriangleMesh::Icosphere(4, 0.025, Point3(0.0), true);
    auto glass = std::make_shared<Dielectric>(Real(1.5), opt.glassRoughness);
    if (opt.sigmaTRgb || opt.albedoRgb) {   // a chromatic medium: the unset one is gray with the scalar option's value
        const Vector3 sigmaT = opt.sigmaTRgb ? *opt.sigmaTRgb : Vector3(opt.sigmaT);
        const Vector3 albedo = opt.albedoRgb ? *opt.albedoRgb : Vector3(opt.mediumAlbedo);
        glass->setMedium(std::make_shared<ChromaticMedium>(sigmaT, albedo, PhaseFunction(opt.phaseG, opt.phaseG2, opt.phaseW)));
        std::cout << "[Blob] chromatic medium: sigma_t (" << sigmaT.x << ", " << sigmaT.y << ", " << sigmaT.z << ") /m, albedo ("
            << albedo.x << ", " << albedo.y << ", " << albedo.z << ")\n";
    }
    else glass->setMedium(std::make_shared<HomogeneousMedium>(opt.sigmaT, opt.mediumAlbedo, PhaseFunction(opt.phaseG, opt.phaseG2, opt.phaseW)));
    blobMaterial = glass;
    std::cout << "[Blob] " << blob.numTriangles() << " triangles, volume " << blob.signedVolume() * 1e6 << " cm^3 (must be > 0: closed, outward), sigma_t "
        << opt.sigmaT << " /m, albedo " << opt.mediumAlbedo << ", g " << opt.phaseG << (opt.phaseW < 1 ? " (weight " + std::to_string(opt.phaseW) + ") + g2 " + std::to_string(opt.phaseG2) : std::string()) << ", glass roughness " << opt.glassRoughness << "\n";
    worldList->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(blob)), glass, opt.accel));
    out.keepAlive.push_back(blobMaterial);
    {
        out.camera.from = Point3(0.0, 0.0, 1.0);
        out.camera.at = Point3(0.0, 0.0, 0.0);
        out.camera.fovDeg = 3.44;
    }
    if (opt.bgGray >= 0) out.background = Vector3(opt.bgGray);
    out.afterRender = [cap = opt.volumeCap](const RenderResult& r) {
        std::cout << "[Blob] paths cut by the limit of " << cap << " collisions: " << r.integrator.volumeCapHits()
            << (r.integrator.volumeCapHits() ? "  (these paths are lost: raise --volume-cap)" : "  (none: no cutoff bias)") << "\n";
    };
    return std::nullopt;
}

// --scene manifold
std::optional<int> BuildManifold(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    auto& studioLights = out.lights;
    auto& manifoldBodies = out.manifoldBodies;
    Point3& causticCenter = out.causticTarget.center;
    Real& causticTargetRadius = out.causticTarget.radius;

    // Metres. A dispersive prism and a point light on a gray disc (BuildPrismWorld, shared with the manifoldtest).
    PrismWorld pw = BuildPrismWorld(opt, opt.accel);
    worldList = pw.world;
    manifoldBodies.push_back({ pw.mesh, pw.glass });
    studioLights.push_back(std::make_shared<PointLight>(pw.lamp, RgbColor(1.0, 0.96, 0.90), Real(0.6)));
    causticCenter = Point3(0.0, 0.035, 0.0);        // only used by --caustic-method photons
    causticTargetRadius = Real(0.06);
    // The manifold NEE alone is the default; photons are added to it only if --caustics asks (they then carry the other chains).
    // With --caustic-method photons the photon map is the only source of the caustic, so it always gets photons.
    out.caustics.photonMap = true;
    out.caustics.defaultMethod = "manifold";
    out.caustics.defaultPhotons = 16000000;
    out.caustics.defaultPhotonsWithManifold = 0;
    out.caustics.photonMethodNeedsPhotons = true;
    {
        out.camera.from = Point3(0.16, 0.26, 0.30);
        out.camera.at = Point3(0.12, 0.0, -0.06);
        out.camera.fovDeg = 38.0;
    }
    return std::nullopt;
}

// --scene fogbeam
// A point light behind a narrow slit sends a thin sheet of light through a dispersive prism, in a ball of fog (a closed surface with no
// optical effect around a homogeneous medium). The sheet is lit directly by next event estimation inside the medium; the rainbow fan behind the
// prism (L S S M) is found by light tracing of the photons, which go on through the fog and are connected to the camera at every collision.
std::optional<int> BuildFogBeam(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    const std::string& accel = opt.accel;

    // Metres. The prism stands at the origin, the light travels along +x.
    auto glass = std::make_shared<Dielectric>(Real(1.55), Real(0.0), Real(0.0), opt.abbe);
    std::vector<Point2> tri;
    for (int k = 0; k < 3; ++k) {
        const Real ang = glm::radians(Real(90 + 120 * k - opt.fogPrismDeg));
        tri.push_back(Point2(Real(0.035) * std::cos(ang), Real(0.035) * std::sin(ang)));
    }
    // (the beam hits a face at 30 degrees when the base is parallel to it; about 51 degrees is the minimum deviation of a 60 degree prism of n = 1.55)
    TriangleMesh prism = MakePrismMesh(tri, Real(0.07));
    prism.transform(glm::dmat3(1.0), Vector3(0, -0.035, 0));
    worldList->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(prism)), glass, accel));

    // The fog: a ball of homogeneous medium inside a null interface (ior 1, smooth).
    auto boundary = std::make_shared<Dielectric>(Real(1.0), Real(0.0));
    auto fog = std::make_shared<HomogeneousMedium>(Real(opt.fogSigma), Real(0.95), PhaseFunction(Real(0.3)));
    boundary->setMedium(fog);
    worldList->add(std::make_shared<Sphere>(Point3(0.0), Real(0.32), boundary));
    glass->setBoundaryMedia(nullptr, fog.get());   // the prism stands in the fog: no fog inside the glass, fog again outside
    out.keepAlive.push_back(fog);

    // The aperture: a black plate at x = -0.42 with a slit 4 mm wide and 2.4 cm high (the sheet is then as tall as the prism where it arrives, 7 cm).
    auto plate = std::make_shared<Lambertian>(RgbColor(0.0));
    const Real xp = Real(-0.42), hw = Real(0.002), y0 = Real(-0.012), y1 = Real(0.012), E = Real(0.6);
    worldList->add(std::make_shared<Quad>(Point3(xp, -E, -E), Vector3(0, 2 * E, 0), Vector3(0, 0, E - hw), plate));            // z < -hw
    worldList->add(std::make_shared<Quad>(Point3(xp, -E, hw), Vector3(0, 2 * E, 0), Vector3(0, 0, E - hw), plate));            // z > hw
    worldList->add(std::make_shared<Quad>(Point3(xp, -E, -hw), Vector3(0, y0 + E, 0), Vector3(0, 0, 2 * hw), plate));          // below the slit
    worldList->add(std::make_shared<Quad>(Point3(xp, y1, -hw), Vector3(0, E - y1, 0), Vector3(0, 0, 2 * hw), plate));          // above the slit
    out.keepAlive.push_back(plate);

    out.lights.push_back(std::make_shared<PointLight>(Point3(-0.62, 0.0, 0.0), RgbColor(1.0, 0.97, 0.92), Real(opt.fogLight)));
    // A faint fill light from above lights the fog a little, so that the prism shows as a dark shape in it. It sends no photons: its caustics
    // (a point light, so nothing can find them by chance) are left out.
    if (opt.fogFill > 0) {
        auto fill = std::make_shared<PointLight>(Point3(0.25, 0.45, 0.35), RgbColor(1.0, 0.97, 0.92), Real(opt.fogFill));
        fill->setPhotonCaustics(false);
        out.lights.push_back(fill);
    }
    out.causticTarget.center = Point3(0.0);
    out.causticTarget.radius = Real(0.056);   // the prism's bounding sphere (circumradius 0.035, half height 0.035)
    out.caustics.photonMap = true;
    out.caustics.defaultPhotons = 20000000;
    out.caustics.lightTraceOnly = true;       // the fan in the fog exists only as light-traced photons
    out.caustics.throughMedia = true;
    out.mediaNee = true;                      // the fog is lit by next event estimation, its boundary does not block shadow rays

    // From above, so that the horizontal fan is seen undistorted (+x to the right, +z down the picture; the fan leaves toward -z).
    out.camera.from = Point3(0.08, 0.60, -0.068);
    out.camera.at = Point3(0.08, 0.0, -0.07);
    out.camera.up = Vector3(0, 0, -1);
    out.camera.fovDeg = 26.0;
    return std::nullopt;
}

// --scene fogspot
// A narrow spot light (a cone) crosses a ball of fog: the beam is lit directly, and the fog around it by multiple scattering. Seen from above.
// The camera path finds the thin beam rarely; light paths from the lamp do (--lvc).
std::optional<int> BuildFogSpot(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto boundary = std::make_shared<Dielectric>(Real(1.0), Real(0.0));
    auto fog = std::make_shared<HomogeneousMedium>(Real(opt.fogSigma), Real(0.9), PhaseFunction(Real(0.5)));
    boundary->setMedium(fog);
    out.world->add(std::make_shared<Sphere>(Point3(0.0), Real(opt.spotPrism ? 0.45 : 0.3), boundary));   // (the prism picture is framed inside the ball)
    auto lamp = std::make_shared<PointLight>(Point3(-0.6, 0.0, 0.0), RgbColor(1.0, 0.97, 0.92), Real(opt.fogLight));
    lamp->setCone(Vector3(1, 0, 0), glm::radians(Real(opt.spotDeg)));
    if (opt.spotFloor) {   // a floor inside the ball, lit by the spot: the fog is lit again by the light that the floor sends back
        out.world->add(std::make_shared<Quad>(Point3(-0.5, -0.12, -0.5), Vector3(0, 0, 1.0), Vector3(1.0, 0, 0), std::make_shared<Lambertian>(RgbColor(out.bptAlbedo))));
        lamp = std::make_shared<PointLight>(Point3(-0.6, 0.12, 0.0), RgbColor(1.0, 0.97, 0.92), Real(opt.fogLight));
        lamp->setCone(Vector3(0.7, -0.24, 0.0), glm::radians(Real(opt.spotDeg)));
    }
    if (opt.spotGlass) {   // a glass ball in the fog: it focuses the spot into a caustic cone (light tracing finds it; the photon map does as well)
        auto ball = std::make_shared<Dielectric>(Real(1.5), Real(0.0));
        ball->setBoundaryMedia(nullptr, fog.get());
        out.world->add(std::make_shared<Sphere>(Point3(0.0), Real(0.04), ball));
        out.keepAlive.push_back(ball);
        out.causticTarget.center = Point3(0.0);
        out.causticTarget.radius = Real(0.042);
    }
    if (opt.spotPrism) {   // the dispersive prism of the fogbeam scene, without the slit: the spot lights all of it
        auto glass = std::make_shared<Dielectric>(Real(1.55), Real(0.0), Real(0.0), opt.abbe);
        std::vector<Point2> tri;
        for (int k = 0; k < 3; ++k) {
            const Real ang = glm::radians(Real(90 + 120 * k - opt.fogPrismDeg));
            tri.push_back(Point2(Real(0.035) * std::cos(ang), Real(0.035) * std::sin(ang)));
        }
        TriangleMesh prism = MakePrismMesh(tri, Real(0.07));
        prism.transform(glm::dmat3(1.0), Vector3(0, -0.035, 0));
        glass->setBoundaryMedia(nullptr, fog.get());
        out.world->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(prism)), glass, opt.accel));
        out.keepAlive.push_back(glass);
        out.causticTarget.center = Point3(0.0);
        out.causticTarget.radius = Real(0.056);
        // the photon map aims at the whole prism and assumes an even spot there, so path tracing + photons needs the wide cone; --lvc takes any angle
        lamp->setCone(Vector3(1, 0, 0), glm::radians(opt.lvc ? Real(opt.spotDeg) : std::max(Real(opt.spotDeg), Real(5.5))));
    }
    out.lights.push_back(lamp);
    out.keepAlive.push_back(fog);
    out.bptLight = lamp.get();
    out.bptMedium = fog.get();
    // With a glass object in the fog its caustics exist only as light-traced photons; path tracing's reference needs them (--lvc finds them itself).
    const bool glassInFog = opt.spotGlass || opt.spotPrism;
    out.caustics.photonMap = true;
    out.caustics.defaultPhotons = (glassInFog && !opt.lvc) ? 20000000 : -1;
    out.caustics.lightTraceOnly = glassInFog;
    out.caustics.throughMedia = true;
    out.mediaNee = true;
    if (opt.spotPrism) {   // from above, the fan leaves toward -z (as in the fogbeam scene)
        out.camera.from = Point3(0.115, 0.60, -0.108);
        out.camera.at = Point3(0.115, 0.0, -0.11);
        out.camera.up = Vector3(0, 0, -1);
        out.camera.fovDeg = 32.0;
    }
    else if (opt.spotFloor) {
        out.camera.from = Point3(0.0, 0.5, 0.42);
        out.camera.at = Point3(0.03, -0.06, 0.0);
        out.camera.up = Vector3(0, 1, 0);
        out.camera.fovDeg = 32.0;
    }
    else {
        out.camera.from = Point3(0.0, 0.62, 0.001);
        out.camera.at = Point3(0.0, 0.0, 0.0);
        out.camera.up = Vector3(0, 0, -1);
        out.camera.fovDeg = 28.0;
    }
    return std::nullopt;
}

// --scene caustics
// The gallery: one glass object (--gal-object ball|prism|column|tumbler|gem) on a gray diffuse floor in front of a gray wall, lit by a narrow
// spot from far away, in a faint haze. The glass bends the light into rainbow-edged caustics on the floor and the wall; the picture is made by
// VolumeBPT (--lvc), which finds them with the light paths: the spot's cone covers the object only, so almost every light path is useful.
// Metres. The object stands at the origin; the lamp and the camera are outside the haze ball (radius 1.2).
std::optional<int> BuildCausticGallery(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& world = out.world;
    const Real y0 = Real(0.0004);   // glass rests a hair above the floor (no coplanar surfaces)
    const std::string& name = opt.galObject;
    out.bptAlbedo = Real(opt.galAlbedo);

    // what each object needs: refractive index, the radius of its bounding sphere, where its lamp and the camera are (azimuth from +z toward +x,
    // elevation; degrees), the height of the point that both look at
    Real ior = Real(1.55), bound = Real(0.05), lampAz = Real(-45), lampEl = Real(30), camAz = Real(18), camEl = Real(20), lookY = Real(0.03);
    std::shared_ptr<Dielectric> glass;
    TriangleMesh mesh;
    bool haveMesh = true;
    if (name == "ball") { ior = 1.6; bound = 0.032; lampAz = -50; lampEl = 28; lookY = 0.02; haveMesh = false; }
    else if (name == "prism") { ior = 1.55; bound = 0.056; lampAz = -90; lampEl = 9; camAz = 35; camEl = 28; lookY = 0.035; }
    else if (name == "column") { ior = 1.9; bound = 0.062; lampAz = -55; lampEl = 32; lookY = 0.045; }
    else if (name == "tumbler") { ior = 1.55; bound = 0.052; lampAz = -50; lampEl = 38; lookY = 0.035; }
    else if (name == "gem") { ior = 2.4; bound = 0.048; lampAz = -45; lampEl = 40; lookY = 0.03; }
    else { std::cerr << "[Caustics] unknown --gal-object " << name << " (ball, prism, column, tumbler, gem)\n"; return 1; }
    const Real sc = Real(opt.galScale);   // the object grows by sc (and so do its bounding sphere and the height that the camera looks at)
    bound *= sc; lookY *= sc;
    if (opt.galAz < 900) lampAz = opt.galAz;
    if (opt.galEl < 900) lampEl = opt.galEl;
    glass = std::make_shared<Dielectric>(ior, Real(0.0), Real(0.0), opt.abbe);

    if (name == "ball") {
        world->add(std::make_shared<Sphere>(Point3(0.0, y0 + 0.025 * sc, 0.0), Real(0.025) * sc, glass));
    }
    else if (name == "prism") {
        std::vector<Point2> tri;
        for (int k = 0; k < 3; ++k) {
            const Real ang = glm::radians(Real(90 + 120 * k - opt.fogPrismDeg));
            tri.push_back(Point2(Real(0.035) * std::cos(ang), Real(0.035) * std::sin(ang)));
        }
        mesh = MakePrismMesh(tri, Real(0.07));
    }
    else if (name == "column") mesh = MakeFacetedLathe({ Point2(0, 0), Point2(0.018, 0), Point2(0.018, 0.07), Point2(0, 0.098) }, 6);
    else if (name == "tumbler") mesh = MakeFacetedLathe({ Point2(0, 0), Point2(0.026, 0), Point2(0.032, 0.07), Point2(0.029, 0.07), Point2(0.0235, 0.008), Point2(0, 0.008) }, 8);
    else if (name == "gem") {
        mesh = MakeFacetedLathe({ Point2(0, 0), Point2(0.024, 0.020), Point2(0.024, 0.024), Point2(0.013, 0.033), Point2(0, 0.033) }, 10);
        mesh.transform(glm::dmat3(1.6 * sc), Vector3(0, y0, 0));   // scaled up: a gem is small
        haveMesh = false;
        world->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(mesh)), glass, opt.accel));
    }
    if (haveMesh && name != "ball") {
        mesh.transform(glm::dmat3(sc), Vector3(0, y0, 0));
        world->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(mesh)), glass, opt.accel));
    }

    // the haze: a ball of fog inside a null interface, around everything; the glass stands in it
    auto fog = std::make_shared<HomogeneousMedium>(Real(opt.galFog), Real(0.9), PhaseFunction(Real(0.3)));
    auto boundary = std::make_shared<Dielectric>(Real(1.0), Real(0.0));
    boundary->setMedium(fog);
    const bool floorOn = opt.galFloor != 0;
    const Real dist = opt.galDist > 0 ? Real(opt.galDist) : (floorOn ? Real(1.32) : Real(0.5));   // lamp and camera, from the object
    const Real fogR = (floorOn ? Real(0.9) : Real(0.8)) * dist;   // the haze ball around everything; the lamp and the camera are outside it
    world->add(std::make_shared<Sphere>(Point3(0.0, 0.08, 0.0), fogR, boundary));
    glass->setBoundaryMedia(nullptr, fog.get());
    out.keepAlive.push_back(fog);

    // the floor and the wall behind (gray, albedo bptAlbedo)
    auto gray = std::make_shared<Lambertian>(RgbColor(out.bptAlbedo));
    if (floorOn) {
        world->add(std::make_shared<Quad>(Point3(-1.5, 0.0, -1.5), Vector3(0, 0, 3.0), Vector3(3.0, 0, 0), gray));             // the floor, normal up
        world->add(std::make_shared<Quad>(Point3(-1.5, 0.0, -0.22 * sc), Vector3(3.0, 0, 0), Vector3(0, 1.5, 0), gray));            // the wall at z = -22 cm (times the scale)
    }
    out.keepAlive.push_back(gray);

    // the lamp: far away
    const Point3 target(0.0, lookY, 0.0);
    auto dirOf = [](Real azDeg, Real elDeg) {
        const Real a = glm::radians(azDeg), e = glm::radians(elDeg);
        return Vector3(std::sin(a) * std::cos(e), std::sin(e), std::cos(a) * std::cos(e));
    };
    const Point3 lampPos = target + dirOf(lampAz, lampEl) * dist;
    auto lamp = std::make_shared<PointLight>(lampPos, RgbColor(1.0, 0.97, 0.92), Real(opt.galLight));
    // the lamp lights the floor and the wall around the object (a cone of --gal-spread degrees); the bidirectional estimators send half of their
    // light paths into the small cone that covers the object, where the caustics are made
    lamp->setCone(target - lampPos, glm::radians(Real(opt.galSpread)));
    lamp->setEmissionTarget(target - lampPos, std::atan(bound * Real(1.1) / dist), Real(0.5));
    out.lights.push_back(lamp);
    out.bptLight = lamp.get();
    out.bptMedium = fog.get();

    out.camera.from = target + dirOf(camAz, camEl) * dist;
    out.camera.at = (opt.galLook.y > Real(-900)) ? Point3(opt.galLook) : Point3(0.0, lookY * Real(0.6), -0.04 * sc);
    out.camera.up = Vector3(0, 1, 0);
    out.camera.fovDeg = Real(opt.galFov);
    return std::nullopt;
}

// --scene crystal
std::optional<int> BuildCrystal(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    auto& studioLights = out.lights;
    auto& manifoldBodies = out.manifoldBodies;
    Point3& causticCenter = out.causticTarget.center;
    Real& causticTargetRadius = out.causticTarget.radius;
    // The photon map is the default. --caustic-method manifold solves the T T chains of the flat-faced pieces exactly and leaves the
    // other chains to fewer photons; with water in the tumbler the chains through its wall go on into the water (four refractions,
    // not handled), so it needs --crystal-inside.
    out.caustics.photonMap = true;
    out.caustics.defaultPhotons = opt.crystalInside ? 24000000 : 16000000;
    out.caustics.defaultPhotonsWithManifold = 8000000;
    out.caustics.defaultRadius = Real(opt.crystalInside ? 0.0012 : 0.0006);
    if (!opt.crystalInside) out.caustics.manifoldUnavailable = "--caustic-method manifold in the crystal scene needs --crystal-inside (no water in the tumbler)";

    // Metres. A 40 cm wide disc in an empty black space; faceted glass objects on it (all inside a 15 cm sphere around
    // (0, 4 cm, 0): the photon targeting sphere), one small warm point light above. The flat facets bend light in one
    // direction each, so the caustics are polygons with sharp edges; dispersion splits their rims into colors.
    auto glass = std::make_shared<Dielectric>(Real(1.9), Real(0.0), Real(0.0), opt.abbe);   // heavy flint: strong bending, total internal reflection
    const Real y0 = Real(0.0004);   // glass rests a hair above the disc (no coplanar surfaces)
    auto placeMat = [&](TriangleMesh mesh, const std::shared_ptr<Material>& material, Real x, Real z, Real yawDeg) {
        mesh.transform(glm::dmat3(glm::rotate(glm::dmat4(1.0), glm::radians(yawDeg), Vector3(0, 1, 0))), Vector3(x, y0, z));
        auto placed = std::make_shared<const TriangleMesh>(std::move(mesh));
        worldList->add(MakeMeshHittable(placed, material, opt.accel));
        if (material.get() == glass.get()) manifoldBodies.push_back({ placed, glass });   // a flat-faced piece the manifold NEE can solve
    };
    auto place = [&](TriangleMesh mesh, Real x, Real z, Real yawDeg) { placeMat(std::move(mesh), glass, x, z, yawDeg); };
    // hexagonal crystal column with a pyramid top
    place(MakeFacetedLathe({ Point2(0, 0), Point2(0.018, 0), Point2(0.018, 0.07), Point2(0, 0.098) }, 6), -0.03, -0.02, 10);
    // cut-glass tumbler: 8 flat facets outside and inside, 7 cm tall
    place(MakeFacetedLathe({ Point2(0, 0), Point2(0.026, 0), Point2(0.032, 0.07), Point2(0.029, 0.07), Point2(0.0235, 0.008), Point2(0, 0.008) }, 8),
        0.065, 0.025, 0);
    // ... filled with water to 5.5 cm: its own faceted body 0.4 mm inside the cavity (every interface is air | dielectric, so no
    // nested-dielectric bookkeeping), a flat free surface on top
    if (!opt.crystalInside) {
        auto rInner = [](Real y) { return Real(0.0235) + Real(0.0055) * (y - Real(0.008)) / Real(0.062); };   // vertex radius of the cavity wall
        const Real gap = Real(0.0004), yb = Real(0.008) + gap, yTop = Real(0.055);
        placeMat(MakeFacetedLathe({ Point2(0, yb), Point2(rInner(yb) - gap, yb), Point2(rInner(yTop) - gap, yTop), Point2(0, yTop) }, 8),
            std::make_shared<Dielectric>(Real(1.333), Real(0.0)), 0.065, 0.025, 0);
    }
    // low octagonal block
    place(MakeFacetedLathe({ Point2(0, 0), Point2(0.022, 0), Point2(0.022, 0.028), Point2(0, 0.028) }, 8), -0.075, 0.055, 22);
    // four-sided pyramid
    place(MakeFacetedLathe({ Point2(0, 0), Point2(0.028, 0), Point2(0, 0.045) }, 4), 0.03, -0.075, 20);
    // cut gem (10 facets, pavilion tip down): total internal reflection throws beams and rainbows
    place(MakeFacetedLathe({ Point2(0, 0), Point2(0.024, 0.020), Point2(0.024, 0.024), Point2(0.013, 0.033), Point2(0, 0.033) }, 10), -0.11, -0.03, 5);
    // ball lens: a smooth sphere focuses the light into a hot spot behind it
    // (its own material instance: the manifold NEE solves the flat-faced pieces only, so the photon map must keep the ball's chains)
    worldList->add(std::make_shared<Sphere>(Point3(0.0, y0 + 0.02, 0.07), Real(0.02),
        std::make_shared<Dielectric>(Real(1.9), Real(0.0), Real(0.0), opt.abbe)));

    // the disc (a closed lathe, 1 cm thick): a dark wood-plank floor under a glossy lacquer coat (it mirrors the glass and the lamp),
    // texture coordinates by planar projection of the 40 cm disc
    TriangleMesh discMesh = MakeLatheMesh({ Point2(0, -0.01), Point2(0.2, -0.01), Point2(0.2, 0), Point2(0, 0) }, 96);
    discMesh.uv.resize(discMesh.p.size());
    for (size_t i = 0; i < discMesh.p.size(); ++i) discMesh.uv[i] = Point2(discMesh.p[i].x / Real(0.4) + Real(0.5), discMesh.p[i].z / Real(0.4) + Real(0.5));
    auto disc = std::make_shared<CoatedDiffuse>(RgbColor(1.0), Real(0.18), Real(1.0), Real(0.0), MakeWoodTexture(1024, 5, opt.seaSeed));
    worldList->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(discMesh)), disc, opt.accel));

    // A low, strong light over a dark floor: the light concentrated by the glass (10-50 times brighter) burns out to white.
    if (opt.crystalInside) {
        // the lamp inside the empty tumbler (center of the cavity, 3.5 cm up); it is the only light, and every ray of it leaves
        // through the glass: NEE is blocked by the walls, so all the light on the floor comes from the photon map
        studioLights.push_back(std::make_shared<PointLight>(Point3(0.065, y0 + 0.035, 0.025), RgbColor(1.0, 0.96, 0.88), Real(0.18)));
        std::cout << "[Crystal] lamp inside the tumbler: every photon leaves through the glass\n";
    }
    else
    studioLights.push_back(std::make_shared<PointLight>(Point3(-0.20, 0.22, 0.14), RgbColor(1.0, 0.96, 0.88), Real(1.0)));   // irradiance about 6 on the floor
    causticCenter = Point3(0.0, 0.04, 0.0);
    causticTargetRadius = Real(0.15);
    {
        out.camera.from = Point3(0.0, 0.26, 0.31);
        out.camera.at = Point3(0.01, 0.025, 0.0);
        out.camera.fovDeg = 32.0;
    }
    return std::nullopt;
}

} // namespace rayt::cli
