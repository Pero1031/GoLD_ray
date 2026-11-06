// Scenes about glass and media: glass spheres, the white furnace, a glass blob, a prism, faceted crystal.
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
    glass->setMedium(std::make_shared<HomogeneousMedium>(opt.sigmaT, opt.mediumAlbedo, opt.phaseG));
    blobMaterial = glass;
    std::cout << "[Blob] " << blob.numTriangles() << " triangles, volume " << blob.signedVolume() * 1e6 << " cm^3 (must be > 0: closed, outward), sigma_t "
        << opt.sigmaT << " /m, albedo " << opt.mediumAlbedo << ", g " << opt.phaseG << ", glass roughness " << opt.glassRoughness << "\n";
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
    {
        out.camera.from = Point3(0.16, 0.26, 0.30);
        out.camera.at = Point3(0.12, 0.0, -0.06);
        out.camera.fovDeg = 38.0;
    }
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
