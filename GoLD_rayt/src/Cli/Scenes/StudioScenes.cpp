// Studio-style scenes: lit objects, triangle meshes, skin.
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
#include "Lights/AreaLight.hpp"
#include "Materials/CoatedDiffuse.hpp"
#include "Materials/DiffuseLight.hpp"
#include "Materials/LayeredSkin.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Scene/Primitive.hpp"

namespace rayt::cli {

// --scene studio
std::optional<int> BuildStudio(SceneContext& ctx, SceneSetup& out) {
    auto& worldList = out.world;
    const auto& matGoldMedium = ctx.mats.goldMedium;
    const auto& matGlassDisp = ctx.mats.glassDisp;
    auto& studioLights = out.lights;

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
    return std::nullopt;
}

// --scene mesh
std::optional<int> BuildMesh(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    const auto& matGoldMedium = ctx.mats.goldMedium;
    auto& studioLights = out.lights;
    out.checksPassed = RunMeshChecks(opt);   // the mesh checks run first

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
    return std::nullopt;
}

// --scene skindemo
std::optional<int> BuildSkinDemo(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    auto& studioLights = out.lights;
    std::shared_ptr<Material> skinMaterial;                       // keeps the shared skin material alive

    // Everything is in metres: a 3 cm sphere and a 2.5 cm ring (OBJ torus scaled by 0.05) on a floor, so the
    // 1-4 mm scatter distances of skin are a visible fraction of the objects. One material (hence one SSS
    // profile) for both, because the exit point of the light is only searched on surfaces with the same profile.
    RgbColor skinColor(0.80, 0.52, 0.42);
    std::shared_ptr<const SubsurfaceProfile> skinProfile;
    if (!opt.skinTable.empty()) {   // a layered-skin table: the table, or Burley's profile with the same albedo and mean radius
        std::string tableErr;
        const auto table = skin::loadTableCsv(opt.skinTable, &tableErr);
        if (!table) { std::cerr << "[SkinDemo] " << tableErr << "\n"; return 1; }
        const Vector3 A(table->rgb[0]->albedo(), table->rgb[1]->albedo(), table->rgb[2]->albedo());
        skinColor = RgbColor(A);
        if (opt.skinMatched) {
            const Vector3 ell(skin::matchedScatterDistance(*table->rgb[0]), skin::matchedScatterDistance(*table->rgb[1]),
                skin::matchedScatterDistance(*table->rgb[2]));
            skinProfile = std::make_shared<SubsurfaceProfile>(ell, Real(1.4));
            std::cout << "[SkinDemo] Burley profile matched to " << opt.skinTable << ": albedo (" << A.x << ", " << A.y << ", " << A.z
                << "), scatter distance (mm) (" << ell.x * 1000 << ", " << ell.y * 1000 << ", " << ell.z * 1000 << ")\n";
        }
        else {
            skinProfile = std::make_shared<SubsurfaceProfile>(table, Real(1.4));
            std::cout << "[SkinDemo] layered skin profile from " << opt.skinTable << ": albedo (" << A.x << ", " << A.y << ", " << A.z << ")\n";
        }
    }
    else if (opt.sss) {
        skinProfile = std::make_shared<SubsurfaceProfile>(Vector3(3.67, 1.37, 0.68) * Real(0.001) * opt.sssScale, Real(1.4));
        std::cout << "[SkinDemo] subsurface scattering on (scale " << opt.sssScale << ")\n";
    }
    else std::cout << "[SkinDemo] no subsurface scattering (plain coated diffuse)\n";
    auto skin = std::make_shared<CoatedDiffuse>(skinColor, Real(0.50), Real(0.4), Real(0));
    if (skinProfile) skin->setSubsurface(skinProfile);
    skinMaterial = skin;

    worldList->add(std::make_shared<Sphere>(Point3(-0.065, 0.03, 0.0), 0.03, skin));

    TriangleMesh ring;
    std::string objErr;
    if (!io::loadObj("assets/Meshes/torus.obj", ring, &objErr)) { std::cerr << "[SkinDemo] " << objErr << "\n"; return 1; }
    ring.transform(glm::dmat3(glm::rotate(glm::dmat4(1.0), glm::radians(Real(65)), Vector3(1, 0, 0))) * Real(0.05),
        Vector3(0.065, 0.036, 0.0));
    worldList->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(ring)), skin, opt.accel));

    worldList->add(std::make_shared<Quad>(Point3(-0.6, 0, -0.6), Vector3(0, 0, 1.2), Vector3(1.2, 0, 0),
        std::make_shared<Lambertian>(RgbColor(0.30, 0.30, 0.32))));                              // floor (normal +y)
    worldList->add(std::make_shared<Quad>(Point3(-0.6, 0, -0.30), Vector3(1.2, 0, 0), Vector3(0, 0.8, 0),
        std::make_shared<Lambertian>(RgbColor(0.10, 0.10, 0.12))));                              // dark backdrop (normal +z)

    const Point3 target(0.0, 0.03, 0.0);
    auto softbox = [&](const Point3& center, Real w, Real h, const EmissionSpectrum& e, Real radiance) {
        const Vector3 n = glm::normalize(target - center);
        const Vector3 e1 = glm::normalize(glm::cross(Vector3(0, 1, 0), n)) * w;
        const Vector3 e2 = glm::normalize(glm::cross(n, e1)) * h;
        studioLights.push_back(AddQuadLight(*worldList, nullptr, center - e1 * 0.5 - e2 * 0.5, e1, e2, e, radiance, true));
    };
    softbox(Point3(-0.34, 0.18, 0.14), 0.25, 0.25, EmissionSpectrum(RgbColor(1.0)), 8.0);             // key: from the left, so the terminator is mid-picture
    softbox(Point3(0.35, 0.10, 0.30), 0.30, 0.20, EmissionSpectrum(RgbColor(1.0)), 1.0);              // fill: front-right, dim
    softbox(Point3(0.10, 0.32, -0.16), 0.25, 0.06, EmissionSpectrum(RgbColor(1.0, 0.93, 0.85)), 60.0); // back light: strip above-behind, out of frame
    out.keepAlive.push_back(skinMaterial);
    {
        out.camera.from = Point3(0.0, 0.075, 0.30);
        out.camera.at = Point3(0.0, 0.032, 0.0);
        out.camera.fovDeg = 30.0;
    }
    return std::nullopt;
}

} // namespace rayt::cli
