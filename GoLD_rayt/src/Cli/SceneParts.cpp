// Building blocks shared by the demo scenes and the self-checks.
#include "pch.h"
#include "Cli/SceneParts.hpp"
#include "Geometry/ProceduralMesh.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

void AddSphereLight(Aggregate& world, Scene& scene, const Point3& center, Real radius,
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

std::shared_ptr<Light> AddQuadLight(Aggregate& world, Scene* scene, const Point3& corner,
    const Vector3& e1, const Vector3& e2, const EmissionSpectrum& color, Real intensity,
    bool registerAsLight, bool twoSided)
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

std::shared_ptr<Hittable> MakeMeshHittable(const std::shared_ptr<const TriangleMesh>& mesh,
    const std::shared_ptr<Material>& material, const std::string& accel)
{
    if (accel == "embree") return std::make_shared<EmbreeMesh>(mesh, material);
    return BuildMeshBVH(mesh, material);
}

bool LoadWaterAbsorption(color::SpectralCurve& sigmaA) {
    color::SpectralCurve perCm;
    std::string err;
    if (!io::loadSpectralCurveCsv("assets/Water/pope_fry_1997_absorption.csv", perCm, &err)) {
        std::cerr << "[Water] " << err << "\n";
        return false;
    }
    sigmaA = color::SpectralCurve();
    for (const auto& s : perCm.samples()) sigmaA.add(s.lambda_nm, s.value * Real(100));
    sigmaA.finalize();
    return true;
}

std::shared_ptr<Texture2D> MakeWoodTexture(int size, int planks, unsigned seed) {
    auto hash = [](int x, int y, unsigned s) {
        unsigned h = static_cast<unsigned>(x) * 374761393u + static_cast<unsigned>(y) * 668265263u + s * 2147483647u;
        h = (h ^ (h >> 13)) * 1274126177u;
        h ^= h >> 16;
        return double(h & 0xFFFFFFu) / double(0x1000000u);
    };
    auto smooth = [](double t) { return t * t * (3.0 - 2.0 * t); };
    auto noise = [&](double x, double y, unsigned s) {   // value noise in [0, 1]
        const int ix = static_cast<int>(std::floor(x)), iy = static_cast<int>(std::floor(y));
        const double u = smooth(x - ix), v = smooth(y - iy);
        const double a = hash(ix, iy, s), b = hash(ix + 1, iy, s), c = hash(ix, iy + 1, s), d = hash(ix + 1, iy + 1, s);
        return (a * (1 - u) + b * u) * (1 - v) + (c * (1 - u) + d * u) * v;
    };
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::vector<double> tint(static_cast<size_t>(planks)), shift(static_cast<size_t>(planks));
    for (int p = 0; p < planks; ++p) { tint[static_cast<size_t>(p)] = 0.85 + 0.30 * U(rng); shift[static_cast<size_t>(p)] = 100.0 * U(rng); }

    std::vector<unsigned char> rgba(static_cast<size_t>(size) * static_cast<size_t>(size) * 4);
    const double plankW = double(size) / double(planks);
    for (int y = 0; y < size; ++y) {
        for (int x = 0; x < size; ++x) {
            const int p = std::min(static_cast<int>(x / plankW), planks - 1);
            const double sh = shift[static_cast<size_t>(p)];
            const double cx = (x - p * plankW) / plankW;   // 0..1 across the board
            const double along = double(y) / double(size);
            const double warp = 0.5 * noise(along * 3.0 + sh, cx * 2.0, seed + 1) + 0.25 * noise(along * 9.0, cx * 5.0 + sh, seed + 2);
            const double rings = 0.5 + 0.5 * std::sin(2.0 * constants::PI * (cx * 6.0 + warp * 3.0 + sh));
            const double fibers = noise(x * 0.8, y * 0.05 + sh, seed + 3);
            const double tone = tint[static_cast<size_t>(p)] * (0.70 + 0.22 * rings + 0.12 * fibers) * ((cx < 0.012 || cx > 0.988) ? 0.45 : 1.0);
            const size_t i = (static_cast<size_t>(y) * static_cast<size_t>(size) + static_cast<size_t>(x)) * 4;
            rgba[i + 0] = static_cast<unsigned char>(std::clamp(0.62 * tone, 0.0, 1.0) * 255.0);
            rgba[i + 1] = static_cast<unsigned char>(std::clamp(0.42 * tone, 0.0, 1.0) * 255.0);
            rgba[i + 2] = static_cast<unsigned char>(std::clamp(0.26 * tone, 0.0, 1.0) * 255.0);
            rgba[i + 3] = 255;
        }
    }
    return Texture2D::fromRgba8(size, size, std::move(rgba));
}

PrismWorld BuildPrismWorld(const CliOptions& opt, const std::string& accel) {
    PrismWorld w;
    w.world = std::make_shared<Aggregate>();
    w.glass = std::make_shared<Dielectric>(Real(1.55), Real(0.0), Real(0.0), opt.abbe);
    std::vector<Point2> tri;
    for (int k = 0; k < 3; ++k) {
        const Real ang = glm::radians(Real(90 + 120 * k));
        tri.push_back(Point2(Real(0.035) * std::cos(ang), Real(0.035) * std::sin(ang)));
    }
    TriangleMesh prism = MakePrismMesh(tri, Real(0.07));
    prism.transform(glm::dmat3(1.0), Vector3(0, 0.0004, 0));   // a hair above the floor: no coplanar surfaces
    w.mesh = std::make_shared<const TriangleMesh>(std::move(prism));
    w.world->add(MakeMeshHittable(w.mesh, w.glass, accel));
    w.world->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(MakeLatheMesh(
        { Point2(0, -0.01), Point2(0.35, -0.01), Point2(0.35, 0), Point2(0, 0) }, 96)),
        std::make_shared<Lambertian>(RgbColor(0.55)), accel));
    return w;
}

} // namespace rayt::cli
