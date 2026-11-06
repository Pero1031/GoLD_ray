// Caustic photon map check.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

int RunCausticTests(const CliOptions& opt) {
    auto sky = std::make_shared<EnvLight>(std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0)))), 1.0);
    auto world = std::make_shared<Aggregate>();
    world->add(std::make_shared<Quad>(Point3(-5, 0, -5), Vector3(0, 0, 10), Vector3(10, 0, 0), std::make_shared<Lambertian>(RgbColor(0.5))));   // normal +y
    const Real R = 0.05;
    const Point3 center(0, R + 0.0005, 0);
    world->add(std::make_shared<Sphere>(center, R, std::make_shared<Dielectric>(Real(1.5), Real(0.0))));
    auto scene = std::make_shared<Scene>(world);
    scene->setEnvLight(sky);
    scene->finalize("spatial");

    CausticPhotonParams P;
    P.numPhotons = static_cast<long long>(opt.spp) * 8000;
    P.radius = opt.causticRadius > 0 ? opt.causticRadius : Real(0.004);
    P.targetCenter = center;
    P.targetRadius = R * Real(1.05);
    CausticPhotonMap map;
    map.build(*scene, P);

    const int n = opt.spp * 4096;
    const int depth = 64;
    bool allPassed = true;
    std::cout << "[Caustic] " << n << " samples per estimate, floor albedo 0.5, uniform sky 1: no caustic at all would give 0.5\n";
    const Point3 targets[] = { Point3(0.0, 0, 0.0), Point3(0.0, 0, 0.08), Point3(0.10, 0, 0.0), Point3(0.04, 0, -0.03), Point3(-0.07, 0, -0.07) };
    for (const Point3& T : targets) {
        const Point3 from = T + Vector3(0.0, 0.15, 0.30);
        const Ray ray(from, glm::normalize(Vector3(T - from)));
        Vector3 meanA, seA, meanB, seB; int badA = 0, badB = 0;
        EstimateRadiance(*scene, ray, n, depth, meanA, seA, badA, nullptr);
        EstimateRadiance(*scene, ray, n, depth, meanB, seB, badB, &map);
        bool ok = (badA == 0 && badB == 0);
        for (int c = 0; c < 3; ++c) ok = ok && std::abs(meanA[c] - meanB[c]) <= 4.0 * std::sqrt(seA[c] * seA[c] + seB[c] * seB[c]) + 0.03 * meanA[c];
        allPassed = allPassed && ok;
        std::cout << "[Caustic]   floor point (" << T.x << ", " << T.z << "):  path tracing " << meanA.y << " +- " << seA.y
            << "   with photons " << meanB.y << " +- " << seB.y << "   ratio R/G/B " << meanB.x / meanA.x << " / " << meanB.y / meanA.y
            << " / " << meanB.z / meanA.z << (ok ? "   OK" : "   FAIL") << "\n";
    }
    std::cout << "[Caustic] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
