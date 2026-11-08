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

    // ---- light tracing (--caustic-trace): the caustic on a floor point the camera sees directly, splatted from the photon paths, plus the rest
    // (direct light and the other paths, the photon gather skipped there) against plain path tracing ----
    {
        const Point3 from(0.0, 0.28, 0.38);
        const int W = 160, H = 160;
        Camera cam(from, Point3(0, 0, 0), Vector3(0, 1, 0), 32.0, 1.0, 0.0, glm::length(from));
        CausticPhotonParams P2 = P;
        P2.numPhotons = static_cast<long long>(opt.spp) * 32000;
        P2.lightTraceCamera = &cam;
        P2.traceWidth = W;
        P2.traceHeight = H;
        CausticPhotonMap map2;
        map2.build(*scene, P2);
        const auto& splat = map2.splat();
        std::cout << "[Caustic] light tracing: a " << W << "x" << H << " pixel camera, the caustic splatted at the receiver; total = splat (5x5 pixels) + paths without the gather\n";
        for (const Point3& T : targets) {
            Real s, t, cosC, dist;
            if (!cam.project(T, s, t, cosC, dist)) { std::cout << "[Caustic]   floor point (" << T.x << ", " << T.z << ") is outside the view   FAIL\n"; allPassed = false; continue; }
            const int px = static_cast<int>(s * W), py = H - 1 - static_cast<int>(t * H);
            Vector3 splatMean(0.0);
            int cnt = 0;
            for (int dy = -2; dy <= 2; ++dy) for (int dx = -2; dx <= 2; ++dx) {
                const int x = px + dx, y = py + dy;
                if (x < 0 || y < 0 || x >= W || y >= H) continue;
                splatMean += splat[static_cast<size_t>(y) * W + x];
                ++cnt;
            }
            splatMean /= Real(std::max(cnt, 1));
            const Ray ray = cam.getRay((Real(px) + Real(0.5)) / Real(W), Real(1) - (Real(py) + Real(0.5)) / Real(H), Point2(0.5, 0.5));
            Vector3 meanA, seA, meanB, seB; int badA = 0, badB = 0;
            EstimateRadiance(*scene, ray, n, depth, meanA, seA, badA, nullptr);
            EstimateRadiance(*scene, ray, n, depth, meanB, seB, badB, &map2, nullptr, true);
            const Vector3 total = meanB + splatMean;
            bool ok = (badA == 0 && badB == 0);
            for (int c = 0; c < 3; ++c) ok = ok && std::abs(meanA[c] - total[c]) <= 4.0 * std::sqrt(seA[c] * seA[c] + seB[c] * seB[c]) + 0.06 * meanA[c];
            allPassed = allPassed && ok;
            std::cout << "[Caustic]   floor point (" << T.x << ", " << T.z << "):  path tracing " << meanA.y << " +- " << seA.y
                << "   light tracing " << total.y << " (caustic part " << splatMean.y << ")   ratio R/G/B " << total.x / meanA.x << " / " << total.y / meanA.y
                << " / " << total.z / meanA.z << (ok ? "   OK" : "   FAIL") << "\n";
        }
    }

    // ---- volumetric caustics: a glass ball in a ball of fog under a uniform sky. Path tracing alone finds the light that the glass focuses into
    // the fog (the sky is large); light tracing with a photon map that goes through the fog finds the same light from the other end ----
    {
        auto fogWorld = std::make_shared<Aggregate>();
        auto boundary = std::make_shared<Dielectric>(Real(1.0), Real(0.0));
        auto fog = std::make_shared<HomogeneousMedium>(Real(8), Real(0.9), PhaseFunction(Real(0.3)));
        boundary->setMedium(fog);
        fogWorld->add(std::make_shared<Sphere>(Point3(0.0), Real(0.25), boundary));
        const Real Rg = 0.04;
        auto ball = std::make_shared<Dielectric>(Real(1.5), Real(0.0));
        ball->setBoundaryMedia(nullptr, fog.get());   // the ball stands in the fog
        fogWorld->add(std::make_shared<Sphere>(Point3(0.0), Rg, ball));
        auto fogScene = std::make_shared<Scene>(fogWorld);
        fogScene->setEnvLight(sky);
        fogScene->finalize("spatial");

        const int nFog = n * 8;   // the caustic is a small part of the fog radiance: more samples
        const Point3 from(0.0, 0.12, 0.7);
        const int W = 160, H = 160;
        Camera cam(from, Point3(0, 0, 0), Vector3(0, 1, 0), 22.0, 1.0, 0.0, glm::length(from));
        CausticPhotonParams P3;
        P3.numPhotons = static_cast<long long>(opt.spp) * 32000;
        P3.radius = Real(0.004);
        P3.targetCenter = Point3(0.0);
        P3.targetRadius = Rg * Real(1.05);
        P3.lightTraceCamera = &cam;
        P3.traceWidth = W;
        P3.traceHeight = H;
        P3.traceMedia = true;
        CausticPhotonMap map3;
        map3.build(*fogScene, P3);
        const auto& splat = map3.splat();
        std::cout << "[Caustic] volumetric light tracing: a glass ball in fog (sigma_t 8, albedo 0.9, g 0.3) under a uniform sky; path tracing alone vs. path tracing + the splat (5x5 pixels)\n";
        const Point3 fogTargets[] = { Point3(0.0, 0.0, -0.06), Point3(0.07, 0.0, 0.0), Point3(-0.06, 0.03, -0.12), Point3(0.02, -0.07, -0.05), Point3(0.0, 0.1, 0.05) };
        for (const Point3& T : fogTargets) {
            Real s, t, cosC, dist;
            if (!cam.project(T, s, t, cosC, dist)) { std::cout << "[Caustic]   fog point is outside the view   FAIL\n"; allPassed = false; continue; }
            const int px = static_cast<int>(s * W), py = H - 1 - static_cast<int>(t * H);
            Vector3 splatMean(0.0);
            int cnt = 0;
            for (int dy = -2; dy <= 2; ++dy) for (int dx = -2; dx <= 2; ++dx) {
                const int x = px + dx, y = py + dy;
                if (x < 0 || y < 0 || x >= W || y >= H) continue;
                splatMean += splat[static_cast<size_t>(y) * W + x];
                ++cnt;
            }
            splatMean /= Real(std::max(cnt, 1));
            const Ray ray = cam.getRay((Real(px) + Real(0.5)) / Real(W), Real(1) - (Real(py) + Real(0.5)) / Real(H), Point2(0.5, 0.5));
            Vector3 meanA, seA, meanB, seB; int badA = 0, badB = 0;
            EstimateRadiance(*fogScene, ray, nFog, depth, meanA, seA, badA, nullptr, nullptr, false, true);
            EstimateRadiance(*fogScene, ray, nFog, depth, meanB, seB, badB, &map3, nullptr, true, true);
            const Vector3 total = meanB + splatMean;
            bool ok = (badA == 0 && badB == 0);
            for (int c = 0; c < 3; ++c) ok = ok && std::abs(meanA[c] - total[c]) <= 4.0 * std::sqrt(seA[c] * seA[c] + seB[c] * seB[c]) + 0.003 * meanA[c];   // the caustic is about 1 % of the fog radiance: a tight bound
            allPassed = allPassed && ok;
            std::cout << "[Caustic]   fog point (" << T.x << ", " << T.y << ", " << T.z << "):  path tracing " << meanA.y << " +- " << seA.y
                << "   with light tracing " << total.y << " (splat " << splatMean.y << ", paths dropped " << meanA.y - meanB.y << ")   ratio R/G/B " << total.x / meanA.x << " / " << total.y / meanA.y
                << " / " << total.z / meanA.z << (ok ? "   OK" : "   FAIL") << "\n";
        }
    }

    std::cout << "[Caustic] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

// ---- VolumeBPT (prototype of LVC-BPT) against path tracing ----
int RunLvcTests(const CliOptions& opt) {
    const int W = 48, H = 27;
    auto boundary = std::make_shared<Dielectric>(Real(1.0), Real(0.0));
    auto fog = std::make_shared<HomogeneousMedium>(Real(3), Real(0.9), PhaseFunction(Real(0.5)));
    boundary->setMedium(fog);
    auto world = std::make_shared<Aggregate>();
    world->add(std::make_shared<Sphere>(Point3(0.0), Real(0.3), boundary));
    auto lamp = std::make_shared<PointLight>(Point3(-0.6, 0.0, 0.0), RgbColor(1.0, 1.0, 1.0), Real(60));
    lamp->setCone(Vector3(1, 0, 0), glm::radians(Real(8)));
    auto scene = std::make_shared<Scene>(world);
    scene->addLight(lamp);
    scene->finalize("uniform");
    const Point3 from(0.0, 0.62, 0.001);
    auto cam = std::make_shared<Camera>(from, Point3(0, 0, 0), Vector3(0, 0, -1), 28.0, double(W) / H, 0.0, glm::length(from));

    // the reference: path tracing with next event estimation in the medium and MIS with the phase-sampled rays
    Film ref(W, H);
    PathIntegrator pt(cam, 32, opt.spp * 64);
    pt.setMediaNee(true);
    pt.render(*scene, ref);

    auto blockMeans = [&](const Vector3* d, std::vector<double>& out, double& total) {
        const int bx = 6, by = 3;   // 6 x 3 blocks of 8 x 9 pixels
        out.assign(bx * by, 0.0);
        total = 0;
        for (int y = 0; y < H; ++y) for (int x = 0; x < W; ++x) {
            const double v = d[static_cast<size_t>(y) * W + x].y;
            out[(y / (H / by)) * bx + x / (W / bx)] += v / double((W / bx) * (H / by));
            total += v / double(W * H);
        }
    };
    std::vector<double> refB; double refTotal;
    blockMeans(ref.getData(), refB, refTotal);
    double refMax = 0; for (double v : refB) refMax = std::max(refMax, v);
    std::cout << "[Lvc] reference (path tracing, " << opt.spp * 64 << " spp): mean G " << refTotal << ", brightest block " << refMax << "\n";

    bool allPassed = true;
    struct Case { const char* name; int mask; };
    for (const Case& c : { Case{ "next event estimation only (s = 0)", 1 }, Case{ "light tracing only (s = n)", 4 }, Case{ "next event estimation + connections (s < n)", 3 }, Case{ "connections + light tracing (s > 0)", 6 }, Case{ "all strategies, balance heuristic", 7 } }) {
        VolumeBptParams bp;
        bp.width = W; bp.height = H; bp.spp = opt.spp * 8; bp.lightPaths = static_cast<long long>(opt.spp) * 100000; bp.maxVertices = 8;
        bp.useNee = (c.mask & 1) != 0; bp.useConnections = (c.mask & 2) != 0; bp.useLightTracing = (c.mask & 4) != 0;
        Film f(W, H);
        std::cout << "[Lvc] " << c.name << "\n";
        VolumeBPT(*scene, *cam, *lamp, *fog, bp).render(f);
        std::vector<double> b; double total;
        blockMeans(f.getData(), b, total);
        double worst = 0; int bad = 0;
        for (size_t i = 0; i < b.size(); ++i) {
            if (refB[i] < 0.05 * refMax) continue;   // dark blocks: mostly noise
            const double e = std::abs(b[i] - refB[i]) / refB[i];
            worst = std::max(worst, e);
            if (e > 0.20) ++bad;   // noise of a few thousand samples per block; a wrong weight shifts the total by more
        }
        const double ratio = total / refTotal;
        const bool ok = std::abs(ratio - 1.0) < 0.03 && bad == 0;
        allPassed = allPassed && ok;
        std::cout << "[Lvc]   mean G " << total << "  ratio to the reference " << ratio << "  worst bright block off by " << 100 * worst << " %" << (ok ? "   OK" : "   FAIL") << "\n";
    }
    std::cout << "[Lvc] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
