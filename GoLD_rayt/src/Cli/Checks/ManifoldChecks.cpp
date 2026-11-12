// Specular manifold NEE checks.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/SceneParts.hpp"
#include "Geometry/ProceduralMesh.hpp"
#include "Materials/Lambertian.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Integrator.hpp"
#include "Render/SpecularNee.hpp"

namespace rayt::cli {

int RunManifoldTests(const CliOptions& opt, bool cup) {
    // The two cases: the prism with a lamp beside it, or a hollow faceted tumbler (8 flat facets outside and inside, as in the crystal scene)
    // with a lamp inside it: every ray of the lamp that reaches the floor crosses the wall (T T), and a chord across the cavity must be rejected.
    std::shared_ptr<Aggregate> world;
    std::vector<SpecularNee::Body> bodies;
    Point3 lampPos, photonCenter;
    Real lampScale = 0.6;
    if (!cup) {
        PrismWorld pw = BuildPrismWorld(opt, "own");
        world = pw.world;
        bodies.push_back({ pw.mesh, pw.glass });
        lampPos = pw.lamp;
        photonCenter = Point3(0, 0.035, 0);
    }
    else {
        world = std::make_shared<Aggregate>();
        auto glass = std::make_shared<Dielectric>(Real(1.9), Real(0.0), Real(0.0), opt.abbe);
        TriangleMesh tumbler = MakeFacetedLathe({ Point2(0, 0), Point2(0.026, 0), Point2(0.032, 0.07), Point2(0.029, 0.07),
            Point2(0.0235, 0.008), Point2(0, 0.008) }, 8);
        tumbler.transform(glm::dmat3(1.0), Vector3(0, 0.0004, 0));
        auto mesh = std::make_shared<const TriangleMesh>(std::move(tumbler));
        world->add(MakeMeshHittable(mesh, glass, "own"));
        world->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(MakeLatheMesh(
            { Point2(0, -0.01), Point2(0.35, -0.01), Point2(0.35, 0), Point2(0, 0) }, 96)),
            std::make_shared<Lambertian>(RgbColor(0.55)), "own"));
        bodies.push_back({ mesh, glass });
        lampPos = Point3(0.004, 0.0354, 0.002);   // inside the cavity, a little off the axis (the exact axis would be a degenerate focus)
        lampScale = 0.12;
        photonCenter = Point3(0, 0.04, 0);
    }
    auto lamp = std::make_shared<PointLight>(lampPos, RgbColor(1.0, 0.96, 0.90), lampScale);
    auto scene = std::make_shared<Scene>(world);
    scene->addLight(lamp);
    scene->finalize("uniform");
    SpecularNee nee(bodies, *scene);
    std::cout << "[Manifold] " << (cup ? "tumbler" : "prism") << ": " << nee.numFacets() << " facets, " << nee.numLights() << " point light(s)\n";

    // Lit floor points (x >= 5 cm: a camera ray from above and behind cannot clip the prism) found with the manifold estimator itself.
    auto floorMat = std::make_shared<Lambertian>(RgbColor(0.55));
    SurfaceInteraction rec;
    rec.n = rec.gn = Vector3(0, 1, 0);
    rec.matPtr = floorMat.get();
    const Vector3 wo = glm::normalize(Vector3(0.3, 1.0, 0.5));
    struct Cand { Point3 p; Real lum; };
    std::vector<Cand> cands;
    Real maxLum = 0;
    auto lumAt = [&](const Point3& p, int samples = 24) {   // mean luminance of the manifold estimate at a floor point (the wavelengths vary)
        rec.p = p;
        Real lum = 0;
        for (int s = 0; s < samples; ++s) {
#if defined(RAYT_SPECTRAL_MODE)
            SampledWavelengths lambda = SampledWavelengths::sampleVisible(sampling::Random());
#else
            SampledWavelengths lambda;
#endif
            const Vector3 c = SpectrumToRGB(nee.sample(*scene, rec, wo, lambda), lambda);
            lum += Real(0.2126) * c.x + Real(0.7152) * c.y + Real(0.0722) * c.z;
        }
        return lum / Real(samples);
    };
    std::vector<Point3> grid;
    if (cup) {   // a ring around the tumbler (radius 5 - 15 cm): its beams land there
        for (int iz = -15; iz <= 15; ++iz)
            for (int ix = -15; ix <= 15; ++ix) {
                const Real r = std::hypot(Real(ix), Real(iz)) * Real(0.01);
                if (r >= Real(0.035) && r <= Real(0.15)) grid.push_back(Point3(0.01 * ix, 0, 0.01 * iz));   // from the foot of the tumbler (r 3.2 cm) outward
            }
    }
    else {
        for (int iz = -20; iz <= 5; ++iz)
            for (int ix = 0; ix <= 27; ++ix) grid.push_back(Point3(0.05 + 0.01 * ix, 0, 0.01 * iz));
    }
    for (const Point3& p : grid) {
        const Real lum = lumAt(p);
        cands.push_back({ p, lum });
        maxLum = std::max(maxLum, lum);
    }
    std::cout << "[Manifold] scanned " << cands.size() << " floor points, brightest manifold radiance " << maxLum << "\n";
    std::vector<Point3> points;
    for (const Cand& c : cands) {
        if (c.lum < Real(0.1) * maxLum || c.lum > Real(0.5) * maxLum) continue;   // mid-level, away from the brightest focus lines
        // Smooth neighborhood: the photon map blurs over its gather radius (3 mm), the manifold estimate does not, so near an edge of the
        // pattern (where the focusing diverges, e.g. with a single IOR) the two cannot agree. Keep points whose value varies by < 15 % over 3 mm
        // (128 samples each: with dispersion one sample carries one wavelength, so the estimate itself fluctuates by several percent).
        bool smooth = true;
        const Real centre = lumAt(c.p, 128);
        for (const Vector3& d : { Vector3(0.003, 0, 0), Vector3(-0.003, 0, 0), Vector3(0, 0, 0.003), Vector3(0, 0, -0.003) })
            if (std::abs(lumAt(c.p + d, 128) - centre) > Real(0.15) * centre) { smooth = false; break; }
        if (!smooth) continue;
        bool spaced = true;   // (not "far": a macro in the Windows headers)
        for (const Point3& q : points) if (glm::length(Vector3(q - c.p)) < Real(0.04)) { spaced = false; break; }
        if (spaced) points.push_back(c.p);
        if (points.size() >= 6) break;
    }
    if (points.empty()) { std::cout << "[Manifold] no lit floor point found: FAIL\n"; return 2; }

    CausticPhotonParams P;
    P.numPhotons = 16000000;
    P.radius = Real(0.003);
    P.targetCenter = photonCenter;
    P.targetRadius = Real(0.06);
    P.onlyTransmissions = 2;   // the manifold NEE handles the chain T T only (other chains: reflections, total internal reflection, ... are not in it)
    CausticPhotonMap map;
    map.build(*scene, P);

    const int n = opt.spp * 256;
    bool allPassed = true;

    // Coverage over the whole grid. The six points above are chosen where the manifold solver finds something, so they cannot show a place where
    // it finds nothing although light arrives (in the hybrid render the T T photons are dropped, so such a place would be a dark hole).
    {
        auto lum = [](const Vector3& c) { return Real(0.2126) * c.x + Real(0.7152) * c.y + Real(0.0722) * c.z; };
        struct Row { Point3 p; Real dM, dP, se; };
        std::vector<Row> rows;
        Real sumM = 0, sumP = 0, maxP = 0;
        for (const Point3& p : grid) {
            const Vector3 outward = glm::normalize(Vector3(p.x, 0, p.z));
            const Point3 from = p + outward * Real(0.12) + Vector3(0.0, 0.20, 0.0);
            const Ray ray(from, glm::normalize(Vector3(p - from)));
            Vector3 m0, s0, mM, sM, mP, sP; int b0 = 0, bM = 0, bP = 0;
            EstimateRadiance(*scene, ray, 512, 32, m0, s0, b0, nullptr, nullptr);
            EstimateRadiance(*scene, ray, 512, 32, mM, sM, bM, nullptr, &nee);
            EstimateRadiance(*scene, ray, 512, 32, mP, sP, bP, &map, nullptr);
            Row r{ p, lum(mM) - lum(m0), lum(mP) - lum(m0), std::sqrt(lum(sM) * lum(sM) + lum(sP) * lum(sP)) };
            rows.push_back(r);
            sumM += r.dM; sumP += r.dP; maxP = std::max(maxP, r.dP);
        }
        int lit = 0, missed = 0, excess = 0;
        Real missedEnergy = 0;
        const Row* worst = nullptr;
        for (const Row& r : rows) {
            if (r.dP < Real(0.05) * maxP || r.dP < 5 * r.se) continue;   // only places where the photon map shows a clear caustic
            ++lit;
            if (r.dM < Real(0.5) * r.dP) { ++missed; missedEnergy += r.dP - r.dM; if (!worst || r.dP - r.dM > worst->dP - worst->dM) worst = &r; }
            if (r.dM > Real(1.5) * r.dP + 4 * r.se) ++excess;
        }
        std::cout << "[Manifold] coverage over " << rows.size() << " floor points: sum manifold " << sumM << ", sum photon map " << sumP
            << " (ratio " << (sumP != 0 ? sumM / sumP : Real(0)) << ")\n";
        std::cout << "[Manifold]   points with a clear photon-map caustic: " << lit << ", of which the manifold gets less than half: " << missed
            << " (missing " << (sumP != 0 ? Real(100) * missedEnergy / sumP : Real(0)) << " % of the photon-map energy), clearly more than the photon map: "
            << excess << "\n";
        {   // where the misses are: by distance from the axis (cup) or from the first grid column (prism), and the first few points
            Real bandM[8] = {}, bandP[8] = {}; int bandMiss[8] = {}, bandLit[8] = {};
            for (const Row& r : rows) {
                const Real dist = cup ? std::hypot(r.p.x, r.p.z) : (r.p.x - Real(0.05));
                const int b = std::clamp(static_cast<int>(dist / Real(0.02)), 0, 7);
                bandM[b] += r.dM; bandP[b] += r.dP;
                if (r.dP >= Real(0.05) * maxP && r.dP >= 5 * r.se) { ++bandLit[b]; if (r.dM < Real(0.5) * r.dP) ++bandMiss[b]; }
            }
            for (int b = 0; b < 8; ++b)
                if (bandLit[b] > 0 || bandP[b] != 0)
                    std::cout << "[Manifold]   band " << b * 2 << "-" << b * 2 + 2 << " cm: manifold/photon " << (bandP[b] != 0 ? bandM[b] / bandP[b] : Real(0))
                        << ", lit " << bandLit[b] << ", missed " << bandMiss[b] << "\n";
            int shown = 0;
            for (const Row& r : rows)
                if (shown < 12 && r.dP >= Real(0.05) * maxP && r.dP >= 5 * r.se && r.dM < Real(0.5) * r.dP) {
                    std::cout << "[Manifold]   miss at (" << r.p.x << ", " << r.p.z << "): manifold " << r.dM << " photon map " << r.dP << "\n";
                    if (shown < 3) { rec.p = r.p; nee.diagnose(*scene, rec, std::cout); }
                    ++shown;
                }
        }
        if (worst) std::cout << "[Manifold]   worst miss at (" << worst->p.x << ", " << worst->p.z << "): manifold " << worst->dM << " photon map " << worst->dP << "\n";
        if (sumP != 0 && std::abs(sumM / sumP - Real(1)) > Real(0.10)) allPassed = false;
        if (lit > 0 && missed > lit / 20) allPassed = false;
    }

    for (const Point3& p : points) {
        const Vector3 outward = glm::normalize(Vector3(p.x, 0, p.z));   // the camera ray comes from outside and above: it cannot clip the glass
        const Point3 from = p + outward * Real(0.12) + Vector3(0.0, 0.20, 0.0);
        const Ray ray(from, glm::normalize(Vector3(p - from)));
        Vector3 m0, s0, mM, sM, mP, sP; int b0 = 0, bM = 0, bP = 0;
        EstimateRadiance(*scene, ray, n, 32, m0, s0, b0, nullptr, nullptr);
        EstimateRadiance(*scene, ray, n, 32, mM, sM, bM, nullptr, &nee);
        EstimateRadiance(*scene, ray, n, 32, mP, sP, bP, &map, nullptr);
        auto lum = [](const Vector3& c) { return Real(0.2126) * c.x + Real(0.7152) * c.y + Real(0.0722) * c.z; };
        const Real dM = lum(mM) - lum(m0), dP = lum(mP) - lum(m0);
        const Real se = std::sqrt(lum(sM) * lum(sM) + lum(sP) * lum(sP));
        // a real comparison needs a caustic signal in the photon map (a pass between two zeros proves nothing) and a manifold solution
        const bool signal = std::abs(dP) > 5 * lum(sP);
        const bool ok = (b0 + bM + bP == 0) && signal && std::abs(dM - dP) <= 4 * se + Real(0.10) * std::abs(dP);
        allPassed = allPassed && ok;
        std::cout << "[Manifold]   floor point (" << p.x << ", " << p.z << "):  direct only " << lum(m0) << "   caustic: manifold " << dM
            << "   photon map " << dP << "   ratio " << (dP != 0 ? dM / dP : Real(0)) << (ok ? "   OK" : "   FAIL") << "\n";
    }
    const auto& c = nee.counters();
    std::cout << "[Manifold] solver: " << nee.solutions() << " solutions in " << nee.queries() << " receiver evaluations\n"
        << "[Manifold]   pairs tried " << c.pairs << " -> Newton converged " << c.converged << " -> inside the facets " << c.insideFacets
        << " -> real refraction " << c.refraction << " -> chord in glass " << c.chordFree << " -> unblocked " << c.unblocked
        << " -> focusing ok " << c.focusOk << "\n";
    std::cout << "[Manifold]   Snell-seed landing on B (B radii <1,<2,<4,<8,<16,more,none): solutions/pairs";
    for (int i = 0; i < 7; ++i) std::cout << "  " << c.landOk[i] << "/" << c.landTried[i];
    std::cout << "\n";
    if (nee.solutions() == 0) allPassed = false;
    std::cout << "[Manifold] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
