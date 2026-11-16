// Manifold connections in VolumeBPT (docs/LVC_BPT.md, "Manifold connections"): checks against the photon map and path tracing.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/SceneParts.hpp"
#include "Geometry/ProceduralMesh.hpp"
#include "Materials/Lambertian.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Integrator.hpp"
#include "Render/SpecularNee.hpp"
#include "Render/VolumeBPT.hpp"

namespace rayt::cli {

namespace {

struct LvcPrismScene {
    PrismWorld pw;
    std::shared_ptr<Scene> scene;
    std::shared_ptr<PointLight> lamp;
    std::shared_ptr<HomogeneousMedium> fog;
    std::shared_ptr<Camera> cam;
    std::vector<SpecularNee::Body> bodies;
};

// The prism of the manifold check (IOR 1.55, no dispersion) on a gray disc, a lamp whose cone covers only the prism, and no fog (the medium of
// VolumeBPT is never entered). The camera looks at the floor beside the prism (the caustic region) or through the prism.
LvcPrismScene BuildLvcPrismScene(const CliOptions& opt, int W, int H, bool throughPrism) {
    LvcPrismScene s;
    CliOptions o = opt;
    o.abbe = 0;
    s.pw = BuildPrismWorld(o, "own");
    s.bodies.push_back({ s.pw.mesh, s.pw.glass });
    s.lamp = std::make_shared<PointLight>(s.pw.lamp, RgbColor(1.0, 1.0, 1.0), Real(0.6));
    const Point3 centre(0, 0.035, 0);
    s.lamp->setCone(Vector3(centre - s.pw.lamp), glm::radians(Real(17)));
    s.scene = std::make_shared<Scene>(s.pw.world);
    s.scene->addLight(s.lamp);
    s.scene->finalize("uniform");
    s.fog = std::make_shared<HomogeneousMedium>(Real(0.01), Real(0.9), PhaseFunction(Real(0.5)));
    const Point3 from = throughPrism ? Point3(0.16, 0.12, 0.10) : Point3(0.30, 0.25, 0.30);
    const Point3 at = throughPrism ? Point3(0.0, 0.01, 0.0) : Point3(0.16, 0.0, -0.06);
    s.cam = std::make_shared<Camera>(from, at, Vector3(0, 1, 0), throughPrism ? 22.0 : 30.0, double(W) / H, 0.0, glm::length(from - at));
    return s;
}

struct BlockStats { std::vector<double> b; double total = 0; };

BlockStats Blocks(const Vector3* d, int W, int H) {
    const int bx = 6, by = 3;
    BlockStats r;
    r.b.assign(bx * by, 0.0);
    for (int y = 0; y < H; ++y) for (int x = 0; x < W; ++x) {
        const double v = d[static_cast<size_t>(y) * W + x].y;
        r.b[(y / (H / by)) * bx + x / (W / bx)] += v / double((W / bx) * (H / by));
        r.total += v / double(W * H);
    }
    return r;
}

// Compares an image with a reference block by block and prints the ratios (the reference value in brackets). True if the mean is within meanTol and
// every block with at least 5 % of the brightest reference block is within blockTol.
bool CompareBlocks(const char* what, const BlockStats& g, const BlockStats& ref, double meanTol, double blockTol) {
    double mx = 0;
    for (double v : ref.b) mx = std::max(mx, v);
    double worst = 0;
    int bad = 0;
    for (size_t i = 0; i < g.b.size(); ++i) {
        if (ref.b[i] < 0.05 * mx) continue;
        const double e = std::abs(g.b[i] - ref.b[i]) / ref.b[i];
        worst = std::max(worst, e);
        if (e > blockTol) ++bad;
    }
    for (int by = 0; by < 3; ++by) {
        std::cout << "[LvcManifold]     ";
        for (int bx = 0; bx < 6; ++bx) std::cout << g.b[by * 6 + bx] / std::max(ref.b[by * 6 + bx], 1e-12) << " (" << ref.b[by * 6 + bx] << ")  ";
        std::cout << "\n";
    }
    const double ratio = g.total / std::max(ref.total, 1e-12);
    const bool ok = std::abs(ratio - 1.0) < meanTol && bad == 0;
    std::cout << "[LvcManifold]   " << what << ": mean G " << g.total << " (ratio " << ratio << "), worst bright block off by " << 100 * worst << " %"
        << (ok ? "   OK" : "   FAIL") << "\n";
    return ok;
}

VolumeBptParams MakeParams(const CliOptions& opt, int W, int H, int mask) {
    VolumeBptParams bp;
    bp.width = W; bp.height = H; bp.spp = opt.spp * 8; bp.lightPaths = static_cast<long long>(opt.spp) * 40000; bp.maxVertices = 6;
    bp.surfaceAlbedo = Real(0.55);   // the albedo of the floor of the prism scene
    bp.useNee = (mask & 1) != 0; bp.useConnections = (mask & 2) != 0; bp.useLightTracing = (mask & 4) != 0;
    bp.mergeRadius = (mask & 8) ? Real(0.003) : Real(0); bp.useMerging = (mask & 8) != 0;
    bp.useManifoldNee = (mask & 16) != 0; bp.useManifoldConn = (mask & 32) != 0;
    return bp;
}

} // namespace

int RunLvcManifoldTests(const CliOptions& optIn) {
    CliOptions opt = optIn;
    opt.spp = std::min(opt.spp, 32);   // the merging of a dense caustic grows faster than the sample count (2.5 minutes at 32; the cases are capped)
    const int W = 48, H = 27;
    bool allPassed = true;

    // ---- U. the API: SpecularNee::chains against SpecularNee::sample (the glue: frame, Fresnel, d omega / d area) ----
    {
        LvcPrismScene ps = BuildLvcPrismScene(opt, W, H, false);
        SpecularNee nee(ps.bodies, *ps.scene);
        auto floorMat = std::make_shared<Lambertian>(RgbColor(0.55));
        SurfaceInteraction rec;
        rec.n = rec.gn = Vector3(0, 1, 0);
        rec.matPtr = floorMat.get();
        const Vector3 wo = glm::normalize(Vector3(0.3, 1.0, 0.5));
        SampledWavelengths lam;
        std::vector<SpecularNee::Chain> ch;
        int compared = 0, bad = 0;
        double worstRel = 0;
        for (int iz = -20; iz <= 5; ++iz) for (int ix = 0; ix <= 27; ++ix) {
            rec.p = Point3(0.05 + 0.01 * ix, 0, 0.01 * iz);
            const long long before = nee.solutions();
            const Spectrum S = nee.sample(*ps.scene, rec, wo, lam);
            const long long nSol = nee.solutions() - before;
            nee.chains(*ps.scene, rec.p, rec.gn, false, ps.lamp->position(), Real(550), ch);
            if (nSol == 0 || static_cast<long long>(ch.size()) != nSol) continue;   // the chains that the outside segments block are not in sample()
            Spectrum C(0.0);
            const Spectrum I = ps.lamp->intensity(lam);
            for (const SpecularNee::Chain& c : ch) {
                const Vector3 wi = -glm::normalize(Vector3(rec.p - c.vA));
                C += rec.matPtr->eval(rec, wo, wi, lam) * (I * (c.fresnel * c.dwdA));
            }
            const double a = SpectrumToRGB(S, lam).y, b = SpectrumToRGB(C, lam).y;
            const double rel = std::abs(a - b) / std::max(std::abs(a), 1e-30);
            worstRel = std::max(worstRel, rel);
            ++compared;
            if (rel > 1e-6) ++bad;
        }
        const bool ok = compared >= 30 && bad == 0;
        std::cout << "[LvcManifold] chains() against sample(): " << compared << " floor points compared, worst relative difference " << worstRel
            << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = allPassed && ok;
    }

    LvcPrismScene ps = BuildLvcPrismScene(opt, W, H, false);
    SpecularNee nee(ps.bodies, *ps.scene);
    CausticPhotonParams P;
    P.numPhotons = static_cast<long long>(opt.spp) * 60000;
    P.radius = Real(0.003);
    P.targetCenter = Point3(0, 0.035, 0);
    P.targetRadius = Real(0.06);
    P.lightTraceCamera = ps.cam.get();
    P.traceWidth = W;
    P.traceHeight = H;
    CausticPhotonParams PTT = P;
    PTT.onlyTransmissions = 2;
    CausticPhotonMap mapTT, mapAll;
    mapTT.build(*ps.scene, PTT);
    mapAll.build(*ps.scene, P);

    // ---- A. the manifold NEE alone against the photon map of the T T chains: a camera path of ONE vertex (the floor), the light seen across the prism ----
    {
        const BlockStats ref = Blocks(mapTT.splat().data(), W, H);
        VolumeBptParams bp = MakeParams(opt, W, H, 16);
        bp.lightPaths = 1000; bp.maxVertices = 1;
        bp.manifoldOnlyBlocked = false;   // the photon map counts the caustic where the straight segment is open as well
        Film f(W, H);
        std::cout << "[LvcManifold] A. manifold next event estimation only (camera path of one vertex) against the photon map of the T T chains: reference mean G "
            << ref.total << "\n";
        VolumeBPT v(*ps.scene, *ps.cam, *ps.lamp, *ps.fog, bp);
        v.setManifold(&nee);
        v.render(f);
        // the tolerances of the manifold check itself (ManifoldChecks: the sum within 10 % of the photon map, which blurs the edges of the pattern)
        allPassed = CompareBlocks("manifold NEE", Blocks(f.getData(), W, H), ref, 0.10, 0.30) && allPassed;
    }

    // ---- B. adding the manifold strategies must not change the expectation: the whole picture, against path tracing + the photon map ----
    {
        Film refF(W, H);
        PathIntegrator pt(ps.cam, 32, opt.spp * 64);
        pt.setCausticMap(&mapAll);
        pt.setSkipDirectCaustics(true);   // the light-traced image counts the caustics on surfaces the camera sees directly
        pt.render(*ps.scene, refF);
        std::vector<Vector3> refImg(static_cast<size_t>(W) * H);
        for (size_t i = 0; i < refImg.size(); ++i) refImg[i] = refF.getData()[i] + mapAll.splat()[i];
        const BlockStats ref = Blocks(refImg.data(), W, H);
        std::cout << "[LvcManifold] B. whole picture, reference (path tracing + photon map): mean G " << ref.total << "\n";
        // The ordinary strategies are 3 to 4 % below the reference everywhere (the same in every case: not the business of the manifold strategies), and
        // the picture behind the prism needs a strategy for the specular-diffuse-specular path: so the manifold strategies are judged against the
        // ordinary ones (adding a strategy must not change the expectation), and the whole against the reference.
        struct Case { const char* name; int mask; };
        BlockStats base;
        for (const Case& c : { Case{ "connections + light tracing (7)", 7 }, Case{ "7 + manifold NEE", 7 | 16 },
                               Case{ "7 + manifold NEE + manifold connections", 7 | 48 }, Case{ "all with merging and manifold (63)", 63 } }) {
            VolumeBptParams bp = MakeParams(opt, W, H, c.mask);
            Film f(W, H);
            std::cout << "[LvcManifold] " << c.name << "\n";
            VolumeBPT v(*ps.scene, *ps.cam, *ps.lamp, *ps.fog, bp);
            if (c.mask & 48) v.setManifold(&nee);
            v.render(f);
            const BlockStats g = Blocks(f.getData(), W, H);
            if (c.mask == 7) { base = g; std::cout << "[LvcManifold]   (ratio to the reference " << g.total / ref.total << ")\n"; continue; }
            if (c.mask == 63) allPassed = CompareBlocks(c.name, g, ref, 0.06, 0.20) && allPassed;
            else allPassed = CompareBlocks(c.name, g, base, 0.015, 0.12) && allPassed;
        }
        // The point of the strategy: less noise. The noise is read off two independent renders of the same configuration: sigma = rms(a - b) / sqrt(2),
        // relative to the mean (all pixels, so the dark parts count as well).
        auto noise = [&](int mask) {
            std::vector<Vector3> img[2];
            for (int k = 0; k < 2; ++k) {
                VolumeBptParams bp = MakeParams(opt, W, H, mask);
                Film f(W, H);
                VolumeBPT v(*ps.scene, *ps.cam, *ps.lamp, *ps.fog, bp);
                if (mask & 48) v.setManifold(&nee);
                v.render(f);
                img[k].assign(f.getData(), f.getData() + static_cast<size_t>(W) * H);
            }
            double s2 = 0, m = 0;
            for (size_t i = 0; i < img[0].size(); ++i) { const double d = img[0][i].y - img[1][i].y; s2 += d * d; m += 0.5 * (img[0][i].y + img[1][i].y); }
            const double n = double(img[0].size());
            return std::sqrt(s2 / (2.0 * n)) / (m / n);
        };
        const double n7 = noise(7), n55 = noise(7 | 48);
        std::cout << "[LvcManifold]   relative noise of the picture: connections + light tracing " << n7 << ", with the manifold strategies " << n55
            << " (ratio " << n55 / n7 << ")\n";
    }

    // ---- C. the camera looks THROUGH the prism: a chain on the camera side, and on the light side as well (E S S D S S L) ----
    {
        LvcPrismScene pc = BuildLvcPrismScene(opt, W, H, true);
        SpecularNee neeC(pc.bodies, *pc.scene);
        CausticPhotonParams Pc = P;
        Pc.lightTraceCamera = pc.cam.get();
        CausticPhotonMap mapC;
        mapC.build(*pc.scene, Pc);
        Film refF(W, H);
        PathIntegrator pt(pc.cam, 32, opt.spp * 64);
        pt.setCausticMap(&mapC);
        pt.setSkipDirectCaustics(true);
        pt.render(*pc.scene, refF);
        std::vector<Vector3> refImg(static_cast<size_t>(W) * H);
        for (size_t i = 0; i < refImg.size(); ++i) refImg[i] = refF.getData()[i] + mapC.splat()[i];
        const BlockStats ref = Blocks(refImg.data(), W, H);
        std::cout << "[LvcManifold] C. looking through the prism, reference (path tracing + photon map): mean G " << ref.total << "\n";
        struct Case { const char* name; int mask; };
        for (const Case& c : { Case{ "connections + light tracing (7), no strategy for the path through the prism", 7 }, Case{ "7 + merging (15)", 15 },
                               Case{ "7 + manifold (55)", 55 }, Case{ "all (63)", 63 } }) {
            VolumeBptParams bp = MakeParams(opt, W, H, c.mask);
            Film f(W, H);
            std::cout << "[LvcManifold] " << c.name << "\n";
            VolumeBPT v(*pc.scene, *pc.cam, *pc.lamp, *pc.fog, bp);
            if (c.mask & 48) v.setManifold(&neeC);
            v.render(f);
            const bool ok = CompareBlocks(c.name, Blocks(f.getData(), W, H), ref, 0.08, 0.35);
            if (c.mask >= 15) allPassed = ok && allPassed;   // the first case is allowed to be dark behind the prism
        }
    }
    std::cout << "[LvcManifold] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
