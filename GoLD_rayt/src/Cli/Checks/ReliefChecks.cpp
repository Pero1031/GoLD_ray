// Checks of the skin micro-relief (--scene relief): the height field, the bent normals, and the energy of a skin sphere in a white furnace.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/EngineHeaders.hpp"
#include "Materials/SkinRelief.hpp"

namespace rayt::cli {

int RunReliefTests(const CliOptions& opt) {
    bool allPassed = true;
    const int n = std::max(opt.spp, 64) * 200;   // points per estimate (default --spp 256: about 51k)
    std::cout << "[Relief] " << n << " points per estimate\n";

    SkinRelief::Params rp;   // the defaults: cell 0.6 mm, furrows 40 um, pores 80 um
    const SkinRelief relief(rp);

    // ---- 1. The height field: zero or negative, deepest in a pore on a furrow, the same value every time ---------------------------
    {
        double hMin = 0, hMax = -1e9, hSum = 0;
        for (int i = 0; i < n; ++i) {
            const Point3 p(sampling::Random() * 0.01, sampling::Random() * 0.01, sampling::Random() * 0.01);
            const Real h = relief.height(p);
            hMin = std::min<double>(hMin, h); hMax = std::max<double>(hMax, h); hSum += h;
        }
        const double mean = hSum / n;
        const Point3 q(0.00123, 0.00456, 0.00789);
        const bool ok = hMax <= 1e-15 && hMin >= -(rp.furrowDepth * (1 + rp.coarseWeight) + rp.poreDepth) * 1.0001 && hMin < -0.5 * rp.poreDepth
            && mean < 0 && mean > -0.5 * (rp.furrowDepth * (1 + rp.coarseWeight) + rp.poreDepth) && relief.height(q) == relief.height(q);
        std::cout << "[Relief] 1. height field over 1 cm^3: min " << hMin * 1e6 << " um, max " << hMax * 1e6 << " um, mean " << mean * 1e6 << " um"
            << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    // ---- 2. The bent normal: unit length, and on average the smooth normal (the gradient of a bounded field averages to 0) ----------------
    {
        const Vector3 nrm(0, 0, 1), wo(0, 0, 1);
        double sx = 0, sy = 0, sz = 0, qx = 0, qy = 0, tilt = 0, maxTilt = 0;
        bool unit = true;
        for (int i = 0; i < n; ++i) {
            const Point3 p(sampling::Random() * 0.01, sampling::Random() * 0.01, 0.0);
            const Vector3 b = relief.perturb(p, nrm, wo);
            unit = unit && std::abs(glm::length(b) - 1.0) < 1e-9;
            sx += b.x; sy += b.y; sz += b.z; qx += b.x * b.x; qy += b.y * b.y;
            const double t = std::acos(std::clamp<double>(b.z, -1, 1)) * 180.0 / constants::PI;
            tilt += t; maxTilt = std::max(maxTilt, t);
        }
        const double mx = sx / n, my = sy / n;
        const double sex = std::sqrt(qx / n / n), sey = std::sqrt(qy / n / n);
        const bool ok = unit && std::abs(mx) <= 4 * sex + 1e-3 && std::abs(my) <= 4 * sey + 1e-3 && tilt / n > 2.0 && maxTilt < 85.0;
        std::cout << "[Relief] 2. normals on a plane: mean (" << mx << ", " << my << ", " << sz / n << ") +- (" << sex << ", " << sey << "), mean tilt " << tilt / n
            << " deg, largest " << maxTilt << " deg" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    // ---- 3. The same normal from two step sizes of the finite difference: the field is smooth enough to differentiate -------------------
    {
        const Real eps = Real(0.04) * std::min(rp.poreRadius, rp.furrowWidth * rp.cell);   // as in SkinRelief::perturb
        int bad = 0, total = 0;
        for (int i = 0; i < 2000; ++i) {
            const Point3 p(sampling::Random() * 0.005, sampling::Random() * 0.005, 0.0);
            const Real h2 = relief.height(p + Vector3(eps, 0, 0)) - relief.height(p - Vector3(eps, 0, 0));
            const Real h1 = relief.height(p + Vector3(eps / 2, 0, 0)) - relief.height(p - Vector3(eps / 2, 0, 0));
            const double d2 = h2 / (2 * eps), d1 = h1 / eps;
            ++total;
            if (std::abs(d1 - d2) > 0.02 + 0.1 * std::abs(d1)) ++bad;   // near a cell boundary (a kink of F1) the two may differ
        }
        const bool ok = bad <= total / 20;
        std::cout << "[Relief] 3. finite-difference slopes with eps and eps / 2 agree at " << total - bad << " of " << total << " points" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    // ---- 4. Energy: a skin sphere with the relief returns what the smooth one does in a white furnace ----------------------------------
    {
        auto whiteSky = std::make_shared<EnvLight>(std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0)))), 1.0);
        auto build = [&](bool withRelief) {
            auto skin = std::make_shared<CoatedDiffuse>(RgbColor(0.80, 0.52, 0.42), Real(0.5), Real(0.4), Real(0));
            if (withRelief) skin->setRelief(std::make_shared<SkinRelief>(rp));
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Sphere>(Point3(0.0), Real(0.03), skin));
            auto scene = std::make_shared<Scene>(world);
            scene->setEnvLight(whiteSky);
            scene->finalize("spatial");
            return scene;
        };
        const auto smooth = build(false), bumped = build(true);
        const int rays = 200, perRay = std::max(opt.spp, 64) * 2;
        Vector3 sumS(0.0), sumB(0.0);
        for (int i = 0; i < rays; ++i) {   // rays parallel to z, spread over the middle of the sphere
            const Real r = Real(0.02) * std::sqrt(sampling::Random()), a = Real(2) * constants::PI * sampling::Random();
            const Ray ray(Point3(r * std::cos(a), r * std::sin(a), 0.2), Vector3(0, 0, -1), Real(0), nullptr);
            Vector3 mS, sS, mB, sB; int badS = 0, badB = 0;
            EstimateRadiance(*smooth, ray, perRay, 8, mS, sS, badS);
            EstimateRadiance(*bumped, ray, perRay, 8, mB, sB, badB);
            sumS += mS; sumB += mB;
        }
        const Vector3 ms = sumS / Real(rays), mb = sumB / Real(rays);
        bool ok = true;
        for (int c = 0; c < 3; ++c) ok = ok && std::abs(mb[c] - ms[c]) <= 0.03 * ms[c];
        std::cout << "[Relief] 4. skin sphere in a white furnace, mean over " << rays << " rays: smooth (" << ms.x << ", " << ms.y << ", " << ms.z
            << "), with relief (" << mb.x << ", " << mb.y << ", " << mb.z << ")" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    std::cout << "[Relief] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
