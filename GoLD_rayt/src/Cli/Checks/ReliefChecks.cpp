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

    // ---- 5. The pigment field: a factor in (0, 1], the same value every time, spots cover a part of the skin, red is absorbed least -----------
    const SkinPigment pigment{ SkinPigment::Params() };
    {
        double lo = 1.0, sum[3] = { 0, 0, 0 };
        long spotted = 0;
        bool inRange = true;
        for (int i = 0; i < n; ++i) {
            const Point3 p(sampling::Random() * 0.04, sampling::Random() * 0.04, sampling::Random() * 0.04);
            const Vector3 f = pigment.factor(p);
            for (int c = 0; c < 3; ++c) { inRange = inRange && f[c] > 0 && f[c] <= 1.0 + 1e-12; sum[c] += f[c]; lo = std::min<double>(lo, f[c]); }
            if (f.z < 0.9) ++spotted;
        }
        const Point3 q(0.0123, 0.0045, 0.0078);
        const double cover = double(spotted) / n;
        const bool ok = inRange && pigment.factor(q) == pigment.factor(q) && cover > 0.005 && cover < 0.3 && sum[0] > sum[1] && sum[0] > sum[2] && lo < 0.75;
        std::cout << "[Relief] 5. pigment over 4 cm^3: mean factor (" << sum[0] / n << ", " << sum[1] / n << ", " << sum[2] / n << "), darkest " << lo
            << ", " << 100.0 * cover << " % of the points darker than 0.9 in blue" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    // ---- 6. A skin sphere with pigment (and with subsurface scattering) in a white furnace: darker than without it, never brighter -------------
    {
        auto whiteSky = std::make_shared<EnvLight>(std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0)))), 1.0);
        auto profile = std::make_shared<SubsurfaceProfile>(Vector3(3.67, 1.37, 0.68) * Real(0.001), Real(1.4));
        auto build = [&](bool sss, bool withPigment) {
            auto skin = std::make_shared<CoatedDiffuse>(RgbColor(0.80, 0.52, 0.42), Real(0.5), Real(0.4), Real(0));
            if (sss) skin->setSubsurface(profile);
            if (withPigment) skin->setPigment(std::make_shared<SkinPigment>(SkinPigment::Params()), Real(0));
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Sphere>(Point3(0.0), Real(0.03), skin));
            auto scene = std::make_shared<Scene>(world);
            scene->setEnvLight(whiteSky);
            scene->finalize("spatial");
            return scene;
        };
        for (bool sss : { false, true }) {
            const auto plain = build(sss, false), spotted = build(sss, true);
            const int rays = 200, perRay = std::max(opt.spp, 64) * 2;
            Vector3 sumP(0.0), sumS(0.0);
            int bad = 0;
            for (int i = 0; i < rays; ++i) {
                const Real r = Real(0.02) * std::sqrt(sampling::Random()), a = Real(2) * constants::PI * sampling::Random();
                const Ray ray(Point3(r * std::cos(a), r * std::sin(a), 0.2), Vector3(0, 0, -1), Real(0), nullptr);
                Vector3 mP, sP, mS, sS; int bP = 0, bS = 0;
                EstimateRadiance(*plain, ray, perRay, 8, mP, sP, bP);
                EstimateRadiance(*spotted, ray, perRay, 8, mS, sS, bS);
                sumP += mP; sumS += mS; bad += bP + bS;
            }
            const Vector3 mp = sumP / Real(rays), ms = sumS / Real(rays);
            bool ok = bad == 0;
            for (int c = 0; c < 3; ++c) ok = ok && ms[c] <= mp[c] * 1.01 && ms[c] >= mp[c] * 0.75;
            std::cout << "[Relief] 6. furnace" << (sss ? " with subsurface scattering" : "") << ": uniform skin (" << mp.x << ", " << mp.y << ", " << mp.z << "), with pigment ("
                << ms.x << ", " << ms.y << ", " << ms.z << ")" << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }
    }

    // ---- 7. Pores: the mask is in [0, 1], pores differ in size, and their density follows the noise ---------------------------------------
    {
        double top = 0, covered = 0;
        bool inRange = true;
        std::vector<double> block(8, 0.0);   // pore coverage in eight 2 cm blocks along x
        for (int i = 0; i < n; ++i) {
            const Point3 p(sampling::Random() * 0.16, sampling::Random() * 0.01, sampling::Random() * 0.01);
            const Real m = relief.poreMask(p);
            inRange = inRange && m >= 0 && m <= 1.0 + 1e-12;
            top = std::max<double>(top, m);
            if (m > 0.5) { covered += 1; block[std::min(7, int(p.x / 0.02))] += 1; }
        }
        const double lo = *std::min_element(block.begin(), block.end()), hi = *std::max_element(block.begin(), block.end());
        const bool ok = inRange && top > 0.8 && covered > 0 && hi > 1.3 * std::max(lo, 1.0);
        std::cout << "[Relief] 7. pores: mask in [0, 1], largest " << top << ", pore coverage per 2 cm block from " << lo << " to " << hi << " of " << n / 8 << " points" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    // ---- 8. The wrinkle graph simulation (Materials/WrinkleField.hpp): a tile that repeats, the same every time, with wrinkles along the primary direction ---------
    {
        WrinkleField::Params wp;
        wp.tileMm = Real(12); wp.resolution = 512; wp.steps = 120; wp.seed = 3;
        WrinkleField::Level lv;
        lv.spacingMm = Real(0.8); lv.angleDeg = Real(0); lv.base = Real(0.05); lv.continuation = Real(0.6); lv.similarity = Real(0.5); lv.crossing = Real(0.7);
        lv.deposit = Real(0.35); lv.depthUm = Real(30); lv.width = Real(0.15);
        wp.levels = { lv };
        const WrinkleField a(wp), b(wp);
        wp.seed = 4;
        const WrinkleField c(wp);
        const auto& H = a.heightMap();
        const int R = a.resolution();
        const bool same = a.heightMap() == b.heightMap() && a.heightMap() != c.heightMap();
        double lo = 0, hi = 0, mean = 0;
        for (float v : H) { lo = std::min<double>(lo, v); hi = std::max<double>(hi, v); mean += v; }
        mean /= double(H.size());
        // orientation: the energy of the slope across the wrinkles (d/dy, along x they run) against along them (d/dx)
        double ex = 0, ey = 0;
        for (int y = 0; y < R; ++y)
            for (int x = 0; x < R; ++x) {
                const double dx = H[size_t(y) * R + (x + 1) % R] - H[size_t(y) * R + x], dy = H[size_t((y + 1) % R) * R + x] - H[size_t(y) * R + x];
                ex += dx * dx; ey += dy * dy;
            }
        // periodic: a point and the point one tile away have the same height (looked up through the world position)
        const Vector3 up(0, 1, 0);
        bool periodic = true;
        for (int i = 0; i < 200 && periodic; ++i) {
            const Point3 p(sampling::Random() * 0.012, 0.0, sampling::Random() * 0.012);
            periodic = std::abs(a.height(p, up) - a.height(p + Vector3(0.012, 0, 0), up)) < 1e-12 && std::abs(a.height(p, up) - a.height(p + Vector3(0, 0, 0.012), up)) < 1e-12;
        }
        bool poresOk = true;
        for (float v : a.poreMap()) poresOk = poresOk && v >= 0 && v <= 1.0001f;
        const bool ok = same && periodic && poresOk && lo >= -30e-6 * 1.001 && lo < -10e-6 && hi < 5e-6 && mean < 0 && ey > 1.4 * ex;
        std::cout << "[Relief] 8. wrinkle graph (one level, primary direction along x): deterministic " << (same ? "yes" : "NO") << ", periodic " << (periodic ? "yes" : "NO")
            << ", heights " << lo * 1e6 << " .. " << hi * 1e6 << " um, mean " << mean * 1e6 << " um, slope energy across / along the wrinkles " << ey / std::max(ex, 1e-30)
            << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;

        // the bent normals of a relief that uses it: unit length, and on average the smooth normal
        SkinRelief::Params rp0; rp0.furrowDepth = 0; rp0.coarseWeight = 0; rp0.poreDepth = 0;
        SkinRelief wr(rp0);
        wr.setWrinkles(std::make_shared<WrinkleField>(wp));
        const Vector3 nrm(0, 1, 0), wo(0, 1, 0);
        double sx = 0, sz = 0, qx = 0, qz = 0, tilt = 0;
        bool unit = true;
        for (int i = 0; i < n; ++i) {
            const Point3 p(sampling::Random() * 0.012, 0.0, sampling::Random() * 0.012);
            const Vector3 bn = wr.perturb(p, nrm, wo);
            unit = unit && std::abs(glm::length(bn) - 1.0) < 1e-9;
            sx += bn.x; sz += bn.z; qx += bn.x * bn.x; qz += bn.z * bn.z;
            tilt += std::acos(std::clamp<double>(bn.y, -1, 1)) * 180.0 / constants::PI;
        }
        const double mx = sx / n, mz = sz / n, sex = std::sqrt(qx / n / n), sez = std::sqrt(qz / n / n);
        const bool ok2 = unit && std::abs(mx) <= 4 * sex + 1e-3 && std::abs(mz) <= 4 * sez + 1e-3 && tilt / n > 0.3;
        std::cout << "[Relief] 8b. normals from the wrinkle field: mean (" << mx << ", " << mz << ") +- (" << sex << ", " << sez << "), mean tilt " << tilt / n << " deg" << (ok2 ? "   OK" : "   FAIL") << "\n";
        allPassed = ok2 && allPassed;
    }

    // ---- 9. (information) How strongly the two kinds of relief bend the normal, with the settings of scenes/skin_macro.json ------------------------------
    {
        SkinRelief::Params worley;   // the Worley version of the macro scene
        worley.cell = Real(0.0005); worley.furrowDepth = Real(0.000025); worley.furrowWidth = Real(0.1); worley.stretch = Real(1.5); worley.coarseScale = Real(3);
        worley.coarseWeight = Real(0.7); worley.poreDepth = Real(0.0001); worley.poreRadius = Real(0.00006); worley.poreFraction = Real(0.3);
        const SkinRelief a(worley);
        WrinkleField::Params wp; wp.depthScale = Real(1.6); wp.steps = 300;
        SkinRelief::Params none; none.furrowDepth = 0; none.coarseWeight = 0; none.poreDepth = 0;
        SkinRelief b(none);
        b.setWrinkles(std::make_shared<WrinkleField>(wp));
        const Vector3 nrm(0, 1, 0), wo(0, 1, 0);
        for (int k = 0; k < 2; ++k) {
            const SkinRelief& rel = k == 0 ? a : b;
            std::vector<double> tilts;
            tilts.reserve(size_t(n));
            for (int i = 0; i < n; ++i) {
                const Point3 p(sampling::Random() * 0.02, 0.0, sampling::Random() * 0.02);
                tilts.push_back(std::acos(std::clamp<double>(rel.perturb(p, nrm, wo).y, -1, 1)) * 180.0 / constants::PI);
            }
            std::sort(tilts.begin(), tilts.end());
            double mean = 0; for (double t : tilts) mean += t;
            std::cout << "[Relief] 9. (information) " << (k == 0 ? "Worley furrows and pores" : "wrinkle graph") << ": normal tilt mean " << mean / n << " deg, median "
                << tilts[tilts.size() / 2] << ", 90th percentile " << tilts[tilts.size() * 9 / 10] << ", share above 10 deg " << 100.0 * double(tilts.end() - std::upper_bound(tilts.begin(), tilts.end(), 10.0)) / n << " %\n";
        }
    }

    std::cout << "[Relief] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
