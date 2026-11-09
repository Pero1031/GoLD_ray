// Checks of the skin micro-relief (--scene relief): the height field, the bent normals, and the energy of a skin sphere in a white furnace.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/EngineHeaders.hpp"
#include "Materials/SkinRelief.hpp"
#include "Materials/SkinFilm.hpp"
#include "Materials/SkinScatter.hpp"
#include "Materials/SkinLotion.hpp"

namespace rayt::cli {

namespace {

#if defined(RAYT_SPECTRAL_MODE)
constexpr int kN = kSpectralSamples;
#else
constexpr int kN = 3;
#endif

Real average(const Spectrum& s) {
    Real m = 0;
    for (int i = 0; i < kN; ++i) m += s[i];
    return m / Real(kN);
}

struct Estimate {
    double mean = 0, se = 0;
    void add(double v) { ++n; sum += v; sum2 += v * v; }
    void finish() {
        mean = sum / double(n);
        se = std::sqrt(std::max(sum2 / double(n) - mean * mean, 0.0) / double(n));
    }
    long n = 0;
    double sum = 0, sum2 = 0;
};

} // namespace

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

    // ---- 10. A liquid film that pools in the furrows (Materials/SkinFilm.hpp, docs/SKIN_MOISTURE.md) -----------------------------------------
    {
        auto reliefP = std::make_shared<SkinRelief>(rp);

        // (a) The level that fills the furrows gives the mean thickness that was asked for.
        for (double um : { 1.0, 5.0, 20.0 }) {
            SkinFilm::Params fp; fp.thickness = Real(um * 1e-6);
            const SkinFilm film(fp, reliefP);
            const double mean = double(film.meanThickness(300000, 99u)) * 1e6;
            const bool ok = std::abs(mean - um) <= 0.03 * um;
            std::cout << "[Relief] 10a. film of " << um << " um on the Worley relief: mean thickness " << mean << " um (level " << double(film.level()) * 1e6 << " um)"
                << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }

        // (b) The top of the film is flatter the more liquid there is: thin = the skin's normal, thick = almost flat.
        {
            const Vector3 nrm(0, 0, 1), wo(0, 0, 1);
            const int m = std::min(n, 20000);
            bool unit = true;
            auto meanTilt = [&](const std::function<Vector3(const Point3&)>& f) {
                double t = 0;
                for (int i = 0; i < m; ++i) {
                    const Vector3 b = f(Point3(sampling::Random() * 0.01, sampling::Random() * 0.01, 0.0));
                    unit = unit && std::abs(glm::length(b) - 1.0) < 1e-9;
                    t += std::acos(std::clamp<double>(b.z, -1, 1)) * 180.0 / constants::PI;
                }
                return t / m;
            };
            const double skinTilt = meanTilt([&](const Point3& p) { return reliefP->perturb(p, nrm, wo); });
            double prev = 1e9, first = 0, last = 0;
            bool ordered = true;
            std::cout << "[Relief] 10b. mean tilt of the normal: skin " << skinTilt << " deg; film";
            for (double um : { 0.5, 5.0, 40.0, 200.0 }) {
                SkinFilm::Params fp; fp.thickness = Real(um * 1e-6);
                const SkinFilm film(fp, reliefP);
                const double t = meanTilt([&](const Point3& p) { return film.at(p, nrm, wo).normal; });
                ordered = ordered && t <= prev * 1.05;
                prev = t;
                if (um == 0.5) first = t;
                last = t;
                std::cout << " " << um << " um: " << t << " deg,";
            }
            const bool ok = unit && ordered && first > 0.5 * skinTilt && last < 0.15 * skinTilt;
            std::cout << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }

        // (c) The reflectance of the film on skin: water 0.0257 at normal incidence, the skin's own 0.0467 when the liquid has the skin's index, rising to 1 at grazing.
        {
            SkinFilm::Params w; w.thickness = Real(5e-6);
            const SkinFilm water(w, reliefP);
            SkinFilm::Params s = w; s.ior = Real(1.55);
            const SkinFilm same(s, reliefP);
            const double fw = water.fresnelReflectance(Real(1)), fs = same.fresnelReflectance(Real(1)), fg = water.fresnelReflectance(Real(0.01));
            bool rising = true;
            double prev = water.fresnelReflectance(Real(1));
            for (int i = 1; i <= 50; ++i) { const double c = 1.0 - 0.0198 * i; const double f = water.fresnelReflectance(Real(c)); rising = rising && f >= prev - 1e-9; prev = f; }
            const bool ok = std::abs(fw - 0.0257) < 0.0005 && std::abs(fs - 0.0467) < 0.0005 && fg > 0.8 && rising;
            std::cout << "[Relief] 10c. film reflectance at normal incidence: water " << fw << ", index 1.55 " << fs << ", water at 89.4 deg " << fg << ", rising with angle " << (rising ? "yes" : "NO")
                << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }

        // (d) eval() and sample() are two estimators of the same albedo (a film over the relief: the lobe about the film's normal, the diffuse one about the skin's).
        for (int variant = 0; variant < 4; ++variant) {   // without / with subsurface scattering, without / with the white scatter of flaking skin (11 below)
            const bool sss = (variant & 1) != 0, scat = (variant & 2) != 0;
            SkinFilm::Params fp; fp.thickness = Real(15e-6); fp.roughness = Real(0.3);
            auto skin = std::make_shared<CoatedDiffuse>(RgbColor(0.80, 0.52, 0.42), Real(0.3), Real(0.6), Real(0));
            skin->setRelief(reliefP);
            skin->setFilm(std::make_shared<SkinFilm>(fp, reliefP));
            if (scat) { SkinScatter::Params sp; sp.rho0 = Real(0.4); sp.threshold = Real(0.5); skin->setScatter(std::make_shared<SkinScatter>(sp)); skin->setEmulsion(Real(1.0)); }   // flakes and the emulsion layer together
            if (sss) skin->setSubsurface(std::make_shared<SubsurfaceProfile>(Vector3(3.67, 1.37, 0.68) * Real(0.001), Real(1.4)));
#if defined(RAYT_SPECTRAL_MODE)
            const SampledWavelengths lambda = SampledWavelengths::sampleVisible(Real(0.37));
#else
            const SampledWavelengths lambda;
#endif
            const double deg = 40.0 * constants::PI / 180.0;
            const Vector3 wo(std::sin(deg), 0, std::cos(deg));
            double sumS = 0, sumE = 0, varS = 0, varE = 0, cover = 0;
            const int points = 24, per = 6000;
            for (int k = 0; k < points; ++k) {
                SurfaceInteraction rec;
                rec.p = Point3(sampling::Random() * 0.005, sampling::Random() * 0.005, 0.0);
                rec.n = Vector3(0, 0, 1); rec.gn = Vector3(0, 0, 1); rec.wo = wo; rec.matPtr = skin.get();
                skin->perturbShading(rec, wo);
                cover += rec.filmCover;
                Estimate bySample, byEval;
                for (int i = 0; i < per; ++i) {
                    const auto bs = skin->sample(rec, wo, Point2(sampling::Random(), sampling::Random()), lambda, TransportMode::Radiance);
                    // the subsurface event is not a reflection: only the local lobes count here
                    bySample.add((bs && !bs->isSubsurface()) ? double(average(bs->f)) * std::abs(glm::dot(rec.n, bs->wi)) / bs->pdf : 0.0);
                    const Real r = std::sqrt(sampling::Random()), phi = Real(2) * constants::PI * sampling::Random();
                    const Vector3 wi(r * std::cos(phi), r * std::sin(phi), std::sqrt(std::max(Real(0), Real(1) - r * r)));
                    byEval.add(wi.z > Real(1e-3) ? double(average(skin->eval(rec, wo, wi, lambda, TransportMode::Radiance))) * std::abs(glm::dot(rec.n, wi)) * constants::PI / wi.z : 0.0);
                }
                bySample.finish(); byEval.finish();
                sumS += bySample.mean; sumE += byEval.mean; varS += bySample.se * bySample.se; varE += byEval.se * byEval.se;
            }
            const double mS = sumS / points, mE = sumE / points, se = std::sqrt(varS + varE) / points;
            const bool ok = std::abs(mS - mE) <= 4 * se + 0.02 * mE && mE < 1.0;
            std::cout << "[Relief] 10d. film over the relief" << (sss ? " with subsurface scattering" : "") << (scat ? " and the white scatter" : "") << ": albedo from sample() " << mS << ", from eval() " << mE << " +- " << se
                << " (mean film cover " << cover / points << ")" << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }

        // (e) A skin sphere with a film in a white furnace returns about what the dry one does (the film only changes the reflectance a little) and never more than 1.
        {
            auto whiteSky = std::make_shared<EnvLight>(std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0)))), 1.0);
            auto profile = std::make_shared<SubsurfaceProfile>(Vector3(3.67, 1.37, 0.68) * Real(0.001), Real(1.4));
            for (bool sss : { false, true }) {
                auto build = [&](bool withFilm) {
                    auto skin = std::make_shared<CoatedDiffuse>(RgbColor(0.80, 0.52, 0.42), Real(0.5), Real(0.4), Real(0));
                    skin->setRelief(reliefP);
                    if (sss) skin->setSubsurface(profile);
                    if (withFilm) { SkinFilm::Params fp; fp.thickness = Real(10e-6); skin->setFilm(std::make_shared<SkinFilm>(fp, reliefP)); }
                    auto world = std::make_shared<Aggregate>();
                    world->add(std::make_shared<Sphere>(Point3(0.0), Real(0.03), skin));
                    auto scene = std::make_shared<Scene>(world);
                    scene->setEnvLight(whiteSky);
                    scene->finalize("spatial");
                    return scene;
                };
                const auto dry = build(false), wet = build(true);
                const int rays = 200, perRay = std::max(opt.spp, 64) * 2;
                Vector3 sumD(0.0), sumW(0.0);
                int bad = 0;
                for (int i = 0; i < rays; ++i) {
                    const Real r = Real(0.02) * std::sqrt(sampling::Random()), a = Real(2) * constants::PI * sampling::Random();
                    const Ray ray(Point3(r * std::cos(a), r * std::sin(a), 0.2), Vector3(0, 0, -1), Real(0), nullptr);
                    Vector3 mD, sD, mW, sW; int bD = 0, bW = 0;
                    EstimateRadiance(*dry, ray, perRay, 8, mD, sD, bD);
                    EstimateRadiance(*wet, ray, perRay, 8, mW, sW, bW);
                    sumD += mD; sumW += mW; bad += bD + bW;
                }
                const Vector3 md = sumD / Real(rays), mw = sumW / Real(rays);
                bool ok = bad == 0;
                for (int c = 0; c < 3; ++c) ok = ok && mw[c] <= 1.0 && std::abs(mw[c] - md[c]) <= 0.08 * md[c];
                std::cout << "[Relief] 10e. furnace" << (sss ? " with subsurface scattering" : "") << ": dry (" << md.x << ", " << md.y << ", " << md.z << "), with a 10 um film ("
                    << mw.x << ", " << mw.y << ", " << mw.z << ")" << (ok ? "   OK" : "   FAIL") << "\n";
                allPassed = ok && allPassed;
            }
        }
    }

    // ---- 11. The white scatter of flaking skin (Materials/SkinScatter.hpp, docs/SKIN_MOISTURE.md 2.3) ------------------------------------------
    {
        // (a) The index-matching factor.
        {
            const double fa = SkinScatter::match(Real(1.55), Real(1)), fs = SkinScatter::match(Real(1.55), Real(1.55)), fw = SkinScatter::match(Real(1.55), Real(1.33)),
                fg = SkinScatter::match(Real(1.55), Real(1.47));
            const bool ok = fa == 1.0 && fs < 1e-12 && std::abs(fw - 0.125) < 0.002 && std::abs(fg - 0.0151) < 0.001;
            std::cout << "[Relief] 11a. index matching of the flakes (1.55) by a liquid: air " << fa << ", index 1.55 " << fs << ", water " << fw << ", glycerin " << fg
                << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }

        // (b) The field: between 0 and rho0, patchy, deterministic; a film and the hydration lower it by the factors above.
        {
            SkinScatter::Params sp; sp.rho0 = Real(0.1); sp.kappaH = Real(0.4);
            const SkinScatter dry(sp);
            SkinScatter::Params sh = sp; sh.hydration = Real(1);
            const SkinScatter hyd(sh);
            bool inRange = true;
            double sumDry = 0, sumWater = 0, sumHyd = 0, sumGly = 0;
            long flaked = 0;
            const int m = std::min(n, 40000);
            for (int i = 0; i < m; ++i) {
                const Point3 p(sampling::Random() * 0.1, sampling::Random() * 0.1, 0.0);
                const Real r = dry.reflectance(p, Real(1), Real(0));
                inRange = inRange && r >= 0 && r <= sp.rho0 * Real(1.0001);
                sumDry += r; sumHyd += hyd.reflectance(p, Real(1), Real(0));
                sumWater += dry.reflectance(p, Real(1.33), Real(10e-6)); sumGly += dry.reflectance(p, Real(1.47), Real(10e-6));
                if (r > 0.5 * sp.rho0) ++flaked;
            }
            const Point3 q(0.0123, 0.0045, 0.0);
            const double cover = double(flaked) / m;
            const bool ok = inRange && dry.reflectance(q, Real(1), Real(0)) == dry.reflectance(q, Real(1), Real(0)) && cover > 0.05 && cover < 0.7
                && std::abs(sumHyd / sumDry - 0.6) < 0.01 && std::abs(sumWater / sumDry - 0.125) < 0.01 && std::abs(sumGly / sumDry - 0.0151) < 0.005;
            std::cout << "[Relief] 11b. flakes: mean rho " << sumDry / m << " (rho0 " << sp.rho0 << "), " << 100.0 * cover << " % of the points above half of it; hydrated / dry "
                << sumHyd / sumDry << ", under 10 um of water / dry " << sumWater / sumDry << ", of glycerin " << sumGly / sumDry << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }

        // (c) A skin sphere in a white furnace: the white scatter makes it brighter and never above 1; a film cuts it back.
        // (d) The direction of a wet surface: under a water film the skin is darker and more saturated (CIELAB of the furnace colour; Sawayama et al. 2017 saw this in photographs).
        {
            auto whiteSky = std::make_shared<EnvLight>(std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0)))), 1.0);
            auto profile = std::make_shared<SubsurfaceProfile>(Vector3(3.67, 1.37, 0.68) * Real(0.001), Real(1.4));
            auto reliefP = std::make_shared<SkinRelief>(rp);
            auto lab = [](const Vector3& c) {   // linear sRGB (white = 1) -> L*, C*
                const double X = 0.4124 * c.x + 0.3576 * c.y + 0.1805 * c.z, Y = 0.2126 * c.x + 0.7152 * c.y + 0.0722 * c.z, Z = 0.0193 * c.x + 0.1192 * c.y + 0.9505 * c.z;
                auto f = [](double t) { return t > 216.0 / 24389.0 ? std::cbrt(t) : (24389.0 / 27.0 * t + 16.0) / 116.0; };
                const double fx = f(X / 0.9505), fy = f(Y), fz = f(Z / 1.089);
                return std::array<double, 2>{ 116.0 * fy - 16.0, std::hypot(500.0 * (fx - fy), 200.0 * (fy - fz)) };
            };
            for (bool sss : { false, true }) {
                auto build = [&](bool scatter, bool film) {
                    auto skin = std::make_shared<CoatedDiffuse>(RgbColor(0.80, 0.52, 0.42), Real(0.5), Real(0.4), Real(0));
                    skin->setRelief(reliefP);
                    if (sss) skin->setSubsurface(profile);
                    if (film) { SkinFilm::Params fp; fp.thickness = Real(10e-6); skin->setFilm(std::make_shared<SkinFilm>(fp, reliefP)); }
                    if (scatter) { SkinScatter::Params sp; sp.rho0 = Real(0.3); sp.threshold = Real(0.5); sp.softness = Real(0.05); skin->setScatter(std::make_shared<SkinScatter>(sp)); }
                    auto world = std::make_shared<Aggregate>();
                    world->add(std::make_shared<Sphere>(Point3(0.0), Real(0.03), skin));
                    auto scene = std::make_shared<Scene>(world);
                    scene->setEnvLight(whiteSky);
                    scene->finalize("spatial");
                    return scene;
                };
                const auto plain = build(false, false), flaky = build(true, false), wet = build(true, true), wetPlain = build(false, true);
                const int rays = 200, perRay = std::max(opt.spp, 64) * 2;
                Vector3 sP(0.0), sF(0.0), sW(0.0), sWP(0.0);
                int bad = 0;
                for (int i = 0; i < rays; ++i) {
                    const Real r = Real(0.02) * std::sqrt(sampling::Random()), a = Real(2) * constants::PI * sampling::Random();
                    const Ray ray(Point3(r * std::cos(a), r * std::sin(a), 0.2), Vector3(0, 0, -1), Real(0), nullptr);
                    Vector3 m1, s1, m2, s2, m3, s3, m4, s4; int b1 = 0, b2 = 0, b3 = 0, b4 = 0;
                    EstimateRadiance(*plain, ray, perRay, 8, m1, s1, b1);
                    EstimateRadiance(*flaky, ray, perRay, 8, m2, s2, b2);
                    EstimateRadiance(*wet, ray, perRay, 8, m3, s3, b3);
                    EstimateRadiance(*wetPlain, ray, perRay, 8, m4, s4, b4);
                    sP += m1; sF += m2; sW += m3; sWP += m4; bad += b1 + b2 + b3 + b4;
                }
                const Vector3 mp = sP / Real(rays), mf = sF / Real(rays), mw = sW / Real(rays), mwp = sWP / Real(rays);
                bool okC = bad == 0;
                for (int c = 0; c < 3; ++c) okC = okC && mf[c] >= mp[c] * 0.99 && mf[c] <= 1.0 && mw[c] <= mf[c] * 1.01;
                std::cout << "[Relief] 11c. furnace" << (sss ? " with subsurface scattering" : "") << ": no flakes (" << mp.x << ", " << mp.y << ", " << mp.z << "), flakes (" << mf.x << ", " << mf.y << ", " << mf.z
                    << "), flakes under a 10 um water film (" << mw.x << ", " << mw.y << ", " << mw.z << ")" << (okC ? "   OK" : "   FAIL") << "\n";
                allPassed = okC && allPassed;
                const auto d = lab(mf), w = lab(mw), wp = lab(mwp);
                const bool okD = w[0] < d[0] && w[1] > d[1];
                std::cout << "[Relief] 11d. wetting the flaking skin" << (sss ? " (with subsurface scattering)" : "") << ": L* " << d[0] << " -> " << w[0] << ", chroma C* " << d[1] << " -> " << w[1]
                    << " (a skin without flakes: L* " << wp[0] << ", C* " << wp[1] << " under the film)" << (okD ? "   OK (darker, more saturated)" : "   FAIL") << "\n";
                allPassed = okD && allPassed;
            }
        }
    }

    // ---- 12. The time course of a lotion (Materials/SkinLotion.hpp) ------------------------------------------------------------------------
    {
        const SkinLotion::Params lp;
        const bool valid = std::string(SkinLotion::validate(lp)).empty();
        SkinLotion::Params bad = lp; bad.hydrationDecay = bad.hydrationRise;
        const bool rejects = !std::string(SkinLotion::validate(bad)).empty();

        const SkinLotion::State s0 = SkinLotion::at(lp, Real(0)), sTau = SkinLotion::at(lp, lp.waterTau), sEnd = SkinLotion::at(lp, Real(36000));
        const double iorStart = (1.33 * 0.95 + 1.47 * 0.05);
        const bool start = std::abs(double(s0.thickness - lp.dose)) < 1e-12 && std::abs(double(s0.ior) - iorStart) < 1e-6 && s0.hydration == Real(0);
        const bool tauOk = std::abs(double(sTau.water) - double(lp.dose) * 0.95 / std::exp(1.0)) < 1e-12;
        const bool end = std::abs(double(sEnd.thickness) - double(lp.dose) * 0.05) < 1e-9 && std::abs(double(sEnd.ior) - 1.47) < 1e-6 && std::abs(double(sEnd.roughness) - 0.08) < 1e-6 && sEnd.hydration >= 0 && sEnd.hydration < 1e-6;
        bool monotone = true;
        double prevT = 1e9, prevN = 0, peak = 0, peakT = 0;
        for (int i = 0; i <= 3600; ++i) {
            const double t = double(i) * 2.0;   // 0 .. 2 h
            const SkinLotion::State s = SkinLotion::at(lp, Real(t));
            monotone = monotone && double(s.thickness) <= prevT + 1e-15 && double(s.ior) >= prevN - 1e-12 && s.hydration >= 0 && s.hydration <= lp.hydrationPeak * Real(1.0001);
            prevT = double(s.thickness); prevN = double(s.ior);
            if (s.hydration > peak) { peak = double(s.hydration); peakT = t; }
        }
        const double tp = double(SkinLotion::hydrationPeakTime(lp));
        const bool peakOk = std::abs(peak - double(lp.hydrationPeak)) < 0.005 * double(lp.hydrationPeak) && std::abs(peakT - tp) <= 2.0;
        const bool ok = valid && rejects && start && tauOk && end && monotone && peakOk;
        std::cout << "[Relief] 12. lotion over time (default: 10 um, 5 % humectant, water tau " << double(lp.waterTau) << " s): at 0 s film " << double(s0.thickness) * 1e6 << " um, index " << double(s0.ior)
            << "; at one tau " << double(sTau.water) * 1e6 << " um of water; at 10 h " << double(sEnd.thickness) * 1e6 << " um of humectant (index " << double(sEnd.ior) << "); hydration peaks at " << peak << " at "
            << peakT << " s (formula " << tp << " s); thickness falls and index rises with time: " << (monotone ? "yes" : "NO") << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    // ---- 13. An emulsion: oil droplets in the lotion (Materials/SkinLotion.hpp, docs/SKIN_MOISTURE.md 5.5) ------------------------------------------
    {
        // (a) The extinction efficiency of a droplet (anomalous diffraction): 0 for an index of 1, a first maximum of about 3.2 at rho = 2 x (m - 1) = 4.09, continuous at the series.
        {
            double best = 0, bestRho = 0;
            for (int i = 1; i <= 2000; ++i) {
                const double rho = 0.005 * i;   // x = rho / 2 for m = 2
                const double q = double(SkinLotion::dropletEfficiency(Real(2), Real(rho / 2)));
                if (q > best) { best = q; bestRho = rho; }
            }
            const double qa = double(SkinLotion::dropletEfficiency(Real(1.0245), Real(1))), qb = double(SkinLotion::dropletEfficiency(Real(1.0255), Real(1)));   // rho = 0.049 and 0.051: no step between the series and the formula
            const bool ok = SkinLotion::dropletEfficiency(Real(1), Real(7)) == Real(0) && best > 3.1 && best < 3.3 && std::abs(bestRho - 4.09) < 0.3
                && std::abs(qa / (0.049 * 0.049 / 2.0) - 1.0) < 0.01 && std::abs(qb / (0.051 * 0.051 / 2.0) - 1.0) < 0.01;   // Q = rho^2 / 2 for small rho, on both sides of the switch
            std::cout << "[Relief] 13a. droplet efficiency (anomalous diffraction): m = 1 -> 0; first maximum " << best << " at rho = " << bestRho << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }

        // (b) The time course of an emulsion: oil 8 % of 10 um stays, the droplets scatter (transport optical depth about 0.09) and the scatter fades as they merge and as the water dries.
        {
            SkinLotion::Params ep; ep.oil = Real(0.08);
            const bool valid = std::string(SkinLotion::validate(ep)).empty();
            SkinLotion::Params bad = ep; bad.oil = Real(0.95);
            const bool rejects = !std::string(SkinLotion::validate(bad)).empty();
            const SkinLotion::State s0 = SkinLotion::at(ep, Real(0)), sEnd = SkinLotion::at(ep, Real(36000));
            const SkinLotion::State plain = SkinLotion::at(SkinLotion::Params(), Real(0));
            // The scatter never exceeds its start: the water dries (the droplets' index gets close to the medium's, and goes through 1 when the medium is as dense as the oil, so the
            // efficiency need not fall monotonically after that, but it stays tiny) and the droplets merge.
            bool monotone = true;
            for (int i = 0; i <= 600; ++i) monotone = monotone && double(SkinLotion::at(ep, Real(10.0 * i)).emulsionTau) <= double(s0.emulsionTau) * 1.0000001;
            const bool ok = valid && rejects && plain.emulsionTau == Real(0) && std::abs(double(s0.oil) - 0.8e-6) < 1e-12 && s0.emulsionTau > 0.085 && s0.emulsionTau < 0.096
                && sEnd.emulsionTau < 1e-6 && std::abs(double(sEnd.thickness) - 1.3e-6) < 1e-9 && std::abs(double(sEnd.ior) - (1.45 * 0.8 + 1.47 * 0.5) / 1.3) < 1e-6 && monotone;
            std::cout << "[Relief] 13b. emulsion (8 % oil, 1 um droplets, 10 um): transport optical depth " << double(s0.emulsionTau) << " at 0 s, " << double(SkinLotion::at(ep, Real(600)).emulsionTau)
                << " at 600 s, " << double(sEnd.emulsionTau) << " at 10 h; the film at 10 h is " << double(sEnd.thickness) * 1e6 << " um (oil + humectant), index " << double(sEnd.ior)
                << "; a plain lotion has none; never above its start: " << (monotone ? "yes" : "NO") << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }

        // (c) A skin sphere in a white furnace with a film: the emulsion adds a white layer, which brightens the dark channels and never gives more than 1.
        {
            auto whiteSky = std::make_shared<EnvLight>(std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0)))), 1.0);
            auto profile = std::make_shared<SubsurfaceProfile>(Vector3(3.67, 1.37, 0.68) * Real(0.001), Real(1.4));
            auto reliefP = std::make_shared<SkinRelief>(rp);
            for (bool sss : { false, true }) {
                auto build = [&](bool emulsion) {
                    auto skin = std::make_shared<CoatedDiffuse>(RgbColor(0.80, 0.52, 0.42), Real(0.5), Real(0.4), Real(0));
                    skin->setRelief(reliefP);
                    if (sss) skin->setSubsurface(profile);
                    SkinFilm::Params fp; fp.thickness = Real(10e-6);
                    skin->setFilm(std::make_shared<SkinFilm>(fp, reliefP));
                    if (emulsion) skin->setEmulsion(Real(0.5));
                    auto world = std::make_shared<Aggregate>();
                    world->add(std::make_shared<Sphere>(Point3(0.0), Real(0.03), skin));
                    auto scene = std::make_shared<Scene>(world);
                    scene->setEnvLight(whiteSky);
                    scene->finalize("spatial");
                    return scene;
                };
                const auto clear = build(false), milky = build(true);
                const int rays = 200, perRay = std::max(opt.spp, 64) * 2;
                Vector3 sC(0.0), sM(0.0);
                int bad = 0;
                for (int i = 0; i < rays; ++i) {
                    const Real r = Real(0.02) * std::sqrt(sampling::Random()), a = Real(2) * constants::PI * sampling::Random();
                    const Ray ray(Point3(r * std::cos(a), r * std::sin(a), 0.2), Vector3(0, 0, -1), Real(0), nullptr);
                    Vector3 mC, eC, mM, eM; int bC = 0, bM = 0;
                    EstimateRadiance(*clear, ray, perRay, 8, mC, eC, bC);
                    EstimateRadiance(*milky, ray, perRay, 8, mM, eM, bM);
                    sC += mC; sM += mM; bad += bC + bM;
                }
                const Vector3 mc = sC / Real(rays), mm = sM / Real(rays);
                bool ok = bad == 0 && mm.z > mc.z * 1.02;
                for (int c = 0; c < 3; ++c) ok = ok && mm[c] >= mc[c] * 0.97 && mm[c] <= 1.0;   // the red channel (the brightest, 0.80) hardly changes: a white layer over it gives +1 % at most
                std::cout << "[Relief] 13c. furnace" << (sss ? " with subsurface scattering" : "") << ", a 10 um film: clear (" << mc.x << ", " << mc.y << ", " << mc.z << "), with the emulsion layer (tau 0.5) ("
                    << mm.x << ", " << mm.y << ", " << mm.z << ")" << (ok ? "   OK" : "   FAIL") << "\n";
                allPassed = ok && allPassed;
            }
        }
    }

    std::cout << "[Relief] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
