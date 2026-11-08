// Checks of hair rendering (--scene hair): the fiber scattering function (energy, eval vs sample vs pdf, where the lobes point),
// the cylinder-chain intersector against brute force, and a white furnace around a lock of non-absorbing hair.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/EngineHeaders.hpp"
#include "Geometry/HairCurves.hpp"
#include "Geometry/HairGenerator.hpp"
#include "Materials/Hair.hpp"

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

// Uniform direction on the sphere.
Vector3 uniformSphere() {
    const Real z = Real(1) - Real(2) * sampling::Random();
    const Real r = std::sqrt(std::max(Real(0), Real(1) - z * z)), phi = Real(2) * constants::PI * sampling::Random();
    return Vector3(r * std::cos(phi), r * std::sin(phi), z);
}

// Outgoing direction in the hair frame (x along the fiber, z = the surface normal): inclination theta along the fiber and
// azimuth phi in the cross-section measured from the y axis; phi in (0, pi) is the visible side.
Vector3 outgoing(double thetaDeg, double phiDeg) {
    const double t = thetaDeg * constants::PI / 180.0, p = phiDeg * constants::PI / 180.0;
    return Vector3(std::sin(t), std::cos(t) * std::cos(p), std::cos(t) * std::sin(p));
}

} // namespace

int RunHairTests(const CliOptions& opt) {
    bool allPassed = true;
    const int n = std::max(opt.spp, 64) * 4000;   // samples per estimate (default --spp 256: about 1M)
    std::cout << "[Hair] " << n << " samples per estimate\n";

    SurfaceInteraction rec;
    rec.p = Point3(0.0); rec.n = Vector3(0, 0, 1); rec.gn = Vector3(0, 0, 1); rec.dpdu = Vector3(1, 0, 0);   // hair frame = world frame
#if defined(RAYT_SPECTRAL_MODE)
    const SampledWavelengths lambda = SampledWavelengths::sampleVisible(Real(0.37));
#else
    const SampledWavelengths lambda;
#endif

    auto makeHair = [](Vector3 sigma, double betaM, double betaN, double alpha = 2.0) {
        Hair::Params p; p.sigmaA = sigma; p.betaM = Real(betaM); p.betaN = Real(betaN); p.alphaDeg = Real(alpha);
        return Hair(p);
    };

    // ---- 1. Energy: a fiber without absorption scatters everything it receives (sum of A_p = 1) ----------------------
    // The integral of f |cos| over the sphere, by uniform sampling of directions (independent of sample() and pdf()).
    std::cout << "[Hair] 1. energy of a non-absorbing fiber: integral of f |wi.n| over the sphere = 1\n";
    for (double betaM : { 0.3, 0.6 }) {
        const Hair m = makeHair(Vector3(0.0), betaM, betaM);
        for (auto ang : { std::pair<double, double>{ 0, 90 }, { 25, 60 }, { -40, 120 } }) {
            const Vector3 wo = outgoing(ang.first, ang.second);
            Estimate e;
            for (int i = 0; i < n; ++i) {
                const Vector3 wi = uniformSphere();
                e.add(average(m.eval(rec, wo, wi, lambda)) * std::abs(wi.z) * 4.0 * constants::PI);
            }
            e.finish();
            const bool ok = std::abs(e.mean - 1.0) <= 4 * e.se + 0.02;
            std::cout << "[Hair]   beta " << betaM << ", theta_o " << ang.first << " deg, phi_o " << ang.second << " deg: albedo " << e.mean << " +- " << e.se
                << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }
    }

    // ---- 2. pdf integrates to 1; sample() draws from pdf(); the weights f cos / pdf give the albedo -----------------
    std::cout << "[Hair] 2. sample() against pdf() and eval()\n";
    for (double betaM : { 0.25, 0.5 }) {
        const Hair m = makeHair(Vector3(0.0), betaM, betaM);   // (no absorption: albedo 1)
        const Vector3 wo = outgoing(20, 70);
        Estimate pdfInt;
        double gU[3] = { 0, 0, 0 }, gS[3] = { 0, 0, 0 };
        Estimate alb;
        for (int i = 0; i < n; ++i) {
            const Vector3 w = uniformSphere();
            const double pd = m.pdf(rec, wo, w, lambda);
            pdfInt.add(pd * 4.0 * constants::PI);
            for (int c = 0; c < 3; ++c) gU[c] += w[c] * pd * 4.0 * constants::PI / double(n);
            const auto bs = m.sample(rec, wo, Point2(sampling::Random(), sampling::Random()), lambda);
            if (bs && bs->pdf > 0) {
                for (int c = 0; c < 3; ++c) gS[c] += bs->wi[c] / double(n);
                alb.add(average(bs->f) * std::abs(bs->wi.z) / bs->pdf);
            }
            else alb.add(0.0);
        }
        pdfInt.finish(); alb.finish();
        bool ok = std::abs(pdfInt.mean - 1.0) <= 4 * pdfInt.se + 0.02;
        std::cout << "[Hair]   beta " << betaM << ": integral of pdf " << pdfInt.mean << " +- " << pdfInt.se << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
        const double tol = 0.02 + 6.0 / std::sqrt(double(n));
        ok = true;
        for (int c = 0; c < 3; ++c) ok = ok && std::abs(gU[c] - gS[c]) <= tol;
        std::cout << "[Hair]   beta " << betaM << ": mean of sampled wi (" << gS[0] << ", " << gS[1] << ", " << gS[2] << ") vs integral of wi pdf ("
            << gU[0] << ", " << gU[1] << ", " << gU[2] << ")" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
        ok = std::abs(alb.mean - 1.0) <= 4 * alb.se + 0.02;
        std::cout << "[Hair]   beta " << betaM << ": albedo by the weights f cos / pdf " << alb.mean << " +- " << alb.se << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    // ---- 3. Absorption: the albedo of a pigmented fiber is the sum of the lobe attenuations A_p ------------------------
    std::cout << "[Hair] 3. pigmented fiber: albedo = sum of A_p\n";
    {
        const Hair m = makeHair(Hair::sigmaFromMelanin(Real(1.3), Real(0.2)), 0.35, 0.35);
        for (auto ang : { std::pair<double, double>{ 0, 90 }, { 30, 50 } }) {
            const Vector3 wo = outgoing(ang.first, ang.second);
            const auto ap = m.attenuation(wo, lambda);
            Real expected = 0;
            for (int p = 0; p < 4; ++p) expected += average(ap[size_t(p)]);
            Estimate e;
            for (int i = 0; i < n; ++i) {
                const auto bs = m.sample(rec, wo, Point2(sampling::Random(), sampling::Random()), lambda);
                e.add(bs && bs->pdf > 0 ? average(bs->f) * std::abs(bs->wi.z) / bs->pdf : 0.0);
            }
            e.finish();
            const bool ok = std::abs(e.mean - double(expected)) <= 4 * e.se + 0.01;
            std::cout << "[Hair]   theta_o " << ang.first << ", phi_o " << ang.second << ": expected " << double(expected) << "   by sample() " << e.mean << " +- " << e.se
                << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }
    }

    // ---- 4. Where the lobes point (the geometry convention between the shape and the material) -------------------------
    std::cout << "[Hair] 4. lobe directions in the cross-section (smooth fiber, no cuticle tilt)\n";
    {
        // Black fiber: only R survives. Its peak is the mirror direction about the fiber's normal: (-x, -y, z).
        const Hair dark = makeHair(Vector3(60.0), 0.08, 0.1, 0.0);
        const Vector3 wo = outgoing(0, 60);
        const Vector3 mirror(-wo.x, -wo.y, wo.z);
        const Vector3 off = outgoing(0, 100);   // the same inclination, 40 degrees away in azimuth
        const double fm = average(dark.eval(rec, wo, mirror, lambda)), fo = average(dark.eval(rec, wo, Vector3(off.x, off.y, off.z), lambda));
        const bool ok = fm > 20.0 * fo;
        std::cout << "[Hair]   R lobe: f at the mirror direction " << fm << " vs 40 deg away " << fo << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
        // Clear fiber looked at head-on (h = 0): TT goes straight through, wi = -wo.
        const Hair clear = makeHair(Vector3(0.0), 0.08, 0.1, 0.0);
        const Vector3 wz(0, 0, 1);
        const double ft = average(clear.eval(rec, wz, Vector3(0, 0, -1), lambda)), fs = average(clear.eval(rec, wz, Vector3(0, 1, 0), lambda));
        const bool ok2 = ft > 20.0 * fs;
        std::cout << "[Hair]   TT lobe: f straight through " << ft << " vs sideways " << fs << (ok2 ? "   OK" : "   FAIL") << "\n";
        allPassed = ok2 && allPassed;
    }

    // ---- 5. The intersector: one cylinder, the entry only, and the BVH against brute force -------------------------
    std::cout << "[Hair] 5. cylinder chains\n";
    {
        auto mat = std::make_shared<Hair>(Hair::Params());
        HairCurves one(mat);
        one.addStrand({ Point3(0, -1, 0), Point3(0, 1, 0) }, Real(0.05), Real(0.05));   // a vertical cylinder of radius 0.05
        one.build();
        SurfaceInteraction h;
        const Ray fromSide(Point3(-2, 0.3, 0), Vector3(1, 0, 0), Real(0), nullptr);
        bool ok = one.hit(fromSide, h) && std::abs(h.t - 1.95) < 1e-9 && glm::length(h.n - Vector3(-1, 0, 0)) < 1e-9 && std::abs(glm::dot(h.dpdu, Vector3(0, 1, 0))) > 0.999999;
        std::cout << "[Hair]   a ray from the side: t = " << h.t << " (expected 1.95), normal (" << h.n.x << ", " << h.n.y << ", " << h.n.z << ")" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
        const Ray fromInside(Point3(-0.049, 0.0, 0), Vector3(1, 0, 0), Real(0), nullptr);
        SurfaceInteraction h2;
        ok = !one.hit(fromInside, h2) && !one.occluded(fromInside);
        std::cout << "[Hair]   a ray that starts inside does not hit the far wall" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
        const Ray miss(Point3(-2, 0.3, 0.06), Vector3(1, 0, 0), Real(0), nullptr);
        ok = !one.hit(miss, h2);
        std::cout << "[Hair]   a ray passing 0.01 beside the cylinder misses" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }
    {
        auto mat = std::make_shared<Hair>(Hair::Params());
        std::mt19937 rng(5);
        std::uniform_real_distribution<Real> U(0, 1);
        HairCurves all(mat);
        std::vector<std::shared_ptr<HairCurves>> singles;
        for (int s = 0; s < 150; ++s) {
            std::vector<Point3> pts;
            Point3 p(U(rng) - 0.5, U(rng) - 0.5, U(rng) - 0.5);
            Vector3 d = glm::normalize(Vector3(U(rng) - 0.5, U(rng) - 0.5, U(rng) - 0.5));
            for (int k = 0; k < 9; ++k) { pts.push_back(p); d = glm::normalize(d + Vector3(U(rng) - 0.5, U(rng) - 0.5, U(rng) - 0.5) * Real(0.4)); p += d * Real(0.1); }
            const Real r = Real(0.004) + Real(0.004) * U(rng);
            all.addStrand(pts, r, r);
            auto one = std::make_shared<HairCurves>(mat);
            one->addStrand(pts, r, r); one->build();
            singles.push_back(one);
        }
        all.build();
        int agree = 0, total = 0, hits = 0, occBad = 0;
        for (int i = 0; i < 4000; ++i) {
            const Vector3 d = uniformSphere();
            const Point3 o = Point3(U(rng) - 0.5, U(rng) - 0.5, U(rng) - 0.5) * Real(0.5) - d * Real(1.0);
            Ray r(o, d, Real(0), nullptr);
            Real bestT = std::numeric_limits<Real>::infinity();
            for (const auto& s : singles) { Real t; if (s->intersectT(r, t)) bestT = std::min(bestT, t); }
            Real t = 0;
            const bool h = all.intersectT(r, t);
            ++total;
            if (h == std::isfinite(bestT) && (!h || std::abs(t - bestT) < 1e-9)) ++agree;
            hits += h ? 1 : 0;
            if (all.occluded(r) != h) ++occBad;
        }
        const bool ok = agree == total && occBad == 0 && hits > 200;
        std::cout << "[Hair]   BVH over " << all.numSegments() << " segments vs brute force: " << agree << " / " << total << " rays agree (" << hits << " hits), occluded() mismatches " << occBad
            << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    // ---- 5b. One non-absorbing fiber in a white furnace returns the sky (any offset across the width, any tilt) ---------------
    std::cout << "[Hair] 5b. one fiber in a white furnace" << std::endl;
    {
        auto whiteSky = std::make_shared<EnvLight>(std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0)))), 1.0);
        Hair::Params hp; hp.sigmaA = Vector3(0.0); hp.betaM = Real(0.3); hp.betaN = Real(0.3);
        auto mat = std::make_shared<Hair>(hp);
        auto curves = std::make_shared<HairCurves>(mat);
        curves->addStrand({ Point3(0, 0.5, 0), Point3(0, 0.0, 0), Point3(0, -0.5, 0) }, Real(0.01), Real(0.01));
        curves->build();
        auto world = std::make_shared<Aggregate>();
        world->add(curves);
        auto scene = std::make_shared<Scene>(world);
        scene->setEnvLight(whiteSky);
        scene->finalize("spatial");
        const int nn = std::max(opt.spp, 64) * 100;
        for (auto aim : { std::pair<Point3, Vector3>{ Point3(0.0, 0.1, 0.0), Vector3(0, 0, -1) }, { Point3(0.007, 0.1, 0.0), Vector3(0, 0, -1) },
                          { Point3(0.0, 0.1, 0.0), glm::normalize(Vector3(0.3, 0.5, -0.8)) }, { Point3(0.002, 0.1, 0.0), glm::normalize(Vector3(-0.9, -0.1, -0.2)) } }) {
            const Ray ray(aim.first - aim.second * Real(0.5), aim.second, Real(0), nullptr);
            SurfaceInteraction h;
            const bool hits = scene->hit(ray, h);
            Vector3 mean, se; int bad = 0;
            EstimateRadiance(*scene, ray, nn, 64, mean, se, bad);
            const bool ok = hits && ReportLightCase("one fiber, white furnace, ray dir (" + std::to_string(aim.second.x).substr(0, 5) + ", " + std::to_string(aim.second.y).substr(0, 5) + ", " + std::to_string(aim.second.z).substr(0, 5) + ")" + (hits ? "" : "  [misses]"), Vector3(1.0), mean, se, bad);
            allPassed = ok && allPassed;
        }
    }

    // ---- 6. White furnace: a lock of non-absorbing hair in a uniform white sky returns the sky ---------------------------
    std::cout << "[Hair] 6. white furnace around a lock of non-absorbing hair\n";
    {
        auto whiteSky = std::make_shared<EnvLight>(std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0)))), 1.0);
        Hair::Params hp; hp.sigmaA = Vector3(0.0); hp.betaM = Real(0.3); hp.betaN = Real(0.3);
        auto mat = std::make_shared<Hair>(hp);
        auto curves = std::make_shared<HairCurves>(mat);
        hairgen::LockParams lp;
        // A lock this sparse (about 5 fibers across a ray) keeps the paths short: with 600 fibers in the same width (2.7 times over-packed) a
        // non-absorbing lock needs hundreds of bounces and the Russian roulette weights have a heavy tail (the mean comes out low, the error bar too small).
        lp.strands = 150; lp.length = Real(0.10); lp.width = Real(0.012); lp.tipWidth = Real(0.5); lp.segments = 12; lp.waveAmp = Real(0.006); lp.waveLength = Real(0.07);
        lp.rootRadius = Real(0.0004); lp.tipRadius = Real(0.0003); lp.frizz = Real(0.0003);
        hairgen::MakeLock(*curves, lp);
        curves->build();
        auto world = std::make_shared<Aggregate>();
        world->add(curves);
        auto scene = std::make_shared<Scene>(world);
        scene->setEnvLight(whiteSky);
        scene->finalize("spatial");
        const int nn = std::max(opt.spp, 64) * 80;
        // rays through the middle of the lock, starting 0.5 m before it. The centerline of MakeLock sways, so aim at it at arc length s.
        auto centerAt = [&](Real s) {
            const Real k = Real(2) * constants::PI / lp.waveLength;
            return Point3(lp.waveAmp * std::sin(k * s), -s, lp.waveAmp * Real(0.6) * std::sin(k * s * Real(0.7) + Real(1.0)));
        };
        for (auto aim : { std::pair<Point3, Vector3>{ centerAt(Real(0.05)), Vector3(0, 0, -1) }, { centerAt(Real(0.04)), Vector3(1, 0, 0) },
                          { centerAt(Real(0.06)), glm::normalize(Vector3(0.6, 0.2, -0.8)) } }) {
            const Ray ray(aim.first - aim.second * Real(0.5), aim.second, Real(0), nullptr);
            SurfaceInteraction h;
            const bool hits = scene->hit(ray, h);   // the test is empty if the ray misses the lock
            Vector3 mean, se; int bad = 0;
            EstimateRadiance(*scene, ray, nn, 256, mean, se, bad);
            const Vector3 expected(1.0);
            const bool ok = hits && ReportLightCase("hair lock, white furnace, ray (" + std::to_string(aim.second.x).substr(0, 5) + ", " + std::to_string(aim.second.y).substr(0, 5) + ", " + std::to_string(aim.second.z).substr(0, 5) + ")" + (hits ? "" : "  [misses the lock]"), expected, mean, se, bad);
            allPassed = ok && allPassed;
        }
    }

    std::cout << "[Hair] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
