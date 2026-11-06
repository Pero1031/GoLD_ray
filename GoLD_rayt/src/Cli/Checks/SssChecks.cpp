// Subsurface scattering checks.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

int RunSssTests(const CliOptions& opt) {
    bool allPassed = true;
    const int n = opt.spp * 4096;   // samples per test (default --spp 256 -> about 1M)
    const int depth = 4;
    const Real albedo = 0.8, scatterDistance = 0.01, eta = 1.4;
    const Real d = scatterDistance / sss::shape(albedo);   // length scale of the profile (m)
    std::cout << "[Sss] " << n << " samples per test, albedo " << albedo << ", scatter distance " << scatterDistance * 1000
        << " mm (profile length scale " << d * 1000 << " mm)\n";

    // ---- 1. the radial profile itself ---------------------------------------------------------
    {
        const double R = 40.0 * d;
        const int N = 400000;
        double integral = 0;   // Integral of profile(r) 2 pi r dr over [0, R] must equal cdf(R)
        for (int i = 0; i < N; ++i) {
            const double r = (i + 0.5) * R / N;
            integral += sss::profile(r, d) * 2.0 * constants::PI * r * (R / N);
        }
        const bool okNorm = std::abs(integral - sss::cdf(R, d)) < 1e-4;

        bool okInverse = true;
        for (double u : { 0.001, 0.1, 0.5, 0.9, 0.999 }) okInverse = okInverse && std::abs(sss::cdf(sss::sampleRadius(u, d), d) - u) < 1e-9;

        double mean = 0, sq = 0;   // the radial density is (e^-x + e^-x/3) / 4 d: mean 2.5 d, variance 7.75 d^2
        const int M = 1000000;
        for (int i = 0; i < M; ++i) { const double r = sss::sampleRadius(sampling::Random(), d); mean += r; sq += r * r; }
        mean /= M; sq /= M;
        const double se = std::sqrt(7.75) * d / std::sqrt(double(M));
        const bool okMean = std::abs(mean - 2.5 * d) < 4.0 * se;
        std::cout << "[Sss]   profile: Integral 2 pi r P(r) dr = " << integral << " vs cdf " << sss::cdf(R, d) << (okNorm ? " OK" : " FAIL")
            << ";  cdf(sampleRadius(u)) = u" << (okInverse ? " OK" : " FAIL")
            << ";  mean radius " << mean / d << " d vs 2.5 d" << (okMean ? " OK" : " FAIL") << "\n";
        allPassed = allPassed && okNorm && okInverse && okMean;
    }

    auto makeSkin = [&]() {
        auto profile = std::make_shared<SubsurfaceProfile>(Vector3(scatterDistance), eta);
        auto skin = std::make_shared<CoatedDiffuse>(RgbColor(albedo), Real(0.5), Real(0.0), Real(0));   // no coat: pure subsurface
        skin->setSubsurface(profile);
        return skin;
    };
    auto run = [&](const std::string& name, const Scene& scene, const Ray& ray, const Vector3& expected) {
        Vector3 mean, se; int bad = 0;
        EstimateRadiance(scene, ray, n, depth, mean, se, bad);
        std::cout << "[Sss]   " << name << "\n";
        allPassed = ReportLightCase(name, expected, mean, se, bad) && allPassed;
    };

    // ---- 2. white furnace: a flat skin plane under a uniform sky returns exactly its albedo ---------
    {
        auto world = std::make_shared<Aggregate>();
        world->add(std::make_shared<Quad>(Point3(-100, 0, -100), Vector3(0, 0, 200), Vector3(200, 0, 0), makeSkin()));
        auto scene = std::make_shared<Scene>(world);
        auto sky = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        scene->setEnvLight(std::make_shared<EnvLight>(sky, 1.0));
        scene->finalize("spatial");
        run("furnace: skin plane under a uniform sky, expect the albedo",
            *scene, Ray(Point3(-3, 3, 0), glm::normalize(Vector3(1, -1, 0))), Vector3(albedo));
    }

    // ---- 3. knife edge: the profile seen across a shadow edge -----------------------------------
    // A point light 100 m above the edge of a black barrier gives a uniform irradiance with E / pi = 1
    // on one side and a sharp shadow on the other. At signed distance x from the edge the radiance is
    //   albedo * (1 - D(|x|)) lit side,   albedo * D(|x|) shadow side,
    //   D(t) = Integral over phi in (-pi/2, pi/2) of (e^(-s/d) + 3 e^(-s/3d)) / (8 pi),  s = t / cos(phi)
    // (the profile mass beyond a straight line at distance t).
    {
        auto world = std::make_shared<Aggregate>();
        world->add(std::make_shared<Quad>(Point3(-100, 0, -100), Vector3(0, 0, 200), Vector3(200, 0, 0), makeSkin()));
        world->add(std::make_shared<Quad>(Point3(-100, 1, -100), Vector3(0, 0, 200), Vector3(100, 0, 0),
            std::make_shared<Lambertian>(RgbColor(0.0))));   // black barrier over x < 0
        auto scene = std::make_shared<Scene>(world);
        const Real height = 100.0;
        scene->addLight(std::make_shared<PointLight>(Point3(0, height, 0), RgbColor(1.0), constants::PI * height * height));
        scene->finalize("spatial");

        auto dark = [&](double t) {   // fraction of the profile beyond a line at distance t
            const int K = 20000;
            double sum = 0;
            for (int i = 0; i < K; ++i) {
                const double phi = -constants::PI / 2 + (i + 0.5) * constants::PI / K;
                const double s = t / std::cos(phi);
                sum += (std::exp(-s / d) + 3.0 * std::exp(-s / (3.0 * d))) / (8.0 * constants::PI) * (constants::PI / K);
            }
            return sum;
        };
        for (double x : { -3.0, -1.0, -0.4, 0.4, 1.0, 3.0 }) {   // in units of d; negative = shadow side
            const double expected = x < 0 ? albedo * dark(-x * d) : albedo * (1.0 - dark(x * d));
            const Point3 target(x * d, 0, 0);
            const Point3 origin = target + Vector3(0, 0.4, 1.6);
            run("knife edge at x = " + std::to_string(x) + " d  (" + (x < 0 ? "shadow" : "lit") + " side)",
                *scene, Ray(origin, glm::normalize(target - origin)), Vector3(expected));
        }
    }

    std::cout << "[Sss] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
