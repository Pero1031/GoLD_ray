// Subsurface scattering checks.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Core/SpectrumUtils.hpp"
#include "Geometry/Quad.hpp"
#include "Geometry/Sphere.hpp"
#include "Lights/PointLight.hpp"
#include "Materials/CoatedDiffuse.hpp"
#include "Materials/Lambertian.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Integrator.hpp"

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

    // =============================================================================================
    // The random walk (Materials/SubsurfaceWalk.hpp): the same skin, with the light followed through the inside of the object.
    // =============================================================================================
    auto makeWalkSkin = [&](const Vector3& A, const Vector3& distance, Real g, Real coat) {
        auto skin = std::make_shared<CoatedDiffuse>(RgbColor(A), Real(0.5), coat, Real(0));
        skin->setSubsurface(std::make_shared<SubsurfaceProfile>(distance, eta));
        skin->setSubsurfaceWalk(g);
        return skin;
    };

    // ---- 4. the conversion (diffuse color, scatter distance, g) -> coefficients -------------------------------
    // Reference: the albedo of an isotropic half space under diffuse light, exact from Chandrasekhar's H function,
    //   1 / H(mu) = sqrt(1 - a) + a / 2 Integral mu' H(mu') / (mu + mu') dmu',   A = 1 - 2 sqrt(1 - a) Integral H(mu) mu dmu.
    // The engine uses van de Hulst's closed form, an approximation: the check is that its error stays below the documented 0.004.
    {
        auto halfSpaceAlbedo = [](double a) {
            const int N = 200;
            std::vector<double> mu(N), H(N, 1.0);
            for (int i = 0; i < N; ++i) mu[i] = (i + 0.5) / N;
            const double root = std::sqrt(1.0 - a);
            for (int it = 0; it < 20000; ++it) {
                double change = 0;
                for (int i = 0; i < N; ++i) {
                    double sum = 0;
                    for (int j = 0; j < N; ++j) sum += mu[j] * H[j] / (mu[i] + mu[j]);
                    const double h = 1.0 / (root + 0.5 * a * sum / N);
                    change = std::max(change, std::abs(h - H[i]));
                    H[i] = h;
                }
                if (change < 1e-11) break;
            }
            double m1 = 0;
            for (int i = 0; i < N; ++i) m1 += H[i] * mu[i] / N;
            return 1.0 - 2.0 * root * m1;
        };
        double worst = 0, worstRound = 0;
        std::cout << "[Sss]   4. conversion: diffuse color -> single-scattering albedo -> exact half-space albedo:";
        for (double A : { 0.1, 0.3, 0.52, 0.8, 0.95 }) {
            const double s = double(sss::walkSimilarity(Real(A)));
            const double back = halfSpaceAlbedo(1.0 - s * s);
            worst = std::max(worst, std::abs(back - A));
            worstRound = std::max(worstRound, std::abs(double(sss::walkAlbedo(Real(s))) - A));
            std::cout << "  " << A << " -> " << 1.0 - s * s << " -> " << back;
        }
        // g only moves scattering between "real" and "reduced": sigma_a and sigma_s (1 - g) stay.
        const sss::WalkCoefficients k0 = sss::walkCoefficients(Real(0.8), Real(0.002), Real(0)), k8 = sss::walkCoefficients(Real(0.8), Real(0.002), Real(0.8));
        const double reduced0 = double(k0.sigmaT), reduced8 = double(k8.sigmaT - k8.sigmaS + k8.sigmaS * Real(0.2));
        const bool okG = std::abs(reduced0 - 500.0) < 1e-6 && std::abs(reduced8 - 500.0) < 1e-6 && std::abs(double(k0.sigmaT - k0.sigmaS) - double(k8.sigmaT - k8.sigmaS)) < 1e-6;
        const bool ok = worst < 0.004 && worstRound < 1e-9 && okG;
        std::cout << "\n[Sss]      largest error of the closed form " << worst << " (limit 0.004), inverse round trip " << worstRound
            << ", reduced extinction with g = 0 / 0.8: " << reduced0 << " / " << reduced8 << " 1/m" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = allPassed && ok;
    }

    // ---- 5. the lobes of the boundary: sample() against eval() and against the closed forms ------------------------
    {
#if defined(RAYT_SPECTRAL_MODE)
        const SampledWavelengths lambda = SampledWavelengths::sampleVisible(Real(0.37));
#else
        const SampledWavelengths lambda;
#endif
        const Real coat = Real(0.4);
        const auto skin = makeWalkSkin(Vector3(albedo), Vector3(scatterDistance), Real(0), coat);
        SurfaceInteraction rec;
        rec.p = Point3(0.0); rec.gn = Vector3(0, 1, 0); rec.n = rec.gn; rec.frontFace = true; rec.matPtr = skin.get();
        const int M = 400000;

        // (a) from outside: the share that enters is 1 - w F(cos_o), by sampling and by integrating eval() over the inner hemisphere
        const Vector3 wo = glm::normalize(Vector3(0.5, 0.8, 0.2));
        const double expectIn = 1.0 - double(coat) * double(fresnel::fresnelDielectric(wo.y, Real(1), Real(1.5)));
        double sampledIn = 0, sampledAll = 0;
        for (int i = 0; i < M; ++i) {
            const auto s = skin->sample(rec, wo, sampling::Random2D(), lambda, TransportMode::Radiance);
            if (!s || s->pdf <= 0) continue;
            const double w = double(s->f[0]) * std::abs(double(glm::dot(rec.n, s->wi))) / double(s->pdf);
            sampledAll += w;
            if (glm::dot(rec.gn, s->wi) < 0) sampledIn += w;
        }
        sampledIn /= M; sampledAll /= M;
        double integratedIn = 0, pdfIn = 0;
        const int K = 200;
        for (int i = 0; i < K; ++i) for (int j = 0; j < 2 * K; ++j) {   // uniform in (cos, phi) over the hemisphere below the surface
            const double c = (i + 0.5) / K, phi = (j + 0.5) * constants::PI / K, sn = std::sqrt(1.0 - c * c);
            const Vector3 wi(Real(sn * std::cos(phi)), Real(-c), Real(sn * std::sin(phi)));
            const double dw = 2.0 * constants::PI / (K * 2.0 * K);
            integratedIn += double(skin->eval(rec, wo, wi, lambda, TransportMode::Radiance)[0]) * c * dw;
            pdfIn += double(skin->pdf(rec, wo, wi, lambda)) * dw;
        }
        const bool okIn = std::abs(sampledIn - expectIn) < 0.01 * expectIn && std::abs(integratedIn - expectIn) < 0.002 * expectIn && pdfIn > 0.04 && pdfIn < 0.96;
        std::cout << "[Sss]   5a. entering light: by sample() " << sampledIn << ", by eval() " << integratedIn << ", expected 1 - w F = " << expectIn
            << " (probability of the event " << pdfIn << ", all lobes return " << sampledAll << ")" << (okIn ? "   OK" : "   FAIL") << "\n";

        // (b) from inside: every walk leaves, with a mean weight of 1 (the exit lobe is normalized), and only to the outside
        const Vector3 woInside = glm::normalize(Vector3(0.3, -0.9, 0.1));
        double sampledOut = 0; int wrongSide = 0;
        for (int i = 0; i < M; ++i) {
            const auto s = skin->sample(rec, woInside, sampling::Random2D(), lambda, TransportMode::Radiance);
            if (!s || s->pdf <= 0) continue;
            if (glm::dot(rec.gn, s->wi) <= 0) ++wrongSide;
            sampledOut += double(s->f[0]) * std::abs(double(glm::dot(rec.n, s->wi))) / double(s->pdf);
        }
        sampledOut /= M;
        const bool okOut = std::abs(sampledOut - 1.0) < 0.01 && wrongSide == 0 && isBlack(skin->eval(rec, woInside, woInside, lambda, TransportMode::Radiance));
        std::cout << "[Sss]   5b. leaving light: mean weight " << sampledOut << " (expected 1), samples to the inside " << wrongSide << (okOut ? "   OK" : "   FAIL") << "\n";
        allPassed = allPassed && okIn && okOut;
    }

    // ---- 6. a skin sphere in a white furnace against an independent walk -----------------------------------------
    // Without a coat the radiance of the sphere is the fraction of the walks that leave it again (the rest is absorbed):
    // a short walk written here (cosine entry, exponential flights, Henyey-Greenstein) gives the same number without the renderer.
    {
        auto escapeFraction = [](double sigmaT, double sigmaS, double g, double R, int walks) {
            long long out = 0;
            for (int i = 0; i < walks; ++i) {
                const Vector3 c = sampling::CosineSampleHemisphere(sampling::Random2D());
                Vector3 p(0, 0, Real(R)), d(c.x, c.y, -c.z);   // enters at the top, about the inward normal
                for (int step = 0; step < 1000000; ++step) {
                    const double t = -std::log(1.0 - double(sampling::Random())) / sigmaT;
                    const double b = double(glm::dot(p, d)), cc = double(glm::dot(p, p)) - R * R;
                    const double tExit = -b + std::sqrt(std::max(b * b - cc, 0.0));
                    if (t >= tExit) { ++out; break; }
                    if (double(sampling::Random()) * sigmaT >= sigmaS) break;   // absorbed
                    p += d * Real(t);
                    d = sampleHenyeyGreenstein(d, Real(g), sampling::Random2D());
                }
            }
            return double(out) / double(walks);
        };
        auto whiteSky = std::make_shared<EnvLight>(std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0)))), 1.0);
        auto furnace = [&](const std::shared_ptr<CoatedDiffuse>& skin, Real R, int samples, Vector3& mean, Vector3& se, int& bad) {
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Sphere>(Point3(0.0), R, skin));
            auto scene = std::make_shared<Scene>(world);
            scene->setEnvLight(whiteSky);
            scene->finalize("spatial");
            EstimateRadiance(*scene, Ray(Point3(Real(0.3) * R, Real(0.2) * R, Real(10) * R), Vector3(0, 0, -1)), samples, 6, mean, se, bad);
        };
        struct Case { double A, radii, g; };
        for (const Case& k : { Case{ 0.8, 2.0, 0.0 }, Case{ 0.8, 10.0, 0.0 }, Case{ 0.4, 3.0, 0.0 }, Case{ 0.8, 3.0, 0.8 } }) {
            const sss::WalkCoefficients co = sss::walkCoefficients(Real(k.A), scatterDistance, Real(k.g));
            const double R = k.radii * double(scatterDistance);
            const double expected = escapeFraction(double(co.sigmaT), double(co.sigmaS), k.g, R, 400000);
            Vector3 mean, se; int bad = 0;
            furnace(makeWalkSkin(Vector3(Real(k.A)), Vector3(scatterDistance), Real(k.g), Real(0)), Real(R), n / 4, mean, se, bad);
            const std::string name = "6. furnace: sphere of radius " + std::to_string(int(k.radii)) + " scatter distances, color " + std::to_string(k.A).substr(0, 4)
                + ", g " + std::to_string(k.g).substr(0, 3) + ", against an independent walk";
            std::cout << "[Sss]   " << name << "\n";
            allPassed = ReportLightCase(name, Vector3(Real(expected)), mean, se, bad) && allPassed;
        }

        // Colored skin: each channel is its own medium, walked together (one channel chooses the distances, the others are weighted).
        {
            const Vector3 A(0.80, 0.52, 0.42), dist = Vector3(3.67, 1.37, 0.68) * Real(0.001);
            const double R = 0.004;
            Vector3 expected(0.0);
            for (int c = 0; c < 3; ++c) {
                const sss::WalkCoefficients co = sss::walkCoefficients(A[c], dist[c], Real(0));
                expected[c] = Real(escapeFraction(double(co.sigmaT), double(co.sigmaS), 0.0, R, 400000));
            }
            Vector3 mean, se; int bad = 0;
            furnace(makeWalkSkin(A, dist, Real(0), Real(0)), Real(R), n / 4, mean, se, bad);
#if defined(RAYT_SPECTRAL_MODE)
            // The spectral build converts at each wavelength; the RGB reference is only a guide here (no verdict).
            std::cout << "[Sss]   6. furnace: colored skin sphere (4 mm), info: (" << mean.x << ", " << mean.y << ", " << mean.z << "), per-channel walk ("
                << expected.x << ", " << expected.y << ", " << expected.z << ")\n";
#else
            const std::string name = "6. furnace: colored skin sphere (4 mm), each channel against an independent walk";
            std::cout << "[Sss]   " << name << "\n";
            allPassed = ReportLightCase(name, expected, mean, se, bad) && allPassed;
#endif
        }

        // ---- 7. a thick sphere returns the diffuse color (what the conversion is for); g does not change it ------------
        // Not exact: the closed form (check 4), the curvature and, with g, the similarity relation leave a little. The limit is 1 %.
        // The sphere must be large against the diffusion length (6 scatter distances at this color): a radius of 40 scatter distances
        // still returns 3 % more than a flat slab, because light that has wandered sideways finds the convex surface sooner.
        for (double g : { 0.0, 0.8 }) {
            Vector3 mean, se; int bad = 0;
            furnace(makeWalkSkin(Vector3(albedo), Vector3(scatterDistance), Real(g), Real(0)), Real(1000) * scatterDistance, n / 16, mean, se, bad);
            const bool ok = bad == 0 && std::abs(mean.y - albedo) < 0.01 * albedo + 4 * se.y;
            std::cout << "[Sss]   7. thick sphere (1000 scatter distances), g " << g << ": " << mean.y << " +- " << se.y << " against the diffuse color " << albedo
                << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = allPassed && ok;
        }
        {   // the same sphere at 40 scatter distances, for the record (no verdict)
            Vector3 mean, se; int bad = 0;
            furnace(makeWalkSkin(Vector3(albedo), Vector3(scatterDistance), Real(0), Real(0)), Real(40) * scatterDistance, n / 16, mean, se, bad);
            std::cout << "[Sss]   7. info: at 40 scatter distances the sphere returns " << mean.y << " (curvature)\n";
        }

        // ---- 8. information: a thin sphere, the walk against the BSSRDF (which does not know the thickness) ------------------
        for (double radii : { 0.5, 2.0 }) {
            Vector3 mw, sw, mb, sb; int bw = 0, bb = 0;
            furnace(makeWalkSkin(Vector3(albedo), Vector3(scatterDistance), Real(0), Real(0)), Real(radii) * scatterDistance, n / 16, mw, sw, bw);
            furnace(makeSkin(), Real(radii) * scatterDistance, n / 16, mb, sb, bb);
            std::cout << "[Sss]   8. info: sphere of radius " << radii << " scatter distances in the furnace: random walk " << mw.y << ", BSSRDF " << mb.y << "\n";
        }
    }

    std::cout << "[Sss] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
