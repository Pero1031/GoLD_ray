// Volumetric random walk checks.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Core/Sampling.hpp"
#include "Geometry/Quad.hpp"
#include "Geometry/Sphere.hpp"
#include "Materials/Dielectric.hpp"
#include "Materials/LayeredSkin.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Integrator.hpp"

namespace rayt::cli {


// -----------------------------------------------------------------------------
// Volumetric random walk tests (--scene volume)
// -----------------------------------------------------------------------------
static bool ReportVolumeCase(const std::string& name, Real expected, const Vector3& mean, const Vector3& se, int bad, bool informational = false) {
    bool ok = (bad == 0);
    for (int c = 0; c < 3; ++c) ok = ok && std::abs(mean[c] - expected) <= 4.0 * se[c] + 0.005 * std::abs(expected) + 1e-6;
    std::cout << "[Volume]   " << name << "\n"
        << "             expected " << expected << "   got (" << mean.x << ", " << mean.y << ", " << mean.z << ") +- ("
        << se.x << ", " << se.y << ", " << se.z << ")   ratio " << mean.y / expected
        << (bad ? "  NaN samples: " + std::to_string(bad) : std::string())
        << (informational ? "   (information)" : (ok ? "   OK" : "   FAIL")) << "\n";
    return informational || ok;
}

int RunVolumeTests(const CliOptions& opt) {
    bool allPassed = true;
    const int n = opt.spp * 1024;   // samples per test (default --spp 256 -> about 260k)
    const int depth = 512;          // surface bounces; collisions inside the medium are not counted
    std::cout << "[Volume] " << n << " samples per test\n";

    auto whiteSky = [] {
        return std::make_shared<EnvLight>(std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0)))), 1.0);
    };
    // Index-matched slab (z in [0, thickness], infinite sideways) of a medium under a uniform white sky.
    auto slabOf = [&](std::shared_ptr<const Medium> medium, Real thickness) {
        auto boundary = std::make_shared<Dielectric>(Real(1.0), Real(0.0));   // IOR 1: no reflection, no refraction
        boundary->setMedium(std::move(medium));
        auto world = std::make_shared<Aggregate>();
        world->add(std::make_shared<Quad>(Point3(-50, -50, 0), Vector3(0, 100, 0), Vector3(100, 0, 0), boundary));          // normal -z
        world->add(std::make_shared<Quad>(Point3(-50, -50, thickness), Vector3(100, 0, 0), Vector3(0, 100, 0), boundary));  // normal +z
        auto scene = std::make_shared<Scene>(world);
        scene->setEnvLight(whiteSky());
        scene->finalize("spatial");
        return scene;
    };
    auto slab = [&](Real sigmaT, Real albedo, Real thickness, Real g = 0) {
        return slabOf(std::make_shared<HomogeneousMedium>(sigmaT, albedo, g), thickness);
    };
    const Ray alongZ(Point3(0, 0, -1), Vector3(0, 0, 1));
    auto run = [&](const std::string& name, const Scene& scene, const Ray& ray, Real expected, bool informational = false) {
        Vector3 mean, se; int bad = 0;
        EstimateRadiance(scene, ray, n, depth, mean, se, bad);
        allPassed = ReportVolumeCase(name, expected, mean, se, bad, informational) && allPassed;
    };

    // ---- 1. Beer-Lambert: a purely absorbing slab transmits exp(-sigma_t d) of the sky behind it --------------
    for (Real tau : { Real(1), Real(3) }) {
        const auto scene = slab(100, 0, tau / Real(100));
        run("absorbing slab, sigma_t d = " + std::to_string(int(tau)) + ": transmittance exp(-sigma_t d)", *scene, alongZ, std::exp(-tau));
    }

    // ---- 1b. AbsorbingMedium (transmittance weight, no collisions) must agree with Beer-Lambert as well -----------
    {
        color::SpectralCurve flat;   // sigma_a = 100/m at every wavelength
        flat.add(300, 100); flat.add(900, 100); flat.finalize();
        for (Real tau : { Real(1), Real(3) }) {
            const auto scene = slabOf(std::make_shared<AbsorbingMedium>(flat), tau / Real(100));
            run("AbsorbingMedium slab, sigma_a d = " + std::to_string(int(tau)) + ": transmittance exp(-sigma_a d)", *scene, alongZ, std::exp(-tau));
        }
        const auto scene = slabOf(std::make_shared<AbsorbingMedium>(flat, Real(0.5)), Real(0.04));   // scale 0.5: sigma_a d = 2
        run("AbsorbingMedium slab, scale 0.5, sigma_a d = 2: transmittance exp(-2)", *scene, alongZ, std::exp(Real(-2)));
    }

    // ---- 2. white furnace: a medium without absorption inside a lossless boundary returns the sky -----------
    run("scattering slab, albedo 1, sigma_t d = 5, isotropic: white furnace", *slab(100, 1, 0.05), alongZ, 1);
    run("scattering slab, albedo 1, sigma_t d = 5, g = 0.8: white furnace", *slab(100, 1, 0.05, 0.8), alongZ, 1);
    run("scattering slab, albedo 1, sigma_t d = 5, g = -0.5, slanted ray: white furnace", *slab(100, 1, 0.05, -0.5),
        Ray(Point3(0, 0, -1), glm::normalize(Vector3(0.6, 0.0, 0.8))), 1);
    {
        auto glass = std::make_shared<Dielectric>(Real(1.5), Real(0.0));
        glass->setMedium(std::make_shared<HomogeneousMedium>(Real(200), Real(1), Real(0)));
        auto world = std::make_shared<Aggregate>();
        world->add(std::make_shared<Sphere>(Point3(0.0), Real(0.05), glass));
        auto scene = std::make_shared<Scene>(world);
        scene->setEnvLight(whiteSky());
        scene->finalize("spatial");
        run("glass sphere (IOR 1.5) filled with a medium, sigma_t R = 10, albedo 1: white furnace (centre ray)", *scene, alongZ, 1);
        run("glass sphere, same, off-centre ray", *scene, Ray(Point3(0.03, 0, -1), Vector3(0, 0, 1)), 1);

        auto rough = std::make_shared<Dielectric>(Real(1.5), Real(0.224));
        rough->setMedium(std::make_shared<HomogeneousMedium>(Real(200), Real(1), Real(0)));
        auto worldRough = std::make_shared<Aggregate>();
        worldRough->add(std::make_shared<Sphere>(Point3(0.0), Real(0.05), rough));
        auto sceneRough = std::make_shared<Scene>(worldRough);
        sceneRough->setEnvLight(whiteSky());
        sceneRough->finalize("spatial");
        run("rough glass sphere (roughness 0.224): single-scattering GGX loses energy at every boundary hit", *sceneRough, alongZ, 1, true);
    }

    // ---- 3. semi-infinite slab against the standalone Monte Carlo of the layered skin code ------------------
    // Under a uniform sky the radiance along the normal is the directional-hemispherical reflectance at normal
    // incidence (reciprocity), which a pencil beam entering at the origin measures as the total exit weight.
    {
        const double mua = 100, musp = 4000;
        const auto ref = skin::simulateRadialProfile({ skin::Layer{ 0, mua, musp } }, 128, 1e-5, 0.03, 4000000, 99);
        const Real expected = Real(ref.reflected + ref.lostOutside);
        run("semi-infinite slab, mu_a 100/m, mu_s' 4000/m (albedo 0.9756): reflectance vs standalone Monte Carlo " + std::to_string(expected),
            *slab(Real(mua + musp), Real(musp / (mua + musp)), Real(0.5)), alongZ, expected);
    }

    // ---- 4. chromatic media: hero-wavelength sampling with the MIS weights of the other wavelengths ------------
    {
        const Real thickness = 0.03;
        auto chromatic = [&](const Vector3& sigmaT, const Vector3& albedo, Real g = 0) {
            return slabOf(std::make_shared<ChromaticMedium>(sigmaT, albedo, g), thickness);
        };
        // (a) No absorption: whatever the extinction of each wavelength, the sky comes back (the weights must sum correctly).
        run("chromatic slab, sigma_t (40, 150, 600) /m, albedo 1: white furnace", *chromatic(Vector3(40, 150, 600), Vector3(1.0)), alongZ, 1);
        run("chromatic slab, same, g = 0.7, slanted ray: white furnace", *chromatic(Vector3(40, 150, 600), Vector3(1.0), 0.7),
            Ray(Point3(0, 0, -1), glm::normalize(Vector3(0.6, 0.0, 0.8))), 1);

        // (b) A chromatic medium with equal channels is the gray medium.
        Vector3 meanG, seG, meanC, seC; int badG = 0, badC = 0;
        EstimateRadiance(*slab(100, 0.9, thickness), alongZ, n, depth, meanG, seG, badG);
        EstimateRadiance(*chromatic(Vector3(100), Vector3(0.9)), alongZ, n, depth, meanC, seC, badC);
        {
            bool ok = (badG + badC == 0);
            for (int c = 0; c < 3; ++c) ok = ok && std::abs(meanC[c] - meanG[c]) <= 4.0 * std::sqrt(seC[c] * seC[c] + seG[c] * seG[c]) + 0.005 * meanG[c];
            std::cout << "[Volume]   chromatic slab with equal channels against the gray slab (sigma_t 100, albedo 0.9)\n"
                << "             gray (" << meanG.x << ", " << meanG.y << ", " << meanG.z << ")   chromatic (" << meanC.x << ", " << meanC.y << ", " << meanC.z
                << ")   " << (ok ? "OK" : "FAIL") << "\n";
            allPassed = ok && allPassed;
        }

#if !defined(RAYT_SPECTRAL_MODE)
        // (c) RGB build: channel c of a chromatic medium is the gray medium with that channel's coefficients.
        {
            const Vector3 sigmaT(40, 120, 400), albedo(0.95, 0.8, 0.5);
            EstimateRadiance(*chromatic(sigmaT, albedo, 0.3), alongZ, n, depth, meanC, seC, badC);
            for (int c = 0; c < 3; ++c) {
                EstimateRadiance(*slab(sigmaT[c], albedo[c], thickness, 0.3), alongZ, n, depth, meanG, seG, badG);
                const bool ok = (badG + badC == 0) && std::abs(meanC[c] - meanG[c]) <= 4.0 * std::sqrt(seC[c] * seC[c] + seG[c] * seG[c]) + 0.005 * meanG[c];
                std::cout << "[Volume]   chromatic slab, channel " << c << " (sigma_t " << sigmaT[c] << ", albedo " << albedo[c] << ", g 0.3) against the gray slab\n"
                    << "             gray " << meanG[c] << " +- " << seG[c] << "   chromatic " << meanC[c] << " +- " << seC[c] << "   " << (ok ? "OK" : "FAIL") << "\n";
                allPassed = ok && allPassed;
            }
        }
#endif
    }

    // ---- 6. heterogeneous media (a density field, null collisions) ---------------------------------------------
    {
        const Real thickness = 0.06;
        auto hetero = [&](const Vector3& sigmaBar, const Vector3& albedo, std::shared_ptr<const DensityField> field) {
            return slabOf(std::make_shared<ChromaticMedium>(sigmaBar, albedo, PhaseFunction(), std::move(field)), thickness);
        };
        // (a) Without scattering a ray keeps exp(-sigma_bar * integral of q) of the sky: exact, from a numeric integral along the ray.
        {
            const auto field = std::make_shared<NoiseDensity>(Real(40), 3, Real(0.35), Real(0.65), 3u);
            const int steps = 200000;
            double integral = 0, meanQ = 0;
            for (int i = 0; i < steps; ++i) { const double q = (*field)(Point3(0, 0, (i + 0.5) / steps * thickness)); integral += q * thickness / steps; meanQ += q / steps; }
            const Real sigmaBar = 100;
            run("noise density, albedo 0, sigma_bar 100 /m, mean q along the ray " + std::to_string(meanQ) + ": transmittance exp(-sigma_bar integral q)",
                *hetero(Vector3(sigmaBar), Vector3(0.0), field), alongZ, Real(std::exp(-sigmaBar * integral)));
        }
        // (b) Nothing is absorbed: the sky comes back, whatever the field and the extinction of each wavelength.
        {
            const auto field = std::make_shared<NoiseDensity>(Real(60), 4, Real(0.4), Real(0.6), 5u);
            run("noise density, sigma_bar (40, 150, 600), albedo 1: white furnace", *hetero(Vector3(40, 150, 600), Vector3(1.0), field), alongZ, 1);
        }
        // (c) A constant density q = 0.5 on top of sigma_bar 200 is the homogeneous medium of sigma_t 100.
        {
            Vector3 meanG, seG, meanC, seC; int badG = 0, badC = 0;
            EstimateRadiance(*slab(100, 0.9, thickness), alongZ, n, depth, meanG, seG, badG);
            EstimateRadiance(*hetero(Vector3(200), Vector3(0.9), std::make_shared<ConstantDensity>(Real(0.5))), alongZ, n, depth, meanC, seC, badC);
            bool ok = (badG + badC == 0);
            for (int c = 0; c < 3; ++c) ok = ok && std::abs(meanC[c] - meanG[c]) <= 4.0 * std::sqrt(seC[c] * seC[c] + seG[c] * seG[c]) + 0.005 * meanG[c];
            std::cout << "[Volume]   constant density 0.5 on sigma_bar 200 (null collisions) against the homogeneous slab (sigma_t 100, albedo 0.9)\n"
                << "             homogeneous (" << meanG.x << ", " << meanG.y << ", " << meanG.z << ")   null collisions (" << meanC.x << ", " << meanC.y << ", " << meanC.z
                << ")   " << (ok ? "OK" : "FAIL") << "\n";
            allPassed = ok && allPassed;
        }
#if !defined(RAYT_SPECTRAL_MODE)
        // (d) RGB build: with q = 0.5, channel c is the gray medium of half the extinction of that channel.
        {
            const Vector3 sigmaHalf(40, 120, 400), albedo(0.95, 0.8, 0.5);
            Vector3 meanG, seG, meanC, seC; int badG = 0, badC = 0;
            EstimateRadiance(*hetero(sigmaHalf * Real(2), albedo, std::make_shared<ConstantDensity>(Real(0.5))), alongZ, n, depth, meanC, seC, badC);
            for (int c = 0; c < 3; ++c) {
                EstimateRadiance(*slab(sigmaHalf[c], albedo[c], thickness), alongZ, n, depth, meanG, seG, badG);
                const bool ok = (badG + badC == 0) && std::abs(meanC[c] - meanG[c]) <= 4.0 * std::sqrt(seC[c] * seC[c] + seG[c] * seG[c]) + 0.005 * meanG[c];
                std::cout << "[Volume]   chromatic, q = 0.5, channel " << c << " against the gray slab of sigma_t " << sigmaHalf[c] << "\n"
                    << "             gray " << meanG[c] << " +- " << seG[c] << "   null collisions " << meanC[c] << " +- " << seC[c] << "   " << (ok ? "OK" : "FAIL") << "\n";
                allPassed = ok && allPassed;
            }
        }
#endif
    }

    // ---- 5. two-lobe phase function ---------------------------------------------------------------------------
    {
        const PhaseFunction two(Real(0.9), Real(-0.4), Real(0.7));
        // (a) The sampled directions have the mean cosine w g1 + (1 - w) g2 (a statistical check of the density).
        {
            const int m = 2000000;
            const Vector3 dir = glm::normalize(Vector3(0.3, -0.5, 0.8));
            double sum = 0, sum2 = 0;
            for (int i = 0; i < m; ++i) {
                const double c = glm::dot(two.sample(dir, Point2(sampling::Random(), sampling::Random())), dir);
                sum += c; sum2 += c * c;
            }
            const double mean = sum / m, se = std::sqrt((sum2 / m - mean * mean) / m);
            const bool ok = std::abs(mean - two.meanCosine()) <= 4 * se + 1e-4;
            std::cout << "[Volume]   two-lobe phase function (g1 0.9, g2 -0.4, weight 0.7): mean cosine expected " << two.meanCosine()
                << "   got " << mean << " +- " << se << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }
        // (b) Nothing is absorbed: the sky comes back.
        run("two-lobe scattering slab, albedo 1, sigma_t d = 5: white furnace", *slabOf(std::make_shared<HomogeneousMedium>(Real(100), Real(1), two), Real(0.05)), alongZ, 1);
        // (c) Similarity relation: a thick medium with the same sigma_s' = sigma_s (1 - g) and sigma_a looks the same whatever the
        //     phase function (approximately: the reference is the isotropic case of the check above).
        {
            const double mua = 100, musp = 4000, g = two.meanCosine();
            const double mus = musp / (1 - g);
            run("semi-infinite slab, same mu_a and mu_s' with the two-lobe phase function: reflectance close to the isotropic 0.6443",
                *slabOf(std::make_shared<HomogeneousMedium>(Real(mua + mus), Real(mus / (mua + mus)), two), Real(0.5)), alongZ, Real(0.644338), true);
        }
    }

    // ---- 6. tabulated phase function (Mie scattering of cloud droplets) ---------------------------------------
    {
        // (a) A table of Henyey-Greenstein reproduces it, and a table's eval() and sample() describe the same density:
        //     the fraction of samples inside a cone equals the integral of eval() over the cone.
        auto coneFraction = [](const PhaseFunction& ph, double cosMin) {   // integral of eval over cos in [cosMin, 1]
            const int steps = 400000;
            double s = 0;
            for (int i = 0; i < steps; ++i) s += ph.eval(Real(cosMin + (1.0 - cosMin) * (i + 0.5) / steps));
            return s * 2.0 * constants::PI * (1.0 - cosMin) / steps;
        };
        auto consistency = [&](const std::string& label, const PhaseFunction& ph) {
            const int m = 4000000;
            const Vector3 dir = glm::normalize(Vector3(0.3, -0.5, 0.8));
            const double cuts[3] = { std::cos(5.0 * 3.14159265358979 / 180.0), std::cos(60.0 * 3.14159265358979 / 180.0), -0.5 };
            long hits[3] = { 0, 0, 0 };
            double sum = 0;
            for (int i = 0; i < m; ++i) {
                const double c = glm::dot(ph.sample(dir, Point2(sampling::Random(), sampling::Random())), dir);
                sum += c;
                for (int k = 0; k < 3; ++k) hits[k] += (c >= cuts[k]);
            }
            bool ok = std::abs(sum / m - ph.meanCosine()) < 2e-3;
            std::cout << "[Volume]   " << label << ": mean cosine table " << ph.meanCosine() << ", samples " << sum / m << "\n";
            for (int k = 0; k < 3; ++k) {
                const double expected = coneFraction(ph, cuts[k]), got = double(hits[k]) / m;
                const double se = std::sqrt(expected * (1 - expected) / m);
                const bool kOk = std::abs(got - expected) <= 4 * se + 2e-4;
                std::cout << "[Volume]     samples with cos >= " << cuts[k] << ": integral of eval " << expected << ", sampled " << got << (kOk ? "   OK" : "   FAIL") << "\n";
                ok = ok && kOk;
            }
            allPassed = ok && allPassed;
        };
        const double gT = 0.85;
        const auto hgTable = TabulatedPhase::fromFunction([gT](double c) { return double(henyeyGreenstein(Real(c), Real(gT))); }, 3601);
        const PhaseFunction tabulatedHg(hgTable);
        {
            double worst = 0;
            for (double c : { 0.999, 0.99, 0.9, 0.5, 0.0, -0.5, -0.9 })
                worst = std::max(worst, std::abs(tabulatedHg.eval(Real(c)) / henyeyGreenstein(Real(c), Real(gT)) - 1.0));
            const bool ok = worst < 0.01 && std::abs(tabulatedHg.meanCosine() - gT) < 2e-3;
            std::cout << "[Volume]   table of Henyey-Greenstein g = 0.85: worst eval error " << worst * 100 << " %, mean cosine " << tabulatedHg.meanCosine() << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }
        consistency("table of Henyey-Greenstein g = 0.85", tabulatedHg);
        // (b) The same medium with the table and with the formula gives the same slab reflectance (a thick scattering slab).
        run("slab with the table of HG 0.85, albedo 0.97 against the same slab with the formula (reflectance of the formula)",
            *slabOf(std::make_shared<HomogeneousMedium>(Real(200), Real(0.97), tabulatedHg), Real(0.1)), alongZ,
            [&] { Vector3 mean, se; int bad = 0; EstimateRadiance(*slabOf(std::make_shared<HomogeneousMedium>(Real(200), Real(0.97), PhaseFunction(Real(gT))), Real(0.1)), alongZ, n, depth, mean, se, bad); return mean.y; }());
        // (c) The Mie table of cloud droplets (assets/Cloud/mie_reff10.csv, made by tools/studies/mie_phase.py), if it is there.
        std::string err;
        if (auto mie = TabulatedPhase::fromCsv("assets/Cloud/mie_reff10.csv", &err)) {
            const PhaseFunction mieCloud(mie);
            consistency("Mie table of cloud droplets (r_eff 10 um)", mieCloud);
            run("slab of Mie droplets, albedo 1, sigma_t d = 20: white furnace", *slabOf(std::make_shared<HomogeneousMedium>(Real(400), Real(1), mieCloud), Real(0.05)), alongZ, 1);
        }
        else std::cout << "[Volume]   (the Mie table assets/Cloud/mie_reff10.csv is not there: " << err << "; skipped)\n";
    }

    std::cout << "[Volume] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
