// Volumetric random walk checks.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/EngineHeaders.hpp"

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

    std::cout << "[Volume] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
