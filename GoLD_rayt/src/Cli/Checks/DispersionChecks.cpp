// Dispersion checks (--scene dispersion): the Sellmeier presets against their published d-line index and Abbe number, and a diamond in a white furnace.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Geometry/Sphere.hpp"
#include "Materials/Dielectric.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Integrator.hpp"

namespace rayt::cli {

int RunDispersionTests() {
    bool allPassed = true;
    auto report = [&](const std::string& name, bool ok, const std::string& detail = "") {
        std::cout << "[Dispersion]   " << name << (detail.empty() ? "" : "  " + detail) << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };
    std::cout << "[Dispersion] Sellmeier presets of glasses\n";

    // (1) n_d and the Abbe number V_d = (n_d - 1) / (n_F - n_C) against the published values of each glass.
    struct Ref { const char* name; Real nd, vd; };
    const Ref refs[] = { { "diamond", 2.417, 55.3 }, { "sapphire", 1.768, 72.2 }, { "bk7", 1.5168, 64.17 }, { "silica", 1.4585, 67.8 } };
    for (const Ref& r : refs) {
        Sellmeier s;
        if (!sellmeierPreset(r.name, s)) { report(std::string(r.name) + ": preset exists", false); continue; }
        const Real nd = s.n(0.5876), nF = s.n(0.4861), nC = s.n(0.6563);
        const Real vd = (nd - 1) / (nF - nC);
        report(std::string(r.name) + ": n_d and Abbe number", std::abs(nd - r.nd) < 0.001 && std::abs(vd - r.vd) < 0.6,
            "n_d " + std::to_string(nd).substr(0, 6) + " (published " + std::to_string(r.nd).substr(0, 6) + "), V_d " + std::to_string(vd).substr(0, 5) + " (" + std::to_string(r.vd).substr(0, 5) + ")");
        bool decreasing = true;
        for (Real um = 0.38; um < 0.78; um += 0.01) decreasing = decreasing && (s.n(um) > s.n(um + 0.01));
        report(std::string(r.name) + ": the index falls with the wavelength over 380-780 nm", decreasing);
    }

    // (2) The material: the spectral build follows the table, the RGB build is constant at n_d; an Abbe number and Sellmeier terms are exclusive in the scene file.
    {
        Sellmeier s;
        sellmeierPreset("diamond", s);
        Dielectric d(1.5, 0.0);
        d.setSellmeier(s);
#if defined(RAYT_SPECTRAL_MODE)
        const bool ok = d.isDispersive() && std::abs(d.iorAtNm(450.0) - s.n(0.45)) < 1e-9 && d.iorAtNm(450.0) > d.iorAtNm(650.0);
        report("Dielectric with Sellmeier: dispersive, n(450 nm) > n(650 nm), as the table", ok);
#else
        report("Dielectric with Sellmeier: the RGB build uses n_d for every wavelength", !d.isDispersive() && std::abs(d.iorAtNm(450.0) - s.n(0.5876)) < 1e-9);
#endif
    }

    // (3) White furnace: a diamond ball in a uniform white environment shows radiance 1 (a lossless dispersive glass; every wavelength is conserved).
    {
        auto sky = std::make_shared<EnvMap>(Image(32, 16, std::vector<Vector3>(32 * 16, Vector3(1.0))));
        auto world = std::make_shared<Aggregate>();
        auto diamond = std::make_shared<Dielectric>(2.417, 0.0);
        Sellmeier s;
        sellmeierPreset("diamond", s);
        diamond->setSellmeier(s);
        world->add(std::make_shared<Sphere>(Point3(0, 0, 0), 1.0, diamond));
        auto scene = std::make_shared<Scene>(world);
        scene->setEnvLight(std::make_shared<EnvLight>(sky, 1.0));
        scene->finalize("spatial");
        Vector3 mean, se;
        int bad = 0;
        EstimateRadiance(*scene, Ray(Point3(0.3, 0.2, 4), glm::normalize(Vector3(-0.07, -0.05, -1))), 120000, 64, mean, se, bad);
        allPassed = ReportLightCase("a diamond ball (Sellmeier dispersion), white furnace", Vector3(1.0), mean, se, bad) && allPassed;
    }

    std::cout << "[Dispersion] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
