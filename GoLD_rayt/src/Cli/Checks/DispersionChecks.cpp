// Dispersion checks (--scene dispersion): the Sellmeier presets against their published d-line index and Abbe number, and a diamond in a white furnace.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/SceneParts.hpp"
#include "Geometry/Quad.hpp"
#include "Geometry/Sphere.hpp"
#include "IO/EnvMap.hpp"
#include "Lights/EnvLight.hpp"
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

    // (4) Nested dielectrics: a surface carries the index of its outside too (Dielectric::setExteriorIor). A glass slab (n = 1.5) in air and in water (1.333).
    std::cout << "[Dispersion] nested dielectrics: the exterior index of a glass surface\n";
    {
#if defined(RAYT_SPECTRAL_MODE)
        const SampledWavelengths lambda = SampledWavelengths::sampleVisible(0.5);
#else
        const SampledWavelengths lambda;
#endif
        // a. Snell's law with the relative index: the transmitted direction of a smooth surface, sampled directly from the BSDF.
        auto refractedDeg = [&](Real exterior, bool fromInside, Real incidenceDeg, bool& transmitted) {
            Dielectric d(1.5, 0.0);
            d.setExteriorIor(exterior);
            SurfaceInteraction rec;
            rec.p = Point3(0.0); rec.n = Vector3(0, 0, 1); rec.gn = rec.n; rec.frontFace = !fromInside; rec.matPtr = &d;   // rec.n faces the incoming ray
            const Real th = glm::radians(incidenceDeg);
            const Vector3 wo(std::sin(th), 0, std::cos(th));
            transmitted = false;
            Real angle = 0;
            for (int k = 0; k < 64 && !transmitted; ++k) {   // the sample picks reflection or transmission by u: scan u until the transmission shows
                const auto s = d.sample(rec, wo, Point2((k + Real(0.5)) / 64, Real(0.5)), lambda, TransportMode::Radiance);
                if (s && s->isTransmission()) { transmitted = true; angle = glm::degrees(std::acos(std::clamp(-s->wi.z, Real(-1), Real(1)))); }
            }
            return angle;
        };
        bool tAir = false, tWater = false, tirAir = false, tirWater = false;
        const Real aAir = refractedDeg(1.0, false, 45.0, tAir), aWater = refractedDeg(1.333, false, 45.0, tWater);
        const Real wantAir = glm::degrees(std::asin(std::sin(glm::radians(45.0)) / 1.5)), wantWater = glm::degrees(std::asin(std::sin(glm::radians(45.0)) * 1.333 / 1.5));
        report("entering at 45 deg: refracted angle in air " + std::to_string(aAir).substr(0, 5) + " deg (Snell " + std::to_string(wantAir).substr(0, 5)
            + "), with water outside " + std::to_string(aWater).substr(0, 5) + " (" + std::to_string(wantWater).substr(0, 5) + ")",
            tAir && tWater && std::abs(aAir - wantAir) < 0.01 && std::abs(aWater - wantWater) < 0.01);
        // leaving at 50 deg from inside the glass: total internal reflection against air (critical 41.8 deg), transmission against water (critical 62.7 deg)
        const Real lAir = refractedDeg(1.0, true, 50.0, tirAir), lWater = refractedDeg(1.333, true, 50.0, tirWater);
        const Real wantLeave = glm::degrees(std::asin(std::sin(glm::radians(50.0)) * 1.5 / 1.333));
        (void)lAir;
        report("leaving at 50 deg: total internal reflection against air, transmission against water at " + std::to_string(lWater).substr(0, 5) + " deg (" + std::to_string(wantLeave).substr(0, 5) + ")",
            !tirAir && tirWater && std::abs(lWater - wantLeave) < 0.01);

        // b. A slab seen at normal incidence in front of an emitter of radiance 1 (black elsewhere): the transmittance with every internal
        //    bounce is (1 - R) / (1 + R), R the Fresnel reflectance of one interface: 0.9231 against air (R = 0.04), 0.9931 against water.
        for (const Real exterior : { Real(1.0), Real(1.333) }) {
            auto glass = std::make_shared<Dielectric>(1.5, 0.0);
            glass->setExteriorIor(exterior);
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Quad>(Point3(-0.5, -0.5, 0.005), Vector3(1, 0, 0), Vector3(0, 1, 0), glass));    // top, faces +z
            world->add(std::make_shared<Quad>(Point3(-0.5, -0.5, -0.005), Vector3(0, 1, 0), Vector3(1, 0, 0), glass));   // bottom, faces -z
            auto scene = std::make_shared<Scene>(world);
            AddQuadLight(*world, scene.get(), Point3(-2, -2, -1), Vector3(4, 0, 0), Vector3(0, 4, 0), EmissionSpectrum(RgbColor(1.0)), 1.0, true);
            scene->finalize("spatial");
            Vector3 mean, se;
            int bad = 0;
            EstimateRadiance(*scene, Ray(Point3(0, 0, 1), Vector3(0, 0, -1)), 120000, 64, mean, se, bad);
            const Real R = math::sqr((1.5 - exterior) / (1.5 + exterior));
            const Real T = (1 - R) / (1 + R);
            allPassed = ReportLightCase("glass slab at normal incidence, emitter behind, outside index " + std::to_string(exterior).substr(0, 5) + ": transmittance", Vector3(T), mean, se, bad) && allPassed;
        }

        // c. The same slab with water outside in a white furnace: a lossless interface pair keeps the energy whatever the ratio of the indices.
        {
            auto sky = std::make_shared<EnvMap>(Image(32, 16, std::vector<Vector3>(32 * 16, Vector3(1.0))));
            auto glass = std::make_shared<Dielectric>(1.5, 0.0);
            glass->setExteriorIor(1.333);
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Quad>(Point3(-0.5, -0.5, 0.005), Vector3(1, 0, 0), Vector3(0, 1, 0), glass));
            world->add(std::make_shared<Quad>(Point3(-0.5, -0.5, -0.005), Vector3(0, 1, 0), Vector3(1, 0, 0), glass));
            auto scene = std::make_shared<Scene>(world);
            scene->setEnvLight(std::make_shared<EnvLight>(sky, 1.0));
            scene->finalize("spatial");
            Vector3 mean, se;
            int bad = 0;
            EstimateRadiance(*scene, Ray(Point3(0.1, 0.1, 1), glm::normalize(Vector3(-0.3, 0.1, -1))), 60000, 64, mean, se, bad);
            allPassed = ReportLightCase("glass slab with water outside, white furnace", Vector3(1.0), mean, se, bad) && allPassed;
        }
    }

    std::cout << "[Dispersion] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
