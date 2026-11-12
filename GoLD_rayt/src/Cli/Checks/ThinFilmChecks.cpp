// Thin-film interference checks (--scene thinfilm): Airy's formula against closed forms, and the material in a furnace.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Geometry/Sphere.hpp"
#include "Materials/ThinFilm.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Integrator.hpp"
#include "Core/Fresnel.hpp"
#include "Core/Sampling.hpp"

namespace rayt::cli {

int RunThinFilmTests() {
    using thinfilm::Complex;
    constexpr Real kPi = Real(3.14159265358979323846);
    bool allPassed = true;
    auto report = [&](const std::string& name, bool ok, const std::string& detail = "") {
        std::cout << "[ThinFilm]   " << name << (detail.empty() ? "" : "  " + detail) << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };
    std::cout << "[ThinFilm] Airy's formula and the thin-film material\n";

    // (1) No film (d = 0): the plain Fresnel reflectance of the interface, for a dielectric and for a metal substrate.
    {
        Real worstD = 0, worstM = 0;
        std::shared_ptr<color::ComplexIorTable> au;   // gold at 550 nm, read with the engine's own table
        for (const Real cosI : { Real(1), Real(0.8), Real(0.5), Real(0.2), Real(0.05) }) {
            const Real a = thinfilm::reflectance(cosI, 1.0, 1.4, 0.0, Complex(1.5, 0), 550.0);
            worstD = std::max(worstD, std::abs(a - fresnel::fresnelDielectric(cosI, 1.0, 1.5)));
            const Real eta = 0.28, k = 2.9;   // a metal-like complex index
            const Real b = thinfilm::reflectance(cosI, 1.0, 1.4, 0.0, Complex(eta, k), 550.0);
            worstM = std::max(worstM, std::abs(b - fresnel::fresnelConductor(cosI, eta, k)));
        }
        report("d = 0: the plain Fresnel reflectance (dielectric substrate), 5 angles", worstD < 1e-9, "max difference " + std::to_string(worstD));
        report("d = 0: the plain Fresnel reflectance (metal substrate), 5 angles", worstM < 1e-9, "max difference " + std::to_string(worstM));
    }

    // (2) A quarter-wave coating with nF = sqrt(n2) on glass reflects nothing at its design wavelength and normal incidence.
    {
        const Real n2 = 1.5, nF = std::sqrt(n2), lambda0 = 550.0, d = lambda0 / (4.0 * nF);
        const Real r = thinfilm::reflectance(1.0, 1.0, nF, d, Complex(n2, 0), lambda0);
        const Real bare = fresnel::fresnelDielectric(1.0, 1.0, n2);
        report("quarter-wave coating: R = 0 at the design wavelength", r < 1e-12, "R = " + std::to_string(r) + " (bare glass " + std::to_string(bare) + ")");
    }

    // (3) A free film in air at normal incidence: R = F sin^2(beta) / (1 + F sin^2(beta)),  F = 4 R0 / (1 - R0)^2,  R0 = ((n - 1) / (n + 1))^2.
    {
        Real worst = 0;
        const Real n = 1.33;
        const Real R0 = std::pow((n - 1) / (n + 1), Real(2)), F = 4 * R0 / std::pow(Real(1) - R0, Real(2));
        for (const Real d : { Real(80), Real(250), Real(520), Real(1100) })
            for (Real lam = 400; lam <= 700; lam += 25) {
                const Real beta = 2 * kPi * n * d / lam;
                const Real s2 = std::sin(beta) * std::sin(beta);
                worst = std::max(worst, std::abs(thinfilm::reflectance(1.0, 1.0, n, d, Complex(1, 0), lam) - F * s2 / (1 + F * s2)));
            }
        report("free film, normal incidence: Airy's closed form (4 thicknesses x 13 wavelengths)", worst < 1e-9, "max difference " + std::to_string(worst));
    }

    // (4) A thick film averages to the incoherent sum (no interference left): free film in air, R_inc = 2 R0 / (1 + R0).
    {
        const Real n = 1.33, R0 = std::pow((n - 1) / (n + 1), Real(2)), inc = 2 * R0 / (1 + R0);
        Real sum = 0;
        int cnt = 0;
        for (Real lam = 400; lam <= 700; lam += 0.25) { sum += thinfilm::reflectance(1.0, 1.0, n, 40000.0, Complex(1, 0), lam); ++cnt; }
        const Real mean = sum / Real(cnt);
        report("a thick film (40 um): the mean over 400-700 nm is the incoherent sum", std::abs(mean - inc) < 0.01 * inc,
            "mean " + std::to_string(mean) + ", incoherent " + std::to_string(inc));
    }

    // (5) R stays in [0, 1], also for a film on a metal and at grazing angles; a film in glass (incidence from the dense side) has total internal reflection.
    {
        Real lo = 1, hi = 0;
        for (const Real cosI : { Real(1), Real(0.5), Real(0.1), Real(0.01) })
            for (const Real d : { Real(0), Real(100), Real(370) })
                for (Real lam = 400; lam <= 700; lam += 30) {
                    for (const Complex nT : { Complex(1.33, 0), Complex(0.2, 3.4) }) {
                        const Real r = thinfilm::reflectance(cosI, 1.0, 1.45, d, nT, lam);
                        lo = std::min(lo, r); hi = std::max(hi, r);
                    }
                }
        const Real tir = thinfilm::reflectance(0.2, 1.5, 1.33, 200.0, Complex(1.0, 0), 550.0);   // glass -> film -> air at 78 degrees: beyond the critical angle
        report("R in [0, 1] over a grid (dielectric and metal substrates)", lo >= 0 && hi <= 1.0000001, "range " + std::to_string(lo) + " .. " + std::to_string(hi));
        report("total internal reflection from the dense side: R = 1", std::abs(tir - 1.0) < 1e-9, "R = " + std::to_string(tir));
    }

    // (6) The material's sample(): E[weight of the reflection samples] = R and E[weight of the transmission samples] = T, per wavelength.
    {
        FilmThickness th;
        th.base = 330;
        const ThinFilm oil(1.45, th, 1.33);
        SurfaceInteraction rec;
        rec.p = Point3(0, 0, 0); rec.n = Vector3(0, 1, 0); rec.gn = rec.n; rec.frontFace = true;
        const Real cosI = 0.6;
        const Vector3 wo(std::sqrt(1 - cosI * cosI), cosI, 0);
        SampledWavelengths lambda;
#if defined(RAYT_SPECTRAL_MODE)
        const Real wl[4] = { 450, 520, 580, 650 };
        for (int i = 0; i < kSpectralSamples; ++i) lambda.lambda[i] = wl[i];
        constexpr int lanes = kSpectralSamples;
#else
        constexpr int lanes = 3;
#endif
        Spectrum R, T;
        oil.reflectTransmit(rec.p, cosI, true, lambda, R, T);
        const int n = 400000;
        Spectrum sumR(0.0), sumT(0.0);
        for (int s = 0; s < n; ++s) {
            const auto smp = oil.sample(rec, wo, Point2(sampling::Random(), 0.0), lambda, TransportMode::Importance);
            if (!smp) continue;
            if (smp->isReflection()) sumR += smp->f; else sumT += smp->f;
        }
        Real worst = 0;
        std::string line;
        for (int i = 0; i < lanes; ++i) {
            const Real r = sumR[i] / Real(n), t = sumT[i] / Real(n);
            worst = std::max({ worst, std::abs(r - R[i]), std::abs(t - T[i]) });
            line += " R" + std::to_string(i) + " " + std::to_string(R[i]).substr(0, 5) + "/" + std::to_string(r).substr(0, 5);
        }
        report("oil on water: sample() weights average to R and T per wavelength", worst < 0.004, "expected/estimated:" + line);
    }

    // (7) Energy: a free film as a closed sphere in a uniform white environment shows radiance 1 whatever its thickness (R + T = 1 at both walls); so does a coated
    // glass ball (every transmission is refracted back out). Wavelengths oscillate, but the sum over all paths is conserved.
    {
        auto sky = std::make_shared<EnvMap>(Image(32, 16, std::vector<Vector3>(32 * 16, Vector3(1.0))));
        auto furnace = [&](const std::shared_ptr<Material>& m, const std::string& name) {
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Sphere>(Point3(0, 0, 0), 1.0, m));
            auto scene = std::make_shared<Scene>(world);
            scene->setEnvLight(std::make_shared<EnvLight>(sky, 1.0));
            scene->finalize("spatial");
            Vector3 mean, se;
            int bad = 0;
            EstimateRadiance(*scene, Ray(Point3(0.3, 0.2, 4), glm::normalize(Vector3(-0.07, -0.05, -1))), 60000, 24, mean, se, bad);
            allPassed = ReportLightCase(name, Vector3(1.0), mean, se, bad) && allPassed;
        };
        FilmThickness th;
        th.base = 420;
        furnace(std::make_shared<ThinFilm>(1.33, th), "free film (soap bubble), 420 nm, white furnace");
        th.base = 150;
        furnace(std::make_shared<ThinFilm>(1.45, th, 1.33), "oil film on a water ball, 150 nm, white furnace");
    }

    std::cout << "[ThinFilm] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
