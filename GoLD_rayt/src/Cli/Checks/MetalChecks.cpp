// Metal spectral data checks.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {


// -----------------------------------------------------------------------------
// Metal data tests (--scene metals)
// -----------------------------------------------------------------------------

// Largest |difference| of eta and k between two tables on a wavelength grid (relative to max(1, value)).
static Real MaxNkDifference(const color::ComplexIorTable& a, const color::ComplexIorTable& b) {
    Real worst = 0;
    for (int l = 400; l <= 700; l += 25) {
        const auto [ea, ka] = a.evaluateNK(Real(l));
        const auto [eb, kb] = b.evaluateNK(Real(l));
        worst = std::max({ worst, std::abs(ea - eb) / std::max(Real(1), std::abs(ea)),
                                  std::abs(ka - kb) / std::max(Real(1), std::abs(ka)) });
    }
    return worst;
}


// Normal-incidence reflectance curve R(l) of a table, sampled every 5 nm over 300..900 nm.
static color::SpectralCurve ReflectanceCurveOf(const color::ComplexIorTable& t) {
    color::SpectralCurve c;
    for (int l = 300; l <= 900; l += 5) {
        const auto [e, k] = t.evaluateNK(Real(l));
        c.add(Real(l), color::conductorNormalIncidenceReflectance(e, k));
    }
    c.finalize();
    return c;
}

int RunMetalTests(const CliOptions& opt) {
    bool allPassed = true;
    auto report = [&](const std::string& what, Real value, Real limit) {
        const bool ok = value <= limit;
        std::cout << "[Metals]   " << what << ": " << value << " (limit " << limit << ")   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };

    color::ComplexIorTable ref;
    std::string err;
    if (!io::loadComplexIor_RefractiveIndexInfoCsv(AU_IOR_CSV_PATH, ref, &err)) {
        std::cerr << "[Metals] cannot load " << AU_IOR_CSV_PATH << ": " << err << "\n";
        return 1;
    }
    std::cout << "[Metals] reference: " << AU_IOR_CSV_PATH << " (" << ref.size() << " rows)\n";

    // ---- 1. loaders: the same data written in other layouts / units must give the same table ----
    const std::filesystem::path dir = std::filesystem::temp_directory_path();
    auto write = [&](const char* name, const std::function<void(std::ofstream&)>& body) {
        const std::string path = (dir / name).string();
        std::ofstream f(path);
        f.precision(12);
        body(f);
        return path;
    };
    const auto& rows = ref.etaCurve().samples();
    const auto& kRows = ref.kCurve().samples();

    const std::string pathNm = write("rayt_au_nm.csv", [&](std::ofstream& f) {
        f << "# wavelength_nm, n, k\n";
        for (size_t i = 0; i < rows.size(); ++i) f << rows[i].lambda_nm << "," << rows[i].value << "," << kRows[i].value << "\n";
    });
    const std::string pathUm = write("rayt_au_um.tsv", [&](std::ofstream& f) {
        f << "wl\tn\tk\n";
        for (size_t i = 0; i < rows.size(); ++i) f << rows[i].lambda_nm / 1000.0 << "\t" << rows[i].value << "\t" << kRows[i].value << "\n";
    });
    const std::string pathEv = write("rayt_au_ev.txt", [&](std::ofstream& f) {
        f << "eV n k\n";
        for (size_t i = rows.size(); i-- > 0;) f << 1239.841984 / rows[i].lambda_nm << " " << rows[i].value << " " << kRows[i].value << "\n";
    });
    const std::string pathN = write("rayt_au_n.csv", [&](std::ofstream& f) {
        for (const auto& s : rows) f << s.lambda_nm << "," << s.value << "\n";
    });
    const std::string pathK = write("rayt_au_k.csv", [&](std::ofstream& f) {
        for (const auto& s : kRows) f << s.lambda_nm << "," << s.value << "\n";
    });

    auto checkLoad = [&](const char* label, bool loaded, const color::ComplexIorTable& t, const std::string& e) {
        if (!loaded) { std::cout << "[Metals]   " << label << ": load failed (" << e << ")   FAIL\n"; allPassed = false; return; }
        report(std::string("loader ") + label + ", max |d(eta,k)| vs reference", MaxNkDifference(ref, t), 1e-6);
    };
    {
        color::ComplexIorTable t; std::string e;
        bool ok = io::loadComplexIor_ThreeColumn(pathNm, t, io::WavelengthUnit::Nm, &e);   checkLoad("3 columns, nm", ok, t, e);
        ok = io::loadComplexIor_ThreeColumn(pathUm, t, io::WavelengthUnit::Auto, &e);     checkLoad("3 columns, um (tab, auto)", ok, t, e);
        ok = io::loadComplexIor_ThreeColumn(pathEv, t, io::WavelengthUnit::Ev, &e);       checkLoad("3 columns, eV (space)", ok, t, e);
        ok = io::loadComplexIor_SeparateFiles(pathN, pathK, t, io::WavelengthUnit::Nm, &e); checkLoad("separate n / k files", ok, t, e);
        ok = io::loadComplexIor(pathUm, t, io::WavelengthUnit::Auto, &e);                 checkLoad("auto-detect (3 columns)", ok, t, e);
        ok = io::loadComplexIor(AU_IOR_CSV_PATH, t, io::WavelengthUnit::Auto, &e);        checkLoad("auto-detect (RefractiveIndex.info)", ok, t, e);
    }

    // ---- 2. reflectance input: F0 must be reproduced exactly ----
    const color::SpectralCurve R = ReflectanceCurveOf(ref);
    const color::ComplexIorTable fromR = color::complexIorTableFromReflectance(R);
    {
        Real worst = 0;
        for (int l = 360; l <= 830; l += 7) {   // off the 5 nm grid on purpose (interpolation error included)
            const auto [e, k] = fromR.evaluateNK(Real(l));
            worst = std::max(worst, std::abs(color::conductorNormalIncidenceReflectance(e, k) - R.evaluate(Real(l))));
        }
        report("reflectance -> (eta,k): max |F0 - R|, 360..830 nm", worst, 2e-3);

        // percent-valued CSV gives the same table
        const std::string pathR = write("rayt_au_R_percent.csv", [&](std::ofstream& f) {
            f << "wavelength_nm,reflectance_percent\n";
            for (const auto& s : R.samples()) f << s.lambda_nm << "," << s.value * 100.0 << "\n";
        });
        color::SpectralCurve Rloaded;
        if (io::loadReflectanceCsv(pathR, Rloaded, &err)) {
            report("reflectance CSV in percent: max |R - R_loaded|", [&] {
                Real w = 0; for (int l = 400; l <= 700; l += 10) w = std::max(w, std::abs(R.evaluate(Real(l)) - Rloaded.evaluate(Real(l)))); return w; }(), 1e-6);
        }
        else { std::cout << "[Metals]   reflectance CSV load failed: " << err << "   FAIL\n"; allPassed = false; }
        std::error_code ec;
        for (const std::string& p : { pathNm, pathUm, pathEv, pathN, pathK, pathR }) std::filesystem::remove(p, ec);
    }

    // ---- 3. angular error of the reflectance-only model against the true (eta, k) ----
    //   (this is the approximation the edge tint stands for; thresholds are the calibrated values + margin)
    for (const char* metal : { "Au", "Ag", "Cu" }) {
        color::ComplexIorTable truth;
        const std::string path = std::string("assets/Metal_data/Johnson_") + metal + ".csv";
        if (!io::loadComplexIor(path, truth, io::WavelengthUnit::Auto, &err)) {
            std::cout << "[Metals]   cannot load " << path << ": " << err << "   FAIL\n";
            allPassed = false;
            continue;
        }
        const color::ComplexIorTable approx = color::complexIorTableFromReflectance(ReflectanceCurveOf(truth));
        std::cout << "[Metals]   " << metal << " max |F(theta) - F_true(theta)| over 400..700 nm:";
        struct { int deg; Real limit; } checks[] = { {0, 0.003}, {45, 0.03}, {60, 0.06}, {70, 0.10}, {80, 0.20} };
        bool ok = true;
        for (const auto& c : checks) {
            Real worst = 0;
            const Real cosT = std::cos(glm::radians(Real(c.deg)));
            for (int l = 400; l <= 700; l += 10) {
                const auto [e1, k1] = truth.evaluateNK(Real(l));
                const auto [e2, k2] = approx.evaluateNK(Real(l));
                worst = std::max(worst, std::abs(fresnel::fresnelConductor(cosT, e1, k1) - fresnel::fresnelConductor(cosT, e2, k2)));
            }
            std::cout << "  " << c.deg << "deg " << worst;
            ok = ok && worst <= c.limit;
        }
        std::cout << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    }

    // ---- 4. through the renderer: true gold vs gold rebuilt from its reflectance ----
    //   A rough-conductor sphere in a uniform white environment; one ray at the center
    //   (near-normal incidence) and one near the rim (about 70 deg).
    {
        auto envMap = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        auto makeScene = [&](const color::ComplexIorTable& t) {
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Sphere>(Point3(0, 0, 0), 1.0,
                std::make_shared<RoughConductor>(std::make_shared<color::ComplexIorTable>(t), 0.2)));
            auto scene = std::make_shared<Scene>(world);
            scene->setEnvLight(std::make_shared<EnvLight>(envMap, 1.0));
            scene->finalize();
            return scene;
        };
        const auto sceneTrue = makeScene(ref);
        const auto sceneFromR = makeScene(color::complexIorTableFromReflectance(ReflectanceCurveOf(ref)));
        const int n = opt.spp * 1024;
        for (const Real x : { Real(0), Real(0.94) }) {
            const Ray ray(Point3(x, 0, 5), Vector3(0, 0, -1));
            Vector3 m1, s1, m2, s2;
            int bad = 0;
            EstimateRadiance(*sceneTrue, ray, n, 8, m1, s1, bad);
            EstimateRadiance(*sceneFromR, ray, n, 8, m2, s2, bad);
            const Vector3 d = glm::abs(m1 - m2);
            const Real limit = (x == 0) ? Real(0.01) : Real(0.15);   // rim: known edge-tint approximation, regression guard
            const Real worst = std::max({ d.x, d.y, d.z });
            std::cout << "[Metals]   rendered gold, " << (x == 0 ? "center (~0 deg)" : "rim (~70 deg)")
                << ": true (" << m1.x << ", " << m1.y << ", " << m1.z << ")  from R (" << m2.x << ", " << m2.y << ", " << m2.z
                << ")  max diff " << worst << " (limit " << limit << ")   " << (worst <= limit + 4 * std::max({ s1.x, s2.x }) ? "OK" : "FAIL") << "\n";
            allPassed = allPassed && (worst <= limit + 4 * std::max({ s1.x, s2.x }));
        }
    }

    // ---- 5. multiple-scattering compensation (Kulla-Conty) ----
    {
        // 5a. near-perfect conductor in a uniform white environment: must return all the light.
        auto envMap = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        auto perfect = std::make_shared<color::ComplexIorTable>();
        perfect->assign({ { 360.0, 1.0, 100.0 }, { 830.0, 1.0, 100.0 } });
        const int n = opt.spp * 512;
        for (const Real roughness : { Real(0.2), Real(0.5), Real(0.8), Real(1.0) }) {
            std::string line = "[Metals]   perfect conductor, roughness " + std::to_string(roughness).substr(0, 3) + ": albedo (center | rim)";
            bool ok = true;
            for (const bool ms : { false, true }) {
                auto mat = std::make_shared<RoughConductor>(perfect, roughness);
                mat->setMultipleScattering(ms);
                auto world = std::make_shared<Aggregate>();
                world->add(std::make_shared<Sphere>(Point3(0, 0, 0), 1.0, mat));
                auto scene = std::make_shared<Scene>(world);
                scene->setEnvLight(std::make_shared<EnvLight>(envMap, 1.0));
                scene->finalize();
                Vector3 mc, sc, mr, sr;
                int bad = 0;
                EstimateRadiance(*scene, Ray(Point3(0, 0, 5), Vector3(0, 0, -1)), n, 16, mc, sc, bad);
                EstimateRadiance(*scene, Ray(Point3(0.94, 0, 5), Vector3(0, 0, -1)), n, 16, mr, sr, bad);
                line += std::string(ms ? "   with MS: " : "   single: ") + std::to_string(mc.x).substr(0, 5) + " | " + std::to_string(mr.x).substr(0, 5);
                if (ms) ok = std::abs(mc.x - 1.0) <= 4 * sc.x + 0.02 && std::abs(mr.x - 1.0) <= 4 * sr.x + 0.03 && bad == 0;
            }
            std::cout << line << "   " << (ok ? "OK" : "FAIL") << "\n";
            allPassed = allPassed && ok;
        }

        // 5b. gold, roughness 0.6: reciprocity, pdf normalization, and sample() vs uniform integration.
        auto gold = std::make_shared<RoughConductor>(std::make_shared<color::ComplexIorTable>(ref), 0.6);
        SurfaceInteraction rec;
        rec.p = Point3(0.0); rec.n = Vector3(0, 0, 1); rec.gn = Vector3(0, 0, 1);
#if defined(RAYT_SPECTRAL_MODE)
        const SampledWavelengths lambda0 = SampledWavelengths::sampleVisible(Real(0.37));   // fixed wavelengths: a meaningful, repeatable albedo
#else
        const SampledWavelengths lambda0;
#endif
        auto randomDir = [&]() {
            const Real z = sampling::Random(0.05, 1.0), phi = 2 * constants::PI * sampling::Random();
            const Real r = std::sqrt(1 - z * z);
            return Vector3(r * std::cos(phi), r * std::sin(phi), z);
        };
        Real worstRecip = 0;
        for (int i = 0; i < 2000; ++i) {
            const Vector3 a = randomDir(), b = randomDir();
            const Spectrum f1 = gold->eval(rec, a, b, lambda0, TransportMode::Radiance);
            const Spectrum f2 = gold->eval(rec, b, a, lambda0, TransportMode::Radiance);
            const Vector3 v1 = SpectrumToRGB(f1, lambda0), v2 = SpectrumToRGB(f2, lambda0);
            worstRecip = std::max(worstRecip, std::abs(v1.x - v2.x) / std::max(Real(1e-9), std::max(v1.x, v2.x)));
        }
        report("gold r=0.6: reciprocity f(wo,wi) = f(wi,wo), max relative difference", worstRecip, 1e-6);

        for (const Real cosO : { Real(0.9), Real(0.5), Real(0.15) }) {
            const Vector3 wo(std::sqrt(1 - cosO * cosO), 0, cosO);
            double pdfInt = 0, albedoU = 0, albedoS = 0, sqU = 0, sqS = 0;
            const int m = 200000;
            for (int i = 0; i < m; ++i) {
                // uniform hemisphere: integral of pdf, and integral of f cos
                const Real z = sampling::Random(), phi = 2 * constants::PI * sampling::Random(), r = std::sqrt(std::max(Real(0), 1 - z * z));
                const Vector3 wi(r * std::cos(phi), r * std::sin(phi), z);
                pdfInt += gold->pdf(rec, wo, wi, lambda0) * 2 * constants::PI;
                const double u = SpectrumToRGB(gold->eval(rec, wo, wi, lambda0, TransportMode::Radiance), lambda0).x * wi.z * 2 * constants::PI;
                albedoU += u; sqU += u * u;
                // importance sampled through the material's own sample()
                const auto s = gold->sample(rec, wo, Point2(sampling::Random(), sampling::Random()), lambda0, TransportMode::Radiance);
                const double w = s ? SpectrumToRGB(s->f, lambda0).x * std::abs(glm::dot(rec.n, s->wi)) / s->pdf : 0.0;
                albedoS += w; sqS += w * w;
            }
            pdfInt /= m; albedoU /= m; albedoS /= m;
            const double seU = std::sqrt(std::max(0.0, sqU / m - albedoU * albedoU) / m), seS = std::sqrt(std::max(0.0, sqS / m - albedoS * albedoS) / m);
            // The pdf integrates to <= 1: GGX reflections that leave the hemisphere are lost (sample() fails for them).
            const bool ok = pdfInt <= 1.02 && pdfInt > 0.8 && std::abs(albedoU - albedoS) <= 4 * std::hypot(seU, seS) + 1e-3;
            std::cout << "[Metals]   gold r=0.6, cos(theta_o)=" << cosO << ": integral of pdf = " << pdfInt << " (<= 1), albedo (R) by uniform integration "
                << albedoU << " vs by sample() " << albedoS << "   " << (ok ? "OK" : "FAIL") << "\n";
            allPassed = allPassed && ok;
        }
    }

    // ---- 6. CoatedDiffuse (diffuse under a GGX coat, metallic blend) ----
    {
        SurfaceInteraction rec;
        rec.p = Point3(0.0); rec.n = Vector3(0, 0, 1); rec.gn = Vector3(0, 0, 1);
#if defined(RAYT_SPECTRAL_MODE)
        const SampledWavelengths lam = SampledWavelengths::sampleVisible(Real(0.37));
#else
        const SampledWavelengths lam;
#endif
        auto randomDir = [&]() {
            const Real z = sampling::Random(0.05, 1.0), phi = 2 * constants::PI * sampling::Random();
            const Real r = std::sqrt(1 - z * z);
            return Vector3(r * std::cos(phi), r * std::sin(phi), z);
        };
        struct Case { const char* name; Real albedo, roughness, weight, metallic; };
        for (const Case& c : { Case{ "plastic r=0.3", 0.8, 0.3, 1.0, 0.0 }, Case{ "plastic r=0.7", 0.8, 0.7, 1.0, 0.0 },
                               Case{ "half coat, 30% metal r=0.4", 0.8, 0.4, 0.5, 0.3 }, Case{ "metal r=0.4", 0.8, 0.4, 1.0, 1.0 } }) {
            auto mat = std::make_shared<CoatedDiffuse>(RgbColor(c.albedo), c.roughness, c.weight, c.metallic);
            Real worstRecip = 0;
            for (int i = 0; i < 2000; ++i) {
                const Vector3 a = randomDir(), b = randomDir();
                const Vector3 v1 = SpectrumToRGB(mat->eval(rec, a, b, lam, TransportMode::Radiance), lam);
                const Vector3 v2 = SpectrumToRGB(mat->eval(rec, b, a, lam, TransportMode::Radiance), lam);
                worstRecip = std::max(worstRecip, std::abs(v1.x - v2.x) / std::max(Real(1e-9), std::max(v1.x, v2.x)));
            }
            bool ok = worstRecip < 1e-6;
            std::ostringstream line;
            line << "[Metals]   CoatedDiffuse " << c.name << ": reciprocity " << worstRecip;
            for (const Real cosO : { Real(0.8), Real(0.3) }) {
                const Vector3 wo(std::sqrt(1 - cosO * cosO), 0, cosO);
                double albedoU = 0, albedoS = 0, sqU = 0, sqS = 0, pdfInt = 0;
                const int m = 200000;
                for (int i = 0; i < m; ++i) {
                    const Real z = sampling::Random(), phi = 2 * constants::PI * sampling::Random(), r = std::sqrt(std::max(Real(0), 1 - z * z));
                    const Vector3 wi(r * std::cos(phi), r * std::sin(phi), z);
                    pdfInt += mat->pdf(rec, wo, wi, lam) * 2 * constants::PI;
                    const double u = SpectrumToRGB(mat->eval(rec, wo, wi, lam, TransportMode::Radiance), lam).x * wi.z * 2 * constants::PI;
                    albedoU += u; sqU += u * u;
                    const auto s = mat->sample(rec, wo, Point2(sampling::Random(), sampling::Random()), lam, TransportMode::Radiance);
                    const double w = s ? SpectrumToRGB(s->f, lam).x * std::abs(glm::dot(rec.n, s->wi)) / s->pdf : 0.0;
                    albedoS += w; sqS += w * w;
                }
                pdfInt /= m; albedoU /= m; albedoS /= m;
                const double seU = std::sqrt(std::max(0.0, sqU / m - albedoU * albedoU) / m), seS = std::sqrt(std::max(0.0, sqS / m - albedoS * albedoS) / m);
                const bool okHere = pdfInt <= 1.02 && std::abs(albedoU - albedoS) <= 4 * std::hypot(seU, seS) + 1e-3;
                ok = ok && okHere;
                line << ";  cos " << cosO << ": albedo " << albedoU << " vs sample() " << albedoS << ", pdf integral " << pdfInt;
            }
            std::cout << line.str() << "   " << (ok ? "OK" : "FAIL") << "\n";
            allPassed = allPassed && ok;
        }

        // A white (albedo 1) plastic sphere in a uniform white environment can only return <= 1 of the light.
        auto envMap = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        for (const Real rough : { Real(0.2), Real(0.5) }) {
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Sphere>(Point3(0, 0, 0), 1.0, std::make_shared<CoatedDiffuse>(RgbColor(1.0), rough, 1.0, 0.0)));
            auto scene = std::make_shared<Scene>(world);
            scene->setEnvLight(std::make_shared<EnvLight>(envMap, 1.0));
            scene->finalize();
            Vector3 mc, sc;
            int bad = 0;
            EstimateRadiance(*scene, Ray(Point3(0, 0, 5), Vector3(0, 0, -1)), opt.spp * 512, 16, mc, sc, bad);
            const bool ok = mc.x <= 1.0 + 4 * sc.x + 0.005 && mc.x > 0.9 && bad == 0;
            std::cout << "[Metals]   white plastic r=" << rough << " in a white furnace: " << mc.x << " (<= 1)   " << (ok ? "OK" : "FAIL") << "\n";
            allPassed = allPassed && ok;
        }
    }

    // ---- 5. a fading series R(lambda, t): three rounds at 0, 30 and 150 h, curves that are linear in wavelength, so the exact answer is known ----
    {
        const std::filesystem::path seriesDir = dir / "rayt_fading_series";
        std::error_code ec;
        std::filesystem::remove_all(seriesDir, ec);
        const double hours[3] = { 0, 30, 150 };
        auto truth = [](int round, Real nm) { return Real(0.9) - Real(0.2) * Real(round - 1) + Real(0.05) * Real(round) * (nm - Real(560)) / Real(120); };
        std::filesystem::create_directories(seriesDir);
        {
            std::ofstream f((seriesDir / "hours.csv").string());
            f << "round,hours\n";
            for (int r = 0; r < 3; ++r) f << r + 1 << "," << hours[r] << "\n";
        }
        for (int r = 1; r <= 3; ++r) {
            std::filesystem::create_directories(seriesDir / ("r" + std::to_string(r)));
            std::ofstream f((seriesDir / ("r" + std::to_string(r)) / "red.csv").string());
            f.precision(12);
            f << "# test\nwavelength_nm,reflectance\n";
            for (int nm = 440; nm <= 680; nm += 20) f << nm << "," << truth(r, Real(nm)) << "\n";
        }
        // expected curve: weights (w1, w2, w3) of the three rounds
        auto worstAgainst = [&](Real h, Real w1, Real w2, Real w3, bool wantClamped) {
            io::FadingCurve fc; std::string e;
            if (!io::loadFadingCurve(seriesDir.string(), "red", h, fc, &e)) { std::cout << "[Metals]   fading series load failed: " << e << "   FAIL\n"; allPassed = false; return Real(1); }
            if (fc.clamped != wantClamped) { std::cout << "[Metals]   fading series: clamped flag wrong at " << h << " h   FAIL\n"; allPassed = false; }
            Real worst = 0;
            for (int nm = 440; nm <= 680; nm += 10)   // between the stored wavelengths too
                worst = std::max(worst, std::abs(fc.curve.evaluate(Real(nm)) - (w1 * truth(1, Real(nm)) + w2 * truth(2, Real(nm)) + w3 * truth(3, Real(nm)))));
            return worst;
        };
        report("fading series: exactly at the second round (30 h)", worstAgainst(30, 0, 1, 0, false), 1e-9);
        report("fading series: halfway between round 1 and 2 (15 h)", worstAgainst(15, 0.5, 0.5, 0, false), 1e-9);
        report("fading series: a quarter of the way from round 2 to 3 (60 h; the spacing is not even)", worstAgainst(60, 0, 0.75, 0.25, false), 1e-9);
        report("fading series: the first round (0 h)", worstAgainst(0, 1, 0, 0, false), 1e-9);
        report("fading series: before the first round (-5 h) is held, and flagged", worstAgainst(-5, 1, 0, 0, true), 1e-9);
        report("fading series: after the last round (400 h) is held, and flagged", worstAgainst(400, 0, 0, 1, true), 1e-9);
        std::filesystem::remove_all(seriesDir, ec);
    }

    std::cout << "[Metals] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
