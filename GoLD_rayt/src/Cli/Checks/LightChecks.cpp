// Light source checks.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/SceneParts.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

int RunLightTests(const CliOptions& opt) {
    const int n = opt.spp * 4096;   // samples per test (default --spp 256 -> about 1M)
    const int depth = std::max(opt.maxDepth, 4);
    std::cout << "[Lights] " << n << " samples per test, albedo 0.5 ground, probe ray hits the origin\n";
    bool allPassed = true;
    const Real albedo = 0.5;

    const Ray probe(Point3(-3, 3, 0), glm::normalize(Vector3(1, -1, 0)));
    auto run = [&](const std::string& name, const std::shared_ptr<Scene>& scene, const Vector3& expected) {
        Vector3 mean, se;
        int bad = 0;
        EstimateRadiance(*scene, probe, n, depth, mean, se, bad);
        allPassed = ReportLightCase(name, expected, mean, se, bad) && allPassed;
    };

    // Overhead light at distance d = 2 from the probe point.
    // Point light:   Lo = albedo / pi * I / d^2
    run("point light I=10, d=2",
        MakeGroundScene(albedo, [](Aggregate&, Scene& s) {
            s.addLight(std::make_shared<PointLight>(Point3(0, 2, 0), RgbColor(1.0), 10.0)); }),
        Vector3(albedo / constants::PI * 10.0 / 4.0));

    // Sphere light (radius R, radiance L, center distance d, angle theta to the normal):
    //   E = pi L (R/d)^2 cos(theta)   ->   Lo = albedo L (R/d)^2 cos(theta)
    const Real overhead = albedo * 10.0 * (0.5 * 0.5) / 4.0;                       // theta = 0, d = 2
    const Real oblique = albedo * 10.0 * (0.25 / 4.5) * (1.5 / std::sqrt(4.5));    // center (1.5, 1.5, 0)
    for (int variant = 0; variant < 2; ++variant) {
        const bool nee = (variant == 0);
        const std::string tag = nee ? "NEE + MIS" : "BSDF sampling only";
        run("sphere light overhead, R=0.5 L=10  [" + tag + "]",
            MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                AddSphereLight(w, s, Point3(0, 2, 0), 0.5, RgbColor(1.0), 10.0, nee); }),
            Vector3(overhead));
        run("sphere light oblique,  R=0.5 L=10  [" + tag + "]",
            MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                AddSphereLight(w, s, Point3(1.5, 1.5, 0), 0.5, RgbColor(1.0), 10.0, nee); }),
            Vector3(oblique));
        run("sphere light overhead, big R=1.2 L=10  [" + tag + "]",
            MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                AddSphereLight(w, s, Point3(0, 2.5, 0), 1.2, RgbColor(1.0), 10.0, nee); }),
            Vector3(albedo * 10.0 * (1.2 * 1.2) / (2.5 * 2.5)));
        run("colored sphere light (1, .4, .1)*10  [" + tag + "]",
            MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                AddSphereLight(w, s, Point3(0, 2, 0), 0.5, RgbColor(1.0, 0.4, 0.1), 10.0, nee); }),
            Vector3(overhead) * Vector3(1.0, 0.4, 0.1));
    }

    // Rectangular (quad) lights: reference = numerical integration over the quad.
    //   (a) facing the ground, (b) tilted toward it, (c) facing away (must give 0), (d) two-sided, facing away.
    {
        struct QuadCase { const char* name; Point3 corner; Vector3 e1, e2; bool twoSided; };
        const QuadCase cases[] = {
            { "quad light overhead 1x1 at height 2, facing down", Point3(-0.5, 2, -0.5), Vector3(1, 0, 0), Vector3(0, 0, 1), false },
            { "quad light tilted toward the receiver",             Point3(1.0, 1.0, -0.7), Vector3(0, 0, 1.4), Vector3(-0.5, 0.7, 0), false },
            { "quad light CLOSE strip 0.6x1.8 at height 0.3",      Point3(0.1, 0.3, -0.9), Vector3(0.6, 0, 0), Vector3(0, 0, 1.8), false },
            { "quad light facing AWAY (one-sided, expect 0)",      Point3(-0.5, 2, -0.5), Vector3(0, 0, 1), Vector3(1, 0, 0), false },
            { "quad light facing away but two-sided",              Point3(-0.5, 2, -0.5), Vector3(0, 0, 1), Vector3(1, 0, 0), true },
        };
        for (const QuadCase& c : cases) {
            const Real expected = QuadGroundRadiance(c.corner, c.e1, c.e2, c.twoSided, 10.0, albedo);
            for (int variant = 0; variant < 2; ++variant) {
                const bool nee = (variant == 0);
                run(std::string(c.name) + (nee ? "  [NEE + MIS]" : "  [BSDF sampling only]"),
                    MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                        AddQuadLight(w, &s, c.corner, c.e1, c.e2, EmissionSpectrum(RgbColor(1.0)), 10.0, nee, c.twoSided); }),
                    Vector3(expected));
            }
        }
    }

    // Several lights of very different brightness at once. Every selection strategy must give the
    // same (sum of the single-light closed forms) value; the standard error shows how efficient it is.
    {
        const Real dA = 6.0;
        const Real sphereFar = albedo * 400.0 * (0.5 * 0.5) / (dA * dA);                         // bright far sphere
        const Real pointNear = albedo / constants::PI * 2.0 * (0.9 / 1.5) / 2.25;               // dim near point (1.2, 0.9, 0)
        const Point3 qc(0.4, 1.2, -0.4); const Vector3 qe1(1, 0, 0), qe2(0, 0, 0.8);
        const Real quadNear = QuadGroundRadiance(qc, qe1, qe2, false, 3.0, albedo);              // dim near quad
        const Real envLevel = albedo * 0.1;                                                      // uniform sky 0.1
        const Real expected = sphereFar + pointNear + quadNear + envLevel;
        auto dimEnv = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(0.1))));
        std::cout << "[Lights]   (multi-light scene: far sphere " << sphereFar << " + point " << pointNear
            << " + quad " << quadNear << " + sky " << envLevel << ")\n";
        for (const bool withSky : { true, false }) {
            for (const char* sampler : { "uniform", "power", "spatial", "all" }) {
                run(std::string(withSky ? "4 lights (with sky)" : "3 lights (no sky)") + ", light sampler = " + sampler,
                    MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                        AddSphereLight(w, s, Point3(0, dA, 0), 0.5, EmissionSpectrum(RgbColor(1.0)), 400.0, true);
                        s.addLight(std::make_shared<PointLight>(Point3(1.2, 0.9, 0), RgbColor(1.0), 2.0));
                        AddQuadLight(w, &s, qc, qe1, qe2, EmissionSpectrum(RgbColor(1.0)), 3.0, true);
                        if (withSky) s.setEnvLight(std::make_shared<EnvLight>(dimEnv, 1.0)); }, sampler),
                    Vector3(withSky ? expected : expected - envLevel));
            }
        }
    }

    // Light source spectra (EmissionSpectrum). Unit-luminance lights: Y of the result equals the
    // D65 point-light value above, only the color differs.
    {
        const Real pointY = albedo / constants::PI * 10.0 / 4.0;
        const EmissionSpectrum tungsten = EmissionSpectrum::IlluminantA();
        auto pointScene = [&](const EmissionSpectrum& e) {
            return MakeGroundScene(albedo, [&](Aggregate&, Scene& s) {
                s.addLight(std::make_shared<PointLight>(Point3(0, 2, 0), e, 10.0)); });
        };

        // (a) the SPD itself: CIE illuminant A has chromaticity xy = (0.4476, 0.4074).
        color::CIE_XYZ sum;
        for (int l = 380; l <= 780; ++l) {
            const Real S = color::evaluateBlackbody(Real(l), 2856.0);
            const color::CIE_XYZ c = color::evaluateCIE_XYZ_1931_MultiLobe(Real(l));
            sum.X += c.X * S; sum.Y += c.Y * S; sum.Z += c.Z * S;
        }
        const Real x = sum.X / (sum.X + sum.Y + sum.Z), y = sum.Y / (sum.X + sum.Y + sum.Z);
        const bool okA = std::abs(x - 0.4476) < 0.005 && std::abs(y - 0.4074) < 0.005;
        std::cout << "[Lights]   blackbody 2856 K chromaticity: xy = (" << x << ", " << y
            << ") vs CIE A (0.4476, 0.4074)   " << (okA ? "OK" : "FAIL") << "\n";
        allPassed = okA && allPassed;

        // (b) point light with the SPD: expected = sensor RGB of the SPD (Y = 1) x D65-light value.
        run("point light, illuminant A (tungsten)", pointScene(tungsten), tungsten.rgb() * pointY);
        run("point light, 1850 K candle blackbody", pointScene(EmissionSpectrum::Blackbody(1850.0)),
            EmissionSpectrum::Blackbody(1850.0).rgb() * pointY);

        // (c) CSV round trip: write a 10 nm SPD table (arbitrary scale), load it, same light.
        const std::string csvPath = (std::filesystem::temp_directory_path() / "rayt_test_spd.csv").string();
        {
            std::ofstream f(csvPath);
            f << "wavelength_nm,power\n";
            for (int l = 300; l <= 900; l += 10) f << l << "," << 7.0 * color::evaluateBlackbody(Real(l), 2856.0) << "\n";
        }
        color::SpectralCurve curve;
        std::string err;
        if (io::loadSpectralCurveCsv(csvPath, curve, &err)) {
            const EmissionSpectrum fromCsv = EmissionSpectrum::FromCurve(curve);
            const Vector3 d = glm::abs(fromCsv.rgb() / tungsten.rgb() - Vector3(1.0));
            const bool okC = std::max({ d.x, d.y, d.z }) < 0.005;
            std::cout << "[Lights]   SPD from CSV (" << curve.size() << " rows) vs generated A: max RGB difference "
                << std::max({ d.x, d.y, d.z }) * 100.0 << " %   " << (okC ? "OK" : "FAIL") << "\n";
            allPassed = okC && allPassed;
            run("point light, SPD loaded from CSV", pointScene(fromCsv), fromCsv.rgb() * pointY);
        }
        else {
            std::cout << "[Lights]   SPD CSV load failed: " << err << "   FAIL\n";
            allPassed = false;
        }
        std::error_code ec;
        std::filesystem::remove(csvPath, ec);

        // (d) the same SPD through a DiffuseLight sphere (NEE + MIS): expected = rgb x overhead value.
        run("sphere light, illuminant A  [NEE + MIS]",
            MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
                AddSphereLight(w, s, Point3(0, 2, 0), 0.5, tungsten, 10.0, true); }),
            tungsten.rgb() * overhead);

        // (e) white balance on the tungsten light: the result must be neutral with the same luminance.
        color::setWhiteBalance(tungsten);
        run("point light, illuminant A with white balance (neutral)", pointScene(tungsten), Vector3(pointY));
        color::resetWhiteBalance();
    }

    // Environment maps: reference = dense grid integration of EnvMap::eval.
    Vector3 sunDir(0.0);
    auto sunEnv = MakeSunSkyEnv(sunDir);
    allPassed = CheckEnvSampling(*sunEnv, "synthetic sun + sky", n / 4) && allPassed;
    run("environment: synthetic sun + sky (importance sampled)",
        MakeGroundScene(albedo, [&](Aggregate&, Scene& s) {
            s.setEnvLight(std::make_shared<EnvLight>(sunEnv, 1.0)); }),
        GroundRadianceFromEnv(*sunEnv, albedo));
    run("environment: synthetic sun + sky, scale 2",
        MakeGroundScene(albedo, [&](Aggregate&, Scene& s) {
            s.setEnvLight(std::make_shared<EnvLight>(sunEnv, 2.0)); }),
        GroundRadianceFromEnv(*sunEnv, albedo) * 2.0);

    try {
        auto hdr = std::make_shared<EnvMap>(io::loadHDR(opt.env));
        allPassed = CheckEnvSampling(*hdr, opt.env, n / 4) && allPassed;
        run("environment: " + opt.env,
            MakeGroundScene(albedo, [&](Aggregate&, Scene& s) {
                s.setEnvLight(std::make_shared<EnvLight>(hdr, 1.0)); }),
            GroundRadianceFromEnv(*hdr, albedo));
    }
    catch (const std::exception& e) {
        std::cout << "[Lights]   (skipped HDR test: " << e.what() << ")\n";
    }

    std::cout << "[Lights] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
