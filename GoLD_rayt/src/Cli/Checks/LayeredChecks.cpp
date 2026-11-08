// Checks of the layered coat material (--scene layered).
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/EngineHeaders.hpp"

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

// Cosine-weighted average of the dielectric Fresnel reflectance seen from outside: 2 Integral F(mu) mu dmu.
double averageFresnel(double eta) {
    const int n = 4000;
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        const double mu = (i + 0.5) / n;
        sum += fresnel::fresnelDielectric(Real(mu), Real(1), Real(eta)) * mu;
    }
    return 2.0 * sum / n;
}

} // namespace

int RunLayeredTests(const CliOptions& opt) {
    bool allPassed = true;
    const int n = std::max(opt.spp, 64) * 1500;   // samples per estimate (default --spp 256: about 380k)
    std::cout << "[Layered] " << n << " samples per estimate\n";

    SurfaceInteraction rec;
    rec.p = Point3(0.0); rec.n = Vector3(0, 0, 1); rec.gn = Vector3(0, 0, 1);
#if defined(RAYT_SPECTRAL_MODE)
    const SampledWavelengths lambda = SampledWavelengths::sampleVisible(Real(0.37));
#else
    const SampledWavelengths lambda;
#endif
    auto outgoing = [](double deg) { const double t = deg * constants::PI / 180.0; return Vector3(std::sin(t), 0, std::cos(t)); };

    // Directional albedo by sampling the walk (sample()): the mean of f cos / pdf.
    auto albedoBySample = [&](const LayeredCoat& m, const Vector3& wo) {
        Estimate e;
        for (int i = 0; i < n; ++i) {
            const auto bs = m.sample(rec, wo, Point2(sampling::Random(), sampling::Random()), lambda, TransportMode::Radiance);
            e.add(bs ? average(bs->f * (std::max(bs->wi.z, Real(0)) / bs->pdf)) : 0.0);
        }
        e.finish();
        return e;
    };
    // The same albedo by evaluating f(wo, wi) at cosine-distributed directions: the mean of f pi.
    auto albedoByEval = [&](const LayeredCoat& m, const Vector3& wo) {
        Estimate e;
        for (int i = 0; i < n; ++i) {
            const Real r = std::sqrt(sampling::Random()), phi = Real(2) * constants::PI * sampling::Random();
            const Vector3 wi(r * std::cos(phi), r * std::sin(phi), std::sqrt(std::max(Real(0), Real(1) - r * r)));
            e.add(average(m.eval(rec, wo, wi, lambda, TransportMode::Radiance)) * constants::PI);
        }
        e.finish();
        return e;
    };

    // ---- 1. A coat that is nearly smooth over a white-ish diffuse base: the closed form of the doubled interface ----------
    {
        const double eta = 1.5, R = 0.8;
        const double fdrE = averageFresnel(eta), fdrI = 1.0 - (1.0 - fdrE) / (eta * eta);
        LayeredCoat::Params p; p.ior = Real(eta); p.roughness = Real(0.04);
        const LayeredCoat m(p, RgbColor(Real(R)));
        for (double deg : { 0.0, 40.0, 70.0 }) {
            const double mu = std::cos(deg * constants::PI / 180.0);
            const double F = fresnel::fresnelDielectric(Real(mu), Real(1), Real(eta));
            const double expected = F + (1 - F) * R * (1 - fdrI) / (1 - R * fdrI);
            const Estimate e = albedoBySample(m, outgoing(deg));
            const bool ok = std::abs(e.mean - expected) <= 4 * e.se + 0.02 * expected;
            std::cout << "[Layered]   smooth-ish coat (roughness 0.04) over a diffuse base R = 0.8, " << deg << " deg: albedo expected " << expected
                << "   got " << e.mean << " +- " << e.se << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }
    }

    // ---- 2. eval() and sample() are two independent estimators of the same albedo -----------------------------------------
    auto goldTable = std::make_shared<color::ComplexIorTable>();
    {
        std::string err;
        if (!io::loadComplexIor_RefractiveIndexInfoCsv("assets/Metal_data/Johnson_Au.csv", *goldTable, &err)) {
            std::cout << "[Layered] cannot read the gold table: " << err << "\n";
            return 2;
        }
    }
    struct Case { std::string name; std::shared_ptr<LayeredCoat> m; };
    std::vector<Case> cases;
    auto add = [&](const std::string& name, double coatRough, const Vector3& tint, bool gold, double baseRough) {
        LayeredCoat::Params p; p.ior = Real(1.5); p.roughness = Real(coatRough); p.tint = tint;
        if (gold) {
            auto base = std::make_shared<RoughConductor>(goldTable, Real(baseRough), Real(0), Real(1.5));
            cases.push_back({ name, std::make_shared<LayeredCoat>(p, base) });
        }
        else cases.push_back({ name, std::make_shared<LayeredCoat>(p, RgbColor(Real(0.7), Real(0.1), Real(0.05))) });
    };
    add("coat 0.1 over diffuse", 0.1, Vector3(1.0), false, 0);
    add("coat 0.4 over diffuse", 0.4, Vector3(1.0), false, 0);
    add("coat 0.1 over diffuse, tinted", 0.1, Vector3(0.9, 0.6, 0.3), false, 0);
    add("coat 0.15 over gold 0.3", 0.15, Vector3(1.0), true, 0.3);
    add("coat 0.4 over gold 0.5, tinted", 0.4, Vector3(0.95, 0.8, 0.5), true, 0.5);
    for (const Case& c : cases) {
        for (double deg : { 20.0, 65.0 }) {
            const Estimate a = albedoBySample(*c.m, outgoing(deg)), b = albedoByEval(*c.m, outgoing(deg));
            const double se = std::sqrt(a.se * a.se + b.se * b.se);
            const bool ok = std::abs(a.mean - b.mean) <= 4 * se + 0.01 * a.mean;
            std::cout << "[Layered]   " << c.name << ", " << deg << " deg: albedo by sample() " << a.mean << " +- " << a.se
                << "   by eval() " << b.mean << " +- " << b.se << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }
    }

    // ---- 3. Energy: a lossless coat over a white diffuse base never returns more than it receives -------------------------
    {
        LayeredCoat::Params p; p.ior = Real(1.5); p.roughness = Real(0.3);
        const LayeredCoat m(p, RgbColor(Real(1)));
        bool ok = true;
        double worst = 0;
        for (double deg : { 0.0, 30.0, 60.0, 80.0 }) {
            const Estimate e = albedoBySample(m, outgoing(deg));
            worst = std::max(worst, e.mean);
            ok = ok && e.mean <= 1.0 + 4 * e.se + 0.005;
        }
        std::cout << "[Layered]   white diffuse base under a rough coat (0.3): largest albedo " << worst << (ok ? "   OK (<= 1)" : "   FAIL (> 1)") << "\n";
        allPassed = ok && allPassed;
    }

    // ---- 4. Reciprocity of the mean of eval() ------------------------------------------------------------------------------
    {
        LayeredCoat::Params p; p.ior = Real(1.5); p.roughness = Real(0.2);
        const LayeredCoat m(p, RgbColor(Real(0.6)));
        const Vector3 a = glm::normalize(Vector3(0.5, 0.2, 0.8)), b = glm::normalize(Vector3(-0.3, 0.4, 0.7));
        Estimate ab, ba;
        for (int i = 0; i < n; ++i) {
            ab.add(average(m.eval(rec, a, b, lambda, TransportMode::Radiance)));
            ba.add(average(m.eval(rec, b, a, lambda, TransportMode::Radiance)));
        }
        ab.finish(); ba.finish();
        const double se = std::sqrt(ab.se * ab.se + ba.se * ba.se);
        const bool ok = std::abs(ab.mean - ba.mean) <= 4 * se + 0.01 * ab.mean;
        std::cout << "[Layered]   reciprocity: f(a, b) = " << ab.mean << " +- " << ab.se << "   f(b, a) = " << ba.mean << " +- " << ba.se << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    std::cout << "[Layered] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
