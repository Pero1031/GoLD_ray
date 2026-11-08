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

// -----------------------------------------------------------------------------
// Anisotropy checks (--scene aniso): the tangent frame from the geometry, the anisotropic GGX, the rotation of the roughness axes.
// -----------------------------------------------------------------------------
int RunAnisotropyTests(const CliOptions& opt) {
    bool allPassed = true;
    const int n = std::max(opt.spp, 64) * 1500;
    std::cout << "[Aniso] " << n << " samples per estimate\n";
    const double pi = constants::PI;

    // ---- 1. The tangent that the shapes provide ----------------------------------------------------------------------------
    {
        SurfaceInteraction rec;
        bool ok = true;
        Sphere sphere(Point3(0, 0, 0), Real(1), std::make_shared<Lambertian>(RgbColor(Real(0.5))));
        Ray r(Point3(0.6, 0.3, 3.0), glm::normalize(Vector3(0.0, 0.0, -1.0)));
        ok = ok && sphere.hit(r, rec) && std::abs(glm::length(rec.dpdu) - 1.0) < 1e-9 && std::abs(glm::dot(rec.dpdu, rec.gn)) < 1e-9 && std::abs(rec.dpdu.y) < 1e-9;
        std::cout << "[Aniso]   sphere: dpdu is a unit vector along the circle of latitude (perpendicular to the normal and to y)" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;

        Quad quad(Point3(0, 0, 0), Vector3(2, 0, 0), Vector3(0, 0, -2), std::make_shared<Lambertian>(RgbColor(Real(0.5))));
        SurfaceInteraction qrec;
        const bool qok = quad.hit(Ray(Point3(1, 1, -1), Vector3(0, -1, 0)), qrec) && glm::length(qrec.dpdu - Vector3(1, 0, 0)) < 1e-9;
        std::cout << "[Aniso]   quad: dpdu is the first edge" << (qok ? "   OK" : "   FAIL") << "\n";
        allPassed = qok && allPassed;

        auto mesh = std::make_shared<TriangleMesh>();
        mesh->p = { Point3(0, 0, 0), Point3(1, 0, 0), Point3(0, 1, 0) };
        mesh->indices = { 0, 1, 2 };
        auto mat = std::make_shared<Lambertian>(RgbColor(Real(0.5)));
        Triangle plain(mesh, 0, mat);
        SurfaceInteraction t0;
        const bool noUv = plain.hit(Ray(Point3(0.2, 0.2, 1), Vector3(0, 0, -1)), t0) && glm::length(t0.dpdu - Vector3(1, 0, 0)) < 1e-9;
        mesh->uv = { Point2(0, 0), Point2(0, 1), Point2(1, 0) };   // u grows along +y
        Triangle withUv(mesh, 0, mat);
        SurfaceInteraction t1;
        const bool uv = withUv.hit(Ray(Point3(0.2, 0.2, 1), Vector3(0, 0, -1)), t1) && glm::length(t1.dpdu - Vector3(0, 1, 0)) < 1e-9;
        std::cout << "[Aniso]   triangle: dpdu along the first edge without texture coordinates" << (noUv ? "   OK" : "   FAIL")
            << ", along increasing u with them" << (uv ? "   OK" : "   FAIL") << "\n";
        allPassed = noUv && uv && allPassed;
    }

    // ---- 2. The visible-normal density of the anisotropic GGX integrates to 1 over the normals that face wo ------------------
    {
        bool ok = true;
        for (double ax : { 0.15, 0.6 }) {
            const GGXDistribution g(Real(ax), Real(0.3 * ax));
            for (double az : { 0.0, 60.0, 90.0 }) {
                const double th = 55.0 * pi / 180.0, ph = az * pi / 180.0;
                const Vector3 wo(std::sin(th) * std::cos(ph), std::sin(th) * std::sin(ph), std::cos(th));
                double sum = 0;
                const int nt = 900, np = 1800;
                for (int i = 0; i < nt; ++i) {
                    const double t = (i + 0.5) / nt * pi / 2;
                    for (int j = 0; j < np; ++j) {
                        const double p = (j + 0.5) / np * 2 * pi;
                        const Vector3 wh(std::sin(t) * std::cos(p), std::sin(t) * std::sin(p), std::cos(t));
                        if (glm::dot(wo, wh) > Real(0)) sum += g.pdf(wo, wh) * std::sin(t);   // only normals facing wo are visible (and ever occur as a half vector)
                    }
                }
                sum *= (pi / 2 / nt) * (2 * pi / np);
                const bool o = std::abs(sum - 1.0) < 0.01;
                ok = ok && o;
                std::cout << "[Aniso]   GGX (alpha " << ax << " x " << 0.3 * ax << "), wo azimuth " << az << " deg: integral of the visible normal density " << sum << (o ? "   OK" : "   FAIL") << "\n";
            }
        }
        allPassed = ok && allPassed;
    }

    // ---- 3. Rough anisotropic gold: sample() against eval(), and the rotation of the axes -------------------------------------
    auto goldTable = std::make_shared<color::ComplexIorTable>();
    {
        std::string err;
        if (!io::loadComplexIor_RefractiveIndexInfoCsv("assets/Metal_data/Johnson_Au.csv", *goldTable, &err)) {
            std::cout << "[Aniso] cannot read the gold table: " << err << "\n";
            return 2;
        }
    }
    SurfaceInteraction rec;
    rec.p = Point3(0.0); rec.n = Vector3(0, 0, 1); rec.gn = Vector3(0, 0, 1); rec.dpdu = Vector3(1, 0, 0);
#if defined(RAYT_SPECTRAL_MODE)
    const SampledWavelengths lambda = SampledWavelengths::sampleVisible(Real(0.37));
#else
    const SampledWavelengths lambda;
#endif
    auto outgoing = [&](double polarDeg, double azimuthDeg) {
        const double t = polarDeg * pi / 180.0, p = azimuthDeg * pi / 180.0;
        return Vector3(std::sin(t) * std::cos(p), std::sin(t) * std::sin(p), std::cos(t));
    };
    auto albedoBySample = [&](const Material& m, const Vector3& wo) {
        Estimate e;
        for (int i = 0; i < n; ++i) {
            const auto bs = m.sample(rec, wo, Point2(sampling::Random(), sampling::Random()), lambda, TransportMode::Radiance);
            e.add(bs ? average(bs->f * (std::max(bs->wi.z, Real(0)) / bs->pdf)) : 0.0);
        }
        e.finish();
        return e;
    };
    auto albedoByEval = [&](const Material& m, const Vector3& wo) {
        Estimate e;
        for (int i = 0; i < n; ++i) {
            const Real r = std::sqrt(sampling::Random()), phi = Real(2) * constants::PI * sampling::Random();
            const Vector3 wi(r * std::cos(phi), r * std::sin(phi), std::sqrt(std::max(Real(0), Real(1) - r * r)));
            e.add(average(m.eval(rec, wo, wi, lambda, TransportMode::Radiance)) * constants::PI);
        }
        e.finish();
        return e;
    };
    auto gold = [&](double rotationDeg) {
        auto m = std::make_shared<RoughConductor>(goldTable, Real(0.35), Real(-0.9));   // smoother along the tangent: brushed along it
        m->setTangentRotation(Real(rotationDeg));
        return m;
    };
    {
        const auto m = gold(0);
        for (double az : { 0.0, 45.0, 90.0 }) {
            const Vector3 wo = outgoing(60, az);
            const Estimate a = albedoBySample(*m, wo), b = albedoByEval(*m, wo);
            const double se = std::sqrt(a.se * a.se + b.se * b.se);
            const bool ok = std::abs(a.mean - b.mean) <= 4 * se + 0.01 * a.mean;
            std::cout << "[Aniso]   brushed gold (roughness 0.35, anisotropy -0.9), wo 60 deg polar, " << az << " deg azimuth: albedo by sample() " << a.mean
                << "   by eval() " << b.mean << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }
        const double a0 = albedoBySample(*m, outgoing(60, 0)).mean, a90 = albedoBySample(*m, outgoing(60, 90)).mean;
        const bool differs = std::abs(a0 - a90) > 0.01 * a0;
        std::cout << "[Aniso]   the albedo depends on the azimuth (" << a0 << " along the tangent, " << a90 << " across): " << (differs ? "OK" : "FAIL (no anisotropy?)") << "\n";
        allPassed = differs && allPassed;
    }
    {
        // Rotating the axes by rho about the normal = looking from azimuth - rho at the unrotated surface.
        const struct { double az, rot; } cases[] = { { 90, 90 }, { 30, 60 }, { 150, 45 } };
        for (const auto& c : cases) {
            const Estimate a = albedoBySample(*gold(c.rot), outgoing(60, c.az));
            const Estimate b = albedoBySample(*gold(0), outgoing(60, c.az - c.rot));
            const double se = std::sqrt(a.se * a.se + b.se * b.se);
            const bool ok = std::abs(a.mean - b.mean) <= 4 * se + 0.003 * a.mean;
            std::cout << "[Aniso]   azimuth " << c.az << " deg, axes rotated by " << c.rot << " deg: " << a.mean << "   azimuth " << c.az - c.rot << " deg, not rotated: " << b.mean << (ok ? "   OK" : "   FAIL") << "\n";
            allPassed = ok && allPassed;
        }
    }
    {
        // A lacquer over brushed gold: the base keeps the tangent of the surface.
        LayeredCoat::Params p; p.ior = Real(1.5); p.roughness = Real(0.15);
        auto base = std::make_shared<RoughConductor>(goldTable, Real(0.35), Real(-0.9), Real(1.5));
        base->setTangentRotation(Real(30));
        const LayeredCoat m(p, base);
        const Vector3 wo = outgoing(55, 20);
        const Estimate a = albedoBySample(m, wo), b = albedoByEval(m, wo);
        const double se = std::sqrt(a.se * a.se + b.se * b.se);
        const bool ok = std::abs(a.mean - b.mean) <= 4 * se + 0.01 * a.mean;
        std::cout << "[Aniso]   lacquer over brushed gold: albedo by sample() " << a.mean << "   by eval() " << b.mean << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    std::cout << "[Aniso] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
