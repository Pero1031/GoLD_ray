// Glint checks (--scene glints): the glinty GGX distribution (Microfacet/GlintNDF.hpp, docs/GLINTS.md).
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Color/SpectralConversion.hpp"
#include "Core/Sampling.hpp"
#include "Geometry/Quad.hpp"
#include "IO/ComplexIorLoader.hpp"
#include "Materials/RoughConductor.hpp"
#include "Microfacet/GGX.hpp"
#include "Microfacet/GlintNDF.hpp"
#include "Render/Camera.hpp"
#include "Render/Film.hpp"
#include "Render/Integrator.hpp"
#include "Scene/Aggregate.hpp"
#include "Scene/Scene.hpp"

#include <random>
#include <sstream>

namespace rayt::cli {

namespace {

    bool Report(const std::string& name, bool ok, const std::string& detail) {
        std::cout << "[Glints]   " << name << "\n             " << detail << (ok ? "   OK" : "   FAIL") << "\n";
        return ok;
    }

    std::string Num(double v, int digits = 3) {
        std::ostringstream s;
        s.precision(digits);
        s << v;
        return s.str();
    }

    struct Rng {
        std::mt19937_64 g;
        explicit Rng(std::uint64_t seed) : g(seed) {}
        Real operator()() { return static_cast<Real>(g() >> 11) * Real(1.0 / 9007199254740992.0); }
    };

    /// The filter of a pixel whose standard deviations on the surface are sigma along a direction at `angle` and sigma / aspect across
    /// it, sized so that the distribution's level is `level`.
    GlintFootprint Footprint(const GlintNDF& ndf, Real level, Real aspect = 1, Real angle = 0) {
        const Real sigma = std::ldexp(Real(1), 0) * std::pow(Real(2), level) / (GlintNDF::kCellsPerSigma * std::sqrt(ndf.params().density));
        const Real a = sigma * sigma, b = a / (aspect * aspect), c = std::cos(angle), s = std::sin(angle);
        return GlintFootprint{ a * c * c + b * s * s, (a - b) * c * s, a * s * s + b * c * c };
    }

    Vector3 Direction(Real thetaDeg, Real phiDeg) {
        const Real t = thetaDeg * constants::PI / Real(180), p = phiDeg * constants::PI / Real(180);
        return Vector3(std::sin(t) * std::cos(p), std::sin(t) * std::sin(p), std::cos(t));
    }

    struct Mean {
        double mean = 0, error = 0;
        long long n = 0;
    };

    /// The average of f over random texture points (and whatever else f draws), in batches until its standard error is below `target`
    /// or `maxN` samples are spent.
    template <class F>
    Mean Average(F&& f, std::uint64_t seed, double target, long long maxN) {
        Rng rng(seed);
        double s = 0, s2 = 0;
        long long n = 0;
        const long long batch = 200000;
        Mean m;
        while (n < maxN) {
            for (long long i = 0; i < batch; ++i) {
                const double v = f(rng);
                s += v;
                s2 += v * v;
            }
            n += batch;
            m.mean = s / double(n);
            m.error = std::sqrt(std::max(s2 / double(n) - m.mean * m.mean, 0.0) / double(n));
            m.n = n;
            if (m.error < target) break;
        }
        return m;
    }

    Point2 RandomUv(Rng& rng) { return Point2(rng() * Real(64), rng() * Real(64)); }

    /// A surface that glows with the size of its pixel: what a camera ray tells the hit (SurfaceInteraction::pixelSpacing).
    class PixelSizeProbe final : public Material {
    public:
        Spectrum eval(const SurfaceInteraction&, const Vector3&, const Vector3&, const SampledWavelengths&, TransportMode) const override { return Spectrum(0.0); }
        std::optional<BSDFSample> sample(const SurfaceInteraction&, const Vector3&, const Point2&, const SampledWavelengths&, TransportMode) const override { return std::nullopt; }
        Real pdf(const SurfaceInteraction&, const Vector3&, const Vector3&, const SampledWavelengths&) const override { return 0; }
        Spectrum emitted(const SurfaceInteraction& rec, const Vector3&, const SampledWavelengths& lambda) const override {
            return RgbToIlluminantSpectrum(Vector3(rec.pixelSpacing), lambda);
        }
    };

}   // namespace

int RunGlintTests(const CliOptions&) {
    bool all = true;
    std::cout << "[Glints] the glinty GGX distribution (Kemppinen et al. 2025): the map to the disc, the facets of a level, the average over the surface\n";

    // (a) The map T: there and back, its Jacobian against finite differences, and the smooth distribution against the engine's GGX.
    {
        GlintParams p;
        p.alphaX = Real(0.3);
        p.alphaY = Real(0.12);
        const GlintNDF ndf(p);
        const GGXDistribution ggx(p.alphaX, p.alphaY);
        Rng rng(1);
        double worstBack = 0, worstJ = 0, worstD = 0;
        for (int i = 0; i < 2000; ++i) {
            const Vector3 w = Direction(rng() * Real(88), rng() * Real(360));
            const Point2 u = ndf.toDisc(w);
            worstBack = std::max(worstBack, double(glm::length(ndf.fromDisc(u) - w)));
            const Real h = Real(1e-5);
            const Vector3 dx = (ndf.fromDisc(Point2(u.x + h, u.y)) - ndf.fromDisc(Point2(u.x - h, u.y))) / (2 * h);
            const Vector3 dy = (ndf.fromDisc(Point2(u.x, u.y + h)) - ndf.fromDisc(Point2(u.x, u.y - h))) / (2 * h);
            const double solidAnglePerArea = glm::length(glm::cross(dx, dy));
            if (math::sqr(u.x - Real(0.5)) + math::sqr(u.y - Real(0.5)) < Real(0.2499))
                worstJ = std::max(worstJ, std::fabs(double(ndf.jacobian(w)) * solidAnglePerArea - 1));
            worstD = std::max(worstD, std::fabs(double(ndf.smoothD(w) / ggx.D(w)) - 1));
        }
        all = Report("a  the map to the disc: back to the same normal (1e-9), its Jacobian against finite differences (1e-4), D against the engine's GGX (1e-9)",
            worstBack < 1e-9 && worstJ < 1e-4 && worstD < 1e-9,
            "back " + Num(worstBack) + ", Jacobian off by " + Num(worstJ) + ", D off by " + Num(worstD)) && all;
    }

    // (b) Where every facet of a level is summed (level 0: one angular cell; level 1: 2 x 2), the distribution is the plain sum of its
    //     facets: against the exhaustive sum, which visits six instead of four standard deviations of the filter.
    {
        GlintParams p;
        p.microRoughness = Real(0.08);
        const GlintNDF ndf(p);
        double worst = 0;
        std::string detail;
        for (const Real level : { Real(-3), Real(0), Real(1) }) {
            const GlintFootprint fp = Footprint(ndf, level, level == Real(1) ? Real(5) : Real(1), Real(0.6));
            Rng rng(2);
            double diff = 0, total = 0;
            for (int i = 0; i < 200000; ++i) {
                const Vector3 w = Direction(rng() * Real(40), rng() * Real(360));
                const Point2 uv = RandomUv(rng);
                const double fast = ndf.ratio(w, uv, fp), full = ndf.ratioExhaustive(w, uv, fp, std::max(0, int(level)));
                diff += std::fabs(fast - full);
                total += full;
            }
            worst = std::max(worst, diff / total);
            detail += std::string(detail.empty() ? "" : ", ") + "level " + Num(level) + ": " + Num(diff / total);
        }
        all = Report("b  levels 0 and 1 against the plain sum of their facets (the filter cut at 4 against 6 standard deviations: 0.2 %)",
            worst < 2e-3, "sum of differences / sum: " + detail) && all;
    }

    // (c) The average over the surface is the smooth distribution: D* / D averaged over random texture points, for fixed half vectors
    //     and footprints at several levels (isotropic and stretched). Expected off by the filter's cut (-0.013 %) and, far from the
    //     normal, by the disc's rim.
    {
        struct Case { Real alphaX, alphaY, micro, level, aspect, theta, phi; };
        const Case cases[] = {
            { Real(0.3), Real(0.3), Real(0.08), Real(-2.0), 1, 0, 0 },
            { Real(0.3), Real(0.3), Real(0.08), Real(0.5), 1, 10, 40 },
            { Real(0.3), Real(0.3), Real(0.03), Real(2.3), 1, 15, 200 },
            { Real(0.3), Real(0.3), Real(0.01), Real(5.7), 1, 20, 100 },
            { Real(0.3), Real(0.3), Real(0.01), Real(9.2), 1, 5, 300 },
            { Real(0.3), Real(0.3), Real(0.02), Real(3.4), 8, 12, 20 },
            { Real(0.4), Real(0.1), Real(0.02), Real(4.5), 3, 8, 75 },
            { Real(0.1), Real(0.1), Real(0.005), Real(6.5), 1, 4, 10 },
        };
        bool ok = true;
        std::string detail;
        double worst = 0;
        int index = 0;
        for (const Case& c : cases) {
            GlintParams p;
            p.alphaX = c.alphaX;
            p.alphaY = c.alphaY;
            p.microRoughness = c.micro;
            const GlintNDF ndf(p);
            const GlintFootprint fp = Footprint(ndf, c.level, c.aspect, Real(0.4));
            const Vector3 w = Direction(c.theta, c.phi);
            const Mean m = Average([&](Rng& rng) { return double(ndf.ratio(w, RandomUv(rng), fp)); }, 100 + index, 0.003, 40000000);
            const double off = m.mean - 1;
            ok = ok && std::fabs(off) < 4 * m.error + 0.003 && m.error < 0.01;
            worst = std::max(worst, std::fabs(off));
            detail += std::string(index ? "; " : "") + "level " + Num(c.level) + ": " + Num(off * 100, 2) + " +- " + Num(m.error * 100, 2) + " %";
            ++index;
        }
        all = Report("c  the average over the surface is the smooth GGX (8 cases: levels -2 to 9.2, stretched footprints, anisotropic roughness; 4 standard errors + 0.3 %)",
            ok, detail) && all;
    }

    // (d) Beyond the block summed one by one, the facets are replaced by their expectation: the average of the fast evaluation against
    //     that of the plain sum at level 3 (64 angular cells), and how much less the fast one varies from pixel to pixel.
    {
        GlintParams p;
        p.microRoughness = Real(0.04);
        const GlintNDF ndf(p);
        const GlintFootprint fp = Footprint(ndf, 3);
        const Vector3 w = Direction(12, 130);
        const Mean fast = Average([&](Rng& rng) { return double(ndf.ratio(w, RandomUv(rng), fp)); }, 7, 0.003, 4000000);
        const Mean full = Average([&](Rng& rng) { return double(ndf.ratioExhaustive(w, RandomUv(rng), fp, 3)); }, 8, 0.003, 2000000);
        const double sdFast = fast.error * std::sqrt(double(fast.n)), sdFull = full.error * std::sqrt(double(full.n));
        all = Report("d  level 3: the average with the far facets replaced by their expectation against the plain sum (4 standard errors + 0.3 %)",
            std::fabs(fast.mean - full.mean) < 4 * std::hypot(fast.error, full.error) + 0.003,
            "fast " + Num(fast.mean, 5) + " +- " + Num(fast.error, 2) + ", plain sum " + Num(full.mean, 5) + " +- " + Num(full.error, 2)
                + "; standard deviation over pixels " + Num(sdFast) + " against " + Num(sdFull)) && all;
    }

    // (e) Energy: D* / D averaged over half vectors drawn from D cos (the normals of a hemisphere, stretched) and over the surface: 1
    //     for a distribution that keeps the smooth one's projected area. The model is not exact here (the paper says so): the lobes
    //     of facets near the disc's rim reach outside it, and the expectation is taken over the square that holds the disc.
    {
        std::string detail;
        double worst = 0;
        int index = 0;
        for (const Real alpha : { Real(0.1), Real(0.3), Real(0.6) })
            for (const Real micro : { Real(0.01), Real(0.05) })
                for (const Real level : { Real(1.5), Real(6.5) }) {
                    GlintParams p;
                    p.alphaX = p.alphaY = alpha;
                    p.microRoughness = micro;
                    const GlintNDF ndf(p);
                    const GlintFootprint fp = Footprint(ndf, level);
                    const Mean m = Average([&](Rng& rng) {
                        const Vector3 s = sampling::CosineSampleHemisphere(Point2(rng(), rng()));
                        const Vector3 w = glm::normalize(Vector3(alpha * s.x, alpha * s.y, std::max(s.z, Real(1e-6))));
                        return double(ndf.ratio(w, RandomUv(rng), fp));
                    }, 300 + index, 0.003, 20000000);
                    worst = std::max(worst, std::fabs(m.mean - 1) - 4 * m.error);
                    detail += std::string(index ? "; " : "") + "a " + Num(alpha) + " micro " + Num(micro) + " level " + Num(level) + ": " + Num((m.mean - 1) * 100, 2)
                        + " +- " + Num(m.error * 100, 2) + " %";
                    ++index;
                }
        all = Report("e  energy against the smooth distribution, averaged over the surface (12 cases; within 2 %)", worst < 0.02, detail) && all;
    }

    // (f) A pixel that holds very many facets sees the smooth distribution: the spread of D* / D over pixels at levels 12 and 16, and
    //     exactly 1 beyond the last level.
    {
        GlintParams p;
        const GlintNDF ndf(p);
        std::string detail;
        double worst = 0;
        for (const Real level : { Real(12), Real(16) }) {
            const GlintFootprint fp = Footprint(ndf, level);
            Rng rng(11);
            double s2 = 0;
            const int n = 20000;
            for (int i = 0; i < n; ++i) {
                const Vector3 w = Direction(rng() * Real(30), rng() * Real(360));
                s2 += math::sqr(double(ndf.ratio(w, RandomUv(rng), fp)) - 1);
            }
            const double rms = std::sqrt(s2 / n);
            worst = std::max(worst, rms);
            detail += "level " + Num(level) + ": rms " + Num(rms * 100, 2) + " %; ";
        }
        const Real beyond = ndf.ratio(Direction(10, 0), Point2(3, 4), Footprint(ndf, Real(GlintNDF::kMaxLevel) + Real(0.5)));
        all = Report("f  many facets per pixel: D* / D stays near 1 (rms 3 %), and is 1 beyond the last level", worst < 0.03 && beyond == Real(1),
            detail + "beyond the last level " + Num(beyond)) && all;
    }

    // (g) Information: the paper's softened roulette (epsilon 0.1, for animation) between two levels changes the average a little.
    {
        GlintParams p;
        p.microRoughness = Real(0.03);
        p.rouletteSoftness = Real(0.1);
        const GlintNDF soft(p);
        const GlintFootprint fp = Footprint(soft, Real(2.5));
        const Vector3 w = Direction(10, 50);
        const Mean m = Average([&](Rng& rng) { return double(soft.ratio(w, RandomUv(rng), fp)); }, 21, 0.003, 20000000);
        std::cout << "[Glints]   g  (information) softened roulette, level 2.5: the average is off by " << Num((m.mean - 1) * 100, 2) << " +- " << Num(m.error * 100, 2) << " %\n";
    }

    // (h) The material: a rough gold with glints on a plane tilted against all three coordinate planes (three arrangements of
    //     facets, mixed by a roulette), seen through a pixel of a given size: its BRDF averaged over the plane against the smooth
    //     gold's. And where the pixel's size is unknown, the smooth BRDF itself.
    {
        auto goldTable = std::make_shared<color::ComplexIorTable>();
        std::string err;
        if (!io::loadComplexIor_RefractiveIndexInfoCsv("assets/Metal_data/Johnson_Au.csv", *goldTable, &err))
            all = Report("h  the material (cannot read the gold table: " + err + ")", false, "") && all;
        else {
            RoughConductor smooth(goldTable, Real(0.3)), glinty(goldTable, Real(0.3));
            smooth.setMultipleScattering(false);
            glinty.setMultipleScattering(false);
            GlintParams gp;
            gp.density = Real(1e4);
            gp.microRoughness = Real(0.03);
            glinty.setGlints(gp, Real(0.5));
            const SampledWavelengths lambda = SampledWavelengths::sampleVisible(Real(0.37));
            const Vector3 n = glm::normalize(Vector3(1, 2, 3));
            const Vector3 t1 = glm::normalize(glm::cross(n, Vector3(0, 0, 1))), t2 = glm::cross(n, t1);
            // view and light 25 and 40 deg from the normal, not in one plane with it
            const Vector3 wo = glm::normalize(n * std::cos(Real(0.44)) + t1 * std::sin(Real(0.44)));
            const Vector3 wi = glm::normalize(n * std::cos(Real(0.70)) + (t1 * Real(-0.8) + t2 * Real(0.6)) * std::sin(Real(0.70)));
            SurfaceInteraction rec;
            rec.n = rec.gn = n;
            rec.wo = wo;
            const Real reference = smooth.eval(rec, wo, wi, lambda, TransportMode::Radiance)[0];
            bool ok = reference > Real(0);
            std::string detail;
            int index = 0;
            for (const Real spacing : { Real(0.004), Real(0.04), Real(0.6) }) {
                const Mean m = Average([&](Rng& rng) {
                    SurfaceInteraction r = rec;
                    r.p = Point3(5, -3, 2) + t1 * (rng() * Real(40)) + t2 * (rng() * Real(40));
                    r.pixelSpacing = spacing;
                    return double(glinty.eval(r, wo, wi, lambda, TransportMode::Radiance)[0] / reference);
                }, 500 + index, 0.003, 40000000);
                ok = ok && std::fabs(m.mean - 1) < 4 * m.error + 0.003 && m.error < 0.01;
                detail += "pixel " + Num(spacing) + ": " + Num((m.mean - 1) * 100, 2) + " +- " + Num(m.error * 100, 2) + " %; ";
                ++index;
            }
            rec.p = Point3(1, 2, 3);
            const Real unknown = glinty.eval(rec, wo, wi, lambda, TransportMode::Radiance)[0];
            all = Report("h  the material on a tilted plane: its BRDF averaged over the plane against the smooth gold's (4 standard errors + 0.3 %); the smooth BRDF where the pixel's size is unknown",
                ok && unknown == reference, detail + "unknown pixel: " + (unknown == reference ? "the same" : "different")) && all;
        }
    }

    // (j) Sampling the facets: the half vectors that sample() draws against the density pdf() claims, at one place on the surface.
    //     Six moments over the disc (1, x, y, x^2, y^2, xy) of the drawn points against the same moments of the density integrated
    //     with uniform points of the unit square; the first is the share of draws that give a sample.
    {
        GlintParams p;
        p.microRoughness = Real(0.02);
        const GlintNDF ndf(p);
        bool ok = true;
        std::string detail;
        int index = 0;
        struct Case { Real level, aspect; };
        for (const Case& c : { Case{ Real(-1.5), 1 }, Case{ Real(0.6), 1 }, Case{ Real(2.4), 4 }, Case{ Real(4.3), 1 } }) {
            const GlintFootprint fp = Footprint(ndf, c.level, c.aspect, Real(0.5));
            const Point2 uv(Real(3.21), Real(7.77));
            const Real weight = index == 2 ? Real(0.6) : Real(1);
            const std::uint32_t stream = static_cast<std::uint32_t>(index);
            auto moments = [](const Point2& x, double m[6]) {
                m[0] = 1; m[1] = x.x; m[2] = x.y; m[3] = x.x * x.x; m[4] = x.y * x.y; m[5] = x.x * x.y;
            };
            double s[6] = {}, s2[6] = {}, t[6] = {}, t2[6] = {};
            const int n = 3000000;
            Rng rng(700 + index), rng2(800 + index);
            const std::function<Real()> draw = [&rng] { return rng(); };
            for (int i = 0; i < n; ++i) {
                Vector3 wh;
                if (ndf.sample(uv, fp, weight, stream, draw, wh)) {
                    double m[6];
                    moments(ndf.toDisc(wh), m);
                    for (int k = 0; k < 6; ++k) { s[k] += m[k]; s2[k] += m[k] * m[k]; }
                }
                const Point2 x(rng2(), rng2());
                if (math::sqr(x.x - Real(0.5)) + math::sqr(x.y - Real(0.5)) < Real(0.25)) {
                    const Vector3 w = ndf.fromDisc(x);
                    const double density = double(ndf.pdf(w, uv, fp, weight, stream) / ndf.jacobian(w));   // per area of the disc
                    double m[6];
                    moments(x, m);
                    for (int k = 0; k < 6; ++k) { const double v = density * m[k]; t[k] += v; t2[k] += v * v; }
                }
            }
            double worst = 0;
            for (int k = 0; k < 6; ++k) {
                const double a = s[k] / n, b = t[k] / n;
                const double ea = std::sqrt(std::max(s2[k] / n - a * a, 0.0) / n), eb = std::sqrt(std::max(t2[k] / n - b * b, 0.0) / n);
                const double excess = std::fabs(a - b) - 4 * std::hypot(ea, eb) - 0.002;
                worst = std::max(worst, excess);
                ok = ok && excess < 0;
            }
            detail += std::string(index ? "; " : "") + "level " + Num(c.level) + ": " + Num(s[0] / n * 100, 3) + " % of the draws give a sample, the density integrates to "
                + Num(t[0] / n * 100, 3) + " %";
            ++index;
        }
        all = Report("j  sampling the facets: six moments of the drawn half vectors against the density's (4 levels; 4 standard errors + 0.002)", ok, detail) && all;
    }

    // (k) The material's sample() with and without the facets' sampling: two estimators of the same reflectance under uniform light
    //     (the BRDF x cosine / density, averaged), at one place on the tilted plane; and how much less the first one varies.
    {
        auto goldTable = std::make_shared<color::ComplexIorTable>();
        std::string err;
        if (!io::loadComplexIor_RefractiveIndexInfoCsv("assets/Metal_data/Johnson_Au.csv", *goldTable, &err))
            all = Report("k  the material's sampling (cannot read the gold table: " + err + ")", false, "") && all;
        else {
            const SampledWavelengths lambda = SampledWavelengths::sampleVisible(Real(0.37));
            const Vector3 n = glm::normalize(Vector3(1, 2, 3));
            const Vector3 t1 = glm::normalize(glm::cross(n, Vector3(0, 0, 1)));
            const Vector3 wo = glm::normalize(n * std::cos(Real(0.44)) + t1 * std::sin(Real(0.44)));
            bool ok = true;
            std::string detail;
            int index = 0;
            for (const Real spacing : { Real(0.002), Real(0.004), Real(0.012) }) {
                RoughConductor facets(goldTable, Real(0.3)), smoothOnly(goldTable, Real(0.3));
                GlintParams gp;
                gp.density = Real(1e4);
                gp.microRoughness = Real(0.01);
                facets.setGlints(gp, Real(0.5));
                smoothOnly.setGlints(gp, Real(0.5));
                smoothOnly.setGlintSampling(false);
                SurfaceInteraction rec;
                rec.n = rec.gn = n;
                rec.wo = wo;
                rec.p = Point3(Real(1.3), Real(-0.4), Real(2.2));
                rec.pixelSpacing = spacing;
                auto estimate = [&](const RoughConductor& m, std::uint64_t seed) {
                    return Average([&](Rng& rng) {
                        const auto bs = m.sample(rec, wo, Point2(rng(), rng()), lambda, TransportMode::Radiance);
                        if (!bs) return 0.0;
                        return double(bs->f[0] * std::max(glm::dot(n, bs->wi), Real(0)) / bs->pdf);
                    }, seed, 0.004, 4000000);
                };
                const Mean a = estimate(facets, 900 + index), b = estimate(smoothOnly, 950 + index);
                ok = ok && std::fabs(a.mean - b.mean) < 4 * std::hypot(a.error, b.error) + 0.003 * b.mean;
                const double sdA = a.error * std::sqrt(double(a.n)), sdB = b.error * std::sqrt(double(b.n));
                detail += std::string(index ? "; " : "") + "pixel " + Num(spacing) + ": " + Num(a.mean, 4) + " +- " + Num(a.error, 2) + " against " + Num(b.mean, 4) + " +- " + Num(b.error, 2)
                    + ", standard deviation per sample " + Num(sdA) + " against " + Num(sdB);
                ++index;
            }
            all = Report("k  the material's reflectance under uniform light, sampled with and without the facets' sampling (3 pixel sizes; 4 standard errors + 0.3 %)", ok, detail) && all;
        }
    }

    // (i) The pixel's size that a camera ray hands to the surface: a plane 5 m in front of a pinhole camera (36 deg high, 9 pixels)
    //     sees, at the picture's centre, 5 x 2 tan(18 deg) / 9; in an orthographic view 2 m high, 2 / 9 at any distance.
    {
        auto probe = std::make_shared<PixelSizeProbe>();
        auto world = std::make_shared<Aggregate>();
        world->add(std::make_shared<Quad>(Point3(-50, -50, -5), Vector3(100, 0, 0), Vector3(0, 100, 0), probe));
        auto scene = std::make_shared<Scene>(world);
        scene->finalize("uniform");
        const int n = 9;
        auto centre = [&](const std::shared_ptr<Camera>& camera, bool on) {
            Film film(n, n);
            PathIntegrator integrator(camera, 2, 64);
            integrator.setPixelFootprints(on);
            integrator.render(*scene, film);
            const Vector3 v = film.getData()[static_cast<size_t>(n / 2) * n + n / 2];
            return double(0.2126 * v.x + 0.7152 * v.y + 0.0722 * v.z);
        };
        auto pinhole = std::make_shared<Camera>(Point3(0.0), Point3(0, 0, -1), Vector3(0, 1, 0), 36.0, 1.0, 0.0, 5.0);
        auto ortho = std::make_shared<Camera>(Point3(0.0), Point3(0, 0, -1), Vector3(0, 1, 0), 36.0, 1.0, 0.0, 5.0);
        ortho->setOrthographic(2);
        const double expectP = 5 * 2 * std::tan(18 * 3.14159265358979323846 / 180) / n, expectO = 2.0 / n;
        const double gotP = centre(pinhole, true), gotO = centre(ortho, true), off = centre(pinhole, false);
        all = Report("i  the pixel's size at the hit: pinhole and orthographic cameras against closed forms (1 %), 0 when not asked for",
            std::fabs(gotP / expectP - 1) < 0.01 && std::fabs(gotO / expectO - 1) < 0.01 && off == 0,
            "pinhole " + Num(gotP, 5) + " (expected " + Num(expectP, 5) + "), orthographic " + Num(gotO, 5) + " (expected " + Num(expectO, 5) + "), off " + Num(off)) && all;
    }

    std::cout << "[Glints] RESULT: " << (all ? "PASS" : "FAIL") << "\n";
    return all ? 0 : 2;
}

} // namespace rayt::cli
