// Multiple-bounce microfacet checks (--scene multibounce): Microfacet/MultiBounce.hpp and its use in RoughConductor (docs/MULTIBOUNCE.md).
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Core/Sampling.hpp"
#include "IO/ComplexIorLoader.hpp"
#include "Materials/RoughConductor.hpp"
#include "Microfacet/GGX.hpp"
#include "Microfacet/MultiBounce.hpp"

#include <atomic>
#include <functional>
#include <random>
#include <sstream>
#include <thread>

namespace rayt::cli {

namespace {

    namespace mb = multibounce;

    bool Report(const std::string& name, bool ok, const std::string& detail) {
        std::cout << "[MultiBounce]   " << name << "\n                 " << detail << (ok ? "   OK" : "   FAIL") << "\n";
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

    Vector3 Direction(Real thetaDeg, Real phiDeg) {
        const Real t = thetaDeg * constants::PI / Real(180), p = phiDeg * constants::PI / Real(180);
        return Vector3(std::sin(t) * std::cos(p), std::sin(t) * std::sin(p), std::cos(t));
    }

    struct Mean {
        double mean = 0, error = 0;
        double robust = 0;   ///< median of the batch means (Heitz's estimator has a heavy tail; see MultiBounce.hpp)
    };

    /// Mean and standard error of f(rng) over n samples, in 32 batches spread over the cores (each batch with its own seed).
    template <class F>
    Mean ParallelMean(F&& f, long long n, std::uint64_t seed) {
        constexpr int kBatches = 32;
        const int threads = std::max(1, std::min<int>(kBatches, int(std::thread::hardware_concurrency())));
        std::vector<double> sums(kBatches), sums2(kBatches);
        std::vector<long long> counts(kBatches);
        std::atomic<int> next{ 0 };
        std::vector<std::thread> pool;
        for (int t = 0; t < threads; ++t) {
            pool.emplace_back([&] {
                for (int b; (b = next++) < kBatches;) {
                    const long long count = n / kBatches + (b < n % kBatches ? 1 : 0);
                    Rng rng(seed * 7919 + std::uint64_t(b) * 104729 + 1);
                    double a = 0, a2 = 0;
                    for (long long i = 0; i < count; ++i) { const double v = f(rng); a += v; a2 += v * v; }
                    sums[b] = a; sums2[b] = a2; counts[b] = count;
                }
            });
        }
        for (auto& th : pool) th.join();
        double s = 0, s2 = 0;
        std::vector<double> means;
        for (int b = 0; b < kBatches; ++b) {
            s += sums[b]; s2 += sums2[b];
            if (counts[b] > 0) means.push_back(sums[b] / double(counts[b]));
        }
        Mean m;
        m.mean = s / double(n);
        m.error = std::sqrt(std::max(0.0, s2 / double(n) - m.mean * m.mean) / double(n));
        std::sort(means.begin(), means.end());
        m.robust = means.empty() ? 0 : (means[(means.size() - 1) / 2] + means[means.size() / 2]) / 2;
        return m;
    }

    /// Integral over the upper hemisphere of g(m) dm by the midpoint rule in (theta, phi), with theta finer near the pole.
    template <class G>
    double HemisphereIntegral(G&& g, int nt = 1600, int np = 360) {
        double sum = 0;
        for (int i = 0; i < nt; ++i) {
            // theta = (pi/2) x^2 crowds the nodes near the pole, where a smooth GGX lobe sits
            const double x = (i + 0.5) / nt, t = 0.5 * constants::PI * x * x, dt = constants::PI * x / nt;
            const double st = std::sin(t), ct = std::cos(t);
            for (int j = 0; j < np; ++j) {
                const double p = 2 * constants::PI * (j + 0.5) / np;
                sum += g(Vector3(st * std::cos(p), st * std::sin(p), ct)) * st * dt * (2 * constants::PI / np);
            }
        }
        return sum;
    }

    Spectrum One(Real) { return Spectrum(1.0); }

} // namespace

int RunMultiBounceTests(const CliOptions&) {
    bool all = true;
    std::cout << "[MultiBounce] light bouncing more than once between the facets of a rough conductor (Heitz et al. 2016, Zhang et al. 2026)\n";

    // ---- a: visible normals from either hemisphere ---------------------------------------------------------------------------------
    {
        struct Case { Real ax, ay, theta, phi; };
        const Case cases[] = { { 0.5, 0.5, 30, 0 }, { 1.0, 1.0, 70, 20 }, { 0.5, 0.5, 110, 0 }, { 1.0, 1.0, 160, 40 }, { 0.3, 0.8, 125, 60 }, { 0.8, 0.3, 80, 130 } };
        bool okArea = true, okMoments = true;
        double worstArea = 0, worstMoment = 0;
        for (const Case& c : cases) {
            const GGXDistribution dist(c.ax, c.ay);
            const Vector3 w = Direction(c.theta, c.phi);
            const double area = mb::projectedArea(dist, w);
            const double numeric = HemisphereIntegral([&](const Vector3& m) { return std::max(0.0, double(glm::dot(w, m))) * dist.D(m); });
            worstArea = std::max(worstArea, std::abs(numeric / area - 1));
            okArea = okArea && std::abs(numeric / area - 1) < 2e-3;
            // moments of the sampled normals against the same integrals of pdf(m) g(m)
            auto moment = [&](auto&& g) {
                return HemisphereIntegral([&](const Vector3& m) { return std::max(0.0, double(glm::dot(w, m))) * dist.D(m) / area * g(m); });
            };
            const std::function<double(const Vector3&)> gs[] = {
                [](const Vector3& m) { return double(m.x); }, [](const Vector3& m) { return double(m.y); },
                [](const Vector3& m) { return double(m.z); }, [](const Vector3& m) { return double(m.z * m.z); },
                [](const Vector3& m) { return double(m.x * m.z); } };
            int gi = 0;
            for (const auto& g : gs) {
                const Mean s = ParallelMean([&](Rng& rng) { const Real u1 = rng(), u2 = rng(); return g(mb::sampleVisibleNormal(c.ax, c.ay, w, u1, u2)); }, 2000000, 11 + gi++);
                const double d = std::abs(s.mean - moment(g));
                worstMoment = std::max(worstMoment, d / std::max(s.error, 1e-9));
                okMoments = okMoments && d < 5 * s.error + 1e-4;
            }
        }
        all = Report("a  visible normals from above and from below: the projected area is the integral of max(0, w.m) D(m) (0.2 %)", okArea,
            "worst relative difference " + Num(worstArea, 2)) && all;
        all = Report("a  sampled normals: five moments against the integrals of max(0, w.m) D(m) / area (6 directions, 5 standard errors)", okMoments,
            "worst difference " + Num(worstMoment, 2) + " standard errors") && all;
    }

    // ---- b: the single bounce of both models is the GGX term with the height-correlated G2 ----------------------------------------
    {
        struct Case { Real alpha, ti, to, phi; };
        const Case cases[] = { { 0.5, 30, 50, 180 }, { 1.0, 0, 40, 0 }, { 1.0, 60, 60, 180 }, { 0.5, 70, 20, 90 } };
        bool okZ = true, okH = true;
        std::string detail;
        int index = 0;
        for (const Case& c : cases) {
            const GGXDistribution dist(c.alpha, c.alpha);
            const Vector3 wi = Direction(c.ti, 0), wo = Direction(c.to, c.phi);
            const Vector3 h = glm::normalize(wi + wo);
            const double analytic = dist.D(h) * dist.G(wi, wo) / (4 * wi.z);
            Rng rng(5);
            const double zhang = mb::zhangRho(dist, c.alpha, c.alpha, wi, wo, One, rng, 1, 1)[0];
            okZ = okZ && std::abs(zhang / analytic - 1) < 1e-9;
            const Mean heitz = ParallelMean([&](Rng& r) { return double(mb::heitzRho(dist, c.alpha, c.alpha, wi, wo, One, r, 1, 1)[0]); }, 4000000, 21 + index);
            okH = okH && std::abs(heitz.mean - analytic) < 5 * heitz.error;
            detail += std::string(index++ ? "; " : "") + Num((heitz.mean / analytic - 1) * 100, 2) + " +- " + Num(heitz.error / analytic * 100, 2) + " %";
        }
        all = Report("b  Zhang, one bounce = D G2 / (4 cos_i) (1e-9)", okZ, "4 configurations") && all;
        all = Report("b  Heitz's walk, first hit only, averages to D G2 / (4 cos_i) (5 standard errors)", okH, detail) && all;
    }

    // ---- c: Zhang's segment term is reciprocal ---------------------------------------------------------------------------------------
    {
        Rng rng(31);
        double worst = 0;
        for (int trial = 0; trial < 20000; ++trial) {
            const int n = 1 + int(rng() * mb::kMaxBounces);
            const int turn = 1 + int(rng() * n);   // e[turn..n] go up, e[0..turn-1] down
            Real lam[mb::kMaxBounces + 2], lamR[mb::kMaxBounces + 2];
            bool up[mb::kMaxBounces + 2], upR[mb::kMaxBounces + 2];
            for (int j = 0; j <= n; ++j) { lam[j] = Real(3) * rng() * rng(); up[j] = j >= turn; }
            for (int j = 0; j <= n; ++j) { lamR[j] = lam[n - j]; upR[j] = !up[n - j]; }
            const double a = mb::zhangSegment(lam, up, n), b = mb::zhangSegment(lamR, upR, n);
            worst = std::max(worst, std::abs(a - b) / std::max(a, 1e-300));
        }
        all = Report("c  Zhang's segment term of a path equals that of the reversed path (20000 random paths, 1e-12)", worst < 1e-12,
            "worst relative difference " + Num(worst, 2)) && all;
        // five-bounce paths of each type, against an independent implementation of the supplement's Eq. 5 (Python)
        const Real lam[6] = { 0.3, 0.5, 0.2, 0.7, 0.4, 0.9 };
        const bool ups[4][6] = { { false, true, true, true, true, true }, { false, false, false, false, false, true },
                                 { false, false, true, true, true, true }, { false, false, false, true, true, true } };
        const double expected[4] = { 0.04951475539710835, 0.015081536820667254, 0.12142997272368115, 0.14748213953877812 };
        double worstValue = 0;
        std::string values;
        for (int i = 0; i < 4; ++i) {
            const double v = mb::zhangSegment(lam, ups[i], 5);
            worstValue = std::max(worstValue, std::abs(v / expected[i] - 1));
            values += (i ? ", " : "") + Num(v, 6);
        }
        all = Report("c  Zhang's segment term of four five-bounce paths (types 2, 1, 3, 3) = an independent implementation (1e-12)", worstValue < 1e-12,
            values) && all;
    }

    // ---- d: white furnace (F = 1): the directional albedo -------------------------------------------------------------------------
    {
        bool okH = true, okZ = true;
        std::string detailH, detailZ;
        int index = 0;
        for (const Real alpha : { Real(0.3), Real(0.6), Real(1.0) }) {
            const GGXDistribution dist(alpha, alpha);
            for (const Real ti : { Real(0), Real(45), Real(75) }) {
                const Vector3 wi = Direction(ti, 0);
                // albedo = integral of rho over wo, with wo drawn by the cosine (pdf cos / pi)
                auto albedo = [&](bool zhang, int maxOrder) {
                    return ParallelMean([&](Rng& r) {
                        const Real u1 = r(), u2 = r();
                        const Real rr = std::sqrt(u1), p = 2 * constants::PI * u2;
                        const Vector3 wo(rr * std::cos(p), rr * std::sin(p), std::sqrt(std::max(Real(0), 1 - u1)));
                        if (!(wo.z > Real(1e-6))) return 0.0;
                        const Spectrum rho = zhang ? mb::zhangRho(dist, alpha, alpha, wi, wo, One, r, 1, maxOrder)
                                                   : mb::heitzRho(dist, alpha, alpha, wi, wo, One, r, 1, maxOrder);
                        return double(rho[0]) * constants::PI / wo.z;
                    }, 1500000, 41 + index++);
                };
                const Mean h = albedo(false, 64), h10 = albedo(false, mb::kMaxBounces), z = albedo(true, mb::kMaxBounces);
                okH = okH && std::abs(h.robust - 1) < 0.01;
                okZ = okZ && std::abs(z.mean - h10.robust) < 0.02;
                detailH += (detailH.empty() ? "" : ", ") + Num(h.robust, 4) + " (" + Num(h.mean, 4) + ")";
                detailZ += (detailZ.empty() ? "" : ", ") + Num(z.mean, 4) + "/" + Num(h10.robust, 4);
            }
        }
        all = Report("d  white furnace, Heitz (all orders): albedo 1 within 0.01 (alpha 0.3, 0.6, 1; theta_i 0, 45, 75; median of 32 batch means, plain mean in brackets)", okH, detailH) && all;
        all = Report("d  white furnace, Zhang / Heitz (both up to 10 bounces) within 0.02", okZ, detailZ) && all;
    }

    // ---- e: Zhang against Heitz, the multiple bounces of the lobe (the paper's Fig. 5: GGX alpha = 1, 10 bounces) ------------------
    {
        bool ok = true;
        std::string detail;
        double worstAbs = 0;
        int index = 0;
        const Real alpha = 1;
        const GGXDistribution dist(alpha, alpha);
        for (const Real ti : { Real(0), Real(60) }) {
            for (const Real to : { Real(-60), Real(-30), Real(0), Real(30), Real(60) }) {
                const Vector3 wi = Direction(ti, 0), wo = Direction(std::abs(to), to < 0 ? 180 : 0);
                const Mean h = ParallelMean([&](Rng& r) { return double(mb::heitzRho(dist, alpha, alpha, wi, wo, One, r, 2, mb::kMaxBounces)[0]); }, 6000000, 61 + index);
                const Mean z = ParallelMean([&](Rng& r) { return double(mb::zhangRho(dist, alpha, alpha, wi, wo, One, r, 2, mb::kMaxBounces)[0]); }, 2000000, 81 + index);
                ++index;
                const double diff = z.mean - h.robust, err = std::hypot(z.error, h.error);
                worstAbs = std::max(worstAbs, std::abs(diff));
                ok = ok && std::abs(diff) < 4 * err + 5e-3;
                detail += std::string(detail.empty() ? "" : "; ") + "(" + Num(ti, 2) + "," + Num(to, 2) + ") " + Num(h.robust, 3) + " vs " + Num(z.mean, 3);
            }
        }
        all = Report("e  multiple bounces of the lobe, rho(theta_i, theta_o) in the plane: Zhang - Heitz (median of batch means) < 0.005 + 4 standard errors (paper: bias ~1e-3)", ok,
            "worst |difference| " + Num(worstAbs, 2) + "\n                 Heitz vs Zhang: " + detail) && all;
    }

    // ---- f: RoughConductor: sample() and eval() agree (the mixture pdf with the stochastic f), for gold ---------------------------
    {
        auto goldTable = std::make_shared<color::ComplexIorTable>();
        std::string err;
        if (!io::loadComplexIor_RefractiveIndexInfoCsv(AU_IOR_CSV_PATH, *goldTable, &err))
            all = Report("f  rough gold (cannot read the gold table: " + err + ")", false, "") && all;
        else {
            bool ok = true;
            std::string detail;
            int index = 0;
            for (const auto model : { mb::Model::Zhang, mb::Model::Heitz }) {
                RoughConductor gold(goldTable, Real(0.9));
                gold.setMultipleScatteringModel(model);
                const SampledWavelengths lambda = SampledWavelengths::sampleVisible(Real(0.37));
                SurfaceInteraction rec;
                rec.n = rec.gn = Vector3(0, 0, 1);
                const Vector3 wo = Direction(40, 0);
                rec.wo = wo;
                // directional albedo by the material's own sampling, and by eval over cosine-distributed directions
                const Mean bySample = ParallelMean([&](Rng& r) {
                    const auto s = gold.sample(rec, wo, Point2(r(), r()), lambda, TransportMode::Radiance);
                    if (!s) return 0.0;
                    return double(s->f[0]) * std::abs(s->wi.z) / s->pdf;
                }, 2000000, 101 + index);
                const Mean byEval = ParallelMean([&](Rng& r) {
                    const Real u1 = r(), u2 = r();
                    const Real rr = std::sqrt(u1), p = 2 * constants::PI * u2;
                    const Vector3 wi(rr * std::cos(p), rr * std::sin(p), std::sqrt(std::max(Real(0), 1 - u1)));
                    if (!(wi.z > Real(1e-6))) return 0.0;
                    return double(gold.eval(rec, wo, wi, lambda, TransportMode::Radiance)[0]) * constants::PI;
                }, 2000000, 121 + index);
                ++index;
                const double diff = bySample.mean - byEval.mean, e = std::hypot(bySample.error, byEval.error);
                ok = ok && std::abs(diff) < 4 * e + 2e-3 * byEval.mean;
                detail += std::string(detail.empty() ? "" : "; ") + (model == mb::Model::Zhang ? "Zhang " : "Heitz ") + Num(bySample.mean, 4) + " vs " + Num(byEval.mean, 4);
            }
            all = Report("f  RoughConductor (gold, roughness 0.9): albedo by sample() = by eval() (4 standard errors + 0.2 %)", ok, detail) && all;
        }
    }

    std::cout << "[MultiBounce] RESULT: " << (all ? "PASS" : "FAIL") << "\n";
    return all ? 0 : 2;
}

} // namespace rayt::cli
