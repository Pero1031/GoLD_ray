// Leaf checks (--scene leaf): PROSPECT (Materials/Prospect.hpp) against an independent implementation, and the leaf BSDF (Materials/Leaf.hpp).
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Core/Sampling.hpp"
#include "Geometry/Frame.hpp"
#include "Materials/Leaf.hpp"
#include "Materials/Prospect.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace rayt::cli {

namespace {

    bool Report(const std::string& name, bool ok, const std::string& detail) {
        std::cout << "[Leaf]   " << name << "\n             " << detail << (ok ? "   OK" : "   FAIL") << "\n";
        return ok;
    }

    std::string Fmt(double v, int prec = 4) {
        std::ostringstream s;
        s << std::setprecision(prec) << v;
        return s.str();
    }

    leaf::Pigments Green() { return leaf::Pigments{}; }
    leaf::Pigments Red() {
        leaf::Pigments p; p.N = Real(1.6); p.chl = 5; p.car = 8; p.ant = 15; p.ewt = Real(0.012); p.lma = Real(0.006); return p;
    }
    leaf::Pigments Brown() {
        leaf::Pigments p; p.N = Real(2.2); p.chl = 2; p.car = 3; p.brown = Real(0.6); p.ewt = Real(0.004); p.lma = Real(0.007); return p;
    }

    /// A hit on a leaf in the plane z = 0, seen from above (front face) or from below (back face).
    SurfaceInteraction Hit(bool fromAbove) {
        SurfaceInteraction rec{};
        rec.p = Point3(0.0);
        rec.gn = rec.n = Vector3(0, 0, fromAbove ? 1 : -1);
        rec.frontFace = fromAbove;
        rec.uv = UV(0.5, 0.5);
        return rec;
    }

    /// The path's wavelengths: in the spectral build the given ones; in the RGB build they are ignored.
    SampledWavelengths At(Real a, Real b, Real c, Real d) {
        SampledWavelengths w;
#if defined(RAYT_SPECTRAL_MODE)
        w.lambda = { a, b, c, d };
#else
        (void)a; (void)b; (void)c; (void)d;
#endif
        return w;
    }

    int Channels() {
#if defined(RAYT_SPECTRAL_MODE)
        return kSpectralSamples;
#else
        return 3;
#endif
    }

    Vector3 Cosine(const Point2& u, bool up) {
        Vector3 d = sampling::CosineSampleHemisphere(u);
        if (!up) d.z = -d.z;
        return d;
    }

    // (a) PROSPECT against tools/lib/prospect.py (numpy, scipy's E1) on the same coefficient table
    bool CheckProspect(const leaf::ProspectTable& table) {
        bool ok = true;
        // E1 against scipy.special.exp1
        const double e1x[] = { 0.01, 0.5, 1.0, 2.5, 30.0 };
        const double e1v[] = { 4.037929576538113, 0.5597735947761608, 0.2193839343955205, 0.024914917870269736, 3.021552010688813e-15 };
        double worst = 0;
        for (int i = 0; i < 5; ++i) worst = std::max(worst, std::abs(double(leaf::expintE1(Real(e1x[i]))) / e1v[i] - 1));
        ok &= Report("E1(x) for x = 0.01 ... 30 against scipy", worst < 1e-12, "worst relative difference " + Fmt(worst, 3));
        // tav against the same formula in numpy
        const double tn[] = { 1.33, 1.47, 1.5115 };
        const double t90[] = { 0.9340691507005972, 0.9126526720476862, 0.9065355562992299 };
        const double t40[] = { 0.978835708878099, 0.9622825488401417, 0.9569218055075076 };
        worst = 0;
        for (int i = 0; i < 3; ++i) {
            worst = std::max(worst, std::abs(double(leaf::tav(Real(90), Real(tn[i]))) - t90[i]));
            worst = std::max(worst, std::abs(double(leaf::tav(Real(40), Real(tn[i]))) - t40[i]));
        }
        ok &= Report("tav(90, n), tav(40, n) for n = 1.33, 1.47, 1.51", worst < 1e-13, "worst difference " + Fmt(worst, 3));

        struct Ref { const char* leaf; int alpha; int nm; double R, T; };
        const Ref refs[] = {
            { "green", 90, 450, 0.09136226275588075, 0.0013262606396651877 }, { "green", 90, 550, 0.1950410817797399, 0.14248664678531245 },
            { "green", 90, 650, 0.09405719651091143, 0.02392047883095247 },  { "green", 90, 700, 0.17169864377910582, 0.1282619979774713 },
            { "green", 40, 450, 0.04125106516246376, 0.001399403660530915 }, { "green", 40, 550, 0.15116726533201982, 0.15025279838113756 },
            { "green", 40, 650, 0.045496092528681586, 0.025202684346952907 }, { "green", 40, 700, 0.12738699492991676, 0.13512363181687964 },
            { "red", 90, 450, 0.0942689338313102, 0.006583399564751313 },     { "red", 90, 550, 0.11544030100173767, 0.0397545183721591 },
            { "red", 90, 650, 0.2734898823061181, 0.21491049595366205 },      { "red", 90, 700, 0.40997318967974583, 0.3525591416173057 },
            { "red", 40, 450, 0.04431803859231852, 0.00694647279284137 },     { "red", 40, 550, 0.0672278904147405, 0.041921315214268715 },
            { "red", 40, 650, 0.23454687924499873, 0.22643030813240528 },     { "red", 40, 700, 0.37840851743315623, 0.37142000278164666 },
            { "brown", 90, 450, 0.15478399282225244, 0.02909123204909484 },   { "brown", 90, 550, 0.32355872163594435, 0.1439575362285291 },
            { "brown", 90, 650, 0.35322412417847415, 0.17753909089044156 },   { "brown", 90, 700, 0.43223876424506147, 0.24263565324559738 },
            { "brown", 40, 450, 0.10817049152385345, 0.0306956079380708 },    { "brown", 40, 550, 0.28668968420702445, 0.15180385779574762 },
            { "brown", 40, 650, 0.31855510264863407, 0.18705568975345702 },   { "brown", 40, 700, 0.4018652337419373, 0.25561593606677535 },
        };
        worst = 0;
        for (const Ref& r : refs) {
            const std::string name = r.leaf;
            const leaf::Pigments p = name == "green" ? Green() : name == "red" ? Red() : Brown();
            const leaf::ReflTrans rt = leaf::prospect(table, p, Real(r.nm), Real(r.alpha));
            worst = std::max({ worst, std::abs(double(rt.R) - r.R), std::abs(double(rt.T) - r.T) });
        }
        ok &= Report("R and T of a green, a red and a brown leaf, 450-700 nm, alpha 90 and 40, against numpy", worst < 1e-10,
            "worst difference " + Fmt(worst, 3));

        // a leaf without any absorber reflects and transmits everything
        leaf::Pigments clear; clear.chl = clear.car = clear.ant = clear.brown = clear.ewt = clear.lma = 0;
        worst = 0;
        for (Real N : { Real(1), Real(1.5), Real(3) })
            for (Real nm : { Real(450), Real(600) }) {
                const leaf::ReflTrans rt = leaf::prospect(table, clear, nm);
                worst = std::max(worst, std::abs(double(rt.R + rt.T) - 1));
            }
        ok &= Report("no absorber: R + T = 1 (N = 1, 1.5, 3)", worst < 1e-12, "worst |R + T - 1| " + Fmt(worst, 3));
        return ok;
    }

    // (b) under light that is uniform over all directions, the diffuse lobes reflect Rd and transmit T, exactly as PROSPECT says
    bool CheckUniformLight(const std::shared_ptr<const leaf::ProspectTable>& table) {
        bool ok = true;
        const SampledWavelengths lam = At(Real(450), Real(550), Real(650), Real(700));
        for (const auto& [name, p] : { std::pair{ "green", Green() }, std::pair{ "red", Red() } }) {
            const Leaf leaf(table, p);
            const SurfaceInteraction rec = Hit(true);
            const Leaf::Optics o = leaf.opticsAt(rec, lam);
            // (1/pi) Integral Integral f |cos_o| |cos_i| dwo dwi: cosine-distributed wo and wi make each sample pi f
            const int n = 2'000'000;
            std::vector<double> sumR(4, 0), sumT(4, 0);
            for (int i = 0; i < n; ++i) {
                const Vector3 wo = Cosine(Point2(sampling::Random(), sampling::Random()), true);
                const Vector3 wr = Cosine(Point2(sampling::Random(), sampling::Random()), true);
                const Vector3 wt = Cosine(Point2(sampling::Random(), sampling::Random()), false);
                const Spectrum fr = leaf.evalDiffuse(rec, wo, wr, lam), ft = leaf.evalDiffuse(rec, wo, wt, lam);
                for (int c = 0; c < Channels(); ++c) { sumR[c] += constants::PI * fr[c]; sumT[c] += constants::PI * ft[c]; }
            }
            double worst = 0;
            std::string detail;
            for (int c = 0; c < Channels(); ++c) {
                const double r = sumR[c] / n, t = sumT[c] / n;
                worst = std::max({ worst, std::abs(r / std::max(double(o.Rd[c]), 1e-9) - 1) * (o.Rd[c] > 0.01 ? 1 : 0),
                    std::abs(t / std::max(double(o.T[c]), 1e-9) - 1) * (o.T[c] > 0.01 ? 1 : 0) });
                detail += (c ? ", " : "") + Fmt(r, 4) + "/" + Fmt(o.Rd[c], 4) + " " + Fmt(t, 4) + "/" + Fmt(o.T[c], 4);
            }
            ok &= Report(std::string(name) + " leaf: reflected / Rd, transmitted / T per channel", worst < 0.004, detail + ";  worst " + Fmt(100 * worst, 2) + " %");
        }
        return ok;
    }

    // (c) reciprocity: f(wo, wi) = f(wi, wo); a transmission seen from the other side of the leaf
    bool CheckReciprocity(const std::shared_ptr<const leaf::ProspectTable>& table) {
        const Leaf leaf(table, Red());
        const SampledWavelengths lam = At(Real(480), Real(560), Real(620), Real(690));
        const SurfaceInteraction above = Hit(true), below = Hit(false);
        double worstR = 0, worstT = 0;
        for (int i = 0; i < 20000; ++i) {
            const Vector3 a = Cosine(Point2(sampling::Random(), sampling::Random()), true);
            const Vector3 b = Cosine(Point2(sampling::Random(), sampling::Random()), true);
            const Vector3 c = Cosine(Point2(sampling::Random(), sampling::Random()), false);
            const Spectrum r1 = leaf.eval(above, a, b, lam), r2 = leaf.eval(above, b, a, lam);
            const Spectrum t1 = leaf.eval(above, a, c, lam), t2 = leaf.eval(below, c, a, lam);
            for (int k = 0; k < Channels(); ++k) {
                worstR = std::max(worstR, std::abs(double(r1[k] - r2[k])) / std::max(double(r1[k]), 1e-6));
                worstT = std::max(worstT, std::abs(double(t1[k] - t2[k])) / std::max(double(t1[k]), 1e-6));
            }
        }
        return Report("reciprocity: reflection f(a, b) = f(b, a); transmission from above = from below (20 000 pairs)",
            worstR < 1e-10 && worstT < 1e-10, "worst relative difference " + Fmt(worstR, 3) + ", " + Fmt(worstT, 3));
    }

    // (d) sample() against eval(): the directional albedo two ways, the pdf integrates to at most 1; (e) energy
    bool CheckSampling(const std::shared_ptr<const leaf::ProspectTable>& table) {
        bool ok = true;
        const SampledWavelengths lam = At(Real(450), Real(550), Real(650), Real(700));
        leaf::Pigments clear; clear.chl = clear.car = clear.ant = clear.brown = clear.ewt = clear.lma = 0;
        for (const auto& [name, p] : { std::pair{ "red", Red() }, std::pair{ "clear (no absorber)", clear } }) {
            for (bool fromAbove : { true, false }) {
                const Leaf leaf(table, p, Real(0.55), Real(0.75));
                const SurfaceInteraction rec = Hit(fromAbove);
                for (Real deg : { Real(0), Real(50), Real(80) }) {
                    const Real th = deg * constants::PI / Real(180);
                    const Vector3 wo(std::sin(th), 0, std::cos(th) * (fromAbove ? 1 : -1));
                    const int n = 400'000;
                    double sS = 0, sS2 = 0, sE = 0, sE2 = 0, sP = 0;
                    for (int i = 0; i < n; ++i) {
                        auto s = leaf.sample(rec, wo, Point2(sampling::Random(), sampling::Random()), lam);
                        double v = 0;
                        if (s) for (int c = 0; c < Channels(); ++c) v += double(s->f[c]) * std::abs(glm::dot(rec.n, s->wi)) / s->pdf / Channels();
                        sS += v; sS2 += v * v;
                        const Vector3 w = sampling::UniformSampleSphere(Point2(sampling::Random(), sampling::Random()));
                        const Spectrum f = leaf.eval(rec, wo, w, lam);
                        double e = 0;
                        for (int c = 0; c < Channels(); ++c) e += double(f[c]) * std::abs(glm::dot(rec.n, w)) * 4 * constants::PI / Channels();
                        sE += e; sE2 += e * e;
                        sP += double(leaf.pdf(rec, wo, w, lam)) * 4 * constants::PI;
                    }
                    const double mS = sS / n, mE = sE / n;
                    const double se = std::sqrt((sS2 / n - mS * mS) / n + (sE2 / n - mE * mE) / n);
                    const double pdfInt = sP / n;
                    const double z = std::abs(mS - mE) / std::max(se, 1e-12);
                    const bool good = z < 4 && pdfInt < 1.01 && pdfInt > 0.95 && mS <= 1.0 + 3 * se;
                    ok &= Report(std::string(name) + (fromAbove ? ", upper face" : ", lower face") + ", theta_o " + Fmt(deg, 3) + ": albedo by sample() = by eval(); pdf integral",
                        good, Fmt(mS, 5) + " vs " + Fmt(mE, 5) + " (" + Fmt(z, 2) + " sigma), pdf integral " + Fmt(pdfInt, 4));
                }
            }
        }
        return ok;
    }

} // namespace

int RunLeafTests(const CliOptions& /*opt*/) {
    std::string err;
    auto table = leaf::ProspectTable::standard(&err);
    if (!table) {
        std::cerr << "[Leaf] " << err << "\n";
        return 1;
    }
    sampling::SetRunSeed(20261124);
    bool ok = true;
    std::cout << "[Leaf] (a) PROSPECT\n";
    ok &= CheckProspect(*table);
    std::cout << "[Leaf] (b) the BSDF under uniform light\n";
    ok &= CheckUniformLight(table);
    std::cout << "[Leaf] (c) reciprocity\n";
    ok &= CheckReciprocity(table);
    std::cout << "[Leaf] (d) sampling, (e) energy\n";
    ok &= CheckSampling(table);
    {
        const Leaf red(table, Red()), green(table, Green());
        const Vector3 a = red.rgbReflectance(), b = red.rgbTransmittance(), c = green.rgbReflectance(), d = green.rgbTransmittance();
        std::cout << "[Leaf] linear sRGB under D65: red leaf R (" << a.x << ", " << a.y << ", " << a.z << ") T (" << b.x << ", " << b.y << ", " << b.z
                  << "); green leaf R (" << c.x << ", " << c.y << ", " << c.z << ") T (" << d.x << ", " << d.y << ", " << d.z << ")\n";
    }
    std::cout << (ok ? "[Leaf] ALL PASSED\n" : "[Leaf] SOME CHECKS FAILED\n");
    return ok ? 0 : 1;
}

} // namespace rayt::cli
