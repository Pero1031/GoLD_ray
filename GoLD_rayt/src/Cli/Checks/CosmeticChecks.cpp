// Foundation layer checks (--scene cosmetic): the two particles of Materials/CosmeticLayer.hpp (SGGX platelets, Henyey-Greenstein
// diffusers) against numerical integrals, the random walk through the slab against energy conservation, a closed form for single
// scattering and reciprocity, and a table of what a layer reflects and lets through.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/EngineHeaders.hpp"
#include "Materials/CosmeticLayer.hpp"

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

// Integral of f over the unit sphere (midpoint rule in theta and phi).
template <class F>
Real IntegrateSphere(F f, int nTheta = 720, int nPhi = 1440) {
    Real sum = 0;
    const Real dT = constants::PI / nTheta, dP = constants::TWO_PI / nPhi;
    for (int i = 0; i < nTheta; ++i) {
        const Real th = (i + Real(0.5)) * dT, st = std::sin(th), ct = std::cos(th);
        Real row = 0;
        for (int j = 0; j < nPhi; ++j) {
            const Real ph = (j + Real(0.5)) * dP;
            row += f(Vector3(st * std::cos(ph), st * std::sin(ph), ct));
        }
        sum += row * st;
    }
    return sum * dT * dP;
}

Vector3 Dir(Real polarDeg, Real azimuthDeg, bool down) {
    const Real t = polarDeg * constants::PI / 180, p = azimuthDeg * constants::PI / 180;
    return Vector3(std::sin(t) * std::cos(p), std::sin(t) * std::sin(p), down ? -std::cos(t) : std::cos(t));
}

std::string Fixed(Real v, int digits = 4) {
    std::ostringstream s;
    s << std::fixed << std::setprecision(digits) << v;
    return s.str();
}

// What a layer does with light arriving from the air: reflected, passed straight through, passed after scattering (RGB each).
struct Budget { Vector3 reflected{ 0 }, scattered{ 0 }; Real straight = 0, lost = 0, bounces = 0; };

Budget Measure(const CosmeticLayer& layer, const Vector3& w, int n) {
    Budget b;
    for (int i = 0; i < n; ++i) {
        const CosmeticLayer::Walk k = layer.walk(w);
        b.bounces += k.bounces;
        if (k.lost) { b.lost += 1; continue; }
        if (k.up) b.reflected += k.weight;
        else if (k.bounces == 0) b.straight += 1;
        else b.scattered += k.weight;
    }
    const Real inv = Real(1) / n;
    b.reflected *= inv; b.scattered *= inv; b.straight *= inv; b.lost *= inv; b.bounces *= inv;
    return b;
}

// The four foundations the paper fitted (its table 2): two dewy, two matte.
struct Product { const char* name; CosmeticLayer::Params p; };

std::vector<Product> PaperProducts() {
    auto make = [](Real alpha, Real tilt, Vector3 lp, Real cd, Vector3 ld, Real g1, Real g2, Real wg) {
        CosmeticLayer::Params p;
        p.plateletRoughness = alpha; p.plateletTiltDeg = tilt; p.plateletAlbedo = lp; p.diffuserFraction = cd;
        p.diffuserAlbedo = ld; p.g1 = g1; p.g2 = g2; p.lobeWeight = wg;
        return p;
    };
    return {
        { "dewy 1", make(0.13, 1.47, Vector3(0.99, 0.85, 0.75), 0.90, Vector3(0.99, 0.98, 0.95), 0.55, 0.09, 1.0) },
        { "dewy 2", make(0.002, 5.5, Vector3(0.8), 0.99, Vector3(0.97, 0.93, 0.87), 0.77, -0.23, 0.79) },
        { "matte 1", make(0.15, 2.01, Vector3(0.91, 0.89, 0.85), 0.89, Vector3(0.95, 0.94, 0.91), 0.37, -0.25, 0.62) },
        { "matte 2", make(0.49, 52.2, Vector3(0.93, 0.83, 0.71), 0.29, Vector3(0.9, 0.9, 0.84), 0.24, -0.22, 0.51) },
    };
}

} // namespace

int RunCosmeticTests(const CliOptions& opt) {
    bool allPassed = true;
    auto report = [&](const std::string& name, bool ok, const std::string& detail = "") {
        std::cout << "[Cosmetic]   " << name << (detail.empty() ? "" : "  " + detail) << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };
    const int scale = std::max(1, opt.spp / 64);   // --spp scales the samples (default 256: 4x)
    std::cout << "[Cosmetic] a layer of foundation: spherical diffusers and platelets in a thin slab\n";

    // (1) The platelets. The normals' distribution must give the projected area, the phase function must integrate to 1, and the
    //     sampler must draw from it.
    std::cout << "[Cosmetic] 1. platelets (SGGX)\n";
    for (const auto& [alpha, tilt] : { std::pair<Real, Real>{ 0.15, 0.0 }, { 0.5, 0.0 }, { 0.3, 35.0 } }) {
        CosmeticLayer::Params p;
        p.plateletRoughness = alpha; p.plateletTiltDeg = tilt;
        const CosmeticLayer layer(p);
        Real worstArea = 0, worstNorm = 0, worstSample = 0;
        for (const Vector3& w : { Dir(10, 20, true), Dir(50, 200, true), Dir(80, 90, true), Dir(40, 300, false) }) {
            const Vector3 wi = -w;
            const Real area = IntegrateSphere([&](const Vector3& n) { return layer.plateletD(n) * std::max(Real(0), glm::dot(n, wi)); });
            worstArea = std::max(worstArea, std::abs(area / layer.plateletArea(w) - 1));
            const Real norm = IntegrateSphere([&](const Vector3& wo) { return layer.plateletPhase(w, wo); });
            worstNorm = std::max(worstNorm, std::abs(norm - 1));
            // the mean of a test function under the sampler against its integral with the phase function
            auto test = [](const Vector3& wo) { return (Real(1) + wo.z) * (Real(1) + wo.z) + Real(0.5) * wo.x + Real(0.3) * wo.y * wo.y; };
            const Real expected = IntegrateSphere([&](const Vector3& wo) { return layer.plateletPhase(w, wo) * test(wo); });
            Real mean = 0;
            const int n = 100000 * scale;
            for (int i = 0; i < n; ++i) mean += test(layer.samplePlatelet(w, sampling::Random(), sampling::Random()));
            worstSample = std::max(worstSample, std::abs(mean / n - expected) / std::max(Real(0.05), std::abs(expected)));
        }
        const std::string label = "roughness " + Fixed(alpha, 2) + ", tilt " + Fixed(tilt, 0) + " deg: ";
        report(label + "the normals' distribution gives the projected area", worstArea < 0.01, "largest difference " + Fixed(100 * worstArea, 2) + " %");
        report(label + "the phase function integrates to 1", worstNorm < 0.01, "largest difference " + Fixed(100 * worstNorm, 2) + " %");
        report(label + "the sampler draws from the phase function", worstSample < 0.015, "largest difference " + Fixed(100 * worstSample, 2) + " %");
    }

    // (2) The diffusers: two Henyey-Greenstein lobes.
    std::cout << "[Cosmetic] 2. diffusers (two Henyey-Greenstein lobes)\n";
    {
        CosmeticLayer::Params p;
        p.g1 = 0.7; p.g2 = -0.3; p.lobeWeight = 0.6;
        const CosmeticLayer layer(p);
        const Vector3 w = Dir(35, 50, true);
        const Real norm = IntegrateSphere([&](const Vector3& wo) { return layer.diffuserPhase(w, wo); });
        Real meanCos = 0;
        const int n = 400000 * scale;
        for (int i = 0; i < n; ++i) meanCos += glm::dot(w, layer.sampleDiffuser(w, sampling::Random(), sampling::Random(), sampling::Random()));
        const Real expected = p.lobeWeight * p.g1 + (1 - p.lobeWeight) * p.g2;
        report("the phase function integrates to 1", std::abs(norm - 1) < 0.005, "integral " + Fixed(norm, 5));
        report("the mean cosine of the samples is w g1 + (1 - w) g2", std::abs(meanCos / n - expected) < 0.005, "sampled " + Fixed(meanCos / n) + ", expected " + Fixed(expected));
    }

    // (3) Energy: particles that absorb nothing send every path out again, upward or downward.
    std::cout << "[Cosmetic] 3. energy (albedo 1)\n";
    {
        Real worst = 0, worstLost = 0;
        for (Real t : { 0.25, 1.0, 4.0 })
            for (Real cd : { 1.0, 0.6, 0.0 }) {
                CosmeticLayer::Params p;
                p.thickness = t; p.diffuserFraction = cd; p.diffuserAlbedo = Vector3(1.0); p.plateletAlbedo = Vector3(1.0); p.plateletTiltDeg = 20;
                const Budget b = Measure(CosmeticLayer(p), Dir(40, 0, true), 50000 * scale);
                worst = std::max(worst, std::abs(b.reflected.y + b.scattered.y + b.straight + b.lost - 1));
                worstLost = std::max(worstLost, b.lost);
            }
        report("reflected + transmitted = 1 for thickness 0.25, 1, 4 and three mixtures", worst < 1e-9 && worstLost < 1e-4,
            "largest difference " + Fixed(worst, 12) + ", paths dropped at the bounce limit " + Fixed(100 * worstLost, 4) + " %");
        CosmeticLayer::Params p;
        p.thickness = 0.001; p.diffuserAlbedo = Vector3(1.0);
        const Budget thin = Measure(CosmeticLayer(p), Dir(0, 0, true), 100000);
        report("a layer of thickness 0.001 lets the light straight through", thin.straight > 0.998, "straight " + Fixed(thin.straight, 5));
    }

    // (4) Single scattering in a slab of isotropic diffusers against the closed form
    //       R1 = (albedo / 2) Integral_0^1 mu / (mu0 + mu) (1 - exp(-t (1 / mu0 + 1 / mu))) dmu,   straight = exp(-t / mu0).
    std::cout << "[Cosmetic] 4. single scattering, isotropic diffusers (closed form)\n";
    {
        CosmeticLayer::Params p;
        p.thickness = 1.0; p.diffuserFraction = 1.0; p.g1 = 0; p.g2 = 0; p.diffuserAlbedo = Vector3(0.8);
        const CosmeticLayer layer(p);
        const Real mu0 = std::cos(30 * constants::PI / 180);
        Real r1 = 0;
        const int steps = 200000;
        for (int i = 0; i < steps; ++i) {
            const Real mu = (i + Real(0.5)) / steps;
            r1 += mu / (mu0 + mu) * (1 - std::exp(-p.thickness * (1 / mu0 + 1 / mu)));
        }
        r1 *= Real(0.8) / 2 / steps;
        Real once = 0, straight = 0;
        const int n = 1000000 * scale;
        for (int i = 0; i < n; ++i) {
            const CosmeticLayer::Walk k = layer.walk(Dir(30, 0, true));
            if (k.lost) continue;
            if (k.up && k.bounces == 1) once += k.weight.y;
            if (!k.up && k.bounces == 0) straight += 1;
        }
        once /= n; straight /= n;
        report("light reflected after one collision", std::abs(once / r1 - 1) < 0.01, "walk " + Fixed(once, 5) + ", closed form " + Fixed(r1, 5));
        report("light that passes without a collision", std::abs(straight / std::exp(-p.thickness / mu0) - 1) < 0.01,
            "walk " + Fixed(straight, 5) + ", exp(-t / mu0) " + Fixed(std::exp(-p.thickness / mu0), 5));
    }

    // (5) Reciprocity, and the two estimators against each other.
    std::cout << "[Cosmetic] 5. reciprocity; the connected estimator against the plain walk\n";
    {
        CosmeticLayer::Params p;
        p.thickness = 1.0; p.diffuserFraction = 0.6; p.plateletRoughness = 0.3; p.plateletTiltDeg = 20;
        p.diffuserAlbedo = Vector3(0.95, 0.85, 0.7); p.plateletAlbedo = Vector3(0.7, 0.8, 0.95); p.g1 = 0.5; p.g2 = -0.3; p.lobeWeight = 0.7;
        const CosmeticLayer layer(p);
        Real worst = 0;
        const int n = 300000 * scale;
        for (const auto& [a, b] : { std::pair<Vector3, Vector3>{ Dir(20, 0, true), Dir(50, 140, false) }, { Dir(60, 90, true), Dir(10, 270, false) }, { Dir(35, 200, true), Dir(35, 20, false) } }) {
            Vector3 fwd(0), back(0);
            for (int i = 0; i < n; ++i) { fwd += layer.reflectionSample(a, b); back += layer.reflectionSample(-b, -a); }
            fwd /= Real(n) * b.z;          // f(a -> b)
            back /= Real(n) * (-a.z);      // f(b -> a)
            for (int c = 0; c < 3; ++c) worst = std::max(worst, std::abs(fwd[c] / back[c] - 1));
        }
        report("f(a -> b) = f(b -> a) for three pairs of directions (tilted platelets, coloured particles)", worst < 0.03, "largest difference " + Fixed(100 * worst, 2) + " %");

        // the reflectance: the plain walk, and the connected estimator integrated over the upper half sphere
        const Vector3 w = Dir(40, 30, true);
        const Budget plain = Measure(layer, w, n);
        Vector3 connected(0);
        for (int i = 0; i < n; ++i) {
            const Vector3 wo = sampling::UniformSampleSphere(sampling::Random2D());
            connected += layer.reflectionSample(w, Vector3(wo.x, wo.y, std::abs(wo.z))) * constants::TWO_PI;
        }
        connected /= Real(n);
        Real diff = 0;
        for (int c = 0; c < 3; ++c) diff = std::max(diff, std::abs(connected[c] / plain.reflected[c] - 1));
        report("both estimators give the same reflectance", diff < 0.03,
            "walk (" + Fixed(plain.reflected.x) + ", " + Fixed(plain.reflected.y) + ", " + Fixed(plain.reflected.z) + "), connected (" + Fixed(connected.x) + ", " + Fixed(connected.y) + ", " + Fixed(connected.z) + ")");
    }

    // (6) Information: the four foundations of the paper (its table 2). What a layer of a given thickness reflects, lets straight
    //     through and lets through after scattering (green channel; light arriving 30 degrees from the normal), and how sharp the
    //     highlight of a thick layer is (the paper measured thick layers on black: thickness 16).
    std::cout << "[Cosmetic] 6. the four foundations fitted in the paper (information)\n";
    std::cout << "[Cosmetic]   product   thickness   reflected   straight   scattered through   absorbed   collisions\n";
    for (const Product& prod : PaperProducts()) {
        for (Real t : { 0.25, 0.5, 1.0, 2.0, 4.0, 16.0 }) {
            CosmeticLayer::Params p = prod.p;
            p.thickness = t;
            const Budget b = Measure(CosmeticLayer(p), Dir(30, 90, true), (t > 4 ? 20000 : 100000) * scale);
            std::cout << "[Cosmetic]   " << std::left << std::setw(8) << prod.name << std::right << std::setw(9) << Fixed(t, 2) << std::setw(13) << Fixed(b.reflected.y, 3)
                << std::setw(11) << Fixed(b.straight, 3) << std::setw(16) << Fixed(b.scattered.y, 3) << std::setw(14) << Fixed(1 - b.reflected.y - b.straight - b.scattered.y - b.lost, 3)
                << std::setw(11) << Fixed(b.bounces, 1) << (b.lost > 0.001 ? "   (dropped " + Fixed(100 * b.lost, 1) + " %)" : "") << "\n";
        }
    }
    // The highlight of a thick layer, as the paper's instrument saw it: light arrives 40 degrees from the normal in the plane the
    // platelets are tilted in; the reflected light is collected in cones of 5 degrees half angle around directions in that plane
    // (+: the mirror side, -: back towards the light). Light per steradian, per unit of arriving light, green channel.
    std::cout << "[Cosmetic]   thick layer (16) lit at 40 degrees: reflected light per steradian in the plane of incidence\n";
    std::cout << "[Cosmetic]   product    back -40   normal 0   mirror +40   brightest cone (angle)   brightest / normal\n";
    for (const Product& prod : PaperProducts()) {
        CosmeticLayer::Params p = prod.p;
        p.thickness = 16;
        const CosmeticLayer layer(p);
        const Vector3 w = Dir(40, 90, true);   // travelling towards +y and down
        const int nBins = 17;                  // -80 ... +80 degrees
        const Real cosHalf = std::cos(5 * constants::PI / 180), coneSolidAngle = constants::TWO_PI * (1 - cosHalf);
        std::vector<Real> bins(nBins, 0);
        const int n = 100000 * scale;
        for (int i = 0; i < n; ++i) {
            const CosmeticLayer::Walk k = layer.walk(w);
            if (k.lost || !k.up) continue;
            for (int b = 0; b < nBins; ++b) {
                const Real th = (-80 + 10 * b) * constants::PI / 180;
                if (glm::dot(k.direction, Vector3(0, std::sin(th), std::cos(th))) >= cosHalf) { bins[b] += k.weight.y; break; }
            }
        }
        for (Real& v : bins) v /= Real(n) * coneSolidAngle;
        int best = 0;
        for (int b = 1; b < nBins; ++b) if (bins[b] > bins[best]) best = b;
        std::cout << "[Cosmetic]   " << std::left << std::setw(8) << prod.name << std::right << std::setw(10) << Fixed(bins[4], 4) << std::setw(11) << Fixed(bins[8], 4)
            << std::setw(13) << Fixed(bins[12], 4) << std::setw(14) << Fixed(bins[best], 4) << " (" << std::showpos << (-80 + 10 * best) << std::noshowpos << ")"
            << std::setw(17) << Fixed(bins[best] / std::max(Real(1e-12), bins[8]), 2) << "\n";
    }

    // The same thick layers as the paper plots them (its figures 4 and 6): the reflectance relative to a white diffuser,
    // pi f(theta_i, theta_o), for light arriving 20 ... 70 degrees from the normal and seen -40 ... +70 degrees in the plane of
    // incidence (cones of 5 degrees half angle). One line per product and angle of incidence: r,g,b for every direction.
    // tools/figures can turn these lines into the coloured patches of the paper (1 = as bright as the white reference).
    // The paper does not say which way the platelets tilt; the two products with a visible tilt are printed both ways.
    std::vector<Product> mapped = PaperProducts();
    for (int i : { 1, 3 }) {
        Product reversed = PaperProducts()[static_cast<size_t>(i)];
        reversed.p.plateletTiltDeg = -reversed.p.plateletTiltDeg;
        static const char* names[] = { "", "dewy 2 (tilt towards the light)", "", "matte 2 (tilt towards the light)" };
        reversed.name = names[i];
        mapped.push_back(reversed);
    }
    std::cout << "[Cosmetic]   reflectance relative to white, thick layer (16): theta_o = -40, -30, ..., +70 (r,g,b each)\n";
    for (const Product& prod : mapped) {
        CosmeticLayer::Params p = prod.p;
        p.thickness = 16;
        const CosmeticLayer layer(p);
        const Real cosHalf = std::cos(5 * constants::PI / 180), coneSolidAngle = constants::TWO_PI * (1 - cosHalf);
        for (int thetaI = 20; thetaI <= 70; thetaI += 10) {
            const int nBins = 12;
            std::vector<Vector3> bins(nBins, Vector3(0));
            const int n = 150000;
            const Vector3 w = Dir(Real(thetaI), 90, true);
            for (int i = 0; i < n; ++i) {
                const CosmeticLayer::Walk k = layer.walk(w);
                if (k.lost || !k.up) continue;
                for (int b = 0; b < nBins; ++b) {
                    const Real th = (-40 + 10 * b) * constants::PI / 180;
                    if (glm::dot(k.direction, Vector3(0, std::sin(th), std::cos(th))) >= cosHalf) { bins[b] += k.weight; break; }
                }
            }
            std::cout << "[Cosmetic]   map " << prod.name << " | " << thetaI << " |";
            for (int b = 0; b < nBins; ++b) {
                const Real cosO = std::cos((-40 + 10 * b) * constants::PI / 180);
                const Vector3 v = bins[b] * (constants::PI / (Real(n) * coneSolidAngle * cosO));
                std::cout << " " << Fixed(v.x, 3) << "," << Fixed(v.y, 3) << "," << Fixed(v.z, 3);
            }
            std::cout << "\n";
        }
    }

    // (7) The layer on skin (CoatedDiffuse::setCosmetic).
    std::cout << "[Cosmetic] 7. the layer on skin\n";
#if defined(RAYT_SPECTRAL_MODE)
    const SampledWavelengths lambda = SampledWavelengths::sampleVisible(Real(0.37));
#else
    const SampledWavelengths lambda;
#endif
    {   // (a) the tables of a layer that absorbs nothing: what is reflected and what crosses add up to 1 for every direction
        CosmeticLayer::Params p;
        p.thickness = 1.0; p.diffuserFraction = 0.7; p.diffuserAlbedo = Vector3(1.0); p.plateletAlbedo = Vector3(1.0);
        const CosmeticCoat coat(p);
        Real worst = 0;
        for (Real mu : { 0.15, 0.5, 0.95 }) worst = std::max(worst, std::abs(coat.reflected(mu).y + coat.through(mu).y - 1));
        report("tables of a layer without absorption: reflected + through = 1", worst < 1e-9, "largest difference " + Fixed(worst, 12)
            + "; diffuse light from below sent back " + Fixed(coat.backDiffuse().y, 3) + ", let through " + Fixed(coat.throughDiffuse().y, 3));
    }
    auto makeSkin = [&](bool sss) {
        auto skin = std::make_shared<CoatedDiffuse>(RgbColor(0.80, 0.52, 0.42), Real(0.5), Real(0.4), Real(0));
        if (sss) skin->setSubsurface(std::make_shared<SubsurfaceProfile>(Vector3(3.67, 1.37, 0.68) * Real(0.001), Real(1.4)));
        return skin;
    };
    auto flatHit = [](const Material* m, const Vector3& wo) {
        SurfaceInteraction rec;
        rec.p = Point3(0.0); rec.n = Vector3(0, 0, 1); rec.gn = Vector3(0, 0, 1); rec.wo = wo; rec.matPtr = m;
        return rec;
    };
    {   // (b) a layer of almost no thickness changes nothing
        CosmeticLayer::Params p;
        p.thickness = Real(1e-6);
        auto bare = makeSkin(false), thin = makeSkin(false);
        thin->setCosmetic(std::make_shared<CosmeticCoat>(p, 2000));
        Real worst = 0;
        for (const auto& [o, i] : { std::pair<Vector3, Vector3>{ Dir(30, 0, false), Dir(30, 180, false) }, { Dir(10, 50, false), Dir(60, 200, false) }, { Dir(70, 90, false), Dir(20, 10, false) } }) {
            const SurfaceInteraction rec = flatHit(bare.get(), o);
            const Real a = average(bare->eval(rec, o, i, lambda, TransportMode::Radiance)), b = average(thin->eval(rec, o, i, lambda, TransportMode::Radiance));
            worst = std::max(worst, std::abs(b / a - 1));
            worst = std::max(worst, std::abs(thin->pdf(rec, o, i, lambda) / bare->pdf(rec, o, i, lambda) - 1));
        }
        report("skin under a layer of thickness 1e-6: eval() and pdf() as without a layer", worst < 1e-3, "largest difference " + Fixed(100 * worst, 4) + " %");
    }
    for (int variant = 0; variant < 4; ++variant) {   // (c) sample() and eval() describe the same reflection: dewy / matte, without / with subsurface scattering
        const bool sss = (variant & 1) != 0;
        CosmeticLayer::Params p;
        CosmeticLayer::preset((variant & 2) ? "matte2" : "dewy1", p);
        p.thickness = 1.0;
        auto skin = makeSkin(sss);
        skin->setCosmetic(std::make_shared<CosmeticCoat>(p));
        const Vector3 wo = Dir(40, 0, false);
        const SurfaceInteraction rec = flatHit(skin.get(), wo);
        Estimate bySample, byEval;
        const int n = 60000 * scale;
        for (int i = 0; i < n; ++i) {
            const auto bs = skin->sample(rec, wo, Point2(sampling::Random(), sampling::Random()), lambda, TransportMode::Radiance);
            bySample.add((bs && !bs->isSubsurface()) ? double(average(bs->f)) * std::abs(glm::dot(rec.n, bs->wi)) / bs->pdf : 0.0);   // the local lobes only
            const Real r = std::sqrt(sampling::Random()), phi = Real(2) * constants::PI * sampling::Random();
            const Vector3 wi(r * std::cos(phi), r * std::sin(phi), std::sqrt(std::max(Real(0), Real(1) - r * r)));
            byEval.add(wi.z > Real(1e-3) ? double(average(skin->eval(rec, wo, wi, lambda, TransportMode::Radiance))) * constants::PI : 0.0);
        }
        bySample.finish(); byEval.finish();
        const double se = std::sqrt(bySample.se * bySample.se + byEval.se * byEval.se);
        report(std::string((variant & 2) ? "matte 2" : "dewy 1") + " at thickness 1 on skin" + (sss ? " with subsurface scattering" : "") + ": albedo from sample() and from eval()",
            std::abs(bySample.mean - byEval.mean) <= 4 * se + 0.02 * byEval.mean && byEval.mean < 1.0,
            Fixed(bySample.mean) + " and " + Fixed(byEval.mean) + " +- " + Fixed(se));
    }
    {   // (c2) on a skin with relief: the layer's normal is the smooth one when it levels (1) and the bent one when it does not (0), and
        //      sample() and eval() still describe the same reflection when the two normals differ
        const auto relief = std::make_shared<SkinRelief>(SkinRelief::Params{});
        CosmeticLayer::Params p;
        CosmeticLayer::preset("dewy1", p);
        p.thickness = 1.0; p.diffuserFraction = 0.6;
        const auto coat = std::make_shared<CosmeticCoat>(p);
        const Vector3 wo = Dir(40, 0, false);
        for (Real leveling : { Real(1), Real(0), Real(0.5) }) {
            auto skin = makeSkin(false);
            skin->setRelief(relief);
            skin->setCosmetic(coat, leveling);
            double worstTilt = 0, meanBent = 0, sumS = 0, sumE = 0, var = 0;
            const int points = 16, per = 20000 * scale;
            for (int k = 0; k < points; ++k) {
                SurfaceInteraction rec = flatHit(skin.get(), wo);
                rec.p = Point3(sampling::Random() * 0.005, sampling::Random() * 0.005, 0.0);
                skin->perturbShading(rec, wo);
                const Vector3 nl = glm::dot(rec.nl, rec.nl) > 0 ? rec.nl : rec.n;
                worstTilt = std::max(worstTilt, double(std::acos(std::clamp(nl.z, Real(-1), Real(1)))) * 180 / constants::PI);
                meanBent += double(std::acos(std::clamp(rec.n.z, Real(-1), Real(1)))) * 180 / constants::PI / points;
                Estimate bySample, byEval;
                for (int i = 0; i < per; ++i) {
                    const auto bs = skin->sample(rec, wo, Point2(sampling::Random(), sampling::Random()), lambda, TransportMode::Radiance);
                    bySample.add(bs ? double(average(bs->f)) * std::abs(glm::dot(rec.n, bs->wi)) / bs->pdf : 0.0);
                    const Real r = std::sqrt(sampling::Random()), phi = Real(2) * constants::PI * sampling::Random();
                    const Vector3 wi(r * std::cos(phi), r * std::sin(phi), std::sqrt(std::max(Real(0), Real(1) - r * r)));
                    byEval.add(wi.z > Real(1e-3) ? double(average(skin->eval(rec, wo, wi, lambda, TransportMode::Radiance))) * std::abs(glm::dot(rec.n, wi)) * constants::PI / wi.z : 0.0);
                }
                bySample.finish(); byEval.finish();
                sumS += bySample.mean; sumE += byEval.mean; var += bySample.se * bySample.se + byEval.se * byEval.se;
            }
            const double mS = sumS / points, mE = sumE / points, se = std::sqrt(var) / points;
            const bool normalOk = leveling >= 1 ? worstTilt < 0.01 : (leveling <= 0 ? true : worstTilt > 0.01);
            report("relief, leveling " + Fixed(leveling, 1) + ": the layer's normal, and the albedo from sample() and from eval()", normalOk && std::abs(mS - mE) <= 4 * se + 0.02 * mE && mE < 1.0,
                "largest tilt of the layer's normal " + Fixed(worstTilt, 2) + " deg (skin normal: " + Fixed(meanBent, 2) + " deg on average); " + Fixed(mS) + " and " + Fixed(mE) + " +- " + Fixed(se));
        }
    }
    {   // (d) white furnace: a white ball under a layer that absorbs nothing stays white (the skin's mirror-like lobe is only seen through the light that
        //     crosses the layer without a collision, so a few percent may be missing), and a skin ball under a foundation is never brighter than 1
        auto sky = std::make_shared<EnvMap>(Image(32, 16, std::vector<Vector3>(32 * 16, Vector3(1.0))));
        auto furnace = [&](const std::shared_ptr<CoatedDiffuse>& m, Vector3& mean, Vector3& se, int samples = 40000) {
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Sphere>(Point3(0, 0, 0), 1.0, m));
            auto scene = std::make_shared<Scene>(world);
            scene->setEnvLight(std::make_shared<EnvLight>(sky, 1.0));
            scene->finalize("spatial");
            int bad = 0;
            EstimateRadiance(*scene, Ray(Point3(0.3, 0.2, 4), glm::normalize(Vector3(-0.07, -0.05, -1))), samples * scale, 16, mean, se, bad);
            return bad;
        };
        CosmeticLayer::Params white;
        white.thickness = 1.0; white.diffuserFraction = 0.7; white.diffuserAlbedo = Vector3(1.0); white.plateletAlbedo = Vector3(1.0);
        auto ball = std::make_shared<CoatedDiffuse>(RgbColor(1.0, 1.0, 1.0), Real(0.5), Real(0.4), Real(0));
        Vector3 m0, s0, m1, s1, m2, s2;
        int bad = furnace(ball, m0, s0);
        ball->setCosmetic(std::make_shared<CosmeticCoat>(white));
        bad += furnace(ball, m1, s1);
        report("white furnace: a white ball under a layer without absorption", bad == 0 && m1.y <= m0.y + 4 * (s0.y + s1.y) + 0.01 && m1.y >= 0.93 * m0.y,
            "bare " + Fixed(m0.y) + ", with the layer " + Fixed(m1.y) + " +- " + Fixed(s1.y));
        CosmeticLayer::Params dewy;
        CosmeticLayer::preset("dewy1", dewy);
        dewy.thickness = 1.0;
        auto skin = makeSkin(true);
        Vector3 mb, sb;
        bad = furnace(skin, mb, sb);
        skin->setCosmetic(std::make_shared<CosmeticCoat>(dewy));
        bad += furnace(skin, m2, s2);
        report("white furnace: a skin ball with subsurface scattering under dewy 1 (thickness 1) is not brighter than 1", bad == 0 && std::max({ m2.x, m2.y, m2.z }) <= 1.0 + 4 * s2.x,
            "bare (" + Fixed(mb.x, 3) + ", " + Fixed(mb.y, 3) + ", " + Fixed(mb.z, 3) + "), with the layer (" + Fixed(m2.x, 3) + ", " + Fixed(m2.y, 3) + ", " + Fixed(m2.z, 3) + ")");

        // (e) Information: in a uniform white sky only the albedo of the ball counts, so a layer that absorbs nothing should add the same light
        //     to a skin with subsurface scattering as to one without (both balls have the same albedo). The share of the gain that is missing
        //     on the scattering skin is the energy the assembly of layer and skin loses. Both skins share one layer table, so its error cancels.
        std::cout << "[Cosmetic] 8. information: white furnace, a layer without absorption on a skin ball (0.80, 0.52, 0.42) without and with subsurface scattering\n";
        std::cout << "[Cosmetic]   (diffuser share, thickness) | bare r,g,b without / with | gain r,g,b without | gain r,g,b with | share of the gain missing r,g,b (+- 1 se)\n";
        auto triple = [](const Vector3& v, int digits = 3) { return Fixed(v.x, digits) + "," + Fixed(v.y, digits) + "," + Fixed(v.z, digits); };
        const int many = 200000;
        auto plainSkin = makeSkin(false), scatteringSkin = makeSkin(true);
        Vector3 bareP, bareS, seBareP, seBareS;
        int badE = furnace(plainSkin, bareP, seBareP, many) + furnace(scatteringSkin, bareS, seBareS, many);
        for (Real cd : { Real(1.0), Real(0.7) })
            for (Real t : { Real(0.25), Real(1.0), Real(4.0) }) {
                CosmeticLayer::Params w;
                w.thickness = t; w.diffuserFraction = cd; w.diffuserAlbedo = Vector3(1.0); w.plateletAlbedo = Vector3(1.0);
                const auto coat = std::make_shared<CosmeticCoat>(w);
                auto a = makeSkin(false), b = makeSkin(true);
                a->setCosmetic(coat); b->setCosmetic(coat);
                Vector3 ma, sa, mbb, sbb;
                badE += furnace(a, ma, sa, many) + furnace(b, mbb, sbb, many);
                const Vector3 gainP = ma - bareP, gainS = mbb - bareS;
                Vector3 missing(0), missingSe(0);
                for (int c = 0; c < 3; ++c) {
                    missing[c] = 1 - gainS[c] / gainP[c];
                    const Real seP = std::hypot(sa[c], seBareP[c]), seS = std::hypot(sbb[c], seBareS[c]);
                    missingSe[c] = (gainS[c] / gainP[c]) * std::hypot(seP / gainP[c], seS / gainS[c]);
                }
                std::cout << "[Cosmetic]   (" << Fixed(cd, 1) << ", " << Fixed(t, 2) << ") | " << triple(bareP) << " / " << triple(bareS) << " | " << triple(gainP) << " | " << triple(gainS)
                    << " | " << triple(missing * Real(100), 1) << " % (+- " << triple(missingSe * Real(100), 1) << ")\n";
            }
        if (badE > 0) std::cout << "[Cosmetic]   (" << badE << " samples with a bad value were dropped)\n";
    }

    std::cout << "[Cosmetic] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
