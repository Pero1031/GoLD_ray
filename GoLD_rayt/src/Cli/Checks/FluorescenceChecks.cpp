// Fluorescence checks (--scene fluorescence): the fluorescent material (Materials/Fluorescent.hpp) and the path tracer's wavelength shift.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/LightPresets.hpp"
#include "Cli/SceneParts.hpp"
#include "Core/Sampling.hpp"
#include "Geometry/Sphere.hpp"
#include "Lights/AreaLight.hpp"
#include "Materials/DiffuseLight.hpp"
#include "Materials/Fluorescent.hpp"
#include "Render/Integrator.hpp"
#include "Scene/Primitive.hpp"

#include <filesystem>
#include <sstream>

namespace rayt::cli {

namespace {

    /// A synthetic matrix: a flat reflectance and a fluorescence band (absorption around 360 nm, emission around 450 nm), RIT layout.
    std::shared_ptr<BispectralMatrix> SyntheticMatrix() {
        auto m = std::make_shared<BispectralMatrix>();
        m->title = "synthetic: reflectance 0.8, absorbs near 360 nm, emits near 450 nm";
        m->emFirst = 380; m->emStep = 10; m->nEm = 41;
        m->exFirst = 300; m->exStep = 10; m->nEx = 49;
        m->d.assign(static_cast<size_t>(m->nEm) * static_cast<size_t>(m->nEx), Real(0));
        for (int r = 0; r < m->nEm; ++r)
            for (int c = 0; c < m->nEx; ++c) {
                const Real lo = m->em(r), li = m->ex(c);
                Real v = 0;
                if (c == m->diagonalColumn(r)) v = Real(0.8);
                else if (li < lo) v = Real(0.06) * std::exp(-math::sqr((li - 360) / 30)) * std::exp(-math::sqr((lo - 450) / 30));
                m->d[static_cast<size_t>(r) * static_cast<size_t>(m->nEx) + static_cast<size_t>(c)] = v;
            }
        return m;
    }

    bool Report(const std::string& name, bool ok, const std::string& detail) {
        std::cout << "[Fluo]   " << name << "\n             " << detail << (ok ? "   OK" : "   FAIL") << "\n";
        return ok;
    }

    /// The excitation sampler against the matrix: the share of samples in each column's band.
    bool CheckSampler(const Fluorescent& mat, Real lambdaO, int n) {
        const BispectralMatrix& M = mat.matrix();
        const Real F = mat.fluorescence(lambdaO);
        std::vector<int> hits(static_cast<size_t>(M.nEx), 0);
        int outside = 0;
        for (int i = 0; i < n; ++i) {
            const Real li = mat.sampleExcitation(lambdaO, sampling::Random());
            const int c = static_cast<int>(std::floor((li - M.exFirst) / M.exStep + Real(0.5)));
            if (c < 0 || c >= M.nEx) { ++outside; continue; }
            ++hits[static_cast<size_t>(c)];
        }
        double worst = 0;
        for (int c = 0; c < M.nEx; ++c) {
            const double p = double(mat.excitationDensity(lambdaO, M.ex(c)) * M.exStep / F);
            const double se = std::sqrt(std::max(p * (1 - p), 1e-12) / n);
            worst = std::max(worst, std::abs(double(hits[static_cast<size_t>(c)]) / n - p) / (se + 1e-4));
        }
        std::ostringstream d;
        d << "lambda_o " << lambdaO << " nm, F " << F << ", " << n << " samples: largest |share - D/F| = " << worst << " (standard errors + 1e-4)";
        return Report("excitation sampler against the matrix", worst <= 4 && outside == 0, d.str());
    }

    /// reemit(): the albedo is R + F, and it shifts with probability F / (R + F).
    bool CheckReemit(const Fluorescent& mat, Real lambdaO, int n) {
        const Real R = mat.reflectance(lambdaO), F = mat.fluorescence(lambdaO);
        int shifted = 0;
        bool albedoOk = true;
        for (int i = 0; i < n; ++i) {
            SurfaceInteraction rec;
            SampledWavelengths lambda;
            lambda.lambda[0] = lambdaO;
            mat.reemit(rec, lambda, sampling::Random(), sampling::Random());
            if (lambda.lambda[0] != lambdaO) ++shifted;
            albedoOk = albedoOk && std::abs(rec.reemitAlbedo - (R + F)) < Real(1e-6);
        }
        const double p = double(F / (R + F)), got = double(shifted) / n, se = std::sqrt(p * (1 - p) / n);
        std::ostringstream d;
        d << "lambda_o " << lambdaO << " nm: albedo R + F = " << R + F << (albedoOk ? "" : " (NOT always set)") << ", shifted " << got
          << " of the paths, expected F / (R + F) = " << p << " +- " << se;
        return Report("reemit(): albedo and the share of shifted paths", albedoOk && std::abs(got - p) <= 4 * se + 1e-4, d.str());
    }

    /// The radiance of a convex fluorescent sphere inside an emitting sphere of radiance Le(lambda): every point receives Le from all
    /// directions, so it sends L(lo) = R(lo) Le(lo) + Integral b(lo, li) Le(li) dli; the film's RGB is the integral of L against its sensor.
    bool CheckFurnace(const std::shared_ptr<Fluorescent>& mat, const EmissionSpectrum& light, const std::string& name, int n) {
        auto world = std::make_shared<Aggregate>();
        auto scene = std::make_shared<Scene>(world);
        auto emitter = std::make_shared<DiffuseLight>(light, Real(1), true);   // two-sided: it shines inward
        auto shell = std::make_shared<Sphere>(Point3(0.0), Real(5), emitter);
        auto area = std::make_shared<AreaLight>(shell, emitter);
        world->add(std::make_shared<Primitive>(shell, emitter, area));
        scene->addLight(area);
        world->add(std::make_shared<Sphere>(Point3(0.0), Real(1), mat));
        scene->finalize("spatial");

        Vector3 expected(0.0);
#if defined(RAYT_SPECTRAL_MODE)
        for (int l = 360; l <= 830; ++l) {
            const Real lo = Real(l);
            Real L = mat->reflectance(lo) * light.at(lo);
            for (Real li = Real(295.125); li < Real(785); li += Real(0.25)) L += mat->excitationDensity(lo, li) * light.at(li) * Real(0.25);
            const Real w = (l == 360 || l == 830) ? Real(0.5) : Real(1);
            expected += color::rgbSensorWeights(lo) * (L * w);
        }
#else
        expected = mat->rgbUnderD65() * light.rgb();
#endif
        Vector3 mean, se;
        int bad = 0;
        EstimateRadiance(*scene, Ray(Point3(0, 0, -3), Vector3(0, 0, 1)), n, 64, mean, se, bad);
        bool ok = bad == 0;
        for (int c = 0; c < 3; ++c) ok = ok && std::abs(mean[c] - expected[c]) <= 4 * se[c] + Real(0.003) * std::abs(expected[c]) + Real(1e-6);
        std::ostringstream d;
        d << "expected (" << expected.x << ", " << expected.y << ", " << expected.z << ")   got (" << mean.x << ", " << mean.y << ", " << mean.z
          << ") +- (" << se.x << ", " << se.y << ", " << se.z << ")" << (bad ? "  NaN samples: " + std::to_string(bad) : std::string());
        return Report(name, ok, d.str());
    }

} // namespace

int RunFluorescenceTests(const CliOptions& opt) {
    const int n = std::max(1, opt.spp) * 64 * 1024;   // --spp scales the samples
    bool all = true;
    std::cout << "[Fluo] fluorescent material (Materials/Fluorescent.hpp), " <<
#if defined(RAYT_SPECTRAL_MODE)
        "spectral build"
#else
        "RGB build (no wavelength shift: the surface is a diffuse color under D65)"
#endif
        << "\n";

    EmissionSpectrum d65uv;
    MakeLightEmission("d65uv", d65uv);
    const auto synthetic = SyntheticMatrix();
    const auto mat = std::make_shared<Fluorescent>(synthetic, true);
    const auto off = std::make_shared<Fluorescent>(synthetic, false);

    all = CheckSampler(*mat, Real(455), 400000) && all;
    all = CheckReemit(*mat, Real(452), 400000) && all;
    all = CheckFurnace(mat, d65uv, "synthetic matrix inside an emitting sphere (D65 with ultraviolet): R Le + Integral b Le", n) && all;
    all = CheckFurnace(off, d65uv, "the same with the fluorescence off: R Le", n) && all;

    // The RIT database, if it is on this PC (tools/data/rit_bispectral.py, docs/FLUORESCENCE.md).
    const std::string dir = "assets/Fluorescence/rit_bispectral/Labsphere/";
    for (const char* file : { "CIPLAW8.BFC", "HERPIORA.BFC" }) {
        const std::string path = dir + file;
        if (!std::filesystem::exists(path)) { std::cout << "[Fluo]   " << path << " not found: the measured matrices are not checked (information)\n"; continue; }
        auto m = std::make_shared<BispectralMatrix>();
        BispectralCleaning cleaning;
        std::string err;
        if (!LoadBfcMatrix(path, *m, &cleaning, &err)) { all = Report("read " + path, false, err); continue; }
        const bool shape = m->nEm == 41 && m->nEx == 49 && m->emFirst == 380 && m->exFirst == 300;
        std::ostringstream d;
        d << "\"" << m->title << "\": " << m->nEm << " x " << m->nEx << ", set to 0: " << cleaning.negative << " negative (sum " << cleaning.negativeSum
          << "), " << cleaning.antiStokes << " anti-Stokes (sum " << cleaning.antiStokesSum << "), noise " << cleaning.noise << ": "
          << cleaning.belowNoise << " fluorescence entries below 3 x noise (sum " << cleaning.belowNoiseSum << ")";
        all = Report("read " + path, shape, d.str()) && all;
        all = CheckFurnace(std::make_shared<Fluorescent>(m, true), d65uv, std::string(file) + " inside the emitting sphere", n) && all;
    }

    std::cout << "[Fluo] RESULT: " << (all ? "PASS" : "FAIL") << "\n";
    return all ? 0 : 2;
}

} // namespace rayt::cli
