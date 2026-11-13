// Helpers of the self-checks.
#include "pch.h"
#include "Cli/CheckSupport.hpp"
#include "Cli/SceneParts.hpp"
#include "Geometry/Sphere.hpp"
#include "Materials/Lambertian.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Integrator.hpp"
#include "Render/SpecularNee.hpp"

namespace rayt::cli {

void EstimateRadiance(const Scene& scene, const Ray& ray, int n, int maxDepth,
    Vector3& mean, Vector3& se, int& badSamples, const CausticPhotonMap* caustics, const SpecularNee* nee, bool skipDirectCaustics, bool mediaNee)
{
    PathIntegrator integrator(nullptr, maxDepth, 1);
    integrator.setCausticMap(caustics);
    integrator.setSkipDirectCaustics(skipDirectCaustics);
    integrator.setSpecularNee(nee);
    integrator.setMediaNee(mediaNee);
    double sx = 0, sy = 0, sz = 0, qx = 0, qy = 0, qz = 0;
    int bad = 0;

    // OpenMP 2.0 (MSVC): scalar reductions and a signed int loop counter only.
#pragma omp parallel for reduction(+:sx,sy,sz,qx,qy,qz,bad)
    for (int i = 0; i < n; ++i) {
#if defined(RAYT_SPECTRAL_MODE)
        SampledWavelengths lambda = SampledWavelengths::sampleVisible(sampling::Random());
#else
        SampledWavelengths lambda;
#endif
        const Spectrum L = integrator.Li(ray, scene, lambda);
        const Vector3 c = SpectrumToRGB(L, lambda);
        if (math::hasNonFinite(c)) { ++bad; continue; }
        sx += c.x; sy += c.y; sz += c.z;
        qx += c.x * c.x; qy += c.y * c.y; qz += c.z * c.z;
    }

    const double inv = 1.0 / double(n);
    mean = Vector3(sx, sy, sz) * inv;
    const Vector3 var(std::max(0.0, qx * inv - mean.x * mean.x),
        std::max(0.0, qy * inv - mean.y * mean.y),
        std::max(0.0, qz * inv - mean.z * mean.z));
    se = Vector3(std::sqrt(var.x * inv), std::sqrt(var.y * inv), std::sqrt(var.z * inv));
    badSamples = bad;
}

Vector3 PhotonVarianceAt(const Scene& scene, const CausticPhotonMap& map, const Ray& ray) {
    SurfaceInteraction rec;
    if (!scene.hit(ray, rec)) return Vector3(0.0);
#if defined(RAYT_SPECTRAL_MODE)
    const SampledWavelengths lambda = SampledWavelengths::sampleVisible(Real(0.5));
#else
    const SampledWavelengths lambda;
#endif
    Vector3 variance(0.0);
    map.radiance(rec, -ray.d, lambda, &variance);
    return variance;
}

bool ReportLightCase(const std::string& name, const Vector3& expected,
    const Vector3& mean, const Vector3& se, int bad)
{
    bool ok = (bad == 0);
    for (int c = 0; c < 3; ++c) {
        ok = ok && std::abs(mean[c] - expected[c]) <= 4.0 * se[c] + 0.005 * std::abs(expected[c]) + 1e-6;
    }
    std::cout << "[Lights]   " << name << "\n"
        << "             expected (" << expected.x << ", " << expected.y << ", " << expected.z << ")"
        << "   got (" << mean.x << ", " << mean.y << ", " << mean.z << ") +- ("
        << se.x << ", " << se.y << ", " << se.z << ")"
        << "   ratio (" << mean.x / expected.x << ", " << mean.y / expected.y << ", " << mean.z / expected.z << ")"
        << (bad ? "  NaN samples: " + std::to_string(bad) : std::string())
        << "   " << (ok ? "OK" : "FAIL") << "\n";
    return ok;
}

std::shared_ptr<Scene> MakeGroundScene(Real albedo,
    const std::function<void(Aggregate&, Scene&)>& addLights, const char* lightSampler)
{
    auto world = std::make_shared<Aggregate>();
    world->add(std::make_shared<Sphere>(Point3(0, -1000, 0), 1000.0,
        std::make_shared<Lambertian>(RgbColor(albedo))));
    auto scene = std::make_shared<Scene>(world);
    addLights(*world, *scene);
    scene->finalize(lightSampler);
    return scene;
}

Real QuadGroundRadiance(const Point3& corner, const Vector3& e1, const Vector3& e2,
    bool twoSided, Real L, Real albedo)
{
    const Vector3 nn = glm::cross(e1, e2);
    const Real area = glm::length(nn);
    const Vector3 n = nn / area;
    const int N = 600;
    Real E = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            const Point3 q = corner + e1 * ((Real(i) + 0.5) / N) + e2 * ((Real(j) + 0.5) / N);
            const Vector3 d = q;                         // from the origin to the light point
            const Real r2 = glm::dot(d, d), r = std::sqrt(r2);
            const Real cosR = d.y / r;                   // receiver normal +y
            const Real cosL = glm::dot(n, -d) / r;       // > 0: the emitting side faces the receiver
            if (cosR <= 0 || (cosL <= 0 && !twoSided)) continue;
            E += cosR * std::abs(cosL) / r2;
        }
    }
    E *= L * area / Real(N * N);
    return albedo / constants::PI * E;
}

std::shared_ptr<EnvMap> MakeSunSkyEnv(Vector3& sunDirOut) {
    const int W = 256, H = 128;
    const Real elev = glm::radians(50.0), azim = glm::radians(30.0);
    const Vector3 sun(std::cos(elev) * std::cos(azim), std::sin(elev), std::cos(elev) * std::sin(azim));
    sunDirOut = sun;

    std::vector<Vector3> px(static_cast<size_t>(W) * H);
    for (int y = 0; y < H; ++y) {
        const Real theta = constants::PI * (Real(y) + 0.5) / Real(H);     // row 0 = zenith
        for (int x = 0; x < W; ++x) {
            const Real phi = 2.0 * constants::PI * (Real(x) + 0.5) / Real(W) - constants::PI;
            const Vector3 d(std::sin(theta) * std::cos(phi), std::cos(theta), std::sin(theta) * std::sin(phi));
            Vector3 c = Vector3(0.8, 0.9, 1.0) * (0.2 + 0.3 * std::max(d.y, 0.0));
            if (std::acos(std::clamp(glm::dot(d, sun), -1.0, 1.0)) < 0.05) c += Vector3(500.0, 450.0, 400.0);
            px[static_cast<size_t>(y) * W + x] = c;
        }
    }
    return std::make_shared<EnvMap>(Image(W, H, std::move(px)));
}

Vector3 GroundRadianceFromEnv(const EnvMap& env, Real albedo) {
    const int nTheta = 1024, nPhi = 2048;
    Vector3 E(0.0);
    for (int i = 0; i < nTheta; ++i) {
        const Real theta = 0.5 * constants::PI * (Real(i) + 0.5) / Real(nTheta);
        const Real dTheta = 0.5 * constants::PI / Real(nTheta);
        const Real dPhi = 2.0 * constants::PI / Real(nPhi);
        for (int j = 0; j < nPhi; ++j) {
            const Real phi = 2.0 * constants::PI * (Real(j) + 0.5) / Real(nPhi);
            const Vector3 d(std::sin(theta) * std::cos(phi), std::cos(theta), std::sin(theta) * std::sin(phi));
            E += env.eval(d) * (std::cos(theta) * std::sin(theta) * dTheta * dPhi);
        }
    }
    return E * (albedo / constants::PI);
}

bool CheckEnvSampling(const EnvMap& env, const std::string& label, int n) {
    // Integral of pdf over the sphere on a deterministic (theta, phi) grid (a small bright sun
    // would make a Monte Carlo estimate of this far too noisy). The integral of Le by the same grid at two resolutions: the difference between
    // them bounds the grid's own error (a small sun is not resolved exactly by either).
    double pdfIntegral = 0.0;
    Vector3 gridLe(0.0), gridLeCoarse(0.0);
    for (int level = 0; level < 2; ++level) {
        const int nTheta = level == 0 ? 1024 : 512, nPhi = 2 * nTheta;
        const double dTheta = constants::PI / nTheta, dPhi = 2.0 * constants::PI / nPhi;
        Vector3 acc(0.0);
        for (int i = 0; i < nTheta; ++i) {
            const Real theta = (Real(i) + 0.5) * dTheta;
            for (int j = 0; j < nPhi; ++j) {
                const Real phi = (Real(j) + 0.5) * dPhi;
                const Vector3 d(std::sin(theta) * std::cos(phi), std::cos(theta), std::sin(theta) * std::sin(phi));
                if (level == 0) pdfIntegral += env.pdf(d) * std::sin(theta) * dTheta * dPhi;
                acc += env.eval(d) * (std::sin(theta) * dTheta * dPhi);
            }
        }
        (level == 0 ? gridLe : gridLeCoarse) = acc;
    }
    const double gridErr = std::abs(gridLe.y - gridLeCoarse.y);

    // The importance-sampled integral and its standard error (Le / pdf varies: the pdf is constant over a texel, Le is interpolated).
    double sum = 0, sum2 = 0;
    int mismatches = 0;   // sample() and pdf() disagree (> 0.1 %)
    int isCount = 0;
    for (int i = 0; i < n; ++i) {
        Vector3 wi(0.0);
        Real pdf = 0;
        const Vector3 Le = env.sample(Point2f(float(sampling::Random()), float(sampling::Random())), wi, pdf);
        if (pdf > 0) {
            const double v = Le.y / pdf;
            sum += v; sum2 += v * v;
            if (std::abs(env.pdf(wi) - pdf) > 1e-3 * pdf) ++mismatches;   // only at texel borders (float round-off)
            ++isCount;
        }
    }
    const double isLe = sum / double(n);
    const double isSe = std::sqrt(std::max(sum2 / double(n) - isLe * isLe, 0.0) / double(n));
    const double mismatchRate = double(mismatches) / std::max(isCount, 1);

    const bool okPdf = std::abs(pdfIntegral - 1.0) < 0.01;
    const bool okMatch = mismatchRate < 1e-3;
    const bool okLe = std::abs(isLe - gridLe.y) <= 4 * isSe + 2 * gridErr + 1e-9;   // a fixed 3 % failed a few per cent of the runs
    std::cout << "[Lights]   env sampling (" << label << "): Integral pdf dw = " << pdfIntegral
        << (okPdf ? " OK" : " FAIL") << ";  sample()/pdf() mismatch rate = " << mismatchRate
        << (okMatch ? " OK" : " FAIL") << ";  Integral Le dw (G): grid " << gridLe.y << " (grid error about " << gridErr
        << ") vs importance " << isLe << " +- " << isSe << (okLe ? " OK" : " FAIL") << "\n";
    return okPdf && okMatch && okLe;
}

} // namespace rayt::cli
