// GPU checks (--scene gpu, docs/GPU.md): the kernel code (include/Kernel) on the CPU against the engine's functions,
// and on the GPU against the CPU. Without the GPU module or a CUDA device the GPU part prints SKIP.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Accelerators/BVH.hpp"
#include "Core/Fresnel.hpp"
#include "Core/Sampling.hpp"
#include "Geometry/Frame.hpp"
#include "Geometry/Sphere.hpp"
#include "Cli/KernelScenePack.hpp"
#include "Cli/Scene.hpp"
#include "Cli/SceneFile.hpp"
#include "Gpu/GpuDevice.hpp"
#include "Kernel/Lights.hpp"
#include "Kernel/PathTracer.hpp"
#include "Kernel/SelfTest.hpp"
#include "Lights/EnvLight.hpp"
#include "Materials/Dielectric.hpp"
#include "Materials/LayeredCoat.hpp"
#include "Materials/Medium.hpp"
#include "Materials/RoughConductor.hpp"
#include "Render/Camera.hpp"
#include "Render/Integrator.hpp"
#include "Scene/Scene.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace rayt::cli {

namespace {

    bool Report(const std::string& name, bool ok, const std::string& detail) {
        std::cout << "[GPU]    " << name << "\n             " << detail << (ok ? "   OK" : "   FAIL") << "\n";
        return ok;
    }

    std::string Fmt(double v, int prec = 3) {
        std::ostringstream s;
        s << std::setprecision(prec) << v;
        return s.str();
    }

    /// Cases for evalMathCase: random ones plus edges (u = 0, poles of the basis, -0, total internal reflection, cos < 0).
    std::vector<kernel::MathCase> MakeCases(int n) {
        std::vector<kernel::MathCase> cases;
        kernel::Pcg32 rng;
        rng.seed(20261126, 7);
        auto add = [&](float u1, float u2, float nx, float ny, float nz, float cosT, float etaI, float etaT) {
            kernel::MathCase c{};
            c.seed = 0x5eedull + uint64_t(cases.size()) * 977u;
            c.pixel = uint32_t(cases.size() * 31u);
            c.sample = uint32_t(cases.size() % 1024u);
            c.u1 = u1; c.u2 = u2; c.nx = nx; c.ny = ny; c.nz = nz; c.cosTheta = cosT; c.etaI = etaI; c.etaT = etaT;
            cases.push_back(c);
        };
        add(0, 0, 0, 0, 1, 1, 1, 1.5f);
        add(0.5f, 0.999999f, 0, 0, -1, 0, 1.5f, 1);
        add(0.25f, 0.5f, 0.6f, 0.8f, -0.0f, -0.3f, 1, 1.33f);   // n.z = -0: copysign picks the lower branch
        add(0.75f, 0.25f, 1, 0, 0, 0.2f, 1.33f, 1);
        const float etas[4][2] = { { 1, 1.5f }, { 1.5f, 1 }, { 1, 1.33f }, { 1.33f, 1 } };
        while (int(cases.size()) < n) {
            const float nz = 1 - 2 * rng.nextFloat(), phi = 6.2831853f * rng.nextFloat(), r = std::sqrt(std::max(0.0f, 1 - nz * nz));
            const int e = int(rng.nextUint() % 4u);
            add(rng.nextFloat(), rng.nextFloat(), r * std::cos(phi), r * std::sin(phi), nz, -0.2f + 1.2f * rng.nextFloat(), etas[e][0], etas[e][1]);
        }
        return cases;
    }

    /// The engine's functions for the same case, in double: what the kernel's double values must equal.
    void EngineValues(const kernel::MathCase& c, double* v) {
        const Point2 u(Real(c.u1), Real(c.u2));
        const Vector3 h = sampling::CosineSampleHemisphere(u);
        const Vector3 s = sampling::UniformSampleSphere(u);
        Vector3 n(Real(c.nx), Real(c.ny), Real(c.nz));
        n = n / std::sqrt(n.x * n.x + n.y * n.y + n.z * n.z);   // as the kernel normalizes (glm::normalize rounds differently)
        Vector3 bs, bt;
        frame::makeOrthonormalBasis(n, bs, bt);
        const double vals[] = { h.x, h.y, h.z, s.x, s.y, s.z, bs.x, bs.y, bs.z, bt.x, bt.y, bt.z,
            fresnel::fresnelDielectric(Real(c.cosTheta), Real(c.etaI), Real(c.etaT)) };
        for (int k = kernel::kHemiX; k <= kernel::kFresnel; ++k) v[k] = vals[k - kernel::kHemiX];
    }

    // ---- image checks of the kernel renderer, on scenes/gpu_box.json ----

    /// A test scene file as a kernel scene of w x h pixels; `furnace`: every diffuse albedo 1, no lamp, a uniform environment of
    /// radiance 1 (every pixel then converges to exactly 1).
    bool KernelTestScene(const std::string& path, int w, int h, int spp, int maxDepth, uint64_t seed, bool furnace, KernelSceneHost& host, std::string& why) {
        CliOptions o;
        o.spp = spp;
        o.maxDepth = maxDepth;
        o.seed = static_cast<long long>(seed);
        const auto edit = [furnace](SceneDescription& d) {
            if (!furnace) return;
            d.environment.type = EnvironmentSpec::Type::Uniform;
            d.environment.scale = 1.0;
            for (auto& [name, m] : d.materials) if (m.type == MaterialSpec::Type::Diffuse) m.albedo = Vector3(1.0);
            std::vector<ObjectSpec> kept;
            for (const ObjectSpec& ob : d.objects) if (d.materials.at(ob.material).type != MaterialSpec::Type::Emission) kept.push_back(ob);
            d.objects = kept;
        };
        return PackSceneFile(path, o, w, h, edit, host, why);
    }

    /// The same scene and seed on the CPU kernel renderer and on the GPU: the same bits (fast math: the same mean within the noise), the
    /// denoiser's features (albedo, normal) included; and a render without the features is the same picture (recording them draws no
    /// random number).
    bool PairedCheck(const std::string& path, bool fastMath) {
        KernelSceneHost host;
        std::string why;
        if (!KernelTestScene(path, 128, 96, 8, 16, 11, false, host, why)) return Report(path, false, why);
        std::vector<float> a, b, plain, albedoA, normalA, albedoB, normalB;
        const auto t0 = std::chrono::steady_clock::now();
        const bool cpuOk = RenderKernelCpu(host.view(), a, why, false, &albedoA, &normalA);
        const auto t1 = std::chrono::steady_clock::now();
        const bool gpuOk = cpuOk && RenderKernelGpu(host.view(), b, why, false, nullptr, nullptr, &albedoB, &normalB);
        const auto t2 = std::chrono::steady_clock::now();
        const bool plainOk = gpuOk && RenderKernelCpu(host.view(), plain, why, false);
        if (!cpuOk || !gpuOk || !plainOk) return Report(path + " on the CPU and the GPU", false, why);
        size_t samePixels = 0, sameFeatures = 0, withFeatures = 0;
        double sa = 0, sb = 0, var = 0;
        for (size_t i = 0; i < a.size(); i += 3) {
            samePixels += a[i] == b[i] && a[i + 1] == b[i + 1] && a[i + 2] == b[i + 2];
            bool same = true, any = false;
            for (int c = 0; c < 3; ++c) {
                sa += a[i + c]; sb += b[i + c]; var += double(a[i + c]) * a[i + c];
                same &= albedoA[i + c] == albedoB[i + c] && normalA[i + c] == normalB[i + c];
                any |= albedoA[i + c] != 0.0f;
            }
            sameFeatures += same;
            withFeatures += any;
        }
        const bool unchanged = plain == a;
        const size_t pixels = a.size() / 3;
        const double n = double(a.size()), meanA = sa / n, se = std::sqrt(std::max(var / n - meanA * meanA, 0.0) / n);
        const bool pass = unchanged && (fastMath ? std::abs(sb - sa) / n <= 4 * std::sqrt(2.0) * se : samePixels == pixels && sameFeatures == pixels);
        return Report(std::string("kernel renderer, ") + path + " with the same seed: GPU " + (fastMath ? "(fast math) = CPU within the noise" : "= CPU bit for bit") +
            ", albedo and normal included", pass,
            Fmt(100.0 * double(samePixels) / double(pixels), 5) + " % of the pixels bit-identical, albedo + normal " + Fmt(100.0 * double(sameFeatures) / double(pixels), 5) +
            " % (an albedo in " + Fmt(100.0 * double(withFeatures) / double(pixels), 4) + " %), image mean GPU / CPU " + Fmt(sb / sa, 8) +
            "; without the features " + (unchanged ? "the same picture" : "A DIFFERENT PICTURE") + "; 128 x 96, 8 spp: CPU " +
            Fmt(std::chrono::duration<double>(t1 - t0).count(), 2) + " s, GPU " + Fmt(std::chrono::duration<double>(t2 - t1).count(), 2) + " s (with the GPU's setup)");
    }

    /// The denoiser's features of the kernel renderer against the engine's (PathIntegrator::Li with a PathAOV) on the same camera rays,
    /// one through the centre of each pixel of a w x h grid. Rays whose first surface is specular (smooth glass) are left out: both
    /// renderers choose reflection or refraction there at random. Every other ray records its features at its first surface or at the
    /// environment, so the two must agree up to float rounding (the kernel's geometry is the engine's in float), except where the
    /// tracers (the engine's BVH in double, Embree in float) send a ray along an edge or a silhouette to different primitives. The
    /// environment's color differs more: the kernel maps the direction to the map in float, the engine in double, and the bilinear
    /// lookup turns that difference (about 1e-7 rad, 3e-5 of a texel of a 2048 wide map) into one of the color, times the contrast
    /// between neighbouring texels.
    bool AovEngineCheck(const std::string& path, int w = 64, int h = 48) {
        const int maxDepth = 16;
        KernelSceneHost host;
        std::string why;
        std::vector<kernel::KRay> rays;
        std::vector<PathAOV> engine;
        const auto withEngine = [&](const SceneSetup& setup, const Environment& env, const Camera&) {
            struct Quiet {
                std::ostringstream sink;
                std::streambuf* old = std::cout.rdbuf(sink.rdbuf());
                ~Quiet() { std::cout.rdbuf(old); }
            } quiet;
            Scene scene(std::make_shared<BVH>(setup.world->objects, 4, /*recordInPlace=*/true));
            if (env.light) scene.setEnvLight(env.light);
            for (const auto& l : setup.lights) scene.addLight(l);
            scene.finalize(setup.lightSampler.empty() ? CliOptions{}.lightSampler : setup.lightSampler);
            const kernel::KCamera& c = host.header.camera;
            const kernel::LensArrays<float> lensArrays{ host.lensSurfaces.data(), host.lensPupil.data(), host.lensBladeNormals.data() };
            PathIntegrator integrator(nullptr, maxDepth, 1);
            for (int j = 0; j < h; ++j)
                for (int i = 0; i < w; ++i) {
                    const float u = (float(i) + 0.5f) / float(w), v = (float(j) + 0.5f) / float(h);
                    kernel::V3 o, d;
                    if (c.type == kernel::kCamPinhole) {
                        o = kernel::lens::v3(c.origin);
                        d = kernel::knormalize(kernel::lens::v3(c.lowerLeft) + u * kernel::lens::v3(c.horizontal) + v * kernel::lens::v3(c.vertical) - o);
                    }
                    else {
                        float weight;
                        kernel::lensCameraRay(c, lensArrays, u, v, 0.5f, 0.5f, o, d, weight);
                    }
                    rays.push_back({ { o.x, o.y, o.z }, { d.x, d.y, d.z }, 0.0f, kernel::kFloatMax });
                    SampledWavelengths lambda;
                    PathAOV aov;
                    integrator.Li(Ray(Point3(o.x, o.y, o.z), glm::normalize(Vector3(d.x, d.y, d.z))), scene, lambda, &aov);
                    engine.push_back(aov);
                }
        };
        CliOptions o;
        o.spp = 1;
        o.maxDepth = maxDepth;
        o.seed = 3;
        if (!PackSceneFile(path, o, w, h, {}, host, why, withEngine)) return Report("features " + path, false, why);
        std::vector<kernel::KAov> kernelAovs;
        std::vector<uint32_t> firstMaterial;
        if (!KernelRayAovs(host.view(), rays, kernelAovs, firstMaterial, why)) return Report("features " + path, false, why);
        size_t compared = 0, skipped = 0, differ = 0, valid = 0, atEnvironment = 0;
        double maxAlbedo = 0, maxNormal = 0, maxEnv = 0;
        std::ostringstream examples;
        for (size_t k = 0; k < rays.size(); ++k) {
            if (firstMaterial[k] != kernel::kNoHit && host.materials[firstMaterial[k]].aovSpecular) { ++skipped; continue; }
            ++compared;
            const kernel::KAov& a = kernelAovs[k];
            const PathAOV& e = engine[k];
            double dA = 0, dN = 0;
            for (int c = 0; c < 3; ++c) {
                dA = std::max(dA, std::abs(double(a.albedo[c]) - e.albedo[c]));
                dN = std::max(dN, std::abs(double(a.normal[c]) - e.normal[c]));
            }
            const bool env = e.valid && firstMaterial[k] == kernel::kNoHit;
            const bool same = (a.valid != 0) == e.valid && dA <= (env ? 1e-3 : 1e-5) && dN <= 1e-4;
            if (same) {
                (env ? maxEnv : maxAlbedo) = std::max(env ? maxEnv : maxAlbedo, dA);
                maxNormal = std::max(maxNormal, dN);
                valid += e.valid;
                atEnvironment += env;
            }
            else if (++differ <= 3)
                examples << "; pixel (" << k % size_t(w) << ", " << k / size_t(w) << "): kernel " << (a.valid ? "" : "no ") << "albedo " << Fmt(a.albedo[0]) << " "
                         << Fmt(a.albedo[1]) << " " << Fmt(a.albedo[2]) << " normal " << Fmt(a.normal[0]) << " " << Fmt(a.normal[1]) << " " << Fmt(a.normal[2])
                         << ", engine " << (e.valid ? "" : "no ") << "albedo " << Fmt(e.albedo.x) << " " << Fmt(e.albedo.y) << " " << Fmt(e.albedo.z)
                         << " normal " << Fmt(e.normal.x) << " " << Fmt(e.normal.y) << " " << Fmt(e.normal.z);
        }
        return Report("the kernel renderer's albedo / normal features = the engine's on the same rays, " + path, differ <= compared / 500,
            std::to_string(compared) + " rays (" + std::to_string(skipped) + " left out: smooth glass first), " + std::to_string(differ) +
            " differ (allowed: 1 in 500, edges and silhouettes); the others agree within " + Fmt(maxAlbedo, 2) + " (albedo), " + Fmt(maxNormal, 2) +
            " (normal) and " + Fmt(maxEnv, 2) + " (the environment's color, allowed 1e-3); " + std::to_string(valid) + " with features (" + std::to_string(atEnvironment) +
            " the environment)" + examples.str());
    }

    /// The rays a CPU path trace of a test scene casts (w x h pixels, spp samples, run seed), traced again by the GPU's tracer
    /// (OptiX): where do OptiX and Embree answer differently? The kernel takes only the primitive from them (and a sphere's distance,
    /// and whether a shadow ray is blocked), so these are the only places where the CPU's and the GPU's paths can part.
    bool TracerCheck(const gpu::Module& mod, const std::string& path, int w = 64, int h = 48, int spp = 4, uint64_t seed = 11) {
        KernelSceneHost host;
        std::string why;
        if (!KernelTestScene(path, w, h, spp, 16, seed, false, host, why)) return Report(path, false, why);
        const kernel::KScene s = host.view();
        std::vector<kernel::KRay> rays;
        std::vector<kernel::KHit> cpu;
        std::vector<uint8_t> shadow;
        if (!RecordPathRays(s, rays, cpu, shadow, why)) return Report(path, false, why);
        std::vector<kernel::KHit> gpu(rays.size());
        char err[512] = {};
        if (mod.traceTest(&s, rays.data(), gpu.data(), int(rays.size()), err, int(sizeof(err))) != 0) return Report("tracing " + path + " on the GPU", false, err);
        size_t closest = 0, shadows = 0, hitMiss = 0, otherTriangle = 0, shadowDiffer = 0, sphereHits = 0, sphereDistance = 0;
        std::ostringstream examples;
        int shown = 0;
        for (size_t i = 0; i < rays.size(); ++i) {
            const kernel::KHit& a = cpu[i];
            const kernel::KHit& b = gpu[i];
            if (shadow[i]) {
                ++shadows;
                const bool blockedCpu = a.prim != kernel::kNoHit, blockedGpu = b.prim != kernel::kNoHit;
                shadowDiffer += blockedCpu != blockedGpu;
                if (blockedCpu != blockedGpu && shown < 4) {
                    const kernel::KRay& r = rays[i];
                    examples << "\n               shadow ray " << i << ": Embree " << blockedCpu << ", OptiX " << blockedGpu << " o " << Fmt(r.o.x, 9) << " " << Fmt(r.o.y, 9) << " " << Fmt(r.o.z, 9)
                             << " d " << Fmt(r.d.x, 9) << " " << Fmt(r.d.y, 9) << " " << Fmt(r.d.z, 9) << " tmax " << Fmt(r.tmax, 9);
                    ++shown;
                }
                continue;
            }
            ++closest;
            if (a.prim == b.prim) {
                // a sphere: the kernel takes the distance from the tracer, which both compute with kernel::intersectSphere
                if (a.prim != kernel::kNoHit && a.prim >= s.numTriangles) { ++sphereHits; sphereDistance += a.t != b.t; }
                continue;
            }
            const bool hm = (a.prim == kernel::kNoHit) != (b.prim == kernel::kNoHit);
            (hm ? hitMiss : otherTriangle) += 1;
            if (shown < 4) {
                const kernel::KRay& r = rays[i];
                const uint32_t p = a.prim != kernel::kNoHit ? a.prim : b.prim;
                auto name = [&](const kernel::KHit& x) {
                    if (x.prim == kernel::kNoHit) return std::string("miss");
                    return (x.prim < s.numTriangles ? "triangle " + std::to_string(x.prim) : "sphere " + std::to_string(x.prim - s.numTriangles)) + " at t " + Fmt(x.t, 6);
                };
                examples << "\n               ray " << i << ": Embree " << name(a) << ", OptiX " << name(b);
                if (p < s.numTriangles) {
                    const kernel::KFloat3 n = s.faceNormals[p];
                    examples << "; cos(ray, face) " << Fmt(r.d.x * n.x + r.d.y * n.y + r.d.z * n.z, 3) << ", barycentrics " << Fmt(a.b1, 3) << " " << Fmt(a.b2, 3);
                }
                ++shown;
            }
        }
        const double share = double(hitMiss + otherTriangle + shadowDiffer + sphereDistance) / double(rays.size());
        return Report("tracers on the rays of " + path + ": OptiX = Embree (at most 1e-4 of the rays answered differently)", share <= 1e-4,
            std::to_string(rays.size()) + " rays (" + std::to_string(closest) + " closest-hit, " + std::to_string(shadows) + " shadow): " + std::to_string(otherTriangle) +
            " other primitive, " + std::to_string(hitMiss) + " hit against miss, " + std::to_string(shadowDiffer) + " shadow rays blocked on one side only" +
            (s.numSpheres > 0 ? ", " + std::to_string(sphereDistance) + " of " + std::to_string(sphereHits) + " sphere hits at another distance" : std::string()) + " (" +
            Fmt(100.0 * share, 3) + " %)" + examples.str());
    }

    // ---- materials: the kernel's against the engine's (double), the GPU's against the CPU's (float) ----

    Vector3 RandomUnit(kernel::Pcg32& rng) {
        const Real z = 1 - 2 * rng.nextDouble(), phi = 2 * constants::PI * rng.nextDouble(), r = std::sqrt(std::max(Real(0), 1 - z * z));
        return Vector3(r * std::cos(phi), r * std::sin(phi), z);
    }

    /// Directions at a surface: a geometric normal, a shading normal tilted from it, wo mostly above, wi a third of the time near the
    /// mirror direction of wo, a third near its refraction through an interface of relative index `ior` (where glossy lobes have their
    /// values), and otherwise anywhere; the ray comes from the outside 70 % of the time.
    struct DirCase { Vector3 wo, wi, ns, ng; Real u1, u2, u3; bool front; };
    std::vector<DirCase> MakeDirCases(int n, uint64_t seed, Real ior = Real(1.5)) {
        kernel::Pcg32 rng;
        rng.seed(seed, 3);
        std::vector<DirCase> cases;
        while (int(cases.size()) < n) {
            DirCase c;
            c.ng = RandomUnit(rng);
            c.ns = glm::normalize(c.ng + Real(0.25) * RandomUnit(rng));
            if (glm::dot(c.ns, c.ng) <= Real(0)) c.ns = c.ng;
            c.wo = RandomUnit(rng);
            if (glm::dot(c.wo, c.ng) < Real(0) && rng.nextDouble() < 0.9) c.wo = -c.wo;
            c.front = rng.nextDouble() < 0.7;
            const double pick = rng.nextDouble();
            Vector3 refracted;
            if (pick < 1.0 / 3) c.wi = glm::normalize(2 * glm::dot(c.wo, c.ns) * c.ns - c.wo + Real(0.3) * RandomUnit(rng));
            else if (pick < 2.0 / 3 && math::refractOutward(c.wo, c.ns, c.front ? Real(1) / ior : ior, refracted))
                c.wi = glm::normalize(refracted + Real(0.3) * RandomUnit(rng));
            else c.wi = RandomUnit(rng);
            c.u1 = rng.nextDouble();
            c.u2 = rng.nextDouble();
            c.u3 = rng.nextDouble();
            cases.push_back(c);
        }
        return cases;
    }

    double RelDiff(double a, double b) { return a == b ? 0.0 : std::abs(a - b) / std::max({ std::abs(a), std::abs(b), 1e-300 }); }

    /// The kernel's camera in double against the engine's Camera (built as a render builds it from the scene file, edited by `edit`): the
    /// ray of a random film point and lens sample, its origin, direction and weight; a real lens's blocked samples included (their weight
    /// is 0 and the ray is the pixel's central one).
    bool CameraEngineCheck(const std::string& path, const std::string& label, const std::function<void(SceneDescription&)>& edit = {}) {
        KernelSceneHost host;
        std::string why;
        bool ok = true;
        const auto withEngine = [&](const SceneSetup&, const Environment&, const Camera& camera) {
            kernel::CameraParams<Real> c{};
            std::vector<kernel::LensSurface<Real>> surfaces;
            std::vector<Real> pupil, bladeNormals;
            if (!camera.kernelParams(c, surfaces, pupil, bladeNormals, why)) { ok = Report("camera " + label, false, why); return; }
            const kernel::LensArrays<Real> arrays{ surfaces.data(), pupil.data(), bladeNormals.data() };
            kernel::Pcg32 rng;
            rng.seed(20261128, 5);
            const int n = 8192;
            double worst = 0;
            size_t identical = 0, values = 0, blocked = 0, blockedDiffer = 0;
            for (int i = 0; i < n; ++i) {
                // film points over the whole film and a margin outside it (saturated back onto the edge), lens samples over [0, 1)^2
                const Real s = Real(-0.02) + Real(1.04) * rng.nextDouble(), t = Real(-0.02) + Real(1.04) * rng.nextDouble();
                const Real lu = rng.nextDouble(), lv = rng.nextDouble();
                Real we = Real(1);
                const Ray re = camera.getRay(s, t, Point2(lu, lv), Real(0), kReferenceWavelengthNm, &we);
                Vector3 ok2, dk;
                Real wk = Real(0);
                kernel::lensCameraRay(c, arrays, s, t, lu, lv, ok2, dk, wk);
                const double vals[][2] = { { re.o.x, ok2.x }, { re.o.y, ok2.y }, { re.o.z, ok2.z }, { re.d.x, dk.x }, { re.d.y, dk.y }, { re.d.z, dk.z }, { we, wk } };
                bool same = true;
                for (const auto& v : vals) { worst = std::max(worst, RelDiff(v[0], v[1])); identical += v[0] == v[1]; same &= v[0] == v[1]; ++values; }
                if (we <= Real(0)) { ++blocked; blockedDiffer += !same; }
            }
            const std::string kind = c.type == kernel::kCamLens ? std::to_string(c.numSurfaces) + " surfaces, " + std::to_string(c.numBands) + " pupil bands" +
                (c.stopBlades ? ", " + std::to_string(c.stopBlades) + " blades" : std::string(", a round stop"))
                : c.type == kernel::kCamThinLens ? "a thin lens of radius " + Fmt(c.lensRadius, 3) + (c.blades ? ", " + std::to_string(c.blades) + " blades" : std::string(", round"))
                : std::string("a pinhole");
            ok = Report("kernel camera (double) = engine Camera, " + label + " (" + kind + ")", worst <= 1e-12,
                std::to_string(n) + " rays (" + std::to_string(blocked) + " blocked by the lens): origin, direction and weight " + Fmt(100.0 * double(identical) / double(values), 5) +
                " % bit-identical, largest relative difference " + Fmt(worst) + " (tolerance 1e-12)");
        };
        CliOptions o;
        o.spp = 1;
        if (!PackSceneFile(path, o, 400, 300, edit, host, why, withEngine)) return Report("camera " + label, false, why);
        return ok;
    }

#if !defined(RAYT_SPECTRAL_MODE)
    /// The kernel's conductor in double against the engine's RoughConductor on the same directions: eval, pdf and sample.
    bool ConductorEngineCheck(const std::string& label, const RoughConductor& engine, std::vector<std::pair<std::string, kernel::KMaterial>>& gpuMaterials) {
        kernel::ConductorParams<Real> c{};
        std::string why;
        if (!engine.kernelParams(c, why)) return Report("conductor " + label, false, why);
        const std::vector<DirCase> cases = MakeDirCases(4096, 77);
        SampledWavelengths lambda;
        double worst = 0;
        size_t identical = 0, values = 0, sampleMismatch = 0, nonzero = 0;
        for (const DirCase& dc : cases) {
            SurfaceInteraction rec{};
            rec.p = Point3(0.0);
            rec.n = dc.ns;
            rec.gn = dc.ng;
            rec.frontFace = true;
            const Spectrum fe = engine.eval(rec, dc.wo, dc.wi, lambda, TransportMode::Radiance);
            const Real pe = engine.pdf(rec, dc.wo, dc.wi, lambda);
            const Vector3 fk = kernel::conductorEval(c, dc.wo, dc.wi, dc.ns, dc.ng);
            const Real pk = kernel::conductorPdf(c, dc.wo, dc.wi, dc.ns, dc.ng);
            const double vals[][2] = { { fe.x, fk.x }, { fe.y, fk.y }, { fe.z, fk.z }, { pe, pk } };
            for (const auto& v : vals) { worst = std::max(worst, RelDiff(v[0], v[1])); identical += v[0] == v[1]; ++values; }
            nonzero += pe > 0;
            const std::optional<BSDFSample> se = engine.sample(rec, dc.wo, Point2(dc.u1, dc.u2), lambda, TransportMode::Radiance);
            kernel::BsdfSample<Real> sk{};
            const bool okK = kernel::conductorSample(c, dc.wo, dc.ns, dc.ng, dc.u1, dc.u2, sk);
            if (bool(se) != okK) { ++sampleMismatch; continue; }
            if (!okK) continue;
            const double svals[][2] = { { se->wi.x, sk.wi.x }, { se->wi.y, sk.wi.y }, { se->wi.z, sk.wi.z }, { se->f.x, sk.f.x }, { se->f.y, sk.f.y }, { se->f.z, sk.f.z },
                { se->pdf, sk.pdf } };
            for (const auto& v : svals) { worst = std::max(worst, RelDiff(v[0], v[1])); identical += v[0] == v[1]; ++values; }
        }
        kernel::KMaterial km{};
        km.type = kernel::kMatConductor;
        for (int i = 0; i < 3; ++i) { km.conductor.eta[i] = float(c.eta[i]); km.conductor.k[i] = float(c.k[i]); km.conductor.fms[i] = float(c.fms[i]); }
        km.conductor.alpha = float(c.alpha);
        km.conductor.multipleScattering = c.multipleScattering;
        km.conductor.eAvg = float(c.eAvg);
        for (int i = 0; i < kernel::kAlbedoMu; ++i) km.conductor.eRow[i] = float(c.eRow[i]);
        gpuMaterials.emplace_back("conductor " + label, km);
        return Report("kernel conductor (double) = engine RoughConductor, " + label + " (alpha " + Fmt(c.alpha, 3) + (c.multipleScattering ? ", Kulla-Conty lobe" : ", no multiple-scattering lobe") + ")",
            worst <= 1e-12 && sampleMismatch == 0,
            std::to_string(cases.size()) + " direction pairs (" + std::to_string(nonzero) + " with a nonzero pdf): eval, pdf and sample " + Fmt(100.0 * double(identical) / double(values), 5) +
            " % bit-identical, largest relative difference " + Fmt(worst) + " (tolerance 1e-12), " + std::to_string(sampleMismatch) + " samples accepted by one side only");
    }

    /// The kernel's dielectric in double against the engine's Dielectric on the same directions: eval, pdf and sample (entering and leaving).
    bool DielectricEngineCheck(const std::string& label, const Dielectric& engine, std::vector<std::pair<std::string, kernel::KMaterial>>& gpuMaterials) {
        kernel::DielectricParams<Real> c{};
        std::string why;
        if (!engine.kernelParams(c, why)) return Report("dielectric " + label, false, why);
        const std::vector<DirCase> cases = MakeDirCases(4096, 78, c.ior);
        SampledWavelengths lambda;
        double worst = 0;
        size_t identical = 0, values = 0, sampleMismatch = 0, nonzero = 0, sampled = 0;
        for (size_t i = 0; i < cases.size(); ++i) {
            const DirCase& dc = cases[i];
            SurfaceInteraction rec{};
            rec.p = Point3(0.0);
            rec.n = dc.ns;
            rec.gn = dc.ng;
            rec.frontFace = dc.front;
            const Spectrum fe = engine.eval(rec, dc.wo, dc.wi, lambda, TransportMode::Radiance);
            const Real pe = engine.pdf(rec, dc.wo, dc.wi, lambda);
            const Vector3 fk = kernel::dielectricEval(c, dc.front, dc.wo, dc.wi, dc.ns, dc.ng);
            const Real pk = kernel::dielectricPdf(c, dc.front, dc.wo, dc.wi, dc.ns);
            const double vals[][2] = { { fe.x, fk.x }, { fe.y, fk.y }, { fe.z, fk.z }, { pe, pk } };
            for (const auto& v : vals) { worst = std::max(worst, RelDiff(v[0], v[1])); identical += v[0] == v[1]; ++values; }
            nonzero += pe > 0;
            // A rough surface chooses between reflection and refraction with a number from the engine's own stream: seed the stream,
            // read the number the engine will draw, seed it again, and hand the same number to the kernel.
            sampling::SeedStream(907, i);
            const Real u3 = sampling::Random();
            sampling::SeedStream(907, i);
            const std::optional<BSDFSample> se = engine.sample(rec, dc.wo, Point2(dc.u1, dc.u2), lambda, TransportMode::Radiance);
            kernel::BsdfSample<Real> sk{};
            const bool okK = kernel::dielectricSample(c, dc.front, dc.wo, dc.ns, dc.ng, dc.u1, dc.u2, u3, sk);
            if (bool(se) != okK || (okK && se->isSpecular() != sk.delta)) { ++sampleMismatch; continue; }
            if (!okK) continue;
            ++sampled;
            const double svals[][2] = { { se->wi.x, sk.wi.x }, { se->wi.y, sk.wi.y }, { se->wi.z, sk.wi.z }, { se->f.x, sk.f.x }, { se->f.y, sk.f.y }, { se->f.z, sk.f.z },
                { se->pdf, sk.pdf } };
            for (const auto& v : svals) { worst = std::max(worst, RelDiff(v[0], v[1])); identical += v[0] == v[1]; ++values; }
        }
        kernel::KMaterial km{};
        km.type = kernel::kMatDielectric;
        km.dielectric.ior = float(c.ior);
        km.dielectric.alpha = float(c.alpha);
        km.dielectric.smooth = c.smooth;
        gpuMaterials.emplace_back("dielectric " + label, km);
        return Report("kernel dielectric (double) = engine Dielectric, " + label + " (relative index " + Fmt(c.ior, 4) + (c.smooth ? ", smooth" : ", alpha " + Fmt(c.alpha, 3)) + ")",
            worst <= 1e-12 && sampleMismatch == 0,
            std::to_string(cases.size()) + " direction pairs (" + std::to_string(nonzero) + " with a nonzero pdf, " + std::to_string(sampled) + " samples): eval, pdf and sample " +
            Fmt(100.0 * double(identical) / double(values), 5) + " % bit-identical, largest relative difference " + Fmt(worst) + " (tolerance 1e-12), " +
            std::to_string(sampleMismatch) + " samples accepted by one side only");
    }

    /// The kernel's phase function in double against the engine's PhaseFunction: eval on random cosines (a tenth near the forward and the
    /// backward direction) and sample on random directions and numbers.
    bool PhaseEngineCheck(const std::string& label, const PhaseFunction& engine) {
        kernel::MediumParams<Real> mp{};
        std::vector<Real> tableF, tableCdf;
        engine.kernelParams(mp, tableF, tableCdf);
        const kernel::PhaseTable<Real> tab{ tableF.data(), tableCdf.data() };
        kernel::Pcg32 rng;
        rng.seed(20261128, 9);
        const int n = 4096;
        double worst = 0;
        size_t identical = 0, values = 0;
        for (int i = 0; i < n; ++i) {
            Real c = Real(2) * rng.nextDouble() - Real(1);
            if (i % 10 == 0) c = Real(1) - Real(1e-4) * rng.nextDouble();
            if (i % 10 == 1) c = Real(-1) + Real(1e-4) * rng.nextDouble();
            const Real z = Real(2) * rng.nextDouble() - Real(1), phi = Real(6.283185307179586) * rng.nextDouble(), r = std::sqrt(std::max(Real(0), Real(1) - z * z));
            const Vector3 dir(r * std::cos(phi), r * std::sin(phi), z);
            const Real ux = rng.nextDouble(), uy = rng.nextDouble();
            const Vector3 de = engine.sample(dir, Point2(ux, uy));
            const Vector3 dk = kernel::phaseSample(mp, tab, dir, ux, uy);
            const double vals[][2] = { { engine.eval(c), kernel::phaseEval(mp, tab, c) }, { de.x, dk.x }, { de.y, dk.y }, { de.z, dk.z } };
            for (const auto& v : vals) { worst = std::max(worst, RelDiff(v[0], v[1])); identical += v[0] == v[1]; ++values; }
        }
        return Report("kernel phase function (double) = engine PhaseFunction, " + label, worst <= 1e-12,
            std::to_string(n) + " cosines and samples: " + Fmt(100.0 * double(identical) / double(values), 5) + " % bit-identical, largest relative difference " +
            Fmt(worst) + " (tolerance 1e-12)");
    }

    /// The kernel's density field in double against the engine's NoiseDensity (as ChromaticMedium reads it, clamped to [0, 1]).
    bool NoiseEngineCheck(const std::string& label, const NoiseDensity& engine, Real extent) {
        const ChromaticMedium med(Vector3(100.0), Vector3(0.9), PhaseFunction(), std::make_shared<NoiseDensity>(engine));
        kernel::MediumParams<Real> mp{};
        std::vector<Real> tableF, tableCdf;
        std::string why;
        if (!med.kernelParams(mp, tableF, tableCdf, why)) return Report("density " + label, false, why);
        kernel::Pcg32 rng;
        rng.seed(20261128, 11);
        const int n = 8192;
        size_t identical = 0, inside = 0;
        double worst = 0;
        for (int i = 0; i < n; ++i) {
            const Point3 p(extent * (Real(2) * rng.nextDouble() - Real(1)), extent * (Real(2) * rng.nextDouble() - Real(1)), extent * (Real(2) * rng.nextDouble() - Real(1)));
            const Real qe = med.density(p), qk = kernel::mediumDensity(mp, p);
            worst = std::max(worst, std::abs(qe - qk));
            identical += qe == qk;
            inside += qe > Real(0) && qe < Real(1);
        }
        return Report("kernel density field (double) = engine NoiseDensity, " + label, identical == size_t(n),
            std::to_string(n) + " points (" + std::to_string(inside) + " strictly between 0 and 1): " + Fmt(100.0 * double(identical) / double(n), 5) +
            " % bit-identical, largest difference " + Fmt(worst));
    }

    /// The engine's random number stream as Kernel/Layered.hpp's source: the kernel draws what LayeredCoat would draw next.
    struct EngineRand {
        Real next() { return sampling::Random(); }
        kernel::Vec2<Real> next2() {
            const Point2 p = sampling::Random2D();
            return kernel::Vec2<Real>(p.x, p.y);
        }
    };

    /// The kernel's layered material in double against the engine's LayeredCoat, for each layered material of a scene file (built
    /// by the engine as a render builds it: the base metal in the coat's medium): eval, pdf and sample on the same directions. The
    /// walks draw their random numbers inside LayeredCoat; the engine's stream is seeded before each call, and seeded again before
    /// the kernel's call, which draws from it the same numbers in the same order.
    bool LayeredEngineChecks(const std::string& path, std::vector<std::pair<std::string, kernel::KMaterial>>& gpuMaterials) {
        bool ok = true;
        KernelSceneHost host;
        std::string why;
        std::vector<std::pair<std::string, const LayeredCoat*>> coats;
        std::vector<std::string> results;
        const auto withEngine = [&](const SceneSetup& setup, const Environment&, const Camera&) {
            for (const auto& [name, m] : setup.namedMaterials) {
                const auto* lc = dynamic_cast<const LayeredCoat*>(m);
                if (!lc) continue;
                kernel::LayeredParams<Real> lp{};
                Vector3 albedo(0.0);
                kernel::ConductorParams<Real> c{};
                std::string limit;
                if (!lc->kernelParams(lp, albedo, c, limit)) { ok &= Report("layered " + name, false, limit); continue; }
                const kernel::LayerBase<Real> base{ albedo, lp.conductorBase ? &c : nullptr };
                const std::vector<DirCase> cases = MakeDirCases(2048, 79);
                SampledWavelengths lambda;
                double worst = 0;
                size_t identical = 0, values = 0, sampleMismatch = 0, nonzero = 0, sampled = 0;
                for (size_t i = 0; i < cases.size(); ++i) {
                    const DirCase& dc = cases[i];
                    SurfaceInteraction rec{};   // no tangent: the layers' frame is built from the normal alone, as the kernel's
                    rec.p = Point3(0.0);
                    rec.n = dc.ns;
                    rec.gn = dc.ng;
                    rec.frontFace = true;
                    EngineRand rnd;
                    sampling::SeedStream(911, i);
                    const Spectrum fe = lc->eval(rec, dc.wo, dc.wi, lambda, TransportMode::Radiance);
                    sampling::SeedStream(911, i);
                    const Vector3 fk = kernel::layeredEval(lp, base, dc.wo, dc.wi, dc.ns, dc.ng, rnd);
                    const Real pe = lc->pdf(rec, dc.wo, dc.wi, lambda);
                    const Real pk = kernel::layeredPdf(lp, base, dc.wo, dc.wi, dc.ns, dc.ng);
                    const double vals[][2] = { { fe.x, fk.x }, { fe.y, fk.y }, { fe.z, fk.z }, { pe, pk } };
                    for (const auto& v : vals) { worst = std::max(worst, RelDiff(v[0], v[1])); identical += v[0] == v[1]; ++values; }
                    nonzero += pe > 0;
                    sampling::SeedStream(913, i);
                    const std::optional<BSDFSample> se = lc->sample(rec, dc.wo, Point2(dc.u1, dc.u2), lambda, TransportMode::Radiance);
                    sampling::SeedStream(913, i);
                    kernel::BsdfSample<Real> sk{};
                    const bool okK = kernel::layeredSample(lp, base, dc.wo, dc.ns, dc.ng, kernel::Vec2<Real>(dc.u1, dc.u2), rnd, sk);
                    if (bool(se) != okK) { ++sampleMismatch; continue; }
                    if (!okK) continue;
                    ++sampled;
                    const double svals[][2] = { { se->wi.x, sk.wi.x }, { se->wi.y, sk.wi.y }, { se->wi.z, sk.wi.z }, { se->f.x, sk.f.x }, { se->f.y, sk.f.y },
                        { se->f.z, sk.f.z }, { se->pdf, sk.pdf } };
                    for (const auto& v : svals) { worst = std::max(worst, RelDiff(v[0], v[1])); identical += v[0] == v[1]; ++values; }
                }
                ok &= Report("kernel layered (double) = engine LayeredCoat, " + path + " \"" + name + "\" (" + (lp.conductorBase ? "metal" : "diffuse") + " base, coat alpha " +
                    Fmt(lp.alpha, 3) + (lp.absorbing ? ", tinted" : "") + ", " + std::to_string(lp.samples) + (lp.samples == 1 ? " walk)" : " walks)"),
                    worst <= 1e-12 && sampleMismatch == 0,
                    std::to_string(cases.size()) + " direction pairs (" + std::to_string(nonzero) + " with a nonzero pdf, " + std::to_string(sampled) + " samples): eval, pdf and sample " +
                    Fmt(100.0 * double(identical) / double(values), 5) + " % bit-identical, largest relative difference " + Fmt(worst) + " (tolerance 1e-12), " +
                    std::to_string(sampleMismatch) + " samples accepted by one side only");
                coats.emplace_back(name, lc);
            }
        };
        CliOptions o;
        o.spp = 1;
        if (!PackSceneFile(path, o, 64, 48, {}, host, why, withEngine)) return Report("layered materials of " + path, false, why);
        if (coats.empty()) return Report("layered materials of " + path, false, "none found");
        // the packed materials are in the order of the scene file's names (a std::map), as setup.namedMaterials, which the hook went through
        size_t k = 0;
        for (const kernel::KMaterial& km : host.materials)
            if (km.type == kernel::kMatLayered && k < coats.size()) gpuMaterials.emplace_back("layered " + coats[k++].first, km);
        return ok;
    }
#endif

    /// The kernel's environment (double where the engine uses double) against the engine's EnvMap: eval, sample and pdf.
    bool EnvEngineCheck(const std::string& label, const EnvLight* light) {
        KernelSceneHost host;
        if (!PackEnvironment(light, host)) return Report("environment " + label, false, "no environment");
        const kernel::KScene s = host.view();
        const EnvMap& map = *light->envMap();
        kernel::Pcg32 rng;
        rng.seed(4242, 5);
        double worst = 0;
        size_t identical = 0, values = 0, mismatch = 0;
        const int n = 4096;
        for (int i = 0; i < n; ++i) {
            // a direction (a tenth of them within a degree of a pole, where the mapping is singular) and two sample numbers
            Vector3 dir = RandomUnit(rng);
            if (i % 10 == 0) dir = glm::normalize(Vector3(Real(0.01) * dir.x, dir.y < 0 ? Real(-1) : Real(1), Real(0.01) * dir.z));
            const Point2f u(rng.nextFloat(), rng.nextFloat());
            const Vector3 ee = map.eval(dir);
            const Vector3 ek = kernel::envEval(s, dir);
            const Real pe = map.pdf(dir), pk = kernel::envPdf(s, dir);
            Vector3 wiE, wiK, leK;
            Real pdfE = 0, pdfK = 0;
            const Vector3 leE = map.sample(u, wiE, pdfE);
            kernel::envSample(s, u.x, u.y, wiK, leK, pdfK);
            if ((pdfE > 0) != (pdfK > 0)) { ++mismatch; continue; }
            std::vector<std::array<double, 2>> v = { { ee.x, ek.x }, { ee.y, ek.y }, { ee.z, ek.z }, { pe, pk }, { pdfE, pdfK } };
            if (pdfE > 0) for (int c = 0; c < 3; ++c) { v.push_back({ wiE[c], wiK[c] }); v.push_back({ leE[c], leK[c] }); }
            for (const auto& p : v) { worst = std::max(worst, RelDiff(p[0], p[1])); identical += p[0] == p[1]; ++values; }
        }
        return Report("kernel environment = engine EnvMap, " + label + " (" + std::to_string(s.envWidth) + " x " + std::to_string(s.envHeight) + ")",
            worst <= 1e-12 && mismatch == 0,
            std::to_string(n) + " directions and samples: eval, pdf and sample " + Fmt(100.0 * double(identical) / double(values), 5) + " % bit-identical, largest relative difference " +
            Fmt(worst) + " (tolerance 1e-12), " + std::to_string(mismatch) + " samples valid on one side only");
    }

    /// The kernel's sphere light in double against the engine's Sphere: sampleFrom (with AreaLight::sampleLi's direction) and
    /// pdfSurface, from points outside the sphere (the cone it subtends) and inside it (a uniform point of its surface).
    bool SphereLightEngineCheck() {
        struct Ball { Point3 c; Real r; };
        const Ball balls[] = { { Point3(0.3, 1.2, -0.4), Real(0.25) }, { Point3(5, 20, 3), Real(0.01) }, { Point3(-2, 0.5, 1), Real(3) } };
        kernel::Pcg32 rng;
        rng.seed(5150, 7);
        double worst = 0;
        size_t identical = 0, values = 0, mismatch = 0, inside = 0;
        const int n = 4096;
        for (int i = 0; i < n; ++i) {
            const Ball& ball = balls[i % 3];
            const Sphere engine(ball.c, ball.r, nullptr);
            // a reference point: a fifth inside the sphere, the others at 1.001 to 50 radii from the centre (most of them near it)
            const bool in = i % 5 == 0;
            const Real dist = in ? ball.r * Real(0.95) * rng.nextDouble() : ball.r * (Real(1.001) + Real(49) * rng.nextDouble() * rng.nextDouble());
            const Point3 ref = ball.c + RandomUnit(rng) * dist;
            inside += in;
            const Real u1 = rng.nextDouble(), u2 = rng.nextDouble();
            const auto es = engine.sampleFrom(ref, Point2(u1, u2));
            kernel::LightSample<Real> ks{};
            const bool kok = kernel::sampleSphereLight<Real>(ball.c, ball.r, false, ref, u1, u2, ks);
            if (bool(es) != kok) { ++mismatch; continue; }
            std::vector<std::array<double, 2>> v;
            if (kok) {
                const Vector3 d = es->si.p - ref;
                const Vector3 wiE = d * math::safe_recip(math::safe_sqrt(glm::dot(d, d)));   // AreaLight::sampleLi
                for (int c = 0; c < 3; ++c) { v.push_back({ es->si.p[c], ks.p[c] }); v.push_back({ wiE[c], ks.wi[c] }); }
                v.push_back({ es->pdf, ks.pdf });
                if (es->si.frontFace != ks.emits) ++mismatch;
            }
            // the density toward the sampled point and toward a random direction (inside: where the engine's ray meets the sphere)
            SurfaceInteraction at;
            at.p = ref;
            at.gn = at.n = RandomUnit(rng);
            for (int k = 0; k < 2; ++k) {
                const Vector3 wi = k == 0 && kok ? ks.wi : RandomUnit(rng);
                Point3 hitP = ref;
                if (in) {
                    SurfaceInteraction isect;
                    if (!engine.hit(Ray(at.OffsetRayOrigin(wi), wi), isect)) continue;
                    hitP = isect.p;
                }
                v.push_back({ engine.pdfSurface(at, wi), kernel::sphereLightPdf<Real>(ball.c, ball.r, ref, wi, hitP) });
            }
            for (const auto& q : v) { worst = std::max(worst, RelDiff(q[0], q[1])); identical += q[0] == q[1]; ++values; }
        }
        return Report("kernel sphere light = engine Sphere::sampleFrom / pdfSurface", worst <= 1e-12 && mismatch == 0,
            std::to_string(n) + " samples (" + std::to_string(inside) + " from inside the sphere) and densities: " + Fmt(100.0 * double(identical) / double(values), 5) +
            " % bit-identical, largest relative difference " + Fmt(worst) + " (tolerance 1e-12), " + std::to_string(mismatch) + " samples valid or facing on one side only");
    }

    /// A material's eval, pdf and sample in float on the CPU and on the GPU, on the same directions.
    bool MaterialGpuCheck(const gpu::Module& mod, bool fastMath, const std::string& label, const kernel::KMaterial& m) {
        const std::vector<DirCase> dirs = MakeDirCases(4096, 91);
        std::vector<kernel::MaterialCase> cases;
        auto f3 = [](const Vector3& v) { return kernel::KFloat3{ float(v.x), float(v.y), float(v.z) }; };
        for (const DirCase& d : dirs) {
            kernel::MaterialCase c{};
            // the float directions are renormalized on the CPU, so both devices start from the same float inputs
            const kernel::Vec3<float> ns = kernel::knormalize(kernel::Vec3<float>(float(d.ns.x), float(d.ns.y), float(d.ns.z)));
            c.wo = f3(d.wo); c.wi = f3(d.wi); c.ng = f3(d.ng); c.ns = { ns.x, ns.y, ns.z };
            c.u1 = float(d.u1); c.u2 = float(d.u2); c.u3 = float(d.u3);
            c.front = d.front ? 1u : 0u;
            c.stream = uint32_t(cases.size());
            if (c.u3 >= 1.0f) c.u3 = 0.99999994f;
            if (c.u1 >= 1.0f) c.u1 = 0.99999994f;
            if (c.u2 >= 1.0f) c.u2 = 0.99999994f;
            cases.push_back(c);
        }
        std::vector<kernel::MaterialOut> cpu(cases.size()), gpu(cases.size());
        for (size_t i = 0; i < cases.size(); ++i) kernel::evalMaterialCase(m, cases[i], cpu[i]);
        char err[512] = {};
        if (mod.materialTest(&m, cases.data(), gpu.data(), int(cases.size()), err, int(sizeof(err))) != 0) return Report(label + " on the GPU", false, err);
        size_t identical = 0;
        double worst = 0;
        for (size_t i = 0; i < cases.size(); ++i) {
            const bool same = memcmp(&cpu[i], &gpu[i], sizeof(kernel::MaterialOut)) == 0;
            identical += same;
            if (!same) {
                const float* a = &cpu[i].f[0];
                const float* b = &gpu[i].f[0];
                for (int k = 0; k < 4; ++k) worst = std::max(worst, RelDiff(a[k], b[k]));
            }
        }
        return Report("material " + label + " in float: GPU " + (fastMath ? "(fast math) close to the CPU" : "= CPU bit for bit") + " (eval, pdf, sample)",
            fastMath ? worst <= 1e-3 : identical == cases.size(),
            Fmt(100.0 * double(identical) / double(cases.size()), 5) + " % of " + std::to_string(cases.size()) + " cases bit-identical" +
            (identical == cases.size() ? "" : ", largest relative difference of eval / pdf " + Fmt(worst)));
    }

    /// Mean of all pixel values and its standard error (pixels are independent: each has its own paths).
    void MeanAndError(const std::vector<float>& img, double& mean, double& se) {
        double s = 0, s2 = 0;
        for (float v : img) { s += v; s2 += double(v) * v; }
        const double n = double(img.size());
        mean = s / n;
        se = std::sqrt(std::max(s2 / n - mean * mean, 0.0) / n);
    }

    bool FurnaceCheck(bool onGpu) {
        KernelSceneHost host;
        std::string why;
        if (!KernelTestScene("scenes/gpu_box.json", 96, 72, 16, 1000, 21, true, host, why)) return Report("white furnace", false, why);
        std::vector<float> img;
        const auto t0 = std::chrono::steady_clock::now();
        if (!(onGpu ? RenderKernelGpu(host.view(), img, why, false) : RenderKernelCpu(host.view(), img, why, false)))
            return Report(std::string("white furnace on the ") + (onGpu ? "GPU" : "CPU"), false, why);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        double mean, se;
        MeanAndError(img, mean, se);
        return Report(std::string("kernel renderer, white furnace on the ") + (onGpu ? "GPU" : "CPU") + " (albedo 1 in a uniform environment 1: every pixel 1)",
            std::abs(mean - 1.0) <= 4 * se + 1e-6,
            "mean " + Fmt(mean, 6) + " +- " + Fmt(se, 2) + " (must be 1 within 4 standard errors); 96 x 72, 16 spp, depth 1000, " + Fmt(seconds, 2) + " s");
    }

} // namespace

int RunGpuTests(const CliOptions& /*opt*/) {
    bool ok = true;
    const int n = 4096;
    const std::vector<kernel::MathCase> cases = MakeCases(n);
    // The CPU computes the float values (as the GPU does) and the double reference (which only the CPU computes).
    std::vector<kernel::MathOut> cpu(cases.size());
    std::vector<std::array<double, kernel::kMathValueCount>> ref(cases.size());
    for (size_t i = 0; i < cases.size(); ++i) {
        kernel::evalMathCase(cases[i], cpu[i]);
        kernel::evalMathCaseDouble(cases[i], ref[i].data());
    }

    // 1. The kernel's functions on the CPU in double are the engine's functions (same formulas, same order).
    {
        double maxDiff = 0;
        size_t differing = 0, nans = 0;
        for (size_t i = 0; i < cases.size(); ++i) {
            double e[kernel::kMathValueCount];
            EngineValues(cases[i], e);
            bool same = true;
            for (int k = kernel::kHemiX; k <= kernel::kFresnel; ++k) {   // atan2 / acos have no engine counterpart (the engine calls the library)
                const double d = std::abs(ref[i][k] - e[k]);
                if (std::isnan(d)) ++nans; else maxDiff = std::max(maxDiff, d);
                if (!(ref[i][k] == e[k])) same = false;
            }
            if (!same) ++differing;
        }
        ok &= Report("kernel (CPU, double) = engine functions: cosine hemisphere, uniform sphere, orthonormal basis, dielectric Fresnel", differing == 0,
            std::to_string(cases.size()) + " cases, " + std::to_string(differing) + " not bit-identical (max |diff| " + Fmt(maxDiff) + ", " +
            std::to_string(nans) + " NaN)");
    }

    // 2. The kernel renderer on the CPU (Embree, float): the white furnace.
    ok &= FurnaceCheck(false);

    // 3. The kernel's materials in double against the engine's (the GPU part below then compares them in float, CPU against GPU).
    std::vector<std::pair<std::string, kernel::KMaterial>> gpuMaterials;
    {
        kernel::KMaterial diffuse{};
        diffuse.type = kernel::kMatDiffuse;
        diffuse.albedo[0] = 0.8f; diffuse.albedo[1] = 0.5f; diffuse.albedo[2] = 0.3f;
        gpuMaterials.emplace_back("diffuse", diffuse);
    }
#if defined(RAYT_SPECTRAL_MODE)
    std::cout << "[GPU]    SKIP the conductor and dielectric checks: the kernel models them in RGB (use the RGB build)" << std::endl;
#else
    {
        const Spectrum goldEta(0.143, 0.374, 1.442), goldK(3.983, 2.385, 1.603);   // gold, roughly (RGB)
        RoughConductor rough(goldEta, goldK, Real(0.4)), smooth(goldEta, goldK, Real(0.1)), plain(goldEta, goldK, Real(0.4));
        plain.setMultipleScattering(false);
        ok &= ConductorEngineCheck("rough gold", rough, gpuMaterials);
        ok &= ConductorEngineCheck("smooth gold", smooth, gpuMaterials);
        ok &= ConductorEngineCheck("rough gold without the lobe", plain, gpuMaterials);
        Dielectric glass(Real(1.5), Real(0)), frosted(Real(1.5), Real(0.3)), underWater(Real(1.5), Real(0.2));
        underWater.setExteriorIor(Real(1.333));
        ok &= DielectricEngineCheck("smooth glass", glass, gpuMaterials);
        ok &= DielectricEngineCheck("rough glass", frosted, gpuMaterials);
        ok &= DielectricEngineCheck("rough glass in water", underWater, gpuMaterials);
    }
    {
        // the uniform environment is a white map in the engine too; the HDR map has a sun and a dark lower half
        const Environment uniform = MakeEnvironment(CliOptions{}, EnvKind::Uniform);
        ok &= EnvEngineCheck("the uniform environment", uniform.light.get());
        CliOptions withMap;
        withMap.env = "assets/env/bryanston_park_sunrise_2k.hdr";
        std::ostringstream quiet;
        std::streambuf* old = std::cout.rdbuf(quiet.rdbuf());
        const Environment hdr = MakeEnvironment(withMap, EnvKind::Hdr);
        std::cout.rdbuf(old);
        ok &= EnvEngineCheck("bryanston_park_sunrise_2k.hdr", hdr.light.get());
    }
#endif
    ok &= SphereLightEngineCheck();
#if !defined(RAYT_SPECTRAL_MODE)
    // 3a. The camera (the kernel renders in RGB, where the glass of a real lens is not dispersive): a thin lens (polygonal and round) and
    //     real lenses (the double Gauss with blades, a singlet, a fisheye), kernel in double = engine.
    ok &= CameraEngineCheck("scenes/gpu_dof.json", "gpu_dof.json (six blades)");
    ok &= CameraEngineCheck("scenes/gpu_dof.json", "gpu_dof.json with a round aperture", [](SceneDescription& d) { d.camera.lens.blades = 0; });
    ok &= CameraEngineCheck("scenes/gpu_lens.json", "gpu_lens.json (double Gauss 50 mm f/2.8)");
    ok &= CameraEngineCheck("scenes/gpu_lens.json", "a BK7 singlet 100 mm f/5.6", [](SceneDescription& d) {
        d.camera.lens.lensFile = std::string("assets/Lenses/singlet_bk7_100mm.txt");
        d.camera.lens.fStop = Real(5.6);
        d.camera.lens.blades = 0;
    });
    ok &= CameraEngineCheck("scenes/gpu_lens.json", "a fisheye 10 mm wide open", [](SceneDescription& d) {
        d.camera.lens.lensFile = std::string("assets/Lenses/fisheye_10mm.txt");
        d.camera.lens.fStop.reset();
        d.camera.lens.blades = 0;
    });
    // 3a'. The media: phase functions and the density field, kernel in double = engine.
    ok &= PhaseEngineCheck("isotropic", PhaseFunction(Real(0)));
    ok &= PhaseEngineCheck("Henyey-Greenstein g 0.6", PhaseFunction(Real(0.6)));
    ok &= PhaseEngineCheck("two lobes g 0.97 / -0.5, weight 0.884", PhaseFunction(Real(0.97), Real(-0.5), Real(0.884)));
    {
        std::string err;
        const auto mie = TabulatedPhase::fromCsv("assets/Cloud/mie_reff10.csv", &err);
        if (mie) ok &= PhaseEngineCheck("the Mie table of assets/Cloud/mie_reff10.csv", PhaseFunction(mie));
        else ok &= Report("Mie table", false, err);
    }
    {
        NoiseDensity cloud(Real(22), 5, Real(0.42), Real(0.62), 7u);
        cloud.setSphereMask(Point3(0.0), Real(0.1), Real(0.45));
        ok &= NoiseEngineCheck("the cloud's (22 / m, 5 octaves, a spherical mask)", cloud, Real(0.12));
        ok &= NoiseEngineCheck("the marble's (55 / m, 4 octaves)", NoiseDensity(Real(55), 4, Real(0.30), Real(0.62), 5u), Real(2.0));
    }
    // 3b. The layered materials of the coat scene, kernel in double = engine.
    ok &= LayeredEngineChecks("scenes/gpu_coat.json", gpuMaterials);
    // 3c. The denoiser's features on the same rays as the engine (the kernel renders in RGB only).
    for (const char* path : { "scenes/gpu_box.json", "scenes/gpu_metals.json", "scenes/gpu_glass.json", "scenes/gpu_env.json", "scenes/gpu_spheres.json", "scenes/gpu_coat.json",
                              "scenes/gpu_media.json" })
        ok &= AovEngineCheck(path);
#endif

    // 4. The GPU module.
    std::string why;
    const gpu::Module* mod = gpu::LoadModule(why);
    if (!mod) {
        std::cout << "[GPU]    SKIP the GPU part: " << why << std::endl;
        return ok ? 0 : 2;
    }
    RaytGpuInfo info{};
    if (mod->query(&info) != 0) {
        std::cout << "[GPU]    SKIP the GPU part: " << info.error << " (module " << mod->path << ")" << std::endl;
        return ok ? 0 : 2;
    }
    std::cout << "[GPU]    device 0: " << info.name << ", compute capability " << info.computeMajor << "." << info.computeMinor << ", "
              << info.multiprocessors << " multiprocessors, " << Fmt(double(info.totalMemory) / (1024.0 * 1024.0 * 1024.0), 3) << " GB; driver supports CUDA "
              << info.driverVersion / 1000 << "." << (info.driverVersion % 1000) / 10 << "; module built with CUDA " << info.runtimeVersion / 1000 << "."
              << (info.runtimeVersion % 1000) / 10 << ", float only, "
              << (info.fastMath ? "fast math (--use_fast_math)" : "IEEE division / square root, no implicit FMA") << std::endl;

    std::vector<kernel::MathOut> gpu(cases.size());
    char err[512] = {};
    if (mod->mathTest(cases.data(), gpu.data(), int(cases.size()), err, int(sizeof(err))) != 0)
        return Report("evalMathCase on the GPU", false, err) ? 0 : 2;

    // 5a. Random numbers: integer arithmetic, so the same bits.
    {
        size_t differing = 0;
        for (size_t i = 0; i < cases.size(); ++i) {
            bool same = true;
            for (int k = 0; k < 4; ++k) same &= cpu[i].rng[k] == gpu[i].rng[k];
            for (int k = kernel::kRngFloat0; k <= kernel::kRngFloat1; ++k) same &= cpu[i].f[k] == gpu[i].f[k];
            if (!same) ++differing;
        }
        ok &= Report("random numbers (PCG32 streams keyed by seed, pixel, sample): GPU = CPU bit for bit", differing == 0,
            std::to_string(cases.size()) + " streams, " + std::to_string(differing) + " differ");
    }

    // 5b. Float values. The default module computes the same operations in the same order as the CPU (IEEE division and square root, no
    //    implicit FMA, our own sin / cos), so its floats must be the CPU's bit for bit. A fast-math module (approximate division and
    //    square root, implicit FMAs) cannot be; it must be as close to the double reference as the CPU is. (A fixed tolerance on
    //    |GPU - CPU| would be the wrong test: the Fresnel term just past the critical angle, cos theta_t = sqrt(small), amplifies a
    //    rounding difference in its input to 1.4e-5.)
    {
        struct Group { const char* name; int first, last; };
        const Group groups[] = { { "cosine hemisphere", kernel::kHemiX, kernel::kHemiZ }, { "uniform sphere", kernel::kSphereX, kernel::kSphereZ },
            { "orthonormal basis", kernel::kBasisSX, kernel::kBasisTZ }, { "dielectric Fresnel", kernel::kFresnel, kernel::kFresnel },
            { "atan2, acos, exp and log", kernel::kAtan2, kernel::kLog }, { "the noise density of a medium", kernel::kNoise, kernel::kNoise },
            { "a product that is a subnormal float", kernel::kSubnormal, kernel::kSubnormal } };
        for (const Group& g : groups) {
            double errCpu = 0, errGpu = 0, maxDiff = 0;
            size_t identical = 0, values = 0, nans = 0, worst = 0;
            for (size_t i = 0; i < cases.size(); ++i)
                for (int k = g.first; k <= g.last; ++k) {
                    const double c = cpu[i].f[k], d = gpu[i].f[k];
                    if (std::isnan(c) || std::isnan(d)) { ++nans; continue; }
                    errCpu = std::max(errCpu, std::abs(c - ref[i][k]));
                    errGpu = std::max(errGpu, std::abs(d - ref[i][k]));
                    if (std::abs(c - d) > maxDiff) { maxDiff = std::abs(c - d); worst = i; }
                    identical += c == d;
                    ++values;
                }
            std::ostringstream detail;
            detail << Fmt(100.0 * double(identical) / double(values), 4) << " % of the values bit-identical, max |GPU - CPU| " << Fmt(maxDiff) << " (case " << worst
                   << "); max |float - double|: GPU " << Fmt(errGpu) << ", CPU " << Fmt(errCpu);
            if (nans) detail << ", " << nans << " NaN";
            const bool pass = nans == 0 && (info.fastMath ? errGpu <= 2 * errCpu + 1e-7 : identical == values);
            ok &= Report(std::string("float, ") + g.name + (info.fastMath ? ": GPU (fast math) as accurate as the CPU (<= 2 x CPU + 1e-7)" : ": GPU = CPU bit for bit"),
                pass, detail.str());
        }
    }

    // 6. The materials in float on the GPU against the CPU.
    for (const auto& [label, m] : gpuMaterials) ok &= MaterialGpuCheck(*mod, info.fastMath != 0, label, m);

    // 7. The kernel renderer on the GPU: the white furnace, and the test scenes against the CPU kernel renderer with the same seed. With
    //    the same float operations on both devices and the hit point recomputed by the kernel from the triangle the tracer reports,
    //    the two pictures must be the same bits (fast math: the same mean within the noise).
    ok &= FurnaceCheck(true);
    ok &= PairedCheck("scenes/gpu_box.json", info.fastMath != 0);
#if !defined(RAYT_SPECTRAL_MODE)
    ok &= PairedCheck("scenes/gpu_metals.json", info.fastMath != 0);
    ok &= PairedCheck("scenes/gpu_glass.json", info.fastMath != 0);
    ok &= PairedCheck("scenes/gpu_env.json", info.fastMath != 0);
    ok &= PairedCheck("scenes/gpu_spheres.json", info.fastMath != 0);
    ok &= PairedCheck("scenes/gpu_coat.json", info.fastMath != 0);
    ok &= PairedCheck("scenes/gpu_dof.json", info.fastMath != 0);
    ok &= PairedCheck("scenes/gpu_lens.json", info.fastMath != 0);
    ok &= PairedCheck("scenes/gpu_media.json", info.fastMath != 0);
    ok &= PairedCheck("scenes/gpu_cloud.json", info.fastMath != 0);
#endif

    // 8. The tracers on the rays the CPU's paths cast: the only place left where the two devices' paths could part (coincident surfaces,
    //    hits exactly on an edge).
    ok &= TracerCheck(*mod, "scenes/gpu_box.json");
#if !defined(RAYT_SPECTRAL_MODE)
    ok &= TracerCheck(*mod, "scenes/gpu_metals.json");
    ok &= TracerCheck(*mod, "scenes/gpu_glass.json");
    ok &= TracerCheck(*mod, "scenes/gpu_env.json");
    ok &= TracerCheck(*mod, "scenes/gpu_spheres.json");
    ok &= TracerCheck(*mod, "scenes/gpu_coat.json");
    ok &= TracerCheck(*mod, "scenes/gpu_lens.json");
    ok &= TracerCheck(*mod, "scenes/gpu_media.json");
    ok &= TracerCheck(*mod, "scenes/gpu_cloud.json");
#endif
    return ok ? 0 : 2;
}

} // namespace rayt::cli
