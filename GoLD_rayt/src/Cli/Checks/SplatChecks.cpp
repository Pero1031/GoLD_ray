// Gaussian splat checks (Geometry/GaussianSplats.hpp): the ray-Gaussian response, compositing, the BVH, the 3DGS PLY, the transform
// and the integrator hook.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Core/Sampling.hpp"
#include "Geometry/GaussianSplats.hpp"
#include "Geometry/TriangleMesh.hpp"
#include "Lights/AreaLight.hpp"
#include "Materials/DiffuseLight.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Integrator.hpp"
#include "Scene/Primitive.hpp"

namespace rayt::cli {

namespace {

constexpr Real kC0 = 0.28209479177387814;

// A particle of degree 0 whose color (before linearization) is `rgb`: f_dc = (rgb - 0.5) / C0.
SplatSource MakeSplat(const Point3& mean, const Vector3& scale, const glm::dvec4& rot, Real opacity, const Vector3& rgb) {
    SplatSource s;
    s.mean = mean;
    s.scale = scale;
    s.rot = glm::normalize(rot);
    s.opacity = opacity;
    s.sh = { float((rgb.x - 0.5) / kC0), float((rgb.y - 0.5) / kC0), float((rgb.z - 0.5) / kC0) };
    return s;
}

SplatParams PlainParams() {
    SplatParams p;
    p.linearize = false;
    return p;
}

// The Gaussian's value along a ray, maximized numerically (golden-section search on the convex quadratic form): independent of
// the closed form of the engine.
Real NumericAlpha(const SplatSource& s, const Ray& r, Real& tBest) {
    const glm::dmat3 R = GaussianSplats::quatToMatrix(s.rot);
    const glm::dmat3 S = glm::dmat3(s.scale.x, 0, 0, 0, s.scale.y, 0, 0, 0, s.scale.z);
    const glm::dmat3 sigma = R * S * S * glm::transpose(R);
    const glm::dmat3 inv = glm::inverse(sigma);
    auto q = [&](Real t) { const Vector3 x = Vector3(r.o) + r.d * t - Vector3(s.mean); return glm::dot(x, inv * x); };
    Real a = -50, b = 50;
    const Real phi = (std::sqrt(5.0) - 1) / 2;
    for (int i = 0; i < 200; ++i) {
        const Real c = b - phi * (b - a), d = a + phi * (b - a);
        if (q(c) < q(d)) b = d; else a = c;
    }
    tBest = (a + b) / 2;
    return s.opacity * std::exp(-0.5 * q(tBest));
}

} // namespace

bool RunSplatChecks(const CliOptions& opt) {
    bool allPassed = true;
    auto check = [&](const std::string& what, bool pass, const std::string& detail) {
        std::cout << "[Splats] " << what << ": " << detail << "   " << (pass ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && pass;
    };

    // a. one isotropic particle: through the center alpha = opacity, at a distance b alpha = opacity exp(-b^2 / 2 s^2); the color is composited once
    {
        const Vector3 rgb(1.0, 0.5, 0.25);
        GaussianSplats cloud({ MakeSplat(Point3(0.0), Vector3(0.1), glm::dvec4(1, 0, 0, 0), 0.8, rgb) }, 0, PlainParams());
        Vector3 L;
        Real T;
        cloud.gather(Ray(Point3(0, 0, -1), Vector3(0, 0, 1)), constants::INFINITY_VAL, L, T);
        const bool ok1 = std::abs(T - 0.2) < 1e-6 && glm::length(L - rgb * 0.8) < 1e-6;   // the parameters are stored as float
        cloud.gather(Ray(Point3(0.1, 0, -1), Vector3(0, 0, 1)), constants::INFINITY_VAL, L, T);
        const Real alpha = 0.8 * std::exp(-0.5);
        const bool ok2 = std::abs(T - (1 - alpha)) < 1e-6 && glm::length(L - rgb * alpha) < 1e-6;
        cloud.gather(Ray(Point3(0, 0, -1), Vector3(0, 0, 1)), 0.5, L, T);   // a surface before the particle: nothing
        const bool ok3 = T == 1 && glm::length(L) == 0;
        check("one particle: alpha through the center, at one sigma, and behind a surface", ok1 && ok2 && ok3,
            "T " + std::to_string(T) + ", radiance of the center ray " + std::to_string(0.8) + " x color");
    }

    // b. an anisotropic, rotated particle on random rays: the closed-form response equals the numerical maximum of the Gaussian along the ray
    {
        const SplatSource s = MakeSplat(Point3(0.3, -0.2, 0.5), Vector3(0.05, 0.2, 0.1), glm::dvec4(0.8, 0.3, -0.4, 0.3), 0.9, Vector3(0.7));
        GaussianSplats cloud({ s }, 0, PlainParams());
        int tested = 0, bad = 0;
        Real maxDiff = 0;
        std::vector<std::pair<Real, Real>> hits;
        for (int i = 0; i < 2000; ++i) {
            const Point3 o = Point3(s.mean) + Vector3(sampling::Random(-0.3, 0.3), sampling::Random(-0.3, 0.3), sampling::Random(-0.3, 0.3)) - Vector3(0, 0, 2);
            const Vector3 d = glm::normalize(Vector3(sampling::Random(-0.3, 0.3), sampling::Random(-0.3, 0.3), 1.0));
            const Ray r(o, d);
            Real tNum;
            const Real aNum = NumericAlpha(s, r, tNum);
            cloud.collectBruteForce(r, constants::INFINITY_VAL, hits);
            if (aNum < Real(1) / 255) { if (!hits.empty()) ++bad; continue; }   // below alphaMin: the engine must not report it either
            ++tested;
            if (hits.size() != 1) { ++bad; continue; }
            const Real diff = std::abs(hits[0].second - aNum) + std::abs(hits[0].first - tNum);
            maxDiff = std::max(maxDiff, diff);
            if (diff > 1e-6) ++bad;
        }
        check("rotated anisotropic particle, closed form vs numerical maximum on " + std::to_string(tested) + " rays", bad == 0 && tested > 500,
            "largest difference (alpha + t) " + std::to_string(maxDiff) + ", " + std::to_string(bad) + " bad");
    }

    // c. compositing: five particles along a ray, inserted in a scrambled order, against the front-to-back sum; a surface at 3.5 cuts the rest
    {
        const Real a[5] = { 0.5, 0.3, 0.9, 0.2, 0.6 };
        const Vector3 c[5] = { Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1), Vector3(1, 1, 0), Vector3(0.5) };
        std::vector<SplatSource> src;
        for (int k : { 3, 0, 4, 1, 2 }) src.push_back(MakeSplat(Point3(0, 0, 1 + k), Vector3(0.05), glm::dvec4(1, 0, 0, 0), a[k], c[k]));
        GaussianSplats cloud(src, 0, PlainParams());
        auto expected = [&](int n, Vector3& L, Real& T) {
            L = Vector3(0.0); T = 1;
            for (int k = 0; k < n; ++k) { L += T * a[k] * c[k]; T *= 1 - a[k]; }
        };
        Vector3 L, Le;
        Real T, Te;
        cloud.gather(Ray(Point3(0, 0, 0), Vector3(0, 0, 1)), constants::INFINITY_VAL, L, T);
        expected(5, Le, Te);
        const bool ok1 = glm::length(L - Le) < 1e-6 && std::abs(T - Te) < 1e-6;
        cloud.gather(Ray(Point3(0, 0, 0), Vector3(0, 0, 1)), 3.5, L, T);
        expected(3, Le, Te);
        const bool ok2 = glm::length(L - Le) < 1e-6 && std::abs(T - Te) < 1e-6;
        const Real Ts = cloud.transmittance(Ray(Point3(0, 0, 0), Vector3(0, 0, 1)), 3.5);
        check("compositing of 5 particles (scrambled order) and the cut by a surface", ok1 && ok2 && std::abs(Ts - Te) < 1e-6,
            "T " + std::to_string(T) + " = " + std::to_string(Te) + ", transmittance() agrees");
    }

    // d. BVH vs brute force: a random cloud, random rays, the same (sorted) particles hence the same radiance and transmittance
    {
        std::vector<SplatSource> src;
        for (int i = 0; i < 3000; ++i) {
            const Point3 m(sampling::Random(-1, 1), sampling::Random(-1, 1), sampling::Random(-1, 1));
            const Vector3 sc(sampling::Random(0.01, 0.08), sampling::Random(0.01, 0.08), sampling::Random(0.01, 0.08));
            const glm::dvec4 q(sampling::Random(-1, 1), sampling::Random(-1, 1), sampling::Random(-1, 1), sampling::Random(-1, 1));
            src.push_back(MakeSplat(m, sc, q, sampling::Random(0.05, 0.95), Vector3(sampling::Random(), sampling::Random(), sampling::Random())));
        }
        GaussianSplats cloud(src, 0, PlainParams());
        int mismatches = 0, nonEmpty = 0;
        std::vector<std::pair<Real, Real>> hits;
        for (int i = 0; i < 5000; ++i) {
            const Point3 o(sampling::Random(-1.5, 1.5), sampling::Random(-1.5, 1.5), sampling::Random(-1.5, 1.5));
            const Vector3 d = glm::normalize(Vector3(sampling::Random(-1, 1), sampling::Random(-1, 1), sampling::Random(-1, 1)));
            const Ray r(o, d);
            const Real tMax = sampling::Random() < 0.5 ? constants::INFINITY_VAL : sampling::Random(0.2, 2.0);
            cloud.collectBruteForce(r, tMax, hits);
            std::sort(hits.begin(), hits.end());
            Real Te = 1;
            for (const auto& h : hits) { Te *= 1 - h.second; if (Te < 1e-3) { Te = 0; break; } }
            Vector3 L;
            Real T;
            cloud.gather(r, tMax, L, T);
            if (!hits.empty()) ++nonEmpty;
            if (std::abs(T - Te) > 1e-9) ++mismatches;
        }
        check("BVH (" + std::to_string(cloud.numNodes()) + " nodes) vs brute force on 5000 rays (" + std::to_string(nonEmpty) + " touch the cloud), transmittance mismatches",
            mismatches == 0 && nonEmpty > 1000, std::to_string(mismatches));
    }

    // e. the 3DGS PLY: written with the file's raw values (logit opacity, log scale, unnormalized quaternion, channel-major f_rest), read back
    {
        const std::string path = (std::filesystem::temp_directory_path() / "rayt_test_splats.ply").string();
        const int nRest = 45;
        // particle 0: degree-1 harmonics with known values; particles 1, 2: filler
        std::vector<float> rest(nRest, 0.f);
        const float dc[3] = { 0.9f, -0.3f, 0.2f };
        const float l1[3][3] = { { 0.1f, -0.2f, 0.3f }, { 0.05f, 0.15f, -0.25f }, { -0.3f, 0.1f, 0.2f } };   // [channel][m = -1, 0, 1]
        for (int c = 0; c < 3; ++c) for (int m = 0; m < 3; ++m) rest[static_cast<size_t>(c * 15 + m)] = l1[c][m];
        {
            std::ofstream f(path, std::ios::binary);
            f << "ply\nformat binary_little_endian 1.0\nelement vertex 3\nproperty float x\nproperty float y\nproperty float z\nproperty float nx\nproperty float ny\nproperty float nz\n"
              << "property float f_dc_0\nproperty float f_dc_1\nproperty float f_dc_2\n";
            for (int i = 0; i < nRest; ++i) f << "property float f_rest_" << i << "\n";
            f << "property float opacity\nproperty float scale_0\nproperty float scale_1\nproperty float scale_2\nproperty float rot_0\nproperty float rot_1\nproperty float rot_2\nproperty float rot_3\nend_header\n";
            auto put = [&](float v) { f.write(reinterpret_cast<const char*>(&v), 4); };
            for (int i = 0; i < 3; ++i) {
                put(float(i)); put(0.f); put(0.f); put(0.f); put(0.f); put(0.f);
                for (float v : dc) put(v);
                for (float v : rest) put(i == 0 ? v : 0.f);
                put(1.0f);                                  // opacity logit 1 -> sigmoid = 0.7310586
                put(std::log(0.1f)); put(std::log(0.2f)); put(std::log(0.3f));
                put(2.f); put(0.f); put(0.f); put(0.f);     // quaternion (2, 0, 0, 0) -> normalized to identity
            }
        }
        std::vector<SplatSource> src;
        int degree = -1;
        std::string err;
        const bool loaded = GaussianSplats::loadPly(path, src, degree, &err);
        bool ok = loaded && src.size() == 3 && degree == 3;
        if (ok) {
            const SplatSource& s = src[0];
            ok = std::abs(s.opacity - 1 / (1 + std::exp(-1.0))) < 1e-6 && glm::length(s.scale - Vector3(0.1, 0.2, 0.3)) < 1e-6
                && glm::length(s.rot - glm::dvec4(1, 0, 0, 0)) < 1e-12 && s.sh.size() == 48 && std::abs(s.sh[3 + 15] - l1[1][0]) < 1e-7;
            // the color for a direction, from the degree-1 formula by hand: C0 dc + C1 (-y sh1 + z sh2 - x sh3) + 0.5
            SplatParams p = PlainParams();
            p.degree = 1;
            GaussianSplats cloud(src, degree, p);
            const Vector3 d = glm::normalize(Vector3(0.3, -0.5, 0.8));
            const Real C1 = 0.4886025119029199;
            Vector3 want;
            for (int c = 0; c < 3; ++c) want[c] = std::max(0.0, kC0 * dc[c] + C1 * (-d.y * l1[c][0] + d.z * l1[c][1] - d.x * l1[c][2]) + 0.5);
            const Vector3 got = cloud.colorForDirection(0, d);
            ok = ok && glm::length(got - want) < 1e-6;
            check("3DGS PLY round trip: 3 particles, degree 3, activations, channel-major harmonics, degree-1 color by hand", ok,
                loaded ? "color (" + std::to_string(got.x) + ", " + std::to_string(got.y) + ", " + std::to_string(got.z) + ")" : err);
        }
        else check("3DGS PLY round trip", false, loaded ? "count " + std::to_string(src.size()) + ", degree " + std::to_string(degree) : err);
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }

    // f. the scene transform: scale 2, a quarter turn about y, a translation; the harmonics turn with the cloud
    {
        SplatSource s = MakeSplat(Point3(1, 0, 0), Vector3(0.1, 0.2, 0.3), glm::dvec4(1, 0, 0, 0), 0.75, Vector3(0.6));
        s.sh.resize(12, 0.f);
        s.sh[3 + 2] = 1.0f;   // red channel, the "x" term (m = 1 of degree 1): brighter toward -x in the cloud's frame (the term is -C1 x)
        SplatParams p = PlainParams();
        p.scale = 2;
        p.rotation = glm::dmat3(glm::rotate(glm::dmat4(1.0), glm::radians(90.0), Vector3(0, 1, 0)));
        p.translate = Vector3(0, 1, 0);
        GaussianSplats cloud({ s }, 1, p);
        const Point3 center = Point3(p.rotation * Vector3(2, 0, 0) + p.translate);   // (0, 1, -2)
        Vector3 L;
        Real T;
        cloud.gather(Ray(center + Vector3(0, 0, -5), Vector3(0, 0, 1)), constants::INFINITY_VAL, L, T);
        const bool okCenter = std::abs(T - 0.25) < 1e-9;
        const AABB& b = cloud.bounds();
        const Vector3 bc = (b.min + b.max) * 0.5;
        const bool okBounds = glm::length(bc - Vector3(center)) < 1e-9 && b.max.x - b.min.x > b.max.z - b.min.z;   // the 0.3 axis (z of the cloud) now lies along world x
        // the harmonics: the world direction R x_cloud must give the color of the cloud direction x_cloud, by hand red = C0 dc - C1 + 0.5
        const Vector3 xWorld = p.rotation * Vector3(1, 0, 0);
        const Vector3 cw = cloud.colorForDirection(0, xWorld);
        const Vector3 cl(std::max(0.0, 0.1 - 0.4886025119029199 + 0.5), 0.6, 0.6);
        const Vector3 cOther = cloud.colorForDirection(0, Vector3(1, 0, 0));   // world x is the cloud's z: no x term
        const bool okSh = glm::length(cw - cl) < 1e-6 && glm::length(cOther - Vector3(0.6)) < 1e-6;
        check("scene transform: center, bounds, harmonics turned with the cloud", okCenter && okBounds && okSh,
            "T " + std::to_string(T) + ", center (" + std::to_string(bc.x) + ", " + std::to_string(bc.y) + ", " + std::to_string(bc.z) + ")");
    }

    // g. the integrator: a black slab of alpha 1/2 between a quad light and the ground. The camera ray crosses it once (x 1/2) and every
    //    ray from the ground to the light crosses it once more (shadow rays through shadowTransmittance, BSDF-sampled rays through
    //    gather): the ground is seen at a quarter of its radiance.
    {
        const Point3 corner(-0.5, 2, -0.5);
        const Vector3 e1(1, 0, 0), e2(0, 0, 1);
        auto mesh = std::make_shared<TriangleMesh>();
        mesh->p = { corner, corner + e1, corner + e1 + e2, corner + e2 };
        mesh->indices = { 0, 1, 2, 0, 2, 3 };
        const std::shared_ptr<const TriangleMesh> meshC = mesh;
        const Real albedo = 0.5;
        auto scene = MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
            auto mat = std::make_shared<DiffuseLight>(EmissionSpectrum(RgbColor(1.0)), 10.0);
            for (size_t i = 0; i < meshC->numTriangles(); ++i) {
                auto tri = std::make_shared<Triangle>(meshC, i, mat);
                auto light = std::make_shared<AreaLight>(tri, mat);
                w.add(std::make_shared<Primitive>(tri, mat, light));
                s.addLight(light);
            }
        });
        // a flat, wide particle at y = 1 (standard deviation 1 mm across, 1 km along): alpha = 1/2 for any ray crossing the plane
        SplatSource slab = MakeSplat(Point3(0, 1, 0), Vector3(1000, 0.001, 1000), glm::dvec4(1, 0, 0, 0), 0.5, Vector3(0.0));
        scene->setSplats(std::make_shared<const GaussianSplats>(std::vector<SplatSource>{ slab }, 0, PlainParams()));
        Vector3 m, se;
        int bad = 0;
        EstimateRadiance(*scene, Ray(Point3(-3, 3, 0), glm::normalize(Vector3(1, -1, 0))), opt.spp * 1024, 8, m, se, bad);
        const Real expected = 0.25 * QuadGroundRadiance(corner, e1, e2, false, 10.0, albedo);
        check("integrator: a slab of alpha 1/2 under the quad light, crossed by the camera ray and the light's rays: ground radiance x 1/4 = " + std::to_string(expected),
            std::abs(m.x - expected) <= 4 * se.x + 0.005 * expected && bad == 0, std::to_string(m.x) + " +- " + std::to_string(se.x));
    }

    std::cout << "[Splats] checks: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed;
}

} // namespace rayt::cli
