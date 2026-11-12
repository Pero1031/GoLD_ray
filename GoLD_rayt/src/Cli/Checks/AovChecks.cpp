// Denoiser feature (AOV) checks.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/SceneParts.hpp"
#include "Geometry/Sphere.hpp"
#include "Materials/CoatedDiffuse.hpp"
#include "Materials/Dielectric.hpp"
#include "Materials/Lambertian.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Integrator.hpp"

namespace rayt::cli {

int RunAovTests() {
    bool allPassed = true;
    const Vector3 groundAlbedo(0.3, 0.5, 0.7), coatedAlbedo(0.8, 0.2, 0.1), lightColor(1.0, 0.4, 0.1);
    const int maxDepth = 8;

    auto buildScene = [&](bool withSky) {
        auto world = std::make_shared<Aggregate>();
        world->add(std::make_shared<Sphere>(Point3(0, -1000, 0), 1000.0, std::make_shared<Lambertian>(RgbColor(groundAlbedo))));
        world->add(std::make_shared<Sphere>(Point3(-1.5, 0.5, 0), 0.5, std::make_shared<Dielectric>(1.5, 0.0)));
        world->add(std::make_shared<Sphere>(Point3(1.5, 0.5, 0), 0.5,
            std::make_shared<CoatedDiffuse>(RgbColor(coatedAlbedo), Real(0.4))));
        auto scene = std::make_shared<Scene>(world);
        AddSphereLight(*world, *scene, Point3(0, 3, -3), 0.5, EmissionSpectrum(RgbColor(lightColor)), 10.0, true);
        if (withSky) {
            auto sky = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(0.5))));
            scene->setEnvLight(std::make_shared<EnvLight>(sky, 1.0));
        }
        scene->finalize("spatial");
        return scene;
    };
    const auto sceneSky = buildScene(true);
    const auto sceneNoSky = buildScene(false);

    struct Result { int n = 0, valid = 0; Vector3 albedo = Vector3(0.0), normal = Vector3(0.0); int nearGround = 0; };
    auto trace = [&](const Scene& scene, const Ray& ray, int n) {
        Result res;
        PathIntegrator integrator(nullptr, maxDepth, 1);
        for (int i = 0; i < n; ++i) {
#if defined(RAYT_SPECTRAL_MODE)
            SampledWavelengths lambda = SampledWavelengths::sampleVisible(sampling::Random());
#else
            SampledWavelengths lambda;
#endif
            PathAOV aov;
            integrator.Li(ray, scene, lambda, &aov);
            ++res.n;
            if (!aov.valid) continue;
            ++res.valid;
            res.albedo += aov.albedo;
            res.normal += aov.normal;
            if (glm::length(aov.albedo - groundAlbedo) < 1e-6) ++res.nearGround;
        }
        if (res.valid > 0) { res.albedo /= Real(res.valid); res.normal /= Real(res.valid); }
        return res;
    };
    auto report = [&](const std::string& name, bool ok) {
        std::cout << "[Aov]   " << name << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };
    auto close = [](const Vector3& a, const Vector3& b, Real tol) { return glm::length(a - b) <= tol; };

    std::cout << "[Aov] albedo / normal features of a camera path\n";

    {   // diffuse ground (a huge sphere, so its normal is +y only to within 0.002)
        const Result r = trace(*sceneSky, Ray(Point3(0.5, 2, 2), Vector3(0, -1, 0)), 200);
        report("diffuse ground: albedo = material color, normal = +y",
            r.valid == r.n && close(r.albedo, groundAlbedo, 1e-6) && close(r.normal, Vector3(0, 1, 0), 5e-3));
    }
    {   // coated diffuse sphere, hit head-on (normal points back at the ray origin)
        const Result r = trace(*sceneSky, Ray(Point3(1.5, 0.5, 3), Vector3(0, 0, -1)), 200);
        report("coated diffuse sphere: albedo = base color, normal = +z",
            r.valid == r.n && close(r.albedo, coatedAlbedo, 1e-6) && close(r.normal, Vector3(0, 0, 1), 1e-6));
    }
    {   // emitter: the hue of the emission, scaled to a maximum of 1
        const Result r = trace(*sceneSky, Ray(Point3(0, 3, 0), Vector3(0, 0, -1)), 200);
        report("sphere light: albedo = emission hue (max 1), normal = +z",
            r.valid == r.n && close(r.albedo, lightColor, 1e-6) && close(r.normal, Vector3(0, 0, 1), 1e-6));
    }
    {   // smooth glass sphere over the ground: the ray is followed through it, so nearly every path ends
        // with the ground color (the rest are Fresnel reflections that escape to the sky). If the glass
        // itself were recorded, the albedo would be white (the Material default) and this fraction ~0.
        const Result r = trace(*sceneSky, Ray(Point3(-1.5, 3, 0), Vector3(0, -1, 0)), 2000);
        const double groundFraction = double(r.nearGround) / double(r.n);
        std::cout << "[Aov]     (through glass: " << 100.0 * groundFraction << " % of the paths record the ground color)\n";
        // (about 1 path in 2000 bounces around inside the sphere until the depth limit and records nothing)
        report("smooth glass is skipped: albedo is the ground behind it", r.valid >= 0.99 * r.n && groundFraction > 0.85);
    }
    {   // a ray that leaves the scene
        const Result sky = trace(*sceneSky, Ray(Point3(0, 5, 5), Vector3(0, 1, 0)), 2000);
        report("escaping ray with an environment light: albedo = environment color, normal = 0",
            sky.valid == sky.n && close(sky.albedo, Vector3(0.5), 0.03) && close(sky.normal, Vector3(0.0), 1e-9));
        const Result none = trace(*sceneNoSky, Ray(Point3(0, 5, 5), Vector3(0, 1, 0)), 50);
        report("escaping ray without environment: no feature recorded", none.valid == 0);
    }

    std::cout << "[Aov] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
