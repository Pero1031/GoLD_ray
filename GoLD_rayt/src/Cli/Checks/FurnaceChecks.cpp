// White-furnace check: every object must return its expected mean color.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Render/Camera.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Film.hpp"
#include "Render/Integrator.hpp"
#include "Scene/Scene.hpp"

namespace rayt::cli {

bool CheckFurnace(const Scene& scene, const Camera& camera, const Film& film,
    const std::vector<FurnaceObject>& objects)
{
    struct Acc { Vector3 sum = Vector3(0.0), sum2 = Vector3(0.0); long long n = 0; };
    std::map<const Material*, Acc> acc;   // key nullptr = background (environment)

    const int w = film.width();
    const int h = film.height();
    for (int j = 0; j < h; ++j) {
        for (int i = 0; i < w; ++i) {
            // Same pixel convention as PathIntegrator::render (film row = h - 1 - j).
            const Ray r = camera.getRay((Real(i) + 0.5) / Real(w), (Real(j) + 0.5) / Real(h), Point2(0.5, 0.5));
            SurfaceInteraction rec;
            const Material* m = scene.hit(r, rec) ? rec.matPtr : nullptr;
            Acc& a = acc[m];
            const Vector3 v = film.getData()[static_cast<size_t>(h - 1 - j) * static_cast<size_t>(w) + static_cast<size_t>(i)];
            a.sum += v;
            a.sum2 += v * v;
            a.n++;
        }
    }

    // The mean over an object's pixels, judged against the expected value with 4 standard errors (from the spread of the pixels: their own noise
    // and the true variation over the object together, so this is a bound) plus the object's documented systematic term. A fixed tolerance alone
    // failed at a low --spp (the dispersive glass in the spectral build: a few samples per pixel carry one wavelength each).
    std::cout << "\n[Furnace] mean linear RGB per object vs the expected value (4 standard errors from the pixels + the systematic term, per channel)\n";
    bool allPassed = true;
    for (const FurnaceObject& obj : objects) {
        const auto it = acc.find(obj.material);
        if (it == acc.end() || it->second.n == 0) {
            std::cout << "[Furnace]   " << obj.label << ": no pixels found -> FAIL\n";
            allPassed = false;
            continue;
        }
        const Real n = Real(it->second.n);
        const Vector3 mean = it->second.sum / n;
        const Vector3 se = glm::sqrt(glm::max(it->second.sum2 / n - mean * mean, Vector3(0.0)) / n);
        bool ok = true;
        for (int c = 0; c < 3; ++c) ok = ok && std::abs(mean[c] - obj.expected) <= 4 * se[c] + obj.tolerance;
        allPassed = allPassed && ok;
        std::cout << "[Furnace]   " << obj.label << ": (" << mean.x << ", " << mean.y << ", " << mean.z << ") +- (" << se.x << ", " << se.y << ", " << se.z
            << ")  expected " << obj.expected << " +- " << obj.tolerance << " systematic"
            << "  [" << it->second.n << " px]  " << (ok ? "OK" : "FAIL") << "\n";
    }
    std::cout << "[Furnace] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed;
}

} // namespace rayt::cli
