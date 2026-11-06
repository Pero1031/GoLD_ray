// White-furnace check: every object must return its expected mean color.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

bool CheckFurnace(const Scene& scene, const Camera& camera, const Film& film,
    const std::vector<FurnaceObject>& objects)
{
    struct Acc { Vector3 sum = Vector3(0.0); long long n = 0; };
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
            a.sum += film.getData()[static_cast<size_t>(h - 1 - j) * static_cast<size_t>(w) + static_cast<size_t>(i)];
            a.n++;
        }
    }

    std::cout << "\n[Furnace] mean linear RGB per object vs the expected value (+- tolerance, per channel)\n";
    bool allPassed = true;
    for (const FurnaceObject& obj : objects) {
        const auto it = acc.find(obj.material);
        if (it == acc.end() || it->second.n == 0) {
            std::cout << "[Furnace]   " << obj.label << ": no pixels found -> FAIL\n";
            allPassed = false;
            continue;
        }
        const Vector3 mean = it->second.sum / Real(it->second.n);
        const bool ok = std::abs(mean.x - obj.expected) <= obj.tolerance
            && std::abs(mean.y - obj.expected) <= obj.tolerance
            && std::abs(mean.z - obj.expected) <= obj.tolerance;
        allPassed = allPassed && ok;
        std::cout << "[Furnace]   " << obj.label << ": (" << mean.x << ", " << mean.y << ", " << mean.z
            << ")  expected " << obj.expected << " +- " << obj.tolerance
            << "  [" << it->second.n << " px]  " << (ok ? "OK" : "FAIL") << "\n";
    }
    std::cout << "[Furnace] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed;
}

} // namespace rayt::cli
