// Report of the metalprobe scene: mean film color per angle of incidence.
#include "pch.h"
#include "Cli/Tools.hpp"
#include "Render/Camera.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Film.hpp"
#include "Render/Integrator.hpp"
#include "Scene/Scene.hpp"

namespace rayt::cli {

void ReportMetalProbe(const Scene& scene, const Camera& camera, const Film& film,
    const Material* metal, const std::string& label)
{
    const Real centersDeg[] = { 0, 15, 30, 45, 60, 70, 75, 80, 85 };
    constexpr int kBins = static_cast<int>(sizeof(centersDeg) / sizeof(centersDeg[0]));
    constexpr Real kHalfWidthDeg = 1.5;
    struct Acc { Vector3 sum = Vector3(0.0); Real cosSum = 0; long long n = 0; };
    Acc acc[kBins];

    const int w = film.width();
    const int h = film.height();
    for (int j = 0; j < h; ++j) {
        for (int i = 0; i < w; ++i) {
            const Ray r = camera.getRay((Real(i) + 0.5) / Real(w), (Real(j) + 0.5) / Real(h), Point2(0.5, 0.5));
            SurfaceInteraction rec;
            if (!scene.hit(r, rec) || rec.matPtr != metal) continue;
            const Real cosTheta = std::abs(glm::dot(r.d, rec.n));
            const Real thetaDeg = std::acos(std::min(cosTheta, Real(1))) * Real(180) / constants::PI;
            for (int b = 0; b < kBins; ++b) {
                if (std::abs(thetaDeg - centersDeg[b]) > kHalfWidthDeg) continue;
                Acc& a = acc[b];
                a.sum += film.getData()[static_cast<size_t>(h - 1 - j) * static_cast<size_t>(w) + static_cast<size_t>(i)];
                a.cosSum += cosTheta;
                a.n++;
            }
        }
    }

    std::cout << "\n[MetalProbe] " << label << ": mean linear RGB per angle of incidence\n";
    std::cout << std::fixed << std::setprecision(4);
    for (int b = 0; b < kBins; ++b) {
        if (acc[b].n == 0) continue;
        const Vector3 mean = acc[b].sum / Real(acc[b].n);
        std::cout << "[MetalProbe]   theta " << std::setw(4) << std::setprecision(1) << centersDeg[b]
            << "  mean cos " << std::setprecision(4) << (acc[b].cosSum / Real(acc[b].n))
            << "  RGB (" << mean.x << ", " << mean.y << ", " << mean.z << ")  [" << acc[b].n << " px]\n";
    }
    std::cout.unsetf(std::ios::floatfield);
}

} // namespace rayt::cli
