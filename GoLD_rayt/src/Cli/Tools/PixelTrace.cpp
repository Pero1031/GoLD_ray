// Debug aid: prints a few light paths of one pixel (--trace x,y).
#include "pch.h"
#include "Cli/Tools.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

void TracePixelDebug(const Scene& scene, const Camera& camera, int px, int py, int w, int h,
    int maxDepth, int numPaths,
    const std::vector<std::pair<std::string, const Material*>>& names)
{
    auto label = [&](const Material* m) {
        for (const auto& [n, p] : names) if (p == m) return n;
        return std::string("?");
    };
    const SampledWavelengths lambda;   // RGB-mode default (wavelength-independent materials)
    std::cout << "[Trace] pixel (" << px << ", " << py << ") of " << w << "x" << h << ", " << numPaths << " paths\n";
    for (int k = 0; k < numPaths; ++k) {
        // Same convention as PathIntegrator: image row 0 is the top, camera t = 1 at the top.
        Ray r = camera.getRay((Real(px) + 0.5) / Real(w), Real(1) - (Real(py) + 0.5) / Real(h), Point2(0.5, 0.5));
        std::cout << "  path " << k << ":";
        for (int depth = 0; depth < maxDepth; ++depth) {
            SurfaceInteraction rec;
            if (!scene.hit(r, rec)) {
                std::cout << " -> ENV [ray o=(" << r.o.x << "," << r.o.y << "," << r.o.z << ") d=("
                    << r.d.x << "," << r.d.y << "," << r.d.z << ") tMin=" << r.tMin << " tMax=" << r.tMax << "]";
                break;
            }
            std::cout << " -> " << label(rec.matPtr) << (rec.frontFace ? "[in]" : "[out]")
                << "@(" << rec.p.x << "," << rec.p.y << "," << rec.p.z << ")";
            const auto s = rec.matPtr->sample(rec, -r.d, Point2(sampling::Random(), sampling::Random()), lambda);
            if (!s) { std::cout << " (no sample: " << (rec.frontFace ? "?" : "TIR") << ")"; break; }
            std::cout << (s->isTransmission() ? "(T)" : "(R)");
            r = rec.SpawnRay(s->wi);
        }
        std::cout << "\n";
    }
}

} // namespace rayt::cli
