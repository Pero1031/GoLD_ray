// An object mask and a depth map from camera rays alone (--mask, --depth-map); see Render/MaskDepth.hpp.
#include "pch.h"
#include "Render/MaskDepth.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <iostream>

#include "Core/Interaction.hpp"
#include "Core/Sampling.hpp"
#include "Materials/Material.hpp"
#include "Render/Camera.hpp"
#include "Render/Film.hpp"
#include "Scene/Scene.hpp"

namespace rayt {

    namespace {

        /// The first surface the ray meets. Crossed: invisible boundaries, and translucent layers that let the ray through
        /// (stochastic: with probability 1 - opacity, as the path tracer does; else when the opacity is below 1/2).
        bool FirstSurface(const Scene& scene, Ray r, bool stochastic, SurfaceInteraction& rec) {
            for (int k = 0; k < 64; ++k) {
                if (!scene.hit(r, rec)) return false;
                bool through = rec.matPtr->isNullInterface();
                if (!through) {
                    const Real opacity = rec.matPtr->opacity(rec);
                    if (opacity < Real(1)) through = stochastic ? (sampling::Random() >= opacity) : (opacity < Real(0.5));
                }
                if (!through) return true;
                r = rec.SpawnRay(r.d);
            }
            return true;   // as the path tracer: after 64 layers the ray stops at the last one
        }

    } // namespace

    void RenderMaskAndDepth(const Scene& scene, const Camera& camera, Film& film, const std::vector<const Material*>& maskMaterials, int samplesPerSide) {
        const int width = film.width(), height = film.height();
        const bool wantMask = film.hasMask(), wantDepth = film.hasDepth();
        if (!wantMask && !wantDepth) return;
        const int n = std::max(samplesPerSide, 1);
        const Vector3 axis = camera.viewAxis();
        auto masked = [&](const Material* m) { return std::find(maskMaterials.begin(), maskMaterials.end(), m) != maskMaterials.end(); };

#pragma omp parallel for schedule(dynamic, 1)
        for (int j = 0; j < height; ++j) {
            sampling::SeedStream(3, 0, static_cast<std::uint64_t>(j));   // its own stream: the picture's random numbers are untouched
            for (int i = 0; i < width; ++i) {
                SurfaceInteraction rec;
                if (wantMask) {
                    int covered = 0;
                    for (int a = 0; a < n; ++a)
                        for (int b = 0; b < n; ++b) {
                            const Real u = (Real(i) + (Real(a) + sampling::Random()) / Real(n)) / Real(width);
                            const Real v = (Real(j) + (Real(b) + sampling::Random()) / Real(n)) / Real(height);
                            const Point2 lens = sampling::Random2D();
                            const Real time = camera.hasShutter() ? sampling::Random() : Real(0);
                            const Ray r = camera.getRay(u, v, lens, time);
                            if (FirstSurface(scene, r, true, rec) && masked(rec.matPtr)) ++covered;
                        }
                    film.setMask(i, height - 1 - j, Real(covered) / Real(n * n));
                }
                if (wantDepth) {
                    const Ray r = camera.centerRay((Real(i) + Real(0.5)) / Real(width), (Real(j) + Real(0.5)) / Real(height));
                    const Real z = FirstSurface(scene, r, false, rec) ? glm::dot(rec.p - camera.origin(), axis) : Real(0);
                    film.setDepth(i, height - 1 - j, z);
                }
            }
        }

        // A one-line summary: how much the mask covers, how far the surfaces are.
        if (wantMask) {
            double sum = 0;
            for (const float c : film.mask()) sum += c;
            std::cout << "[MaskDepth] mask: " << n << " x " << n << " samples per pixel, the objects cover "
                << 100.0 * sum / double(std::max<size_t>(film.mask().size(), 1)) << " % of the picture\n";
        }
        if (wantDepth) {
            float lo = 0, hi = 0;
            size_t hits = 0;
            for (const float z : film.depth()) {
                if (!(z > 0)) continue;
                lo = (hits == 0) ? z : std::min(lo, z);
                hi = std::max(hi, z);
                ++hits;
            }
            std::cout << "[MaskDepth] depth: " << (hits ? std::to_string(lo) + " .. " + std::to_string(hi) + " (scene units) along the view axis" : std::string("no surface seen"))
                << ", " << 100.0 * double(film.depth().size() - hits) / double(std::max<size_t>(film.depth().size(), 1)) << " % of the pixels see none (0)\n";
        }
    }

} // namespace rayt
