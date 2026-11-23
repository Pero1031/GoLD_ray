// Experimental cel shading with outlines (--toon); see Render/ToonIntegrator.hpp.
#include "pch.h"
#include "Render/ToonIntegrator.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <vector>

#include "Core/Interaction.hpp"
#include "Core/Sampling.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Color/SpectralConversion.hpp"
#include "Lights/EnvLight.hpp"
#include "Materials/Material.hpp"
#include "Render/Camera.hpp"
#include "Render/Film.hpp"
#include "Scene/Scene.hpp"

namespace rayt {

    namespace {

        /// What one sub-sample saw first (for the outlines and the rim) and its color.
        struct Sub {
            Vector3 color = Vector3(0.0);
            Vector3 albedo = Vector3(0.0);                         // linear, for the rim
            Vector3 p = Vector3(0.0);
            Vector3 n = Vector3(0.0);
            Real depth = std::numeric_limits<Real>::infinity();   // along the view axis; infinite: nothing hit
            Real footprint = 0;                                    // the sub-sample's size at that depth
            const void* id = nullptr;                              // the material (null: nothing hit)
            Real lit = 0;                                          // 1: lit, 0: in shadow (shaded directly only)
            const ToonLook* look = nullptr;                        // the material's own recipe, if any
            Vector3 lineColor = Vector3(0.0);                      // linear
            Real lineWidth = 1;                                    // pixels
        };

        /// Widens a mask by `radius` cells in x and y (a square), separably.
        void Dilate(std::vector<unsigned char>& mask, int w, int h, int radius) {
            if (radius <= 0) return;
            std::vector<unsigned char> tmp(mask.size(), 0);
#pragma omp parallel for schedule(static)
            for (int y = 0; y < h; ++y) {
                int last = -1000000;   // the last marked cell at or left of x + radius
                for (int x = -radius; x < w; ++x) {
                    const int xr = x + radius;
                    if (xr < w && mask[static_cast<size_t>(y) * w + xr]) last = xr;
                    if (x >= 0) tmp[static_cast<size_t>(y) * w + x] = (last >= x - radius) ? 1 : 0;
                }
            }
#pragma omp parallel for schedule(static)
            for (int x = 0; x < w; ++x) {
                int last = -1000000;
                for (int y = -radius; y < h; ++y) {
                    const int yr = y + radius;
                    if (yr < h && tmp[static_cast<size_t>(yr) * w + x]) last = yr;
                    if (y >= 0) mask[static_cast<size_t>(y) * w + x] = (last >= y - radius) ? 1 : 0;
                }
            }
        }

        inline Real Smoothstep(Real a, Real b, Real x) {
            const Real t = std::clamp((x - a) / (b - a), Real(0), Real(1));
            return t * t * (Real(3) - Real(2) * t);
        }

        inline Real Luminance(const Vector3& c) { return Real(0.2126) * c.x + Real(0.7152) * c.y + Real(0.0722) * c.z; }

        inline Real ToSrgb(Real c) {
            c = std::clamp(c, Real(0), Real(1));
            return c <= Real(0.0031308) ? Real(12.92) * c : Real(1.055) * std::pow(c, Real(1) / Real(2.4)) - Real(0.055);
        }
        inline Real ToLinear(Real c) {
            c = std::clamp(c, Real(0), Real(1));
            return c <= Real(0.04045) ? c / Real(12.92) : std::pow((c + Real(0.055)) / Real(1.055), Real(2.4));
        }
        inline Vector3 ToSrgb(const Vector3& c) { return Vector3(ToSrgb(c.x), ToSrgb(c.y), ToSrgb(c.z)); }
        inline Vector3 ToLinear(const Vector3& c) { return Vector3(ToLinear(c.x), ToLinear(c.y), ToLinear(c.z)); }

        /// The first surface the ray meets; invisible boundaries and see-through layers are crossed (as Render/MaskDepth.cpp).
        bool FirstSurface(const Scene& scene, Ray r, SurfaceInteraction& rec) {
            for (int k = 0; k < 64; ++k) {
                if (!scene.hit(r, rec)) return false;
                bool through = rec.matPtr->isNullInterface();
                if (!through) {
                    const Real opacity = rec.matPtr->opacity(rec);
                    if (opacity < Real(1)) through = sampling::Random() >= opacity;
                }
                if (!through) return true;
                r = rec.SpawnRay(r.d);
            }
            return true;
        }

        SampledWavelengths PathWavelengths() {
#if defined(RAYT_SPECTRAL_MODE)
            return SampledWavelengths::sampleVisible(sampling::Random());
#else
            return SampledWavelengths();
#endif
        }

        /// A light's sample point is unblocked (aimed at the sampled point for area lights, along wi up to tMax otherwise).
        bool Visible(const Scene& scene, const SurfaceInteraction& ref, const Light& light, const LightLiSample& ls) {
            if (!std::isfinite(ls.tMax) || ls.tMax <= Real(0)) return !scene.occluded(ref.SpawnRay(ls.wi));
            const Point3 origin = ref.OffsetRayOrigin(ls.wi);
            Vector3 dir = ls.wi;
            Real dist = ls.tMax;
            if (light.type() == LightType::Area) {
                const Vector3 to = ls.pLight - origin;
                dist = glm::length(to);
                if (dist <= Real(0)) return true;
                dir = to / dist;
            }
            Ray shadow(origin, dir, Real(0), ref.medium);
            shadow.time = ref.time;
            shadow.tMax = dist * (Real(1) - Real(1e-4));
            return !scene.occluded(shadow);
        }

    } // namespace

    void ToonIntegrator::render(const Scene& scene, Film& film) {
        const auto t0 = std::chrono::steady_clock::now();
        const int width = film.width(), height = film.height();
        // s x s sub-samples per pixel, from the sample count; at most 12 M sub-samples in all
        int s = std::clamp(static_cast<int>(std::lround(std::sqrt(double(std::max(m_spp, 1))))), 1, 4);
        while (s > 1 && double(width) * height * s * s > 12e6) --s;
        const int sw = width * s, sh = height * s;
        std::vector<Sub> subs(static_cast<size_t>(sw) * static_cast<size_t>(sh));
        const ToonSettings& set = m_set;
        const ToonLooks* looks = set.looks.get();
        std::cout << "[Toon] " << width << "x" << height << ", " << s << " x " << s << " sub-samples per pixel, "
            << scene.numLights() << " light(s)" << (scene.hasLights() ? "" : " (a key light from the environment, or from above the camera)")
            << (looks ? ", " + std::to_string(looks->size()) + " materials with their own recipe" : std::string()) << "\n";

        const Vector3 axis = m_camera->viewAxis();
        const Vector3 camRight = m_camera->right();
        const Vector3 camUp = glm::normalize(glm::cross(camRight, axis));
        // A scene without lights of its own: the key light comes from the brightest part of the environment (the strength-weighted
        // mean direction, e.g. the sun), else from above, a little to the left and from behind the camera.
        Vector3 keyFallback = glm::normalize(Vector3(0, 1, 0) * Real(1.0) - camRight * Real(0.6) - axis * Real(0.6));
        if (!scene.hasLights() && scene.hasEnvLight()) {
            SampledWavelengths lambda = PathWavelengths();
            LightSampleContext ctx{ m_camera->origin(), Vector3(0, 1, 0) };
            Vector3 sum(0.0);
            for (int gy = 0; gy < 64; ++gy)
                for (int gx = 0; gx < 64; ++gx) {
                    const auto ls = scene.envLight()->sampleLi(ctx, Point2((gx + Real(0.5)) / 64, (gy + Real(0.5)) / 64), lambda);
                    if (!ls || ls->pdf <= Real(0)) continue;
                    const Real strength = Luminance(SpectrumToRGB(ls->Li, lambda)) / ls->pdf;
                    if (std::isfinite(strength) && strength > Real(0)) sum += strength * ls->wi;
                }
            if (glm::length(sum) > Real(0) && glm::normalize(sum).y > Real(0.05)) keyFallback = glm::normalize(sum);
            std::cout << "[Toon] key light from the environment: direction (" << keyFallback.x << ", " << keyFallback.y << ", " << keyFallback.z << ")\n";
        }

        // The angle and offset between neighbouring sub-samples' rays (perspective and orthographic views)
        const Ray c0 = m_camera->centerRay(Real(0.5), Real(0.5));
        const Ray c1 = m_camera->centerRay(Real(0.5) + Real(1) / Real(sw), Real(0.5));
        const Real subAngle = glm::length(glm::normalize(c1.d) - glm::normalize(c0.d));
        const Real subOffset = glm::length(c1.o - c0.o);

        // A ToonLook's edge size is MMD's: about 1 / 400 of the picture's height per unit (1 pixel for edge 0.5 at 800 pixels)
        const Real lookLineScale = Real(height) / Real(400);
        auto lookOf = [&](const Material* m) -> const ToonLook* {
            if (!looks) return nullptr;
            const auto it = looks->find(m);
            return it == looks->end() ? nullptr : &it->second;
        };

        // The light's share that reaches a point (see the header), and the direction of the light that lights it most.
        auto lightShare = [&](const SurfaceInteraction& rec, const Vector3& n, const Vector3& wo, const SampledWavelengths& lambda, Vector3& keyDir) -> Real {
            LightSampleContext ctx{ rec.p, rec.gn };
            ctx.medium = rec.medium;
            ctx.time = rec.time;
            Real sum = 0, lit = 0, keyBest = 0;
            keyDir = keyFallback;
            for (const auto& light : scene.lights()) {
                // N.L of an area light: the mean over a 6 x 6 grid of points on it (a large light wraps further round, smoothly);
                // the shadow: one ray to the grid point that lights it most (a hard shadow; the centre of the light can be below the
                // horizon where the light still reaches the point). A point light: its one point.
                const int g = (light->type() == LightType::Area) ? 6 : 1;
                Real strengthSum = 0, cosSum = 0, best = 0;
                std::optional<LightLiSample> bestSample;
                for (int gy = 0; gy < g; ++gy)
                    for (int gx = 0; gx < g; ++gx) {
                        const Point2 u((Real(gx) + Real(0.5)) / Real(g), (Real(gy) + Real(0.5)) / Real(g));
                        const auto ls = light->sampleLi(ctx, u, lambda);
                        if (!ls || ls->pdf <= Real(0)) continue;
                        const Real strength = Luminance(SpectrumToRGB(ls->Li, lambda)) / ls->pdf;
                        if (!(strength > Real(0)) || !std::isfinite(strength)) continue;
                        strengthSum += strength;
                        const Real c = glm::dot(n, ls->wi);
                        if (c > Real(0) && glm::dot(rec.gn, ls->wi) * glm::dot(rec.gn, wo) > Real(0)) {
                            cosSum += strength * c;
                            if (strength * c > best) { best = strength * c; bestSample = ls; }
                        }
                    }
                sum += strengthSum;
                if (strengthSum > keyBest && bestSample) { keyBest = strengthSum; keyDir = bestSample->wi; }
                if (cosSum <= Real(0)) continue;
                if (bestSample && Visible(scene, rec, *light, *bestSample)) lit += cosSum;
            }
            if (sum > Real(0)) return lit / sum;
            const Real c = glm::dot(n, keyFallback);
            return (c > Real(0) && !scene.occluded(rec.SpawnRay(keyFallback))) ? c : Real(0);
        };

        // MMD's recipe in display values; returns the linear color, its alpha and the lit factor.
        auto shadeLook = [&](const ToonLook& L, const SurfaceInteraction& rec, const Vector3& n, Real share, const Vector3& keyDir,
            const Vector3& wo, Vector3& albedoLinear, Real& alpha) -> Vector3 {
            Vector3 c = glm::clamp(L.diffuse * set.lightColor + L.ambient, Vector3(0.0), Vector3(1.0));
            alpha = L.alpha;
            if (L.texture) {
                c *= ToSrgb(L.texture->rgb(rec.uv.x, rec.uv.y));
                alpha *= L.texture->alpha(rec.uv.x, rec.uv.y);
            }
            albedoLinear = ToLinear(c);
            if (L.sphere && (L.sphereMode == 1 || L.sphereMode == 2)) {
                const Real su = glm::dot(n, camRight) * Real(0.5) + Real(0.5);
                const Real sv = Real(0.5) - glm::dot(n, camUp) * Real(0.5);
                const Vector3 sp = ToSrgb(L.sphere->rgb(su, sv));
                c = (L.sphereMode == 1) ? c * sp : c + sp;
            }
            const Vector3 shadow = L.toon ? ToSrgb(L.toon->rgb(Real(0.5), Real(1) - Real(0.5) / Real(L.toon->height()))) : L.toonShadow;
            const Real comp = std::clamp(Real(3) * share, Real(0), Real(1));
            c *= shadow + (Vector3(1.0) - shadow) * comp;
            if (L.shininess > Real(0) && (L.specular.x + L.specular.y + L.specular.z) > Real(0)) {
                const Vector3 h = glm::normalize(keyDir + wo);
                c += L.specular * set.lightColor * std::pow(std::max(Real(0), glm::dot(n, h)), L.shininess) * comp;
            }
            return ToLinear(glm::clamp(c, Vector3(0.0), Vector3(1.0)));
        };

        // Shades one camera ray; fills `sub` with the first opaque surface it met.
        auto shade = [&](Ray r, Sub& sub) -> Vector3 {
            SampledWavelengths lambda = PathWavelengths();
            const Point3 eye = r.o;   // the camera ray's origin (r moves on through glass and overlay layers)
            Vector3 weight(1.0), add(0.0);
            bool recorded = false;
            for (int step = 0; step < 16; ++step) {
                SurfaceInteraction rec;
                if (!FirstSurface(scene, r, rec)) {
                    Vector3 bg(0.0);
                    if (!recorded && scene.hasCameraBackground()) bg = scene.cameraBackground();
                    else if (scene.hasEnvLight()) {
                        const Vector3 L = glm::max(SpectrumToRGB(scene.envLight()->Le(r, lambda), lambda), Vector3(0.0));
                        bg = L / (Vector3(1.0) + L);   // a soft compression: the backdrop stays below white
                    }
                    return add + weight * glm::clamp(bg, Vector3(0.0), Vector3(1.0));
                }
                const Material* m = rec.matPtr;
                const ToonLook* look = lookOf(m);
                Vector3 wo = -glm::normalize(r.d);

                // An eyebrow just behind the hair shows through it
                if (look && look->hair && !recorded && set.hairSeeThrough > Real(0)) {
                    Ray probe = rec.SpawnRay(r.d);
                    probe.tMax = set.hairSeeThrough;
                    SurfaceInteraction behind;
                    for (int k = 0; k < 4 && scene.hit(probe, behind); ++k) {
                        const ToonLook* bl = lookOf(behind.matPtr);
                        if (bl && bl->showThroughHair) { rec = behind; m = behind.matPtr; look = bl; break; }
                        if (!(bl && bl->hair)) break;   // only through hair
                        const Real left = probe.tMax - glm::length(behind.p - probe.o);
                        if (left <= Real(0)) break;
                        probe = behind.SpawnRay(r.d);
                        probe.tMax = left;
                    }
                }

                m->perturbShading(rec, wo);
                Vector3 n = glm::normalize(rec.n);
                if (look && look->faceNormals) n = glm::normalize((rec.p - look->faceCenter) / (look->faceRadii * look->faceRadii));
                else if (glm::dot(n, wo) < Real(0)) n = -n;

                // An overlay layer: alpha blended over what lies behind it
                if (look && look->overlay) {
                    Vector3 keyDir, alb;
                    Real a = 1;
                    const Real share = lightShare(rec, n, wo, lambda, keyDir);
                    const Vector3 c = shadeLook(*look, rec, n, share, keyDir, wo, alb, a);
                    a = std::clamp(a, Real(0), Real(1));
                    add += weight * a * c;
                    weight *= (Real(1) - a);
                    r = rec.SpawnRay(r.d);
                    continue;
                }

                const Vector3 albedo = glm::clamp(m->aovAlbedo(rec), Vector3(0.0), Vector3(1.0));
                if (!recorded) {
                    recorded = true;
                    sub.p = rec.p;
                    sub.n = n;
                    sub.depth = glm::dot(rec.p - m_camera->origin(), axis);
                    sub.footprint = subOffset + subAngle * glm::length(rec.p - eye);
                    sub.id = m;
                    sub.albedo = albedo;
                    sub.look = look;
                    sub.lineColor = albedo * Real(0.15) + Vector3(0.02);
                    sub.lineWidth = set.lineWidth;
                    if (look) { sub.lineColor = ToLinear(look->edgeColor); sub.lineWidth = set.lineWidth * look->edgeSize * lookLineScale; }
                }

                if (!look && m->emittedLuminance() > Real(0)) return add + weight * albedo;   // a light: its hue, unshaded

                if (!look && m->isSpecular()) {
                    // A mirror: the reflected view, tinted by its color. Smooth glass: seen straight through (no refraction, as cel
                    // shading draws glass), tinted, with a little white on top.
                    bool glass = false;
                    std::optional<BSDFSample> mirror;
                    for (int k = 0; k < 16 && !glass; ++k) {
                        // glass or mirror: asked at normal incidence (a ray inside glass may be totally reflected at wo)
                        const auto bs = m->sample(rec, (k < 8) ? n : wo, sampling::Random2D(), lambda);
                        if (!bs || bs->pdf <= Real(0)) continue;
                        if (bs->isTransmission()) glass = true;
                        else if (!mirror && k >= 8) mirror = bs;
                    }
                    if (glass) {
                        add += weight * Real(0.12);
                        weight *= albedo * Real(0.88);
                        r = rec.SpawnRay(-wo);
                        continue;
                    }
                    if (!mirror) return add;
                    weight *= albedo;
                    r = rec.SpawnRay(glm::normalize(mirror->wi));
                    continue;
                }

                Vector3 keyDir;
                const Real share = lightShare(rec, n, wo, lambda, keyDir);
                Vector3 col;
                Real t;
                if (look) {
                    Vector3 alb;
                    Real a = 1;
                    col = shadeLook(*look, rec, n, share, keyDir, wo, alb, a);
                    t = std::clamp(Real(3) * share, Real(0), Real(1));
                    if (sub.id == m) sub.albedo = alb;
                }
                else {
                    t = Smoothstep(set.threshold - set.softness, set.threshold + set.softness, share);
                    col = albedo * (set.shadowTone + (Vector3(1.0) - set.shadowTone) * t);
                }
                if (sub.id == m) sub.lit = t;
                return add + weight * glm::clamp(col, Vector3(0.0), Vector3(1.0));
            }
            return add;
        };

#pragma omp parallel for schedule(dynamic, 1)
        for (int y = 0; y < sh; ++y) {
            sampling::SeedStream(5, 0, static_cast<std::uint64_t>(y));
            for (int x = 0; x < sw; ++x) {
                Sub& sub = subs[static_cast<size_t>(y) * sw + x];
                const Ray r = m_camera->centerRay((Real(x) + Real(0.5)) / Real(sw), (Real(y) + Real(0.5)) / Real(sh));
                sub.color = shade(r, sub);
            }
        }

        // ---- outlines: mark the nearer side of every edge between neighbouring sub-samples ----
        // ---- rim: the lit sub-samples near a silhouette (where the neighbour lies behind, or is the backdrop) ----
        std::vector<unsigned char> edge(subs.size(), 0), sil(subs.size(), 0);
        const Real cosCrease = std::cos(set.creaseDeg * Real(3.14159265358979) / Real(180));
        auto test = [&](size_t a, size_t b) {
            const Sub& A = subs[a];
            const Sub& B = subs[b];
            if (!A.id && !B.id) return;
            const size_t nearer = (A.depth <= B.depth) ? a : b;
            const Sub& N = subs[nearer];
            bool step = false, other = false;
            if (!A.id || !B.id) step = true;
            else {
                const Real d = std::max(std::abs(glm::dot(B.p - A.p, A.n)), std::abs(glm::dot(A.p - B.p, B.n)));
                step = d > Real(3) * std::min(A.footprint, B.footprint);
                other = A.id != B.id || glm::dot(A.n, B.n) < cosCrease;
            }
            const bool line = N.look ? (step && N.look->edge) : (step || other);
            if (line) edge[nearer] = 1;
            if (step) sil[nearer] = 1;
        };
        for (int y = 0; y < sh; ++y)
            for (int x = 0; x < sw; ++x) {
                const size_t i = static_cast<size_t>(y) * sw + x;
                if (x + 1 < sw) test(i, i + 1);
                if (y + 1 < sh) test(i, i + sw);
            }
        // rim first (on the colors), then the lines over it
        if (set.rim > Real(0)) {
            std::vector<unsigned char> nearSil = sil;
            Dilate(nearSil, sw, sh, static_cast<int>(std::lround(Real(s) * Real(2.5))));   // about 2.5 pixels inward
            for (size_t k = 0; k < subs.size(); ++k)
                if (nearSil[k] && subs[k].lit > Real(0) && subs[k].id)
                    subs[k].color = glm::min(subs[k].color + set.rim * subs[k].lit * subs[k].albedo, Vector3(1.0));
        }
        // each marked sub-sample paints a disc of its own line width in its own line color
        std::vector<Vector3> lineColor(subs.size(), Vector3(-1.0));
        if (set.lineWidth > Real(0)) {
            int rMax = 0;
            for (size_t k = 0; k < subs.size(); ++k)
                if (edge[k]) rMax = std::max(rMax, static_cast<int>(std::ceil((Real(s) * subs[k].lineWidth - Real(1)) / Real(2))));
            rMax = std::min(rMax, 4 * s);
#pragma omp parallel for schedule(dynamic, 4)
            for (int y = 0; y < sh; ++y)
                for (int x = 0; x < sw; ++x) {
                    Real bestD = std::numeric_limits<Real>::infinity();
                    for (int dy = -rMax; dy <= rMax; ++dy)
                        for (int dx = -rMax; dx <= rMax; ++dx) {
                            const int xx = x + dx, yy = y + dy;
                            if (xx < 0 || xx >= sw || yy < 0 || yy >= sh) continue;
                            const size_t k = static_cast<size_t>(yy) * sw + xx;
                            if (!edge[k]) continue;
                            const Real r = std::max(Real(0), (Real(s) * subs[k].lineWidth - Real(1)) / Real(2));
                            const Real d2 = Real(dx * dx + dy * dy);
                            if (d2 <= r * r + Real(0.25) && d2 < bestD) { bestD = d2; lineColor[static_cast<size_t>(y) * sw + x] = subs[k].lineColor; }
                        }
                }
        }

        // ---- average the sub-samples to the pixels (the film's rows run from the top) ----
        size_t lineSubs = 0;
        for (const Vector3& c : lineColor) lineSubs += (c.x >= Real(0)) ? 1 : 0;
        for (int j = 0; j < height; ++j)
            for (int i = 0; i < width; ++i) {
                Vector3 c(0.0);
                for (int b = 0; b < s; ++b)
                    for (int a = 0; a < s; ++a) {
                        const size_t k = static_cast<size_t>(j * s + b) * sw + static_cast<size_t>(i * s + a);
                        c += (lineColor[k].x >= Real(0)) ? lineColor[k] : subs[k].color;
                    }
                film.setPixel(i, height - 1 - j, c / Real(s * s));
            }
        std::cout << "[Toon] outlines cover " << 100.0 * double(lineSubs) / double(std::max<size_t>(subs.size(), 1)) << " % of the picture; "
            << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s\n";
    }

} // namespace rayt
