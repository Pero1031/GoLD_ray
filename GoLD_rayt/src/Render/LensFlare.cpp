// Ghosts of a real lens (--ghosts); see Render/LensFlare.hpp.
#include "pch.h"
#include "Render/LensFlare.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <random>

#include "Color/SpectralConversion.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Lights/Light.hpp"
#include "Render/Camera.hpp"
#include "Scene/Scene.hpp"

namespace rayt {

    namespace {

        constexpr int kMaxBands = 32;
        constexpr int kPilot = 48;   // cells per side of the coarse grid that finds where a ghost's light gets through

        /// Parallel light from one light point: its direction (lens frame) and, per wavelength band, the RGB the film records per unit of
        /// the ghost's throughput (the irradiance across the beam through the band's colour weights, over the irradiance factor at the
        /// film centre, so that a pixel value is in the film's radiance unit).
        struct BeamSet {
            Vector3 dir;
            std::vector<Real> nm;
            std::vector<Vector3> rgb;
        };

        struct Setup {
            const LensSystem* lens = nullptr;
            std::vector<std::pair<int, int>> ghosts;   // the pairs of reflecting surfaces
            Real zTop = 0, radius = 0;                 // the plane in front of the lens the beams cross, and the radius of the grid on it
            Real filmW = 0, filmH = 0, pixelArea = 0, gc = 1;
            int width = 0, height = 0;
        };

        Setup MakeSetup(const Camera& camera, int width, int height, std::pair<int, int> only) {
            Setup su;
            su.lens = camera.lensSystem();
            su.width = width;
            su.height = height;
            if (!su.lens) return su;
            const auto& s = su.lens->surfaces();
            const int n = static_cast<int>(s.size());
            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j)
                    if (!s[static_cast<size_t>(i)].isStop() && !s[static_cast<size_t>(j)].isStop() && (only.first < 0 || (only.first == i && only.second == j)))
                        su.ghosts.emplace_back(i, j);
            su.zTop = s.front().z;
            for (const auto& x : s) su.zTop = std::max(su.zTop, x.z);
            // the grid covers what can enter: the first lens surface's clear aperture, and a stop in front of it
            for (const auto& x : s) {
                su.radius = std::max(su.radius, x.apertureRadius);
                if (!x.isStop()) break;
            }
            su.radius *= Real(1.02);
            su.filmW = su.lens->filmWidth();
            su.filmH = su.lens->filmHeight();
            su.pixelArea = su.filmW / Real(width) * (su.filmH / Real(height));
            su.gc = su.lens->centreIrradiance();
            return su;
        }

        /// Calls f for the light of every light point: every light that is not an environment; an area light as lightSamples points
        /// (stratified), merged into one when they lie within 0.25 deg of each other as seen from the lens (a small lamp). The spectral
        /// build gives f the bands one by one (the glass bends each its own way), the RGB build all at once. Returns the number of
        /// light points that reach the lens.
        int ForEachBeam(const Scene& scene, const Camera& camera, const Setup& su, const LensFlareSettings& st, const std::function<void(const BeamSet&)>& f) {
            const int bands = std::clamp(st.wavelengths, 1, kMaxBands);
            const Real lo = constants::LAMBDA_MIN, hi = constants::LAMBDA_MAX, bandWidth = (hi - lo) / Real(bands);
            std::vector<Real> bandNm(static_cast<size_t>(bands));
            for (int k = 0; k < bands; ++k) bandNm[static_cast<size_t>(k)] = lo + (Real(k) + Real(0.5)) * bandWidth;
            // each band's colour weights, integrated over the band (a smooth spectrum keeps its colour with few bands)
            std::vector<Vector3> bandWeight(static_cast<size_t>(bands), Vector3(0.0));
            constexpr int kSub = 16;
#if defined(RAYT_SPECTRAL_MODE)
            for (int k = 0; k < bands; ++k)
                for (int q = 0; q < kSub; ++q) {
                    // the film's weight at one wavelength, as SpectrumToRGB applies it (the CIE observer or a camera, the white balance)
                    SampledWavelengths one;
                    one.lambda.fill(lo + (Real(k) + (Real(q) + Real(0.5)) / Real(kSub)) * bandWidth);
                    one.pdf.fill(Real(0));
                    one.pdf[0] = Real(1) / Real(kSpectralSamples);
                    Spectrum unit(Real(0));
                    unit[0] = Real(1);
                    bandWeight[static_cast<size_t>(k)] += SpectrumToRGB(unit, one) * (bandWidth / Real(kSub));
                }
#else
            // RGB build: each band's share of each channel, the film sensor under D65 (as Materials/ThinFilm colours a film)
            Vector3 total(0.0);
            for (int k = 0; k < bands; ++k)
                for (int q = 0; q < kSub; ++q) {
                    const Real nm = lo + (Real(k) + (Real(q) + Real(0.5)) / Real(kSub)) * bandWidth;
                    const Vector3 w = color::rgbSensorWeights(nm) * color::evaluateD65(nm);
                    bandWeight[static_cast<size_t>(k)] += w;
                    total += w;
                }
            for (Vector3& w : bandWeight) w = Vector3(w.x / total.x, w.y / total.y, w.z / total.z);
#endif
            LightSampleContext ctx;
            ctx.p = camera.lensPointToWorld(Point3(0, 0, su.zTop));
            ctx.gn = Vector3(0.0);
            int points = 0;
            for (const auto& light : scene.allLights()) {
                const LightType type = light->type();
                if (type == LightType::Infinite) continue;
                const bool delta = type == LightType::DeltaPosition || type == LightType::DeltaDirection;
                const int m = delta ? 1 : std::max(1, st.lightSamples);
                const int side = static_cast<int>(std::ceil(std::sqrt(double(m))));
                std::vector<BeamSet> sets;
                for (int sIdx = 0; sIdx < m; ++sIdx) {
                    const Point2 u((Real(sIdx % side) + Real(0.5)) / Real(side), (Real(sIdx / side) + Real(0.5)) / Real(side));
                    BeamSet set;
                    set.nm = bandNm;
                    set.rgb.assign(static_cast<size_t>(bands), Vector3(0.0));
                    bool reaches = true;
                    for (int k0 = 0; k0 < bands && reaches; k0 += kSpectralSamples) {
                        SampledWavelengths lambda;
                        for (int q = 0; q < kSpectralSamples; ++q) {
                            lambda.lambda[static_cast<size_t>(q)] = bandNm[static_cast<size_t>(std::min(k0 + q, bands - 1))];
                            lambda.pdf[static_cast<size_t>(q)] = Real(1) / (hi - lo);
                        }
                        const auto ls = light->sampleLi(ctx, u, lambda);
                        if (!ls || !(ls->pdf > Real(0))) { reaches = false; break; }
                        Ray shadow(ctx.p, ls->wi);
                        if (ls->tMax > Real(0)) shadow.tMax = ls->tMax * (Real(1) - Real(1e-4));
                        const Vector3 dir = glm::normalize(camera.worldDirToLens(-ls->wi));
                        if (scene.occluded(shadow) || !(dir.z < Real(-1e-6))) { reaches = false; break; }   // blocked, or behind the lens
                        set.dir = dir;
                        const Spectrum e = ls->Li / (ls->pdf * Real(m));   // the irradiance across the beam (its share of the light)
#if defined(RAYT_SPECTRAL_MODE)
                        for (int q = 0; q < kSpectralSamples && k0 + q < bands; ++q)
                            set.rgb[static_cast<size_t>(k0 + q)] = bandWeight[static_cast<size_t>(k0 + q)] * (e[q] / su.gc);
#else
                        const Vector3 rgbE = SpectrumToRGB(e, lambda) / su.gc;
                        for (int k = 0; k < bands; ++k) set.rgb[static_cast<size_t>(k)] = rgbE * bandWeight[static_cast<size_t>(k)];
                        break;   // the RGB build's light has no spectrum: one call covers every band
#endif
                    }
                    if (reaches) sets.push_back(std::move(set));
                }
                if (sets.empty()) continue;
                // a small light: its points as one beam
                Real spread = 0;
                for (const BeamSet& a : sets) spread = std::max(spread, std::acos(std::clamp(glm::dot(a.dir, sets.front().dir), Real(-1), Real(1))));
                if (sets.size() > 1 && spread < Real(0.25) * Real(3.14159265358979323846) / Real(180)) {
                    BeamSet merged = sets.front();
                    Vector3 dirSum = merged.dir;
                    for (size_t s = 1; s < sets.size(); ++s) {
                        dirSum += sets[s].dir;
                        for (int k = 0; k < bands; ++k) merged.rgb[static_cast<size_t>(k)] += sets[s].rgb[static_cast<size_t>(k)];
                    }
                    merged.dir = glm::normalize(dirSum);
                    sets.assign(1, merged);
                }
                for (const BeamSet& set : sets) {
                    ++points;
#if defined(RAYT_SPECTRAL_MODE)
                    for (int k = 0; k < bands; ++k) {
                        if (set.rgb[static_cast<size_t>(k)] == Vector3(0.0)) continue;
                        f(BeamSet{ set.dir, { set.nm[static_cast<size_t>(k)] }, { set.rgb[static_cast<size_t>(k)] } });
                    }
#else
                    f(set);
#endif
                }
            }
            return points;
        }

        /// The beam's plane: two unit vectors across it.
        void BeamBasis(const Vector3& d, Vector3& e1, Vector3& e2) {
            const Vector3 helper = std::fabs(d.x) < Real(0.9) ? Vector3(1, 0, 0) : Vector3(0, 1, 0);
            e1 = glm::normalize(glm::cross(d, helper));
            e2 = glm::cross(d, e1);
        }

        inline void AtomicAdd(double* acc, size_t pixel, const Vector3& v) {
            double* p = acc + 3 * pixel;
#pragma omp atomic
            p[0] += double(v.x);
#pragma omp atomic
            p[1] += double(v.y);
#pragma omp atomic
            p[2] += double(v.z);
        }

        /// Film point (lens frame, metres) -> picture position in pixels (x to the right, y down from the top row), as Camera::filmX/Y map it.
        inline Point2 ToPixel(const Setup& su, Real x, Real y) {
            return Point2((Real(0.5) - x / su.filmW) * Real(su.width), (Real(0.5) + y / su.filmH) * Real(su.height));
        }

        /// A ray of a ghost's grid: where it lands, what clips it, and the RGB it carries per unit of film irradiance factor.
        struct Vertex {
            LensSystem::GhostSample g;
            Vector3 w = Vector3(0.0);
            bool ok = false;
            bool clear = false;   // passes the stop on every crossing and every rim
        };

        /// Whether light at the interpolated point of a ghost passes the stop on every crossing and every rim.
        inline bool Passes(const Setup& su, const Vertex& a, const Vertex& b, const Vertex& c, Real wa, Real wb, Real wc) {
            if (wa * a.g.rim + wb * b.g.rim + wc * c.g.rim > Real(1)) return false;
            for (int p = 0; p < std::min(a.g.stopPasses, 3); ++p)
                if (!su.lens->insideStopNow(wa * a.g.stopX[p] + wb * b.g.stopX[p] + wc * c.g.stopX[p], wa * a.g.stopY[p] + wb * b.g.stopY[p] + wc * c.g.stopY[p]))
                    return false;
            return true;
        }

        /// One triangle of a beam's grid on the film: the light of `power` (m^2 of the beam per unit irradiance) spread over its area.
        void DrawTriangle(const Setup& su, const Vertex& a, const Vertex& b, const Vertex& c, Real power, double* acc) {
            if (a.g.rim > Real(1) && b.g.rim > Real(1) && c.g.rim > Real(1)) return;   // the whole cell misses a clear aperture
            const Point2 pa = ToPixel(su, a.g.x, a.g.y), pb = ToPixel(su, b.g.x, b.g.y), pc = ToPixel(su, c.g.x, c.g.y);
            const Real area2 = (pb.x - pa.x) * (pc.y - pa.y) - (pb.y - pa.y) * (pc.x - pa.x);   // twice the signed area (pixels)
            const Real areaPx = std::fabs(area2) * Real(0.5);
            if (!(areaPx > Real(1e-12))) return;
            // a cell torn apart (its rays clipped or beyond an asphere's domain) spans the picture: no image of anything
            if (areaPx > Real(4) * Real(su.width) * Real(su.height) && !(a.clear && b.clear && c.clear)) return;
            const Real third = Real(1) / Real(3);
            if (areaPx < Real(0.5)) {   // smaller than a pixel: its light into the pixel of its centre, the share of 16 points that pass
                Real share = 1;
                if (!(a.clear && b.clear && c.clear)) {
                    int pass = 0;
                    for (int i = 0; i < 4; ++i)
                        for (int j = 0; j < 4; ++j) {
                            Real s = (Real(i) + Real(0.5)) / Real(4), u = (Real(j) + Real(0.5)) / Real(4);
                            if (s + u > Real(1)) { s = Real(1) - s; u = Real(1) - u; }   // the square folded onto the triangle
                            if (Passes(su, a, b, c, Real(1) - s - u, s, u)) ++pass;
                        }
                    share = Real(pass) / Real(16);
                }
                if (share <= Real(0)) return;
                const Point2 m = (pa + pb + pc) * third;
                const int x = static_cast<int>(std::floor(m.x)), y = static_cast<int>(std::floor(m.y));
                if (x < 0 || y < 0 || x >= su.width || y >= su.height) return;
                AtomicAdd(acc, static_cast<size_t>(y) * su.width + x, (a.w + b.w + c.w) * (third * share * power / su.pixelArea));
                return;
            }
            const Real value = power / (areaPx * su.pixelArea);   // the film irradiance per unit of the vertices' weights
            const int x0 = std::max(0, static_cast<int>(std::floor(std::min({ pa.x, pb.x, pc.x }))));
            const int x1 = std::min(su.width - 1, static_cast<int>(std::floor(std::max({ pa.x, pb.x, pc.x }))));
            const int y0 = std::max(0, static_cast<int>(std::floor(std::min({ pa.y, pb.y, pc.y }))));
            const int y1 = std::min(su.height - 1, static_cast<int>(std::floor(std::max({ pa.y, pb.y, pc.y }))));
            const Real inv = Real(1) / area2;
            // inside the stop and the rims at all three corners: inside everywhere between them (the stop is convex, the rim linear), so one
            // point per pixel decides only whether the pixel belongs to this triangle; else 2 x 2 points test the edges of the light
            const bool clear = a.clear && b.clear && c.clear;
            const int subs = clear ? 1 : 4;
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x) {
                    Vector3 sum(0.0);
                    for (int sub = 0; sub < subs; ++sub) {
                        const Real qx = clear ? Real(x) + Real(0.5) : Real(x) + Real(0.25) + Real(0.5) * Real(sub & 1);
                        const Real qy = clear ? Real(y) + Real(0.5) : Real(y) + Real(0.25) + Real(0.5) * Real(sub >> 1);
                        const Real wa = ((pb.x - qx) * (pc.y - qy) - (pb.y - qy) * (pc.x - qx)) * inv;
                        const Real wb = ((pc.x - qx) * (pa.y - qy) - (pc.y - qy) * (pa.x - qx)) * inv;
                        const Real wc = Real(1) - wa - wb;
                        if (wa < Real(0) || wb < Real(0) || wc < Real(0)) continue;
                        if (!clear && !Passes(su, a, b, c, wa, wb, wc)) continue;
                        sum += a.w * wa + b.w * wb + c.w * wc;
                    }
                    if (sum != Vector3(0.0)) AtomicAdd(acc, static_cast<size_t>(y) * su.width + x, sum * (value / Real(subs)));
                }
        }

        /// One beam set through one ghost: a coarse grid over the whole entrance finds the box where light gets through (a stopped-down
        /// stop passes a small part of it), then the fine grid's rays over that box (geometry at the first band's wavelength) and its triangles.
        void DrawGhost(const Setup& su, const LensFlareSettings& st, const BeamSet& beam, const LensSystem::ReflectanceTable& table, int gi, int gj,
            std::vector<Vertex>& v, std::vector<Vertex>& pilot, double* acc, double& rays) {
            const int bands = static_cast<int>(beam.nm.size());
            Vector3 e1, e2;
            BeamBasis(beam.dir, e1, e2);
            const Point3 centre(0, 0, su.zTop + Real(1e-3));
            auto trace = [&](Real u, Real w, Vertex& x) {
                const Point3 p = centre + e1 * u + e2 * w - beam.dir * Real(0.05);
                Real bandT[kMaxBands];
                x.ok = su.lens->traceGhost(gi, gj, p, beam.dir, beam.nm[0], st.coating, x.g, &table, bandT);
                if (!x.ok) return;
                x.w = Vector3(0.0);
                for (int k = 0; k < bands; ++k) x.w += beam.rgb[static_cast<size_t>(k)] * bandT[k];
                x.clear = x.g.rim <= Real(1);
                for (int p2 = 0; p2 < std::min(x.g.stopPasses, 3) && x.clear; ++p2) x.clear = su.lens->insideStopNow(x.g.stopX[p2], x.g.stopY[p2]);
            };
            const Real r = su.radius;
            const int gp = kPilot, np = gp + 1;
            const Real cp = Real(2) * r / Real(gp);
#pragma omp parallel for schedule(static)
            for (int a = 0; a < np; ++a)
                for (int b = 0; b < np; ++b) trace(-r + Real(a) * cp, -r + Real(b) * cp, pilot[static_cast<size_t>(a) * np + b]);
            rays += double(pilot.size());
            int a0 = np, a1 = -1, b0 = np, b1 = -1;
            for (int a = 0; a < np; ++a)
                for (int b = 0; b < np; ++b) {
                    const Vertex& x = pilot[static_cast<size_t>(a) * np + b];
                    if (!x.ok || !x.clear) continue;
                    a0 = std::min(a0, a); a1 = std::max(a1, a); b0 = std::min(b0, b); b1 = std::max(b1, b);
                }
            if (a1 < 0) return;   // no light of this ghost gets through
            a0 = std::max(0, a0 - 1); a1 = std::min(gp, a1 + 1); b0 = std::max(0, b0 - 1); b1 = std::min(gp, b1 + 1);
            const Real u0 = -r + Real(a0) * cp, w0 = -r + Real(b0) * cp;
            const int g = std::max(2, st.grid), n = g + 1;
            const Real du = Real(a1 - a0) * cp / Real(g), dw = Real(b1 - b0) * cp / Real(g);
#pragma omp parallel for schedule(static)
            for (int a = 0; a < n; ++a)
                for (int b = 0; b < n; ++b) trace(u0 + Real(a) * du, w0 + Real(b) * dw, v[static_cast<size_t>(a) * n + b]);
            rays += double(v.size());
            const Real power = Real(0.5) * du * dw;   // each triangle is half a cell of the beam
#pragma omp parallel for schedule(dynamic, 4)
            for (int a = 0; a < g; ++a)
                for (int b = 0; b < g; ++b) {
                    const size_t k00 = static_cast<size_t>(a) * n + b, k10 = k00 + n, k11 = k10 + 1, k01 = k00 + 1;
                    if (v[k00].ok && v[k10].ok && v[k11].ok) DrawTriangle(su, v[k00], v[k10], v[k11], power, acc);
                    if (v[k00].ok && v[k11].ok && v[k01].ok) DrawTriangle(su, v[k00], v[k11], v[k01], power, acc);
                }
        }

        std::vector<Vector3> ToPicture(const std::vector<double>& acc, int width, int height, Real* peak) {
            std::vector<Vector3> out(static_cast<size_t>(width) * height, Vector3(0.0));
            Real top = 0;
            for (size_t i = 0; i < out.size(); ++i) {
                out[i] = Vector3(Real(acc[3 * i]), Real(acc[3 * i + 1]), Real(acc[3 * i + 2]));
                top = std::max(top, std::max(out[i].x, std::max(out[i].y, out[i].z)));
            }
            if (peak) *peak = top;
            return out;
        }

    } // namespace

    std::vector<Vector3> RenderLensFlare(const Scene& scene, const Camera& camera, int width, int height, const LensFlareSettings& settings,
        LensFlareStats* stats) {
        const auto t0 = std::chrono::steady_clock::now();
        const Setup su = MakeSetup(camera, width, height, settings.only);
        std::vector<double> acc(static_cast<size_t>(width) * height * 3, 0.0);
        LensFlareStats st;
        if (su.lens) {
            const int n = std::max(2, settings.grid) + 1;
            std::vector<Vertex> v(static_cast<size_t>(n) * n), pilot(static_cast<size_t>(kPilot + 1) * (kPilot + 1));
            st.ghosts = static_cast<int>(su.ghosts.size());
            st.lights = ForEachBeam(scene, camera, su, settings, [&](const BeamSet& beam) {
                const LensSystem::ReflectanceTable table = su.lens->makeReflectanceTable(static_cast<int>(beam.nm.size()), beam.nm.data(), settings.coating);
                for (const auto& [i, j] : su.ghosts) DrawGhost(su, settings, beam, table, i, j, v, pilot, acc.data(), st.rays);
            });
        }
        std::vector<Vector3> out = ToPicture(acc, width, height, &st.peak);
        st.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (stats) *stats = st;
        return out;
    }

    std::vector<Vector3> RenderLensFlareBySampling(const Scene& scene, const Camera& camera, int width, int height, const LensFlareSettings& settings,
        int rays) {
        const Setup su = MakeSetup(camera, width, height, settings.only);
        std::vector<double> acc(static_cast<size_t>(width) * height * 3, 0.0);
        if (!su.lens || rays <= 0) return ToPicture(acc, width, height, nullptr);
        std::mt19937_64 rng(12345);
        std::uniform_real_distribution<double> uni(-1.0, 1.0);
        ForEachBeam(scene, camera, su, settings, [&](const BeamSet& beam) {
            Vector3 e1, e2;
            BeamBasis(beam.dir, e1, e2);
            const int bands = static_cast<int>(beam.nm.size());
            const int side = std::max(1, static_cast<int>(std::lround(std::sqrt(double(rays)))));   // jittered strata over the square
            const Real r = su.radius, power = Real(4) * r * r / (Real(side) * Real(side));
            const LensSystem::ReflectanceTable table = su.lens->makeReflectanceTable(bands, beam.nm.data(), settings.coating);
            const Point3 centre(0, 0, su.zTop + Real(1e-3));
            for (const auto& [i, j] : su.ghosts)
                for (int k = 0; k < side * side; ++k) {
                    const Real su1 = (Real(k % side) + Real(0.5) + Real(0.5) * Real(uni(rng))) / Real(side) * Real(2) - Real(1);
                    const Real su2 = (Real(k / side) + Real(0.5) + Real(0.5) * Real(uni(rng))) / Real(side) * Real(2) - Real(1);
                    const Point3 p = centre + e1 * (r * su1) + e2 * (r * su2) - beam.dir * Real(0.05);
                    LensSystem::GhostSample g;
                    Real bandT[kMaxBands];
                    if (!su.lens->traceGhost(i, j, p, beam.dir, beam.nm[0], settings.coating, g, &table, bandT) || g.rim > Real(1)) continue;
                    bool inside = true;
                    for (int q = 0; q < std::min(g.stopPasses, 3); ++q) inside = inside && su.lens->insideStopNow(g.stopX[q], g.stopY[q]);
                    if (!inside) continue;
                    const Point2 px = ToPixel(su, g.x, g.y);
                    const int x = static_cast<int>(std::floor(px.x)), y = static_cast<int>(std::floor(px.y));
                    if (x < 0 || y < 0 || x >= width || y >= height) continue;
                    Vector3 w(0.0);
                    for (int b = 0; b < bands; ++b) w += beam.rgb[static_cast<size_t>(b)] * bandT[b];
                    AtomicAdd(acc.data(), static_cast<size_t>(y) * width + x, w * (power / su.pixelArea));
                }
        });
        return ToPicture(acc, width, height, nullptr);
    }

} // namespace rayt
