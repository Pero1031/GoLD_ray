/**
 * @file Render/VolumeBPT.hpp
 * @brief PROTOTYPE: bidirectional path tracing in a gray homogeneous fog, with a cache of light sub-paths (the idea of LVC-BPT).
 *
 * Scope. One point light (isotropic, or a spot: a cone), a ball of gray homogeneous medium bounded by a null interface (the fog), a pinhole
 * camera, both outside the ball. No other surfaces. The medium is convex, so the transmittance between two points inside it is exp(-sigma_t d).
 *
 * The cache. Light sub-paths are traced first, N of them, and their scattering vertices are kept (positions only: every density below is a
 * function of positions). A camera sample then draws one cached light path at random and connects every vertex of its own path to every
 * vertex of that light path. Sharing light paths between camera samples keeps the estimate unbiased (each pair is still an independent draw
 * from the same distribution) and makes the light side cost almost nothing per camera sample.
 *
 * The strategies. A path with n medium vertices x_1..x_n (from the light) can be made in n + 1 ways: s of its vertices come from the light
 * sub-path and n - s from the camera sub-path.
 *   s = 0       next event estimation: the camera sub-path is connected to the light itself
 *   0 < s < n   a light vertex and a camera vertex are connected
 *   s = n       light tracing: the light sub-path is connected to the camera (a splat)
 * Every path is weighted with the balance heuristic over the strategies that can make it (the counts of samples enter: spp per pixel for the
 * strategies that start at the camera, N for light tracing), which gives the contribution  f(x) / sum_s n_s p_s(x)  for any strategy.
 *
 * Densities. A collision at x, reached from a point a with a direction density p_dir (per steradian), has the density
 *   p_dir * sigma_s * T(a, x) / |x - a|^2          (per volume; sigma_s = albedo * sigma_t: a collision happens with sigma_t T and is a scattering
 *                                                   with probability albedo)
 * where p_dir is the emission density at the light (1 / solid angle of the cone), the phase function at a scattering vertex, and 1 / Omega_pixel
 * (the pixel's solid angle) at the camera. The throughput of a path is then exactly f / p for the strategy that made it, and with the
 * balance heuristic the phase and transmittance factors that the strategies share cancel, so only ratios of these densities matter.
 *
 * Not handled: surfaces and glass (a specular vertex cannot be connected), colored or heterogeneous media, a camera or a light inside the fog,
 * a camera with a lens. Spectral mode: the medium is gray, so each sample carries its own wavelengths only for the color of the light.
 */

#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Core/Types.hpp"
#include "Core/Math.hpp"
#include "Core/Ray.hpp"
#include "Core/Sampling.hpp"
#include "Core/Constants.hpp"
#include "Core/Interaction.hpp"
#include "Core/SpectrumUtils.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Scene/Scene.hpp"
#include "Render/Camera.hpp"
#include "Render/Film.hpp"
#include "Render/Transmittance.hpp"
#include "Lights/PointLight.hpp"
#include "Materials/Material.hpp"
#include "Materials/Medium.hpp"
#include "Color/SpectralConversion.hpp"

namespace rayt {

    struct VolumeBptParams {
        int width = 0, height = 0;
        int spp = 16;                      ///< camera samples per pixel
        long long lightPaths = 1000000;    ///< N: light sub-paths in the cache
        int maxVertices = 6;               ///< scattering vertices per sub-path (a path has at most twice as many)
        bool useLightTracing = true;       ///< strategy s = n
        bool useConnections = true;        ///< strategies 0 < s < n
        bool useNee = true;                ///< strategy s = 0
    };

    class VolumeBPT {
    public:
        VolumeBPT(const Scene& scene, const Camera& camera, const PointLight& light, const Medium& medium, const VolumeBptParams& P)
            : m_scene(scene), m_cam(camera), m_light(light), m_med(medium), m_P(P)
        {
            const SampledWavelengths lam = makeLambda();
            m_sigT = m_med.sigmaTSpectrum(lam)[0];
            m_sigS = m_med.sigmaSSpectrum(lam)[0];
            m_albedo = m_sigS / m_sigT;
            m_lightPos = m_light.position();
            m_camPos = m_cam.origin();
        }

        /// Renders the image into the film (final pixel values, not sums).
        void render(Film& film) {
            const auto t0 = std::chrono::steady_clock::now();
            const int W = m_P.width, H = m_P.height;
            const size_t nPix = static_cast<size_t>(W) * static_cast<size_t>(H);
#ifdef _OPENMP
            const int nThreads = omp_get_max_threads();
#else
            const int nThreads = 1;
#endif

            // ---- 1. the cache of light sub-paths (and the light-traced image) ----
            std::vector<std::vector<Point3>> vPer(nThreads);
            std::vector<std::vector<uint32_t>> countPer(nThreads);
            std::vector<std::vector<float>> t0Per(nThreads);
            std::vector<std::vector<Vector3>> splatPer(nThreads, std::vector<Vector3>(nPix, Vector3(0.0)));
            const long long chunk = 4096;
            const int nChunks = static_cast<int>((m_P.lightPaths + chunk - 1) / chunk);
#pragma omp parallel for schedule(dynamic, 1)
            for (int c = 0; c < nChunks; ++c) {
#ifdef _OPENMP
                const int th = omp_get_thread_num();
#else
                const int th = 0;
#endif
                std::vector<Point3> Z;
                const long long first = static_cast<long long>(c) * chunk, last = std::min(m_P.lightPaths, first + chunk);
                for (long long i = first; i < last; ++i) {
                    Real tFirst = 0;
                    traceLightPath(Z, tFirst);
                    vPer[th].insert(vPer[th].end(), Z.begin(), Z.end());
                    countPer[th].push_back(static_cast<uint32_t>(Z.size()));
                    t0Per[th].push_back(static_cast<float>(std::exp(-m_sigT * tFirst)));
                    if (m_P.useLightTracing && !Z.empty()) lightTrace(Z, tFirst, splatPer[th]);
                }
            }
            // concatenate (order does not matter: every light path is an independent sample)
            m_verts.clear(); m_start.clear(); m_t0.clear();
            m_start.push_back(0);
            for (int th = 0; th < nThreads; ++th) {
                size_t at = 0;
                for (size_t k = 0; k < countPer[th].size(); ++k) {
                    for (uint32_t q = 0; q < countPer[th][k]; ++q) m_verts.push_back(vPer[th][at + q]);
                    at += countPer[th][k];
                    m_start.push_back(static_cast<uint32_t>(m_verts.size()));
                    m_t0.push_back(t0Per[th][k]);
                }
                std::vector<Point3>().swap(vPer[th]);
            }
            m_nStored = m_t0.size();
            const double tLight = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

            // ---- 2. camera sub-paths: connections and next event estimation ----
            std::vector<Vector3> img(nPix, Vector3(0.0));
#pragma omp parallel for schedule(dynamic, 1)
            for (int j = 0; j < H; ++j) {
                std::vector<Point3> Y, X;
                for (int i = 0; i < W; ++i) {
                    Vector3 sum(0.0);
                    for (int k = 0; k < m_P.spp; ++k) sum += cameraSample(i, j, Y, X);
                    img[static_cast<size_t>(H - 1 - j) * static_cast<size_t>(W) + static_cast<size_t>(i)] = sum;
                }
            }
            for (int th = 0; th < nThreads; ++th) for (size_t p = 0; p < nPix; ++p) img[p] += splatPer[th][p];
            for (int y = 0; y < H; ++y) for (int x = 0; x < W; ++x) film.setPixel(x, y, img[static_cast<size_t>(y) * static_cast<size_t>(W) + static_cast<size_t>(x)]);

            const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::cout << "[VolumeBPT] " << m_P.lightPaths << " light paths (" << m_verts.size() << " vertices cached, " << tLight << " s), "
                << double(W) * H * m_P.spp << " camera samples, " << sec << " s in all\n";
        }

    private:
        const Scene& m_scene;
        const Camera& m_cam;
        const PointLight& m_light;
        const Medium& m_med;
        VolumeBptParams m_P;
        Real m_sigT = 1, m_sigS = 1, m_albedo = 1;
        Point3 m_lightPos, m_camPos;

        // the cache: path k has the vertices m_verts[m_start[k] .. m_start[k + 1]) and the transmittance m_t0[k] of its first edge (light -> vertex 1)
        std::vector<Point3> m_verts;
        std::vector<uint32_t> m_start;
        std::vector<float> m_t0;
        size_t m_nStored = 0;

        static SampledWavelengths makeLambda() {
#if defined(RAYT_SPECTRAL_MODE)
            return SampledWavelengths::sampleVisible(sampling::Random());
#else
            return SampledWavelengths();
#endif
        }

        /// Follows a ray from outside the fog through it (delta tracking, null interfaces crossed) and appends its scattering vertices.
        /// tFirst: the distance inside the medium to the first vertex (for the transmittance of the first edge).
        void walk(Ray r, std::vector<Point3>& out, Real& tFirst) const {
            out.clear();
            tFirst = 0;
            for (int crossings = 0; crossings < 64 && static_cast<int>(out.size()) < m_P.maxVertices; ++crossings) {
                SurfaceInteraction rec;
                const bool hit = m_scene.hit(r, rec);
                if (r.medium) {
                    const Real tS = hit ? rec.t : constants::INFINITY_VAL;
                    Real t = 0;
                    for (;;) {
                        t -= std::log(Real(1) - sampling::Random()) / m_sigT;
                        if (t >= tS) break;                                   // reaches the surface
                        if (sampling::Random() >= m_albedo) return;           // absorbed
                        const Point3 x = r.o + r.d * t;
                        if (out.empty()) tFirst = t;
                        out.push_back(x);
                        r = Ray(x, m_med.samplePhase(r.d, Point2(sampling::Random(), sampling::Random())), Real(0), &m_med);
                        goto next;
                    }
                }
                if (!hit) return;
                if (!rec.matPtr->isNullInterface()) return;                   // the prototype has no other surfaces
                {
                    const Medium* nm = rec.matPtr->isMediumBoundary() ? (rec.frontFace ? rec.matPtr->mediumInside() : rec.matPtr->mediumOutside()) : r.medium;
                    if (nm) r = Ray(rec.p, r.d, constants::RAY_EPSILON, nm);
                    else r = rec.SpawnRay(r.d);
                }
            next:;
            }
        }

        void traceLightPath(std::vector<Point3>& Z, Real& tFirst) const {
            const Vector3 dir = m_light.sampleEmissionDir(Point2(sampling::Random(), sampling::Random()));
            walk(Ray(m_lightPos, dir, constants::RAY_EPSILON, nullptr), Z, tFirst);
        }

        /// Transmittance (scalar) of the straight segment a -> b, a inside the fog or outside, through null interfaces.
        Real transmittanceBetween(const Point3& a, const Medium* mediumAtA, const Point3& b, const SampledWavelengths& lam) const {
            const Vector3 d = b - a;
            const Real dist = glm::length(d);
            if (dist <= Real(0)) return Real(0);
            const Spectrum T = shadowTransmittance(m_scene, a, mediumAtA, d / dist, dist - Real(2) * constants::RAY_EPSILON, lam);
            return T[0];
        }

        // ------------------------------------------------------------------------------------------------
        // The path: weight of one path by the balance heuristic, and its contribution
        // ------------------------------------------------------------------------------------------------

        /**
         * Contribution f(x) / D(x) of the path x_1..x_n (X[0..n-1], from the light to the camera), without the light's spectrum.
         * Te0 / Ten: transmittance of the first edge (light -> x_1) and the last edge (x_n -> camera); the edges between vertices are inside the fog.
         */
        Real pathContribution(const Point3* X, int n, Real Te0, Real Ten) const {
            constexpr int kMax = 40;
            if (n < 1 || n + 2 > kMax) return Real(0);
            Point3 P[kMax];
            P[0] = m_lightPos;
            for (int i = 0; i < n; ++i) P[i + 1] = X[i];
            P[n + 1] = m_camPos;

            Vector3 dir[kMax];
            Real d2[kMax], T[kMax];
            for (int i = 0; i <= n; ++i) {
                const Vector3 e = P[i + 1] - P[i];
                d2[i] = glm::dot(e, e);
                if (d2[i] <= Real(0)) return Real(0);
                dir[i] = e / std::sqrt(d2[i]);
            }
            T[0] = Te0;
            T[n] = Ten;
            for (int i = 1; i < n; ++i) T[i] = std::exp(-m_sigT * std::sqrt(d2[i]));
            if (!(T[0] > 0) || !(T[n] > 0)) return Real(0);

            // the pixel's density at x_n: 1 / Omega_pixel
            Real s, t, cosCam, dist;
            if (!m_cam.project(P[n], s, t, cosCam, dist)) return Real(0);
            const Real F = m_cam.focusDistance();
            const Real w = F * F * Real(m_P.width) * Real(m_P.height) / (m_cam.focalPlaneArea() * cosCam * cosCam * cosCam);

            if (!m_light.emitsToward(dir[0])) return Real(0);
            const Real pL = m_light.emissionPdf(dir[0]);

            // phase function at x_k (k = 1..n)
            Real ph[kMax];
            for (int k = 1; k <= n; ++k) ph[k] = m_med.phase(glm::dot(dir[k - 1], dir[k]));

            // f / I: the intensity I(dir_0) is a factor of the light's spectrum
            Real f = w;
            for (int i = 0; i <= n; ++i) f *= T[i] / d2[i];
            for (int k = 1; k <= n; ++k) f *= m_sigS * ph[k];

            // densities of the vertices when sampled from the light (F) and from the camera (C)
            Real pdfF[kMax], pdfC[kMax];
            for (int i = 1; i <= n; ++i) {
                const Real dirF = (i == 1) ? pL : ph[i - 1];
                pdfF[i] = dirF * m_sigS * T[i - 1] / d2[i - 1];
                const Real dirC = (i == n) ? w : ph[i + 1];
                pdfC[i] = dirC * m_sigS * T[i] / d2[i];
            }
            // p_s = prod_{i <= s} pdfF[i] * prod_{i > s} pdfC[i]
            Real suffix[kMax];
            suffix[n] = Real(1);
            for (int i = n - 1; i >= 0; --i) suffix[i] = suffix[i + 1] * pdfC[i + 1];
            Real D = 0, prefix = Real(1);
            const int M = m_P.maxVertices;
            for (int sIdx = 0; sIdx <= n; ++sIdx) {
                if (sIdx > 0) prefix *= pdfF[sIdx];
                if (sIdx > M || n - sIdx > M) continue;                         // this strategy cannot make a path of this length
                if (sIdx == 0 && !m_P.useNee) continue;
                if (sIdx == n && !m_P.useLightTracing) continue;
                if (sIdx > 0 && sIdx < n && !m_P.useConnections) continue;
                const Real count = (sIdx == n) ? Real(m_P.lightPaths) : Real(m_P.spp);
                D += count * prefix * suffix[sIdx];
            }
            if (!(D > Real(0)) || !std::isfinite(D)) return Real(0);
            const Real c = f / D;
            return std::isfinite(c) ? c : Real(0);
        }

        // ------------------------------------------------------------------------------------------------
        // Light tracing (s = n): every vertex of a light path is connected to the camera
        // ------------------------------------------------------------------------------------------------
        void lightTrace(const std::vector<Point3>& Z, Real tFirst, std::vector<Vector3>& splat) const {
            const SampledWavelengths lam = makeLambda();
            const Spectrum I = m_light.intensity(lam);
            const Real Te0 = std::exp(-m_sigT * tFirst);
            for (size_t s = 1; s <= Z.size(); ++s) {
                const Point3& x = Z[s - 1];
                Real sx, sy, cosCam, dist;
                if (!m_cam.project(x, sx, sy, cosCam, dist)) continue;
                const Real Ten = transmittanceBetween(x, &m_med, m_camPos, lam);
                if (!(Ten > 0)) continue;
                const Real c = pathContribution(Z.data(), static_cast<int>(s), Te0, Ten);
                if (c <= Real(0)) continue;
                const Vector3 rgb = SpectrumToRGB(I * c, lam);
                if (math::hasNonFinite(rgb)) continue;
                const int px = std::min(m_P.width - 1, static_cast<int>(sx * Real(m_P.width)));
                const int py = m_P.height - 1 - std::min(m_P.height - 1, static_cast<int>(sy * Real(m_P.height)));
                splat[static_cast<size_t>(py) * static_cast<size_t>(m_P.width) + static_cast<size_t>(px)] += rgb;
            }
        }

        // ------------------------------------------------------------------------------------------------
        // One camera sample: next event estimation and connections to one cached light path
        // ------------------------------------------------------------------------------------------------
        Vector3 cameraSample(int i, int j, std::vector<Point3>& Y, std::vector<Point3>& X) const {
            const Real u = (Real(i) + sampling::Random()) / Real(m_P.width);
            const Real v = (Real(j) + sampling::Random()) / Real(m_P.height);
            const Ray ray = m_cam.getRay(u, v, Point2(0.5, 0.5));
            Real tFirstC = 0;
            walk(ray, Y, tFirstC);
            if (Y.empty()) return Vector3(0.0);
            const Real Tec = std::exp(-m_sigT * tFirstC);   // the first edge seen from the camera: camera -> y_1 = x_n

            const SampledWavelengths lam = makeLambda();
            const Spectrum I = m_light.intensity(lam);
            Spectrum acc(0.0);

            // the cached light path of this sample
            size_t lb = 0, le = 0;
            Real Tl0 = 0;
            if (m_P.useConnections && m_nStored > 0) {
                const size_t k = std::min(m_nStored - 1, static_cast<size_t>(sampling::Random() * Real(m_nStored)));
                lb = m_start[k]; le = m_start[k + 1];
                Tl0 = m_t0[k];
            }

            const int T = static_cast<int>(Y.size());
            for (int t = 1; t <= T; ++t) {
                // ---- s = 0: next event estimation from y_t ----
                if (m_P.useNee) {
                    const Point3& y = Y[t - 1];
                    const Vector3 toL = m_lightPos - y;
                    const Real dl = glm::length(toL);
                    if (dl > Real(0) && m_light.emitsToward(-toL / dl)) {
                        const Real Te0 = transmittanceBetween(y, &m_med, m_lightPos, lam);
                        if (Te0 > Real(0)) {
                            X.clear();
                            for (int q = t - 1; q >= 0; --q) X.push_back(Y[q]);
                            acc += I * pathContribution(X.data(), t, Te0, Tec);
                        }
                    }
                }
                // ---- 0 < s < n: a vertex of the cached light path connected to y_t ----
                for (size_t q = lb; q < le; ++q) {
                    const int s = static_cast<int>(q - lb) + 1;
                    X.clear();
                    for (size_t a = lb; a <= q; ++a) X.push_back(m_verts[a]);
                    for (int a = t - 1; a >= 0; --a) X.push_back(Y[a]);
                    acc += I * pathContribution(X.data(), s + t, Tl0, Tec);
                }
            }
            const Vector3 rgb = SpectrumToRGB(acc, lam);
            return math::hasNonFinite(rgb) ? Vector3(0.0) : rgb;
        }
    };

} // namespace rayt
