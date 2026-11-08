/**
 * @file Render/VolumeBPT.hpp
 * @brief PROTOTYPE: bidirectional path tracing in a gray homogeneous fog and on gray diffuse surfaces, with a cache of light sub-paths
 *        (the idea of LVC-BPT).
 *
 * Scope. One point light (isotropic, or a spot: a cone), a ball of gray homogeneous medium bounded by a null interface (the fog), a pinhole
 * camera, both outside the ball, and any number of other surfaces, all Lambertian with the same gray albedo (VolumeBptParams::surfaceAlbedo).
 * Surfaces may lie inside the fog, outside it, or cross its boundary; they occlude. No glass: a specular vertex cannot be connected.
 *
 * The cache. Light sub-paths are traced first, N of them, and their vertices (scattering points in the fog and points on surfaces) are kept.
 * A camera sample then draws one cached light path at random and connects every vertex of its own path to every vertex of that light path.
 * Sharing light paths between camera samples keeps the estimate unbiased (each pair is still an independent draw from the same distribution)
 * and makes the light side cost almost nothing per camera sample.
 *
 * The strategies. A path with n vertices x_1..x_n (from the light) can be made in n + 1 ways: s of its vertices come from the light
 * sub-path and n - s from the camera sub-path.
 *   s = 0       next event estimation: the camera sub-path is connected to the light itself
 *   0 < s < n   a light vertex and a camera vertex are connected
 *   s = n       light tracing: the light sub-path is connected to the camera (a splat)
 * Every path is weighted with the balance heuristic over the strategies that can make it (the numbers of samples enter: spp per pixel for the
 * strategies that start at the camera, N for light tracing), which gives the contribution  f(x) / sum_s n_s p_s(x)  for any strategy.
 * (A path of n vertices is made by a strategy only if both sub-paths are short enough: s, n - s <= maxVertices.)
 *
 * Densities (per unit of the measure of the vertex: volume in the fog, area on a surface). A vertex x reached from a point a through a
 * direction with density p_dir (per steradian):
 *   fog vertex      p_dir * sigma_s T(a, x) / |x - a|^2                    (sigma_s = albedo sigma_t: a collision is a scattering with probability albedo)
 *   surface vertex  p_dir * rho T(a, x) |cos_x| / |x - a|^2               (rho: the probability that the ray is not absorbed)
 * where T is the transmittance of the segment and p_dir is the emission density at the light (1 / solid angle of the cone), the phase function
 * at a fog vertex, |cos| / pi at a surface vertex, and 1 / Omega_pixel (the pixel's solid angle) at the camera. The throughput of a path is
 * f / p for the strategy that made it. f is the product of the vertex factors (sigma_s phase in the fog, rho / pi on a surface; a surface
 * reflects: its two neighbors lie on the side the surface faces) and of  T |cos_a| |cos_b| / d^2  of every edge (cos = 1 at a fog vertex).
 * The transmittance of every edge is known: the edges of a sub-path are recorded when it is traced, an edge that connects two sub-paths (or a
 * sub-path and the light or the camera) is measured with a shadow ray, which is blocked by surfaces and crosses null interfaces.
 *
 * Glass (smooth dielectrics). A specular vertex is a vertex of a sub-path, but it cannot be the end of a connection (the directions would not
 * match: the probability is zero), so a chain of specular vertices behaves like one long edge. Consequences:
 *   - a strategy that connects at a specular vertex (or lets the light or the camera be connected to one) cannot make the path and is left out
 *     of the sum D;
 *   - the throughput of a specular vertex is the weight W that Material::sample returns (1 for a reflection, 1/eta^2 for a refraction; the
 *     eta^2 factors of a chain that enters and leaves a body cancel, so a chain has the same weight whichever end it was traced from). In the
 *     contribution it enters with the cosine that the sub-path that sampled the vertex used, so that f / p of the strategy that made the path
 *     is exactly the product of the sampled weights;
 *   - in the densities a specular vertex scatters with density 1 (the Dirac factors are left out) and keeps its geometric conversion cos / d^2.
 *     The choice between reflection and refraction has the same probability from both ends and is left out as well.
 * The MIS weights n_s p_s / D sum to one over the strategies of a path whatever the p_s are, so an inexact convention for the Dirac factors
 * could only make the estimate noisier, never biased; the checks against independent estimators show whether it is exact enough.
 *
 * Dispersion. The shape of a path through a dispersive glass depends on the wavelength, so two sub-paths can be connected only if they have
 * the same wavelength. Instead of keeping track of wavelengths per sub-path, VolumeBptParams::passes = K renders the picture K times with one
 * wavelength per pass (stratified over the visible spectrum): inside a pass every sub-path, light and camera alike, has that wavelength, so
 * the wavelengths of connected sub-paths always agree. The density of the wavelength is the same factor in every strategy and cancels in the
 * MIS weights.
 *
 * Not handled: rough dielectrics, colored or heterogeneous media, a camera or a light inside the fog, a camera with a lens, a light that
 * sees itself (the lamp is not a surface). Spectral mode: the medium, the surfaces and the glass are gray, so each sample carries its own
 * wavelengths only for the color of the light.
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
#include "Geometry/Frame.hpp"
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
        int maxVertices = 6;               ///< vertices per sub-path (a path has at most twice as many)
        Real surfaceAlbedo = Real(0.5);    ///< albedo of every (non-null) surface
        bool useLightTracing = true;       ///< strategy s = n
        bool useConnections = true;        ///< strategies 0 < s < n
        bool useNee = true;                ///< strategy s = 0
        /// Test switch: light tracing counts only the paths whose first vertex from the light is a specular one (L S ...), which is what the caustic
        /// photon map with traceMedia counts, so that the two can be compared.
        bool lightTraceOnlyAfterGlass = false;
        /// 1: every sub-path draws its own wavelengths. K > 1 (needed with dispersive glass): the picture is rendered K times, each time with ONE
        /// wavelength shared by every sub-path (see VolumeBPT::render); spp and lightPaths are the totals over all passes.
        int passes = 1;
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
            m_rho = m_P.surfaceAlbedo;
            m_P.maxVertices = std::min(m_P.maxVertices, 16);
            m_spp = m_P.spp;
            m_N = m_P.lightPaths;
            m_lightPos = m_light.position();
            m_camPos = m_cam.origin();
        }

        /// Renders the image into the film (final pixel values, not sums).
        void render(Film& film) {
            const auto t0 = std::chrono::steady_clock::now();
            const int W = m_P.width, H = m_P.height;
            const size_t nPix = static_cast<size_t>(W) * static_cast<size_t>(H);
            const int K = std::max(1, m_P.passes);
            m_spp = std::max(1, m_P.spp / K);                 // per pass
            m_N = std::max<long long>(1, m_P.lightPaths / K);
            std::vector<Vector3> total(nPix, Vector3(0.0));
            double tLight = 0;
            size_t cached = 0;
            for (int k = 0; k < K; ++k) {
                m_passLambdaOn = false;
#if defined(RAYT_SPECTRAL_MODE)
                if (K > 1) {   // one wavelength for the whole pass, stratified over the passes
                    m_passLambda = SampledWavelengths::sampleVisible((Real(k) + sampling::Random()) / Real(K));
                    m_passLambda.terminateSecondary();
                    m_passLambdaOn = true;
                }
#endif
                std::vector<Vector3> img;
                tLight += renderPass(img);
                cached = std::max(cached, m_verts.size());
                for (size_t p = 0; p < nPix; ++p) total[p] += img[p];
            }
            m_passLambdaOn = false;
            for (int y = 0; y < H; ++y) for (int x = 0; x < W; ++x) film.setPixel(x, y, total[static_cast<size_t>(y) * static_cast<size_t>(W) + static_cast<size_t>(x)] / Real(K));

            const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::cout << "[VolumeBPT] " << m_N * K << " light paths (up to " << cached << " vertices cached per pass, " << tLight << " s), "
                << double(W) * H * m_spp * K << " camera samples, " << K << " pass" << (K > 1 ? "es" : "") << ", " << sec << " s in all\n";
        }

    private:
        /// A vertex of a sub-path (compact: the cache holds millions).
        struct Vtx {
            float p[3];
            float n[3];        ///< surface: the normal on the side of the ray that arrived (unit); fog: zero
            float rho;         ///< surface: albedo; fog: 0
            float W;           ///< specular vertex: the weight of the sampled scattering (Material::sample); else 0
            float tIn;         ///< transmittance of the edge to the previous vertex of the sub-path (or from the light / the camera)
            uint8_t surface;   ///< 1: a point on a surface (diffuse or specular)
            uint8_t spec;      ///< 1: a specular vertex (smooth glass): it cannot be connected
            uint8_t inFog;     ///< 1: the rays at this vertex travel in the fog
        };
        static Point3 pos(const Vtx& v) { return Point3(v.p[0], v.p[1], v.p[2]); }
        static Vector3 nrm(const Vtx& v) { return Vector3(v.n[0], v.n[1], v.n[2]); }
        const Medium* mediumOf(const Vtx& v) const { return v.inFog ? &m_med : nullptr; }

        const Scene& m_scene;
        const Camera& m_cam;
        const PointLight& m_light;
        const Medium& m_med;
        VolumeBptParams m_P;
        Real m_sigT = 1, m_sigS = 1, m_albedo = 1, m_rho = 0.5;
        Point3 m_lightPos, m_camPos;

        // the cache: path k has the vertices m_verts[m_start[k] .. m_start[k + 1])
        std::vector<Vtx> m_verts;
        std::vector<uint32_t> m_start;
        size_t m_nStored = 0;

        int m_spp = 1;                 // camera samples per pixel in the current pass
        long long m_N = 1;             // light paths in the current pass
        SampledWavelengths m_passLambda;
        bool m_passLambdaOn = false;

        /// One complete estimate of the picture with the current wavelengths: the cache of light sub-paths (and the light-traced image),
        /// then the camera sub-paths with their connections. Returns the time of the light stage.
        double renderPass(std::vector<Vector3>& img) {
            const auto t0 = std::chrono::steady_clock::now();
            const int W = m_P.width, H = m_P.height;
            const size_t nPix = static_cast<size_t>(W) * static_cast<size_t>(H);
#ifdef _OPENMP
            const int nThreads = omp_get_max_threads();
#else
            const int nThreads = 1;
#endif

            // ---- 1. the cache of light sub-paths (and the light-traced image) ----
            std::vector<std::vector<Vtx>> vPer(nThreads);
            std::vector<std::vector<uint8_t>> countPer(nThreads);
            std::vector<std::vector<Vector3>> splatPer(nThreads, std::vector<Vector3>(nPix, Vector3(0.0)));
            const long long chunk = 4096;
            const int nChunks = static_cast<int>((m_N + chunk - 1) / chunk);
#pragma omp parallel for schedule(dynamic, 1)
            for (int c = 0; c < nChunks; ++c) {
#ifdef _OPENMP
                const int th = omp_get_thread_num();
#else
                const int th = 0;
#endif
                std::vector<Vtx> Z;
                const long long first = static_cast<long long>(c) * chunk, last = std::min(m_N, first + chunk);
                for (long long i = first; i < last; ++i) {
                    SampledWavelengths lam = makeLambda();
                    traceLightPath(Z, lam);
                    if (m_P.useConnections) {   // the cache is only needed by the connections
                        vPer[th].insert(vPer[th].end(), Z.begin(), Z.end());
                        countPer[th].push_back(static_cast<uint8_t>(Z.size()));
                    }
                    if (m_P.useLightTracing && !Z.empty()) lightTrace(Z, lam, splatPer[th]);
                }
            }
            // concatenate (order does not matter: every light path is an independent sample)
            m_verts.clear(); m_start.clear();
            m_start.push_back(0);
            for (int th = 0; th < nThreads; ++th) {
                size_t at = 0;
                for (size_t k = 0; k < countPer[th].size(); ++k) {
                    m_verts.insert(m_verts.end(), vPer[th].begin() + at, vPer[th].begin() + at + countPer[th][k]);
                    at += countPer[th][k];
                    m_start.push_back(static_cast<uint32_t>(m_verts.size()));
                }
                std::vector<Vtx>().swap(vPer[th]);
            }
            m_nStored = m_start.size() - 1;
            const double tLight = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

            // ---- 2. camera sub-paths: connections and next event estimation ----
            img.assign(nPix, Vector3(0.0));
#pragma omp parallel for schedule(dynamic, 1)
            for (int j = 0; j < H; ++j) {
                std::vector<Vtx> Y, X;
                for (int i = 0; i < W; ++i) {
                    Vector3 sum(0.0);
                    for (int k = 0; k < m_spp; ++k) sum += cameraSample(i, j, Y, X);
                    img[static_cast<size_t>(H - 1 - j) * static_cast<size_t>(W) + static_cast<size_t>(i)] = sum;
                }
            }
            for (int th = 0; th < nThreads; ++th) for (size_t p = 0; p < nPix; ++p) img[p] += splatPer[th][p];
            return tLight;
        }

        SampledWavelengths makeLambda() const {
            if (m_passLambdaOn) return m_passLambda;   // one wavelength for the whole pass (dispersion)
#if defined(RAYT_SPECTRAL_MODE)
            return SampledWavelengths::sampleVisible(sampling::Random());
#else
            return SampledWavelengths();
#endif
        }

        /// Follows a ray from outside the fog (delta tracking in the fog, null interfaces crossed, diffuse bounces on surfaces) and appends its vertices.
        /// `lam`: the wavelengths of the sub-path (a dispersive material reduces them to one).
        void walk(Ray r, std::vector<Vtx>& out, SampledWavelengths& lam) const {
            out.clear();
            Real lenInMedium = 0;   // distance travelled in the medium since the last vertex
            for (int iter = 0; iter < 256 && static_cast<int>(out.size()) < m_P.maxVertices; ++iter) {
                SurfaceInteraction rec;
                const bool hit = m_scene.hit(r, rec);
                if (r.medium) {
                    if (!hit) return;                                          // a ray that never leaves the fog: a leaky boundary
                    Real t = 0;
                    bool collided = false;
                    for (;;) {
                        t -= std::log(Real(1) - sampling::Random()) / m_sigT;
                        if (t >= rec.t) break;                                 // reaches the surface
                        if (sampling::Random() >= m_albedo) return;            // absorbed
                        const Point3 x = r.o + r.d * t;
                        lenInMedium += t;
                        Vtx v = {};
                        v.p[0] = float(x.x); v.p[1] = float(x.y); v.p[2] = float(x.z);
                        v.tIn = float(std::exp(-m_sigT * lenInMedium));
                        v.inFog = 1;
                        out.push_back(v);
                        lenInMedium = 0;
                        r = Ray(x, m_med.samplePhase(r.d, Point2(sampling::Random(), sampling::Random())), Real(0), &m_med);
                        collided = true;
                        break;
                    }
                    if (collided) continue;
                    lenInMedium += rec.t;
                }
                if (!hit) return;
                rec.medium = r.medium;
                if (rec.matPtr->isNullInterface()) {                           // the fog's boundary: cross it
                    const Medium* nm = rec.matPtr->isMediumBoundary() ? (rec.frontFace ? rec.matPtr->mediumInside() : rec.matPtr->mediumOutside()) : r.medium;
                    rec.medium = nm;   // SpawnRay carries rec.medium into the new ray
                    if (nm) r = Ray(rec.p, r.d, constants::RAY_EPSILON, nm);
                    else r = rec.SpawnRay(r.d);
                    continue;
                }
                if (rec.matPtr->isSpecular()) {
                    // smooth glass: reflection or refraction (chosen by Fresnel), the weight W goes into the vertex
                    rec.matPtr->prepareWavelengths(lam);
                    const std::optional<BSDFSample> bs = rec.matPtr->sample(rec, -r.d, Point2(sampling::Random(), sampling::Random()), lam);
                    if (!bs || !bs->isSpecular()) return;
                    const Real W = bs->f[0];
                    if (!(W > Real(0)) || !std::isfinite(W)) return;
                    const Vector3 wi = bs->wi;
                    const Medium* before = r.medium;
                    if (rec.matPtr->isMediumBoundary()) {   // a body standing in the fog: the medium changes as the ray enters or leaves it
                        const bool wiOnRaySide = glm::dot(wi, rec.gn) > 0;
                        rec.medium = (wiOnRaySide == rec.frontFace) ? rec.matPtr->mediumOutside() : rec.matPtr->mediumInside();
                    }
                    Vtx v = {};
                    v.p[0] = float(rec.p.x); v.p[1] = float(rec.p.y); v.p[2] = float(rec.p.z);
                    v.n[0] = float(rec.gn.x); v.n[1] = float(rec.gn.y); v.n[2] = float(rec.gn.z);
                    v.W = float(W);
                    v.tIn = float(std::exp(-m_sigT * lenInMedium));
                    v.surface = 1;
                    v.spec = 1;
                    v.inFog = before ? 1 : 0;
                    out.push_back(v);
                    lenInMedium = 0;
                    if (rec.medium) r = Ray(rec.p, wi, constants::RAY_EPSILON, rec.medium);
                    else r = rec.SpawnRay(wi);
                    continue;
                }
                // a diffuse surface: absorbed with probability 1 - rho, else a cosine-weighted bounce
                if (sampling::Random() >= m_rho) return;
                Vtx v = {};
                v.p[0] = float(rec.p.x); v.p[1] = float(rec.p.y); v.p[2] = float(rec.p.z);
                v.n[0] = float(rec.gn.x); v.n[1] = float(rec.gn.y); v.n[2] = float(rec.gn.z);
                v.rho = float(m_rho);
                v.tIn = float(std::exp(-m_sigT * lenInMedium));
                v.surface = 1;
                v.inFog = r.medium ? 1 : 0;
                out.push_back(v);
                lenInMedium = 0;
                const frame::Frame fr(rec.gn);
                const Vector3 wi = fr.localToWorld(sampling::CosineSampleHemisphere(Point2(sampling::Random(), sampling::Random())));
                r = rec.SpawnRay(wi);
            }
        }

        void traceLightPath(std::vector<Vtx>& Z, SampledWavelengths& lam) const {
            const Vector3 dir = m_light.sampleEmissionDir(Point2(sampling::Random(), sampling::Random()));
            walk(Ray(m_lightPos, dir, constants::RAY_EPSILON, nullptr), Z, lam);
        }

        /// Transmittance (scalar) of the straight segment from `a` (travelling in `mediumAtA`) to `b`; zero if a surface is in the way.
        Real transmittanceBetween(const Point3& a, const Medium* mediumAtA, const Point3& b, const SampledWavelengths& lam) const {
            const Vector3 d = b - a;
            const Real dist = glm::length(d);
            if (dist <= Real(0)) return Real(0);
            const Spectrum T = shadowTransmittance(m_scene, a, mediumAtA, d / dist, dist - Real(2) * constants::RAY_EPSILON, lam);
            return T[0];
        }

        /// The point of a vertex that a segment toward `toward` leaves from or arrives at (a surface point is moved off the surface, to the side of `toward`).
        Point3 anchor(const Vtx& v, const Point3& toward) const {
            const Point3 p = pos(v);
            if (!v.surface) return p;
            const Vector3 n = nrm(v);
            const Real side = glm::dot(n, toward - p) >= 0 ? Real(1) : Real(-1);
            return p + n * (side * constants::RAY_EPSILON * Real(4));
        }

        /// Transmittance of the connecting edge between two vertices.
        Real connectEdge(const Vtx& a, const Vtx& b, const SampledWavelengths& lam) const {
            const Point3 pa = anchor(a, pos(b)), pb = anchor(b, pos(a));
            return transmittanceBetween(pa, mediumOf(a), pb, lam);
        }

        // ------------------------------------------------------------------------------------------------
        // The path: weight of one path by the balance heuristic, and its contribution
        // ------------------------------------------------------------------------------------------------

        /**
         * Contribution f(x) / D(x) of the path x_1..x_n (X[0..n-1], from the light to the camera), without the light's spectrum.
         * Tedge[0..n]: transmittance of the edges light -> x_1, x_1 -> x_2, ..., x_n -> camera.
         * sActual: the strategy that made the path (x_1..x_s from the light, the rest from the camera): it decides which cosine a specular
         * vertex's weight is paired with.
         */
        Real pathContribution(const Vtx* X, int n, const Real* Tedge, int sActual) const {
            constexpr int kMax = 40;
            if (n < 1 || n + 2 > kMax) return Real(0);
            Point3 P[kMax];
            P[0] = m_lightPos;
            for (int i = 0; i < n; ++i) P[i + 1] = pos(X[i]);
            P[n + 1] = m_camPos;

            Vector3 dir[kMax];
            Real d2[kMax], T[kMax];
            for (int i = 0; i <= n; ++i) {
                const Vector3 e = P[i + 1] - P[i];
                d2[i] = glm::dot(e, e);
                if (d2[i] <= Real(0)) return Real(0);
                dir[i] = e / std::sqrt(d2[i]);
                T[i] = Tedge[i];
                if (!(T[i] > 0)) return Real(0);
            }

            // surfaces reflect; cosines at the vertices (1 in the fog)
            Real cin[kMax], cout_[kMax];
            for (int k = 1; k <= n; ++k) {
                const Vtx& v = X[k - 1];
                if (v.surface) {
                    const Vector3 nk = nrm(v);
                    const Real a = glm::dot(nk, dir[k - 1]), b = glm::dot(nk, dir[k]);
                    if (!v.spec && a * b >= 0) return Real(0);               // a diffuse surface reflects; glass also refracts
                    cin[k] = std::abs(a); cout_[k] = std::abs(b);
                }
                else { cin[k] = Real(1); cout_[k] = Real(1); }
            }

            // the pixel's density at x_n: 1 / Omega_pixel
            Real s, t, cosCam, dist;
            if (!m_cam.project(P[n], s, t, cosCam, dist)) return Real(0);
            const Real F = m_cam.focusDistance();
            const Real w = F * F * Real(m_P.width) * Real(m_P.height) / (m_cam.focalPlaneArea() * cosCam * cosCam * cosCam);

            if (!m_light.emitsToward(dir[0])) return Real(0);
            const Real pL = m_light.emissionPdf(dir[0]);

            // vertex factors and the densities of the direction a vertex scatters into, seen from the light (F) and from the camera (C)
            Real fv[kMax], scatF[kMax], scatC[kMax];
            for (int k = 1; k <= n; ++k) {
                if (X[k - 1].spec) {
                    fv[k] = Real(0);   // unused: the weight W replaces fv cin cout (below)
                    scatF[k] = Real(1);
                    scatC[k] = Real(1);
                }
                else if (X[k - 1].surface) {
                    fv[k] = X[k - 1].rho / constants::PI;
                    scatF[k] = cout_[k] / constants::PI;
                    scatC[k] = cin[k] / constants::PI;
                }
                else {
                    const Real ph = m_med.phase(glm::dot(dir[k - 1], dir[k]));
                    fv[k] = m_sigS * ph;
                    scatF[k] = ph;
                    scatC[k] = ph;
                }
            }
            scatF[0] = pL;
            scatC[n + 1] = w;

            Real f = w;
            for (int i = 0; i <= n; ++i) f *= T[i] / d2[i];
            for (int k = 1; k <= n; ++k) {
                if (X[k - 1].spec) f *= Real(X[k - 1].W) * (k <= sActual ? cin[k] : cout_[k]);   // f / p of the sampling gives W
                else f *= fv[k] * cin[k] * cout_[k];
            }

            // densities of the vertices when sampled from the light (F) and from the camera (C)
            Real pdfF[kMax], pdfC[kMax];
            for (int i = 1; i <= n; ++i) {
                const bool surf = X[i - 1].surface != 0, spec = X[i - 1].spec != 0;
                pdfF[i] = scatF[i - 1] * (spec ? cin[i] : surf ? X[i - 1].rho * cin[i] : m_sigS) * T[i - 1] / d2[i - 1];
                pdfC[i] = scatC[i + 1] * (spec ? cout_[i] : surf ? X[i - 1].rho * cout_[i] : m_sigS) * T[i] / d2[i];
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
                // a connection cannot end at a specular vertex (nor can the light or the camera be connected to one)
                if (sIdx == 0 ? X[0].spec : sIdx == n ? X[n - 1].spec : (X[sIdx - 1].spec || X[sIdx].spec)) continue;
                if (sIdx == 0 && !m_P.useNee) continue;
                if (sIdx == n && !m_P.useLightTracing) continue;
                if (sIdx > 0 && sIdx < n && !m_P.useConnections) continue;
                const Real count = (sIdx == n) ? Real(m_N) : Real(m_spp);
                D += count * prefix * suffix[sIdx];
            }
            if (!(D > Real(0)) || !std::isfinite(D)) return Real(0);
            const Real c = f / D;
            return std::isfinite(c) ? c : Real(0);
        }

        // ------------------------------------------------------------------------------------------------
        // Light tracing (s = n): every vertex of a light path is connected to the camera
        // ------------------------------------------------------------------------------------------------
        void lightTrace(const std::vector<Vtx>& Z, const SampledWavelengths& lam, std::vector<Vector3>& splat) const {
            if (m_P.lightTraceOnlyAfterGlass && !Z[0].spec) return;
            const Spectrum I = m_light.intensity(lam);
            Real T[40];
            for (size_t s = 1; s <= Z.size(); ++s) {
                const Vtx& v = Z[s - 1];
                if (v.spec) continue;                                            // the camera cannot be connected to glass
                Real sx, sy, cosCam, dist;
                if (!m_cam.project(pos(v), sx, sy, cosCam, dist)) continue;
                for (size_t i = 0; i < s; ++i) T[i] = Z[i].tIn;
                const Point3 pa = anchor(v, m_camPos);
                T[s] = transmittanceBetween(pa, mediumOf(v), m_camPos, lam);
                if (!(T[s] > 0)) continue;
                const Real c = pathContribution(Z.data(), static_cast<int>(s), T, static_cast<int>(s));
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
        Vector3 cameraSample(int i, int j, std::vector<Vtx>& Y, std::vector<Vtx>& X) const {
            const Real u = (Real(i) + sampling::Random()) / Real(m_P.width);
            const Real v = (Real(j) + sampling::Random()) / Real(m_P.height);
            const Ray ray = m_cam.getRay(u, v, Point2(0.5, 0.5));
            SampledWavelengths lam = makeLambda();
            walk(ray, Y, lam);
            if (Y.empty()) return Vector3(0.0);

            const Spectrum I = m_light.intensity(lam);
            Spectrum acc(0.0);

            // the cached light path of this sample
            size_t lb = 0, le = 0;
            if (m_P.useConnections && m_nStored > 0) {
                const size_t k = std::min(m_nStored - 1, static_cast<size_t>(sampling::Random() * Real(m_nStored)));
                lb = m_start[k]; le = m_start[k + 1];
            }

            Real T[40];
            const int nY = static_cast<int>(Y.size());
            for (int t = 1; t <= nY; ++t) {
                const Vtx& yt = Y[t - 1];
                if (yt.spec) continue;                                           // nothing can be connected to glass
                // edges of the camera part, from the connection end to the camera: x_{s+j} -- x_{s+j+1} is the edge recorded at y_{t-j+1}
                // ---- s = 0: next event estimation from y_t ----
                if (m_P.useNee) {
                    const Vector3 toL = m_lightPos - pos(yt);
                    const Real dl = glm::length(toL);
                    if (dl > Real(0) && m_light.emitsToward(-toL / dl)) {
                        T[0] = transmittanceBetween(anchor(yt, m_lightPos), mediumOf(yt), m_lightPos, lam);
                        if (T[0] > Real(0)) {
                            X.clear();
                            for (int q = t - 1; q >= 0; --q) X.push_back(Y[q]);
                            for (int jj = 1; jj <= t; ++jj) T[jj] = Y[t - jj].tIn;
                            acc += I * pathContribution(X.data(), t, T, 0);
                        }
                    }
                }
                // ---- 0 < s < n: a vertex of the cached light path connected to y_t ----
                for (size_t q = lb; q < le; ++q) {
                    const int s = static_cast<int>(q - lb) + 1;
                    const Vtx& zs = m_verts[q];
                    if (zs.spec) continue;
                    const Real Tc = connectEdge(zs, yt, lam);
                    if (!(Tc > Real(0))) continue;
                    X.clear();
                    for (size_t a = lb; a <= q; ++a) X.push_back(m_verts[a]);
                    for (int a = t - 1; a >= 0; --a) X.push_back(Y[a]);
                    for (int e = 0; e < s; ++e) T[e] = m_verts[lb + e].tIn;
                    T[s] = Tc;
                    for (int jj = 1; jj <= t; ++jj) T[s + jj] = Y[t - jj].tIn;
                    acc += I * pathContribution(X.data(), s + t, T, s);
                }
            }
            const Vector3 rgb = SpectrumToRGB(acc, lam);
            return math::hasNonFinite(rgb) ? Vector3(0.0) : rgb;
        }
    };

} // namespace rayt
