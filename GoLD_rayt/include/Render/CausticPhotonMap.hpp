/**
 * @file Render/CausticPhotonMap.hpp
 * @brief Caustic photon map: light paths  L S+ D  (light, one or more perfectly specular bounces, then a non-specular surface).
 *
 * Why. The unidirectional path tracer finds such a path only when a BSDF-sampled ray from the receiver D happens to leave
 * through the specular object and hit the light: with a small or concentrated light (the sun of an environment map, a point
 * light, which can never be hit) that is rare or impossible, so caustics are grainy fireflies at best. Photon mapping
 * traces the same paths from the light, where the specular chain is followed deterministically.
 *
 * Emission (photon targeting). Photons are only sent toward a bounding sphere (center, radius) of the specular objects.
 * A light is picked with a PowerLightSampler and sampled as seen from the sphere center (Light::sampleLi), which gives the
 * radiance Li and the solid-angle density pdf at the center for every light type (area, point, environment). The flux
 * through the disc of radius R perpendicular to that direction is  Li * pi R^2 / (pdf * pSelect), divided by the number of
 * emitted photons N. The photon starts at the sampled point of a finite light (at emitDistance along the direction for an
 * infinite light) and flies toward a uniformly sampled point of the disc. This is exact for R << distance to the light.
 * A point light inside the sphere (a lamp inside a glass vessel) is emitted isotropically instead (flux I * 4 pi / N).
 * EVERY specular object of the scene must lie inside the sphere: the integrator drops all L S+ D paths found by BSDF
 * sampling and relies on the photons for them.
 *
 * Tracing. At a specular surface the photon continues (BSDF sample, same as the camera path: the 1/eta^2 radiance scaling
 * of refraction cancels between entering and leaving, so photons that end in air carry the right flux); an absorbing medium
 * multiplies the flux by its transmittance; a scattering medium, a light or a non-opaque hit ends the photon. At the first
 * non-specular surface after at least one specular bounce the photon is stored; without one it is direct light, which next
 * event estimation already handles, and is dropped.
 *
 * Spectral mode. A photon carries its own wavelengths, so its flux is converted to linear RGB at storage (SpectrumToRGB, the
 * film's own conversion). At the receiver the mean color of the gathered photons (unit luminance) is turned back into a
 * light-like spectrum for the camera path's wavelengths and multiplied by the receiver's BSDF at those wavelengths and by the
 * summed photon luminance. The caustic light is therefore right in RGB, while the metameric detail of the light after the
 * specular chain is lost (dispersion and absorption colors are kept; the color is averaged over one gather radius).
 *
 * Density estimation. Cone-filtered fixed radius:  L = 3 / (pi r^2) * sum_p (1 - d_p / r) f(wo, wi_p) Phi_p,  photons on the
 * same surface side only (normal test) and within a thin slab around the receiver's tangent plane. Biased (blurs detail
 * smaller than r) but consistent as photons increase and the radius shrinks. A hashed uniform grid of cell size r finds them.
 *
 * Not handled: caustics seen through media with scattering, receivers inside a medium, subsurface receivers (their exit lobe).
 */

#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
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
#include "Lights/Light.hpp"
#include "Lights/LightSampler.hpp"
#include "Materials/Material.hpp"
#include "Materials/Medium.hpp"
#include "Color/SpectralConversion.hpp"

namespace rayt {

    /// Photons are traced and gathered in air and inside purely absorbing media (their transmittance is a weight); a scattering
    /// medium is not traced, so neither tracing nor gathering happens inside one.
    inline bool photonsAllowedIn(const Medium* medium) { return !medium || medium->attenuatesContinuously(); }

    struct CausticPhotonParams {
        long long numPhotons = 1000000;   ///< photons emitted (the number stored is smaller: misses and direct light are dropped)
        Real radius = Real(0.0015);       ///< gather radius (m)
        Point3 targetCenter = Point3(0.0);///< bounding sphere of ALL specular objects
        Real targetRadius = Real(0);
        Real emitDistance = Real(10);     ///< infinite lights: photons start this far from the sphere (beyond the scene)
        int maxDepth = 32;                ///< specular bounces per photon
        /// >= 0: store only photons whose chain is exactly this many specular bounces, ALL of them transmissions (T^k). For tests that
        /// compare the map with an estimator of one chain type (the specular manifold NEE handles T T only); -1 = every chain.
        int onlyTransmissions = -1;
    };

    class CausticPhotonMap {
    public:
        struct Photon {
            float p[3];     ///< position
            float n[3];     ///< geometric normal of the surface side the photon arrived on
            float wi[3];    ///< direction toward the light (the reversed travel direction)
            float rgb[3];   ///< flux (linear RGB, already divided by the number of emitted photons)
        };

        bool empty() const { return m_photons.empty(); }
        size_t size() const { return m_photons.size(); }
        Real radius() const { return m_radius; }
        size_t memoryBytes() const { return m_photons.capacity() * sizeof(Photon) + m_start.capacity() * sizeof(uint32_t); }

        /// Traces the photons (OpenMP) and builds the grid.
        void build(const Scene& scene, const CausticPhotonParams& P) {
            const auto t0 = std::chrono::steady_clock::now();
            m_radius = P.radius;
            m_photons.clear();
            m_start.clear();
            std::vector<std::shared_ptr<Light>> lights;   // only the lights whose caustics belong to the photons
            for (const auto& l : scene.allLights()) if (l->photonCaustics()) lights.push_back(l);
            if (lights.empty() || P.targetRadius <= Real(0) || P.numPhotons <= 0) return;

            const PowerLightSampler sampler(lights);
#ifdef _OPENMP
            const int nThreads = omp_get_max_threads();
#else
            const int nThreads = 1;
#endif
            std::vector<std::vector<Photon>> perThread(static_cast<size_t>(nThreads));
            const long long chunk = 4096;
            const int nChunks = static_cast<int>((P.numPhotons + chunk - 1) / chunk);
            const Real invN = Real(1) / Real(P.numPhotons);

#pragma omp parallel for schedule(dynamic, 1)
            for (int c = 0; c < nChunks; ++c) {
#ifdef _OPENMP
                auto& mine = perThread[static_cast<size_t>(omp_get_thread_num())];
#else
                auto& mine = perThread[0];
#endif
                const long long first = static_cast<long long>(c) * chunk;
                const long long last = std::min(P.numPhotons, first + chunk);
                for (long long i = first; i < last; ++i) tracePhoton(scene, P, sampler, invN, mine);
            }

            size_t total = 0;
            for (const auto& v : perThread) total += v.size();
            m_photons.reserve(total);
            for (auto& v : perThread) { m_photons.insert(m_photons.end(), v.begin(), v.end()); std::vector<Photon>().swap(v); }
            buildGrid();

            const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::cout << "[Caustics] " << P.numPhotons << " photons emitted, " << m_photons.size() << " stored ("
                << (m_photons.empty() ? 0.0 : 100.0 * double(m_photons.size()) / double(P.numPhotons)) << " %), gather radius "
                << m_radius * 1000 << " mm, " << double(memoryBytes()) / (1024.0 * 1024.0) << " MB, " << sec << " s\n";
        }

        /**
         * @brief Reflected caustic radiance toward `wo` at a non-specular surface point (light-like Spectrum, same units as Li()).
         * @param rec  the hit (position, geometric normal, material)
         * @param lambda the camera path's wavelengths
         */
        Spectrum radiance(const SurfaceInteraction& rec, const Vector3& wo, const SampledWavelengths& lambda) const {
            if (m_photons.empty()) return Spectrum(0.0);
            const Real r = m_radius, r2 = r * r;
            const Real slab = Real(0.25) * r;   // half thickness around the tangent plane
            const int ix = cellIndex(rec.p.x), iy = cellIndex(rec.p.y), iz = cellIndex(rec.p.z);

#if defined(RAYT_SPECTRAL_MODE)
            Spectrum acc(0.0);          // sum of f * weight * luminance of the photon
            Vector3 colorSum(0.0);      // sum of weight * rgb: its chromaticity is the color of the caustic light
#else
            Spectrum acc(0.0);          // RGB mode: sum of f * weight * rgb, directly
#endif
            bool any = false;
            uint32_t seen[27];
            int nSeen = 0;
            for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
                const uint32_t h = hashCell(ix + dx, iy + dy, iz + dz) & m_mask;
                bool dup = false;
                for (int k = 0; k < nSeen; ++k) if (seen[k] == h) { dup = true; break; }
                if (dup) continue;
                seen[nSeen++] = h;
                for (uint32_t i = m_start[h]; i < m_start[h + 1]; ++i) {
                    const Photon& ph = m_photons[i];
                    const Vector3 d(Real(ph.p[0]) - rec.p.x, Real(ph.p[1]) - rec.p.y, Real(ph.p[2]) - rec.p.z);
                    const Real dist2 = glm::dot(d, d);
                    if (dist2 >= r2) continue;
                    if (std::abs(glm::dot(d, rec.gn)) > slab) continue;
                    if (Real(ph.n[0]) * rec.gn.x + Real(ph.n[1]) * rec.gn.y + Real(ph.n[2]) * rec.gn.z < Real(0.7)) continue;
                    const Vector3 wi(ph.wi[0], ph.wi[1], ph.wi[2]);
                    const Spectrum f = rec.matPtr->eval(rec, wo, wi, lambda);
                    if (isBlack(f)) continue;
                    const Real w = Real(1) - std::sqrt(dist2) / r;   // cone filter
                    const Vector3 rgb(ph.rgb[0], ph.rgb[1], ph.rgb[2]);
#if defined(RAYT_SPECTRAL_MODE)
                    acc += f * (w * luminance(rgb));
                    colorSum += rgb * w;
#else
                    acc += f * Spectrum(rgb * w);
#endif
                    any = true;
                }
            }
            if (!any) return Spectrum(0.0);

            const Real norm = Real(3) / (constants::PI * r2);   // cone filter normalization: 1 / (pi r^2 (1 - 2/3))
#if defined(RAYT_SPECTRAL_MODE)
            // One light-like spectrum for the mean color of the photons (unit luminance). Splitting into the unit RGB basis
            // would be linear only without clamping: the saturated basis spectra are clamped and no longer add up to white.
            const Real Y = luminance(colorSum);
            if (Y <= Real(0)) return Spectrum(0.0);
            return acc * RgbToIlluminantSpectrum(colorSum / Y, lambda) * norm;
#else
            return acc * norm;
#endif
        }

    private:
        std::vector<Photon> m_photons;   // sorted by grid bucket
        std::vector<uint32_t> m_start;   // bucket b holds photons [m_start[b], m_start[b + 1])
        uint32_t m_mask = 0;
        Real m_radius = Real(0.0015);

        static Real luminance(const Vector3& c) { return Real(0.2126) * c.x + Real(0.7152) * c.y + Real(0.0722) * c.z; }
        int cellIndex(Real x) const { return static_cast<int>(std::floor(x / m_radius)); }
        static uint32_t hashCell(int x, int y, int z) {
            return (static_cast<uint32_t>(x) * 73856093u) ^ (static_cast<uint32_t>(y) * 19349663u) ^ (static_cast<uint32_t>(z) * 83492791u);
        }

        void buildGrid() {
            uint32_t size = 1024;
            while (size < 2 * m_photons.size()) size <<= 1;
            m_mask = size - 1;
            std::vector<uint32_t> bucket(m_photons.size());
            m_start.assign(static_cast<size_t>(size) + 1, 0);
            for (size_t i = 0; i < m_photons.size(); ++i) {
                const Photon& p = m_photons[i];
                bucket[i] = hashCell(cellIndex(p.p[0]), cellIndex(p.p[1]), cellIndex(p.p[2])) & m_mask;
                ++m_start[bucket[i] + 1];
            }
            for (size_t b = 0; b < size; ++b) m_start[b + 1] += m_start[b];
            std::vector<Photon> sorted(m_photons.size());
            std::vector<uint32_t> cursor(m_start.begin(), m_start.end() - 1);
            for (size_t i = 0; i < m_photons.size(); ++i) sorted[cursor[bucket[i]]++] = m_photons[i];
            m_photons.swap(sorted);
        }

        /// Emits and follows one photon; appends it to `out` when it ends as a caustic photon.
        static void tracePhoton(const Scene& scene, const CausticPhotonParams& P, const LightSampler& sampler, Real invN,
            std::vector<Photon>& out)
        {
#if defined(RAYT_SPECTRAL_MODE)
            SampledWavelengths lambda = SampledWavelengths::sampleVisible(sampling::Random());
#else
            SampledWavelengths lambda;
#endif
            // ---- emission, as seen from the target center ----
            LightSampleContext ctx;
            ctx.p = P.targetCenter;
            ctx.gn = Vector3(0, 1, 0);
            const auto sl = sampler.sample(ctx, sampling::Random());
            if (!sl || sl->p <= Real(0)) return;
            const Point2 uLight(sampling::Random(), sampling::Random());
            const auto ls = sl->light->sampleLi(ctx, uLight, lambda);
            if (!ls || ls->pdf <= Real(0) || isBlack(ls->Li)) return;

            Point3 origin;
            Vector3 dir;
            Spectrum power;   // flux of this photon
            if (ls->isDelta && ls->tMax < P.targetRadius) {
                // A point light INSIDE the target sphere (for example inside a glass vessel): the disc construction below needs the
                // light outside, so emit uniformly over the sphere of directions. Li = I / d^2 with d = ls->tMax, hence I = Li d^2,
                // and the photon carries I * 4 pi / N.
                const Real cosT = Real(1) - Real(2) * sampling::Random();
                const Real sinT = std::sqrt(std::max(Real(0), Real(1) - cosT * cosT));
                const Real ph = Real(2) * constants::PI * sampling::Random();
                origin = ls->pLight;
                dir = Vector3(sinT * std::cos(ph), cosT, sinT * std::sin(ph));
                power = ls->Li * (ls->tMax * ls->tMax * Real(4) * constants::PI * invN / (ls->pdf * sl->p));
            }
            else {
                const frame::Frame fr(ls->wi);
                const Real rd = P.targetRadius * std::sqrt(sampling::Random());
                const Real phi = Real(2) * constants::PI * sampling::Random();
                const Point3 q = P.targetCenter + fr.localToWorld(Vector3(rd * std::cos(phi), rd * std::sin(phi), Real(0)));
                origin = std::isfinite(ls->tMax) ? ls->pLight : q + ls->wi * P.emitDistance;
                dir = glm::normalize(q - origin);
                const Real discArea = constants::PI * P.targetRadius * P.targetRadius;
                power = ls->Li * (discArea * invN / (ls->pdf * sl->p));
            }
            Spectrum beta(1.0);

            // ---- specular chain ----
            Ray r(origin, dir, constants::RAY_EPSILON, nullptr);
            int specularBounces = 0;
            bool onlyTransmissions = true;   // every specular bounce so far was a transmission
            for (int depth = 0; depth < P.maxDepth; ++depth) {
                SurfaceInteraction rec;
                if (!scene.hit(r, rec)) return;

                if (r.medium) {
                    if (!r.medium->attenuatesContinuously()) return;           // a scattering medium: not traced
                    beta *= r.medium->transmittance(rec.t, lambda);
                }
                rec.medium = r.medium;

                const Real opacity = rec.matPtr->opacity(rec);
                if (opacity < Real(1) && sampling::Random() >= opacity) { r = rec.SpawnRay(r.d); continue; }

                if (rec.areaLight) return;                                      // reached a light
                rec.matPtr->prepareWavelengths(lambda);
                if (!isBlack(rec.matPtr->emitted(rec, -r.d, lambda))) return;

                if (!rec.matPtr->isSpecular()) {
                    const bool chainWanted = P.onlyTransmissions < 0 || (onlyTransmissions && specularBounces == P.onlyTransmissions);
                    if (specularBounces > 0 && chainWanted) {   // (a scattering medium ended the photon above; an absorbing one is fine)
                        const Vector3 rgb = SpectrumToRGB(power * beta, lambda);
                        if (!math::hasNonFinite(rgb)) {
                            Photon ph;
                            ph.p[0] = float(rec.p.x); ph.p[1] = float(rec.p.y); ph.p[2] = float(rec.p.z);
                            ph.n[0] = float(rec.gn.x); ph.n[1] = float(rec.gn.y); ph.n[2] = float(rec.gn.z);
                            const Vector3 wi = -r.d;
                            ph.wi[0] = float(wi.x); ph.wi[1] = float(wi.y); ph.wi[2] = float(wi.z);
                            ph.rgb[0] = float(rgb.x); ph.rgb[1] = float(rgb.y); ph.rgb[2] = float(rgb.z);
                            out.push_back(ph);
                        }
                    }
                    return;   // direct light (no specular bounce) is left to next event estimation
                }

                const Point2 u(sampling::Random(), sampling::Random());
                const auto s = rec.matPtr->sample(rec, -r.d, u, lambda);
                if (!s || !s->isSpecular()) return;
                if (!s->isTransmission()) onlyTransmissions = false;
                beta *= s->f;
                if (isBlack(beta)) return;
                ++specularBounces;

                // Russian roulette on the relative throughput
                if (specularBounces >= 4) {
                    const Real rr = std::clamp(math::maxComponent(beta), Real(0.05), Real(0.95));
                    if (sampling::Random() > rr) return;
                    beta /= rr;
                }

                const Vector3 wi = s->wi;
                if (const Medium* inner = rec.matPtr->mediumInside()) {
                    const bool wiOnRaySide = glm::dot(wi, rec.gn) > 0;
                    rec.medium = (wiOnRaySide == rec.frontFace) ? nullptr : inner;
                }
                if (rec.medium) {
                    r = Ray(rec.p, wi, constants::RAY_EPSILON, rec.medium);
                    r.time = rec.time;
                }
                else r = rec.SpawnRay(wi);
            }
        }
    };

} // namespace rayt
