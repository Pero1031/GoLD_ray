// The caustic photon map (Render/CausticPhotonMap.hpp): emission and tracing of the photons, the hashed grid, the gathering, light tracing.
#include "pch.h"
#include "Render/CausticPhotonMap.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Core/Math.hpp"
#include "Core/Ray.hpp"
#include "Core/Sampling.hpp"
#include "Core/Constants.hpp"
#include "Core/SpectrumUtils.hpp"
#include "Geometry/Frame.hpp"
#include "Scene/Scene.hpp"
#include "Render/Camera.hpp"
#include "Lights/Light.hpp"
#include "Lights/LightSampler.hpp"
#include "Materials/Material.hpp"
#include "Render/Transmittance.hpp"
#include "Color/SpectralConversion.hpp"

namespace rayt {

    bool photonsAllowedIn(const Medium* medium, bool mediaTraced) {
        return !medium || medium->attenuatesContinuously() || (mediaTraced && medium->isSpectrallyConstant());
    }

    void CausticPhotonMap::build(const Scene& scene, const CausticPhotonParams& P) {
        const auto t0 = std::chrono::steady_clock::now();
        m_radius = P.radius;
        m_tracesMedia = P.traceMedia && P.lightTraceCamera;
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
        const bool trace = P.lightTraceCamera && P.traceWidth > 0 && P.traceHeight > 0;
        const size_t nPix = static_cast<size_t>(P.traceWidth) * static_cast<size_t>(P.traceHeight);
        std::vector<std::vector<Vector3>> splatPerThread(static_cast<size_t>(nThreads));
        if (trace) for (auto& v : splatPerThread) v.assign(nPix, Vector3(0.0));
        const long long chunk = 4096;
        const int nChunks = static_cast<int>((P.numPhotons + chunk - 1) / chunk);
        const Real invN = Real(1) / Real(P.numPhotons);

#pragma omp parallel for schedule(dynamic, 1)
        for (int c = 0; c < nChunks; ++c) {
#ifdef _OPENMP
            auto& mine = perThread[static_cast<size_t>(omp_get_thread_num())];
            auto& mySplat = splatPerThread[static_cast<size_t>(omp_get_thread_num())];
#else
            auto& mine = perThread[0];
            auto& mySplat = splatPerThread[0];
#endif
            const long long first = static_cast<long long>(c) * chunk;
            const long long last = std::min(P.numPhotons, first + chunk);
            for (long long i = first; i < last; ++i) tracePhoton(scene, P, sampler, invN, mine, mySplat);
        }

        size_t total = 0;
        for (const auto& v : perThread) total += v.size();
        m_splat.clear();
        if (trace) {
            m_splat.assign(nPix, Vector3(0.0));
            for (const auto& v : splatPerThread) for (size_t i = 0; i < nPix; ++i) m_splat[i] += v[i];
        }
        m_photons.reserve(total);
        for (auto& v : perThread) { m_photons.insert(m_photons.end(), v.begin(), v.end()); std::vector<Photon>().swap(v); }
        buildGrid();

        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cout << "[Caustics] " << P.numPhotons << " photons emitted, " << m_photons.size() << " stored ("
            << (m_photons.empty() ? 0.0 : 100.0 * double(m_photons.size()) / double(P.numPhotons)) << " %), gather radius "
            << m_radius * 1000 << " mm, " << double(memoryBytes()) / (1024.0 * 1024.0) << " MB, " << sec << " s\n";
    }

    Spectrum CausticPhotonMap::radiance(const SurfaceInteraction& rec, const Vector3& wo, const SampledWavelengths& lambda, Vector3* varianceRgb) const {
        if (varianceRgb) *varianceRgb = Vector3(0.0);
        if (m_photons.empty()) return Spectrum(0.0);
        const Real r = m_radius, r2 = r * r;
        const Real slab = Real(0.25) * r;   // half thickness around the tangent plane
        const int ix = cellIndex(rec.p.x), iy = cellIndex(rec.p.y), iz = cellIndex(rec.p.z);

#if defined(RAYT_SPECTRAL_MODE)
        Spectrum acc(0.0);          // sum of f * weight * luminance of the photon
        Vector3 colorSum(0.0);      // sum of weight * rgb: its chromaticity is the color of the caustic light
        Real w2 = 0;                // sum of (weight * luminance)^2: the relative variance of the gather, the photons taken as one colour
#else
        Spectrum acc(0.0);          // RGB mode: sum of f * weight * rgb, directly
        Spectrum acc2(0.0);         // sum of the squared contributions: the variance of the gather (the photons are independent)
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
                const Real wl = w * luminance(rgb);
                acc += f * wl;
                colorSum += rgb * w;
                w2 += wl * wl;
#else
                const Spectrum c = f * Spectrum(rgb * w);
                acc += c;
                acc2 += c * c;
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
        const Spectrum L = acc * RgbToIlluminantSpectrum(colorSum / Y, lambda) * norm;
        if (varianceRgb) {   // relative variance sum (w lum)^2 / (sum w lum)^2 of the gather, applied to the mean: the BSDF is taken as the same for every photon
            const Vector3 m = SpectrumToRGB(L, lambda);
            *varianceRgb = m * m * (w2 / (Y * Y));
        }
        return L;
#else
        if (varianceRgb) *varianceRgb = SpectrumToRGB(acc2, lambda) * (norm * norm);
        return acc * norm;
#endif
    }

    int CausticPhotonMap::cellIndex(Real x) const { return static_cast<int>(std::floor(x / m_radius)); }

    uint32_t CausticPhotonMap::hashCell(int x, int y, int z) {
        return (static_cast<uint32_t>(x) * 73856093u) ^ (static_cast<uint32_t>(y) * 19349663u) ^ (static_cast<uint32_t>(z) * 83492791u);
    }

    void CausticPhotonMap::buildGrid() {
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

    /**
     * Light tracing: connects a caustic photon (flux `flux` arriving from the light at the receiver `rec`, travelling along -wiLight) to the pinhole
     * camera. A pixel's average radiance is its footprint's reflected flux per area:  Phi f cos(theta_x) / (Omega_pix d^2),  with the solid angle
     * Omega_pix = A_plane cos^3(theta_c) / (N F^2) of a pixel at the angle theta_c from the view axis (A_plane: the film on the focal plane at distance F,
     * N pixels) and the distance d to the camera. `flux` is already divided by the number of emitted photons, so the sum over the photons is the
     * image itself. The receiver's BSDF is evaluated at the photon's wavelengths: the dispersion colors are kept exactly.
     */
    void CausticPhotonMap::connectToCamera(const Scene& scene, const CausticPhotonParams& P, const SurfaceInteraction& rec, const Vector3& wiLight,
        const Spectrum& flux, const SampledWavelengths& lambda, std::vector<Vector3>& splat)
    {
        const Camera& cam = *P.lightTraceCamera;
        if (cam.medium()) return;                                  // a camera inside a medium: not handled
        Real s, t, cosCam, dist;
        if (!cam.project(rec.p, s, t, cosCam, dist) || dist <= Real(0)) return;
        const Vector3 toCam = (cam.origin() - rec.p) / dist;
        const Real cosX = glm::dot(rec.gn, toCam);
        if (cosX * glm::dot(rec.gn, wiLight) <= Real(0)) return;   // reflection only: the camera is on the side the photon came from
        const Spectrum f = rec.matPtr->eval(rec, toCam, wiLight, lambda);
        if (isBlack(f)) return;
        Spectrum T(1.0);
        if (P.traceMedia) {   // fog between the receiver and the camera: attenuation, and its boundary is crossed
            T = shadowTransmittance(scene, rec.OffsetRayOrigin(toCam), rec.medium, toCam, dist - Real(2) * constants::RAY_EPSILON, lambda, rec.time);
            if (isBlack(T)) return;
        }
        else if (scene.occluded(rec.SpawnRay(toCam, dist))) return;   // something between the receiver and the camera

        const Real F = cam.focusDistance();
        const Real w = F * F * Real(P.traceWidth) * Real(P.traceHeight) / (cam.focalPlaneArea() * cosCam * cosCam * cosCam);
        const Vector3 rgb = SpectrumToRGB(flux * f * T, lambda) * (std::abs(cosX) * w / (dist * dist));
        if (math::hasNonFinite(rgb)) return;
        const int px = std::min(P.traceWidth - 1, static_cast<int>(s * Real(P.traceWidth)));
        const int py = P.traceHeight - 1 - std::min(P.traceHeight - 1, static_cast<int>(t * Real(P.traceHeight)));   // the film's rows run top to bottom
        splat[static_cast<size_t>(py) * static_cast<size_t>(P.traceWidth) + static_cast<size_t>(px)] += rgb;
    }

    /**
     * Light tracing in a medium: connects a photon that scatters at `x` (it arrived travelling along `dirIn`, flux `flux`) to the camera.
     * A scattering event deposits the flux as radiant intensity  flux * phase(cos)  into the direction of the camera; the pixel that sees
     * the point averages the radiance over its footprint, so (as for a surface, with the phase function in place of  f cos(theta_x) )
     * the contribution is  flux * phase * T / (Omega_pix d^2),  T the transmittance to the camera.
     */
    void CausticPhotonMap::connectMediumToCamera(const Scene& scene, const CausticPhotonParams& P, const Point3& x, const Vector3& dirIn, const Medium& med,
        const Spectrum& flux, const SampledWavelengths& lambda, std::vector<Vector3>& splat)
    {
        const Camera& cam = *P.lightTraceCamera;
        if (cam.medium()) return;
        Real s, t, cosCam, dist;
        if (!cam.project(x, s, t, cosCam, dist) || dist <= Real(0)) return;
        const Vector3 toCam = (cam.origin() - x) / dist;
        const Real ph = med.phase(glm::dot(dirIn, toCam));   // travel direction before and after the collision
        if (ph <= Real(0)) return;
        const Spectrum T = shadowTransmittance(scene, x, &med, toCam, dist - Real(2) * constants::RAY_EPSILON, lambda);
        if (isBlack(T)) return;

        const Real F = cam.focusDistance();
        const Real w = F * F * Real(P.traceWidth) * Real(P.traceHeight) / (cam.focalPlaneArea() * cosCam * cosCam * cosCam);
        const Vector3 rgb = SpectrumToRGB(flux * T, lambda) * (ph * w / (dist * dist));
        if (math::hasNonFinite(rgb)) return;
        const int px = std::min(P.traceWidth - 1, static_cast<int>(s * Real(P.traceWidth)));
        const int py = P.traceHeight - 1 - std::min(P.traceHeight - 1, static_cast<int>(t * Real(P.traceHeight)));
        splat[static_cast<size_t>(py) * static_cast<size_t>(P.traceWidth) + static_cast<size_t>(px)] += rgb;
    }

    /// Emits and follows one photon; appends it to `out` when it ends as a caustic photon.
    void CausticPhotonMap::tracePhoton(const Scene& scene, const CausticPhotonParams& P, const LightSampler& sampler, Real invN,
        std::vector<Photon>& out, std::vector<Vector3>& splat)
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
        int mediumEvents = 0, passes = 0;   // scattering events in media; null interfaces crossed (both bounded)
        bool onlyTransmissions = true;   // every specular bounce so far was a transmission
        const Material* firstMaterial = nullptr;    // the materials of the first two specular bounces (for skipTransmissionPairs)
        const Material* secondMaterial = nullptr;
        for (int depth = 0; depth < P.maxDepth; ++depth) {
            SurfaceInteraction rec;
            if (!scene.hit(r, rec)) return;

            if (r.medium) {
                if (!r.medium->attenuatesContinuously()) {
                    if (!P.traceMedia || !r.medium->isSpectrallyConstant()) return;   // a scattering medium that is not traced
                    // Gray scattering medium: delta tracking (null collisions for a density field), as in the camera path's walk.
                    const Medium& med = *r.medium;
                    const Real sT = med.sigmaTSpectrum(lambda)[0], sS = med.sigmaSSpectrum(lambda)[0];
                    Real tCur = 0;
                    bool scattered = false;
                    for (;;) {
                        tCur -= std::log(Real(1) - sampling::Random()) / sT;
                        if (tCur >= rec.t) break;                          // reaches the surface
                        if (++mediumEvents > 4096) return;
                        const Point3 x = r.o + r.d * tCur;
                        const Real q = med.density(x), u = sampling::Random();
                        if (u >= q) continue;                              // null collision
                        if (u * sT >= q * sS) return;                      // absorbed
                        // A real scattering event. Without a specular bounce before it the photon is direct light of the medium (left to
                        // the camera paths), and its direction was aimed at the target: drop it.
                        if (specularBounces == 0) return;
                        if (P.lightTraceCamera && !splat.empty()) connectMediumToCamera(scene, P, x, r.d, med, power * beta, lambda, splat);
                        const Point2 uPhase(sampling::Random(), sampling::Random());
                        r = Ray(x, med.samplePhase(r.d, uPhase), Real(0), &med);
                        scattered = true;
                        break;
                    }
                    if (scattered) { --depth; continue; }                  // follow the new direction (not a specular bounce)
                }
                else beta *= r.medium->transmittance(rec.t, lambda);
            }
            rec.medium = r.medium;

            // A boundary with no optical effect (the surface of a fog): cross it, switching the medium; not a bounce.
            if (rec.matPtr->isNullInterface() && ++passes < 256) {
                if (rec.matPtr->isMediumBoundary()) rec.medium = rec.frontFace ? rec.matPtr->mediumInside() : rec.matPtr->mediumOutside();
                if (rec.medium) { r = Ray(rec.p, r.d, constants::RAY_EPSILON, rec.medium); r.time = rec.time; }
                else r = rec.SpawnRay(r.d);
                --depth;
                continue;
            }

            const Real opacity = rec.matPtr->opacity(rec);
            if (opacity < Real(1) && sampling::Random() >= opacity) { r = rec.SpawnRay(r.d); continue; }

            if (rec.areaLight) return;                                      // reached a light
            rec.matPtr->prepareWavelengths(lambda);
            if (!isBlack(rec.matPtr->emitted(rec, -r.d, lambda))) return;

            if (!rec.matPtr->isSpecular()) {
                bool chainWanted = P.onlyTransmissions < 0 || (onlyTransmissions && specularBounces == P.onlyTransmissions);
                if (chainWanted && onlyTransmissions && specularBounces == 2 && !P.skipTransmissionPairs.empty()) {
                    const auto& skip = P.skipTransmissionPairs;
                    const bool solved = std::find(skip.begin(), skip.end(), firstMaterial) != skip.end()
                                     && std::find(skip.begin(), skip.end(), secondMaterial) != skip.end();
                    if (solved) chainWanted = false;   // the manifold NEE accounts for this T T chain
                }
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
                    if (P.lightTraceCamera && !splat.empty()) connectToCamera(scene, P, rec, -r.d, power * beta, lambda, splat);
                }
                return;   // direct light (no specular bounce) is left to next event estimation
            }

            const Point2 u(sampling::Random(), sampling::Random());
            const auto s = rec.matPtr->sample(rec, -r.d, u, lambda);
            if (!s || !s->isSpecular()) return;
            if (!s->isTransmission()) onlyTransmissions = false;
            if (specularBounces == 0) firstMaterial = rec.matPtr;
            else if (specularBounces == 1) secondMaterial = rec.matPtr;
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
            if (rec.matPtr->isMediumBoundary()) {
                const bool wiOnRaySide = glm::dot(wi, rec.gn) > 0;
                rec.medium = (wiOnRaySide == rec.frontFace) ? rec.matPtr->mediumOutside() : rec.matPtr->mediumInside();
            }
            if (rec.medium) {
                r = Ray(rec.p, wi, constants::RAY_EPSILON, rec.medium);
                r.time = rec.time;
            }
            else r = rec.SpawnRay(wi);
        }
    }

} // namespace rayt
