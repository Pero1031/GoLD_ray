/**
 * @file Render/Integrator.hpp
 * @brief The path tracer: Integrator (interface) and PathIntegrator (unidirectional path tracing with next-event
 *        estimation, MIS, participating media, subsurface scattering and photon caustics).
 */

#pragma once

#include "Materials/Medium.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <memory>
#include <iostream>
#include <algorithm>
#include <atomic>

#include "Core/Forward.hpp"
#include "Core/Types.hpp"
#include "Core/Math.hpp"
#include "Core/Ray.hpp"
#include "Core/Sampling.hpp"
#include "Core/Constants.hpp"
#include "Core/Interaction.hpp"
#include "Core/SpectrumUtils.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Scene/Scene.hpp" 
#include "Render/Film.hpp"
#include "Render/Camera.hpp"
#include "Render/SubsurfaceProbe.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/SpecularNee.hpp"
#include "Lights/Light.hpp"
#include "Lights/EnvLight.hpp"
#include "Lights/PointLight.hpp"
#include "Lights/LightSampler.hpp"
#include "Materials/Material.hpp"
#include "Color/SpectralConversion.hpp"
#include "IO/EnvMap.hpp"

namespace rayt {

    /**
     * @brief Denoiser auxiliary features of one camera path (see Film::enableAOVs).
     *
     * Recorded at the first vertex that is not perfectly specular, so glass and mirrors show what is
     * seen through / in them; a path that escapes to the environment records the environment color.
     */
    struct PathAOV {
        Vector3 albedo = Vector3(0.0);   ///< lighting-free color, [0,1]
        Vector3 normal = Vector3(0.0);   ///< shading normal facing the ray (zero for the environment)
        bool valid = false;
    };

    // Abstract base class for all rendering algorithms.
    class Integrator {
    public:
        virtual ~Integrator() = default;

        // Pure virtual function
        // The types Scene and Film are used here, so the includes above are required
        virtual void render(const Scene& scene, Film& film) = 0;
    };

    /**
     * @brief Path Tracing Integrator with Multiple Importance Sampling.
     *
     * Features:
     * - Supports arbitrary number of lights via LightSampler
     * - Environment/infinite light support
     * - Multiple Importance Sampling (MIS) for variance reduction
     * - Handles both delta and non-delta lights correctly
     */
    class PathIntegrator : public Integrator {
    public:
        PathIntegrator(std::shared_ptr<Camera> camera, int maxDepth, int spp)
            : m_camera(camera), m_maxDepth(maxDepth), m_spp(spp) {}

        /// Safety limit of collisions inside participating media per path (paths that reach it are dropped, which
        /// biases the image; volumeCapHits() counts them so that the bias can be checked).
        void setMaxVolumeEvents(int n) { m_maxVolumeEvents = n; }
        long long volumeCapHits() const { return m_volumeCapHits.load(); }

        /**
         * @brief Caustics from a photon map (null = off: caustics are then found by BSDF sampling only).
         *
         * With a map the integrator adds the photon estimate at every non-specular vertex and drops the light that a
         * BSDF-sampled ray finds through a chain of specular bounces that started at such a vertex (L S+ D paths), so
         * nothing is counted twice. Every specular object must lie inside the map's target sphere. An empty map is ignored.
         */
        void setCausticMap(const CausticPhotonMap* map) { m_caustics = (map && !map->empty()) ? map : nullptr; }

        /// Specular manifold NEE for point lights seen through a flat-faced glass body (null = off). It adds the caustic of
        /// those lights at every non-specular vertex, computed exactly; do not use it together with a caustic photon map.
        void setSpecularNee(const SpecularNee* nee) { m_specularNee = nee; }

        // The render loop
        virtual void render(const Scene& scene, Film& film) override {
            int width = film.width();
            int height = film.height();

#ifdef _OPENMP
            std::cout << "[PathIntegrator] OpenMP enabled with "
                << omp_get_max_threads() << " threads" << std::endl;
#else
            std::cout << "[PathIntegrator] WARNING: OpenMP NOT enabled!" << std::endl;
#endif

            std::cout << "[PathIntegrator] Rendering " << width << "x" << height
                << " (" << m_spp << " spp)" << std::endl;

            std::atomic<int> linesCompleted(0);
            const int totalLines = height;
            const bool wantAOVs = film.hasAOVs();

#pragma omp parallel for schedule(dynamic, 1)
            for (int j = 0; j < height; ++j) {
                // Progress display
                //std::cout << "\rScanlines remaining: " << (height - j) << " " << std::flush;

                for (int i = 0; i < width; ++i) {
                    Vector3 pixelColor(0.0);   // linear sRGB
                    Vector3 pixelAlbedo(0.0), pixelNormal(0.0);

                    for (int s = 0; s < m_spp; ++s) {
                        // Jitter for anti-aliasing
                        Real u = (Real(i) + sampling::Random()) / Real((width));
                        Real v = (Real(j) + sampling::Random()) / Real((height));

                        Point2 lensSample = sampling::Random2D();

                        Ray r = m_camera->getRay(u, v, lensSample);

                        // Each camera path carries its own set of wavelengths (spectral mode).
                        // Li() may drop the secondary wavelengths (dispersion), so convert with
                        // the wavelengths *after* Li() has run (separate statements on purpose).
                        SampledWavelengths lambda = samplePathWavelengths();
                        PathAOV aov;
                        const Spectrum Lpath = Li(r, scene, lambda, wantAOVs ? &aov : nullptr);
                        pixelColor += SpectrumToRGB(Lpath, lambda);
                        if (wantAOVs && aov.valid) {
                            pixelAlbedo += aov.albedo;
                            pixelNormal += aov.normal;
                        }
                    }
                    pixelColor /= Real(m_spp);
                    pixelAlbedo /= Real(m_spp);
                    pixelNormal /= Real(m_spp);

                    // NaN removal (for debugging)
                    if (math::hasNonFinite(pixelColor)) {
                        #pragma omp critical
                        {
                            std::cerr << "NaN detected at " << i << ", " << j << std::endl;
                        }
                        pixelColor = Vector3(0.0);
                    }

                    // Flip vertically and store
                    film.setPixel(i, height - 1 - j, pixelColor);
                    if (wantAOVs) film.setAOV(i, height - 1 - j, pixelAlbedo, pixelNormal);
                }

                // Update the progress when a row is done (thread-safe)
                int completed = ++linesCompleted;
                if (completed % 10 == 0 || completed == totalLines) {
                    #pragma omp critical
                    {
                        std::cout << "\rProgress: " << completed << "/" << totalLines
                            << " (" << (100 * completed / totalLines) << "%)"
                            << std::flush;
                    }
                }
            }

            std::cout << "\n[PathIntegrator] Done." << std::endl;
        }

        /**
         * @brief Compute radiance along a ray using path tracing.
         * @param lambda Wavelengths carried by this path (spectral mode); unused in RGB mode.
         *               In/out: a dispersive material drops the secondary wavelengths (pdf = 0),
         *               so pass the same object to SpectrumToRGB afterwards.
         * @param aov    If not null, receives the denoiser features of this path (albedo, normal).
         */
        Spectrum Li(Ray r, const Scene& scene, SampledWavelengths& lambda, PathAOV* aov = nullptr) const {
            Spectrum L(0.0);        // Accumulated Radiance
            Spectrum beta(1.0);     // Throughput

            Real lastPdf = 0;
            bool lastSpecular = false;
            bool hasLastBsdf = false;
            // True from a non-specular scattering vertex on until the next non-specular event, as long as only perfectly
            // specular bounces follow: light found now would be an L S+ D path, which the caustic photon map accounts for.
            bool fromReceiver = false;

            SurfaceInteraction lastRec;
            int passThrough = 0;
            int volumeEvents = 0;   // collisions inside participating media (not counted in depth)

            // Chromatic scattering media: the hero (the wavelength / channel whose extinction samples the distances) and the
            // per-wavelength path-density ratios rho_i = p_i / p_hero, kept normalized to mean 1 (see ChromaticMedium).
            // beta already includes the division by the mean of the ratios, so nothing else changes in this function.
            constexpr int kComponents =
#if defined(RAYT_SPECTRAL_MODE)
                kSpectralSamples;
#else
                3;
#endif
            int hero = -1;
            Spectrum rho(1.0);

            for (int depth = 0; depth < m_maxDepth; ++depth) {
                SurfaceInteraction rec;

                const bool hit = scene.hit(r, rec);

                // ========== Participating medium (random walk) ==========
                // The next collision is at the sampled free-flight distance; if the surface comes first the ray
                // reaches it (the transmittance cancels with the sampling density: the weight stays 1).
                // At a collision the path scatters with probability = albedo and is absorbed otherwise
                // (an unbiased analog walk), and no light is sampled there (the boundary hides the lights).
                if (r.medium && r.medium->attenuatesContinuously()) {
                    // Purely absorbing medium (clear water): the survivor fraction of the segment is a weight, not an event.
                    if (hit) beta *= r.medium->transmittance(rec.t, lambda);
                    else break;   // a ray inside a closed boundary that escapes: a leaky mesh, nothing sensible to add
                }
                else if (r.medium && r.medium->isChromatic()) {
                    // Hero-wavelength walk with null collisions (delta tracking).
                    //
                    // Distances are sampled with the extinction of the hero at the densest point, sigma_bar_h (the majorant). At
                    // a tentative collision x the local density q = med.density(x) decides what happens: scatter with
                    // probability q albedo_h, absorb with q (1 - albedo_h), and with 1 - q nothing (a null collision: go on in
                    // the same direction). q = 1 everywhere is the plain chromatic walk.
                    //
                    // Each wavelength i gets the weight a_i = (density of the same event if i were the hero) / (density of the
                    // event as sampled), for the three kinds of event of a path segment of length t:
                    //   scatter:       sigma_s,i exp(-sigma_bar_i t) / (sigma_s,h exp(-sigma_bar_h t))
                    //   null:          sigma_bar_i exp(-sigma_bar_i t) / (sigma_bar_h exp(-sigma_bar_h t))
                    //   surface found: exp(-sigma_bar_i t) / exp(-sigma_bar_h t)
                    // (q is the same for all wavelengths and cancels; the phase function is gray and sampled exactly.)
                    // Normalizing a by m = mean(rho a) is the balance heuristic over the choice of the hero.
                    const Medium& med = *r.medium;
                    if (hero < 0) {
#if defined(RAYT_SPECTRAL_MODE)
                        hero = 0;   // the wavelengths are a randomly rotated set: the first one is a uniform choice among them
#else
                        hero = std::min(2, int(sampling::Random() * Real(3)));   // one of the three channels, chosen once per path
#endif
                    }
                    const Spectrum sT = med.sigmaTSpectrum(lambda);
                    const Spectrum sS = med.sigmaSSpectrum(lambda);
                    const Real tSurface = hit ? rec.t : constants::INFINITY_VAL;
                    // Multiplies beta and rho by the weight of one event; false if the path cannot go on.
                    auto weigh = [&](Real t, int kind) {   // kind: 0 surface, 1 null, 2 scatter
                        Real e[kComponents], eMax = -constants::INFINITY_VAL;
                        for (int i = 0; i < kComponents; ++i) { e[i] = -(sT[i] - sT[hero]) * t; eMax = std::max(eMax, e[i]); }
                        Spectrum a(0.0);
                        Real m = 0;
                        for (int i = 0; i < kComponents; ++i) {   // exp(e - eMax): cannot overflow where the surface is far away
                            const Real k = kind == 2 ? sS[i] / sS[hero] : kind == 1 ? sT[i] / sT[hero] : Real(1);
                            a[i] = std::exp(e[i] - eMax) * k;
                            m += rho[i] * a[i];
                        }
                        m /= Real(kComponents);
                        if (!(m > Real(0)) || !std::isfinite(m)) return false;
                        for (int i = 0; i < kComponents; ++i) { beta[i] *= a[i] / m; rho[i] *= a[i] / m; }
                        return true;
                    };
                    enum { kSurface, kScattered, kEnded } outcome = kSurface;
                    Real tCur = 0, tScatter = 0;
                    for (;;) {
                        const Real tNext = tCur - std::log(Real(1) - sampling::Random()) / sT[hero];
                        if (tNext >= tSurface) {
                            if (!hit) { outcome = kEnded; break; }   // a ray that escapes a closed boundary: a leaky mesh
                            if (!weigh(tSurface - tCur, 0)) outcome = kEnded;
                            break;
                        }
                        if (++volumeEvents > m_maxVolumeEvents) { ++m_volumeCapHits; outcome = kEnded; break; }
                        const Real q = med.density(r.o + r.d * tNext);
                        const Real u = sampling::Random();
                        if (u >= q) {   // null collision
                            if (!weigh(tNext - tCur, 1)) { outcome = kEnded; break; }
                            tCur = tNext;
                            continue;
                        }
                        if (u * sT[hero] >= q * sS[hero]) { outcome = kEnded; break; }   // absorbed (probability 1 - albedo of the hero)
                        if (!weigh(tNext - tCur, 2)) { outcome = kEnded; break; }
                        tScatter = tNext;
                        outcome = kScattered;
                        break;
                    }
                    if (outcome == kEnded) break;
                    if (outcome == kScattered) {
                        const Point3 p = r.o + r.d * tScatter;
                        const Vector3 dir = med.samplePhase(r.d, Point2(sampling::Random(), sampling::Random()));
                        r = Ray(p, dir, Real(0), &med);
                        hasLastBsdf = true;
                        lastSpecular = true;
                        fromReceiver = false;
                        --depth;
                        continue;
                    }
                }
                else if (r.medium) {
                    const Medium& med = *r.medium;
                    const Real tSurface = hit ? rec.t : constants::INFINITY_VAL;
                    const Real tCollision = med.sampleDistance(sampling::Random());
                    if (tCollision < tSurface) {
                        if (++volumeEvents > m_maxVolumeEvents) { ++m_volumeCapHits; break; }
                        if (sampling::Random() >= med.albedo()) break;   // absorbed
                        const Point3 p = r.o + r.d * tCollision;
                        const Vector3 dir = med.samplePhase(r.d, Point2(sampling::Random(), sampling::Random()));
                        r = Ray(p, dir, Real(0), &med);   // p is inside the medium: nothing to skip (an epsilon would hide nearby surfaces)
                        hasLastBsdf = true;
                        lastSpecular = true;   // no NEE was done here: a light hit next is not MIS weighted
                        fromReceiver = false;  // photons are not traced through scattering media: what follows is counted here
                        --depth;               // a collision in a medium is not a surface bounce
                        continue;
                    }
                }
                if (hit) rec.medium = r.medium;   // whatever the ray hits is surrounded by the ray's medium

                if (!hit) {  // Intersection test: the ray missed everything

                    if (!hasLastBsdf && scene.hasCameraBackground()) {   // a camera ray that touched nothing: the backdrop
                        L += beta * RgbToIlluminantSpectrum(scene.cameraBackground(), lambda);
                        break;
                    }

                    if (scene.hasEnvLight()) {  // There is an environment light
                        Spectrum envL = scene.envLight()->Le(r, lambda);  // Evaluate the environment light

                        if (aov && !aov->valid) {   // seen directly (or through glass): the background color
                            aov->albedo = glm::clamp(SpectrumToRGB(envL, lambda), Vector3(0.0), Vector3(1.0));
                            aov->valid = true;
                        }

                        if (hasLastBsdf && !lastSpecular) {  // A BSDF was evaluated before and it was not a specular bounce
                            // MIS with BSDF sampling
                            // Density of reaching this light by NEE = P(select it) * pdfLi,
                            // evaluated for the same context the NEE sample at lastRec used.
                            LightSampleContext ctx{ lastRec.p, lastRec.gn };
                            Real pdfEnv = scene.envLight()->pdfLi(ctx, r.d);
                            if (scene.lightSampler()) pdfEnv *= scene.lightSampler()->pmf(ctx, scene.envLight());

                            Real w = 1.0;
                            if (pdfEnv > 0 && lastPdf > 0) {
                                w = math::powerHeuristic(lastPdf, pdfEnv);
                            }
                            L += beta * envL * w;
                        }
                        else if (!(m_caustics && fromReceiver && scene.envLight()->photonCaustics())) {
                            // Camera ray or specular path: no MIS (an L S+ D path is left to the caustic photons)
                            L += beta * envL;
                        }
                    }
                    break;
                }

                // Translucent layer: with probability 1 - opacity the ray continues straight through.
                {
                    const Real opacity = rec.matPtr->opacity(rec);
                    if (opacity < Real(1) && sampling::Random() >= opacity && ++passThrough < 64) {
                        r = rec.SpawnRay(r.d);
                        --depth;   // not a scattering event
                        continue;
                    }
                }

                // The BSDFs below assume unit-length directions (wo = -r.d).
                Assert(math::isNormalized(r.d) && "ray direction must be unit length");

                // Denoiser features: the first vertex that is not a perfect mirror / smooth glass.
                if (aov && !aov->valid && !rec.matPtr->isSpecular()) {
                    aov->albedo = rec.matPtr->aovAlbedo(rec);
                    aov->normal = rec.n;
                    aov->valid = true;
                }

                // Dispersive materials reduce the path to a single wavelength before anything is evaluated.
                rec.matPtr->prepareWavelengths(lambda);

                // ========== Emissive Surface Hit ==========
                // ※ wo = -r.direction
                Spectrum Le = rec.matPtr->emitted(rec, -r.d, lambda);
                if (!isBlack(Le)) {
                    if (hasLastBsdf && !lastSpecular && rec.areaLight) {
                        LightSampleContext ctx{ lastRec.p, lastRec.gn };
                        Real pdfLight = rec.areaLight->pdfLi(ctx, r.d);
                        if (scene.lightSampler()) pdfLight *= scene.lightSampler()->pmf(ctx, rec.areaLight);
                        Real w = 1.0;
                        if (pdfLight > 0 && lastPdf > 0) {
                            w = math::powerHeuristic(lastPdf, pdfLight);
                        }
                        L += beta * Le * w;
                    }
                    else if (!(m_caustics && fromReceiver && rec.areaLight && rec.areaLight->photonCaustics())) {
                        L += beta * Le;
                    }
                }

                // Remember this vertex for the MIS weight of the NEXT hit (emitter or environment).
                // It must happen AFTER the emissive-hit MIS above, which needs the PREVIOUS vertex.
                lastRec = rec;

                // ========== NEE + BSDF sampling ==========
                // A subsurface material may return a "subsurface" event instead of a direction: the path then
                // moves to an exit point (probe rays), and the same two steps run there with the exit lobe.
                std::optional<BSDFSample> bsdfSample;
                for (int stage = 0; ; ++stage) {
                    // ========== Next Event Estimation (NEE) ==========
                    if (!rec.matPtr->isSpecular()) {
                        L += beta * sampleAllLights(scene, rec, -r.d, lambda);
                        // Caustic light (L S+ D) at this vertex from the photon map (not inside a scattering medium: no photons there)
                        if (m_caustics && photonsAllowedIn(r.medium)) L += beta * m_caustics->radiance(rec, -r.d, lambda);
                        // The same caustic of point lights, solved exactly through a flat-faced glass body (manifold NEE)
                        if (m_specularNee && !r.medium) L += beta * m_specularNee->sample(scene, rec, -r.d, lambda);
                    }

                    // ========== BSDF Sampling for Next Path Vertex ==========
                    Point2 u(sampling::Random(), sampling::Random());
                    bsdfSample = rec.matPtr->sample(rec, -r.d, u, lambda);
                    if (!bsdfSample || !bsdfSample->isSubsurface() || stage > 0) break;

                    SurfaceInteraction exitVertex;
                    Spectrum exitWeight;
                    if (!SampleSubsurfaceExit(scene, rec, lambda, exitVertex, exitWeight)) { bsdfSample.reset(); break; }
                    beta *= bsdfSample->f * exitWeight;
                    rec = exitVertex;
                    lastRec = rec;   // later emitter hits are weighted against the exit vertex
                }

                if (!bsdfSample) {
                    break;
                }

                // 4. Update the throughput (beta)
                // Monte Carlo estimator: beta_new = beta_old * (f * cos_theta / pdf)

                Spectrum f = bsdfSample->f;
                Real pdf = bsdfSample->pdf;
                Vector3 wi = bsdfSample->wi; // the new direction
                Assert(math::isNormalized(wi) && "BSDF sample direction must be unit length");

                lastPdf = pdf;
                lastSpecular = bsdfSample->isSpecular();
                hasLastBsdf = true;
                // A non-delta bounce starts a chain that belongs to the photons; a delta lobe of a material that also has
                // diffuse lobes (a coat) is not a photon path (photons only continue through purely specular materials).
                if (!lastSpecular) { if (photonsAllowedIn(r.medium)) fromReceiver = true; }
                else if (!rec.matPtr->isSpecular()) fromReceiver = false;

                // Is this a specular (delta) bounce?
                if (bsdfSample->isSpecular()) {
                    // Specular case
                    // The PDF is conceptually infinite (Dirac delta), so
                    // f is taken to already hold the weight f / pdf, cancelled inside sample()
                    beta *= f;
                }
                else {
                    Real cosTheta = std::abs(glm::dot(rec.n, wi));
                    if (pdf > 1e-8f)
                        beta *= f * cosTheta / pdf;
                    else
                        break;
                }

                // Stop when the throughput has become 0
                if (isBlack(beta)) break;

                // ========== Russian Roulette ==========
                // Early bounces are never terminated:
                //   - they carry most of the energy and are cheap
                //   - terminating them randomly mostly adds variance
                // From rrStartDepth on, the path continues with a probability based on the
                // largest component of the throughput (beta). Survivors are reweighted by
                // 1 / survival probability, so the estimator stays unbiased.
                // This avoids the bias of a hard maxDepth cutoff while still ending
                // low-contribution paths early to save time.
                constexpr int rrStartDepth = 3;
                if (depth >= rrStartDepth) {
                    Real rrProb = std::clamp(math::maxComponent(beta), Real(0.05), Real(0.95));
                    if (sampling::Random() > rrProb) {
                        break; // terminate this path
                    }
                    beta /= rrProb; // compensate for the survival probability
                }

                // 5. Update the ray
                // r = Ray(rec.p + rec.n * constants::RAY_EPSILON, wi);  old
                // A closed boundary with a medium inside: the new direction goes to the inside or the outside.
                // gn faces the incoming ray, frontFace says that this side is the outside of the surface.
                if (const Medium* inner = rec.matPtr->mediumInside()) {
                    const bool wiOnRaySide = glm::dot(wi, rec.gn) > 0;
                    rec.medium = (wiOnRaySide == rec.frontFace) ? nullptr : inner;
                }
                if (rec.medium) {
                    // The new ray travels in a medium: start exactly on the surface (tMin keeps it from hitting that
                    // surface again) instead of offsetting the origin, which would skip the first collisions: the
                    // epsilon (10 um) is a visible fraction of the mean free path of skin-like media (~0.2 mm).
                    r = Ray(rec.p, wi, constants::RAY_EPSILON, rec.medium);
                    r.time = rec.time;
                }
                else {
                    r = rec.SpawnRay(wi);
                }
            }

            return L;
        }

        /**
         * @brief Progressive rendering: render one sample per pixel.
         */
        void renderOnePass(const Scene& scene, Film& film, int sampleIndex) {
            int width = film.width();
            int height = film.height();

            // Parallelization (OpenMP)
            #pragma omp parallel for schedule(dynamic, 1)
            for (int j = 0; j < height; ++j) {
                for (int i = 0; i < width; ++i) {

                    // 1. Take one sample
                    Real u = (Real(i) + sampling::Random()) / Real(width);
                    Real v = (Real(j) + sampling::Random()) / Real(height);

                    Point2 lensSample = sampling::Random2D();
                    Ray r = m_camera->getRay(u, v, lensSample);

                    // Each camera path carries its own set of wavelengths (spectral mode).
                    SampledWavelengths lambda = samplePathWavelengths();
                    const Spectrum Lpath = Li(r, scene, lambda);   // may drop secondary wavelengths
                    Vector3 L = SpectrumToRGB(Lpath, lambda);

                    // NaN check
                    if (math::hasNonFinite(L)) L = Vector3(0.0);

                    // 2. Add the sample to the film
                    // (Film::addPixel only sums; the normalization by 1/spp is done when the film is written.)

                    film.addPixel(i, height - 1 - j, L);
                }
            }
        }

    private:
        std::shared_ptr<Camera> m_camera;
        int m_maxDepth;
        int m_spp;
        int m_maxVolumeEvents = 100000;
        mutable std::atomic<long long> m_volumeCapHits{ 0 };
        const CausticPhotonMap* m_caustics = nullptr;
        const SpecularNee* m_specularNee = nullptr;

        /**
         * @brief Wavelengths for one camera path.
         *
         * Spectral mode: four stratified wavelengths importance-sampled by visual sensitivity.
         * RGB mode: a dummy set (never read, no random numbers consumed).
         */
        static SampledWavelengths samplePathWavelengths() {
#if defined(RAYT_SPECTRAL_MODE)
            return SampledWavelengths::sampleVisible(sampling::Random());
#else
            return SampledWavelengths();
#endif
        }

        /**
         * @brief Next event estimation: pick one light with the scene's LightSampler, sample it, MIS.
         *
         * The light is chosen with probability p (uniform / power / spatial importance, the
         * environment light included). The estimator and both MIS densities use the SAME combined
         * density of "choose this light, then sample this direction":
         *     pdf_light_strategy = p * light.pdfLi(wi)
         * The BSDF-sampling side (Li) uses sampler->pmf(ctx, light) * pdfLi for the same context.
         */
        Spectrum sampleAllLights(const Scene& scene,
                                 const SurfaceInteraction& rec,
                                 const Vector3& wo,
                                 const SampledWavelengths& lambda) const {
            Spectrum Ld(0.0);
            const LightSampler* sampler = scene.lightSampler();
            if (!sampler || !sampler->hasLights()) return Ld;

            LightSampleContext ctx{ rec.p, rec.gn };

            // ========== Light Selection ==========
            // An exhaustive sampler hands over every reachable light (probability 1 each).
            if (sampler->samplesAll()) {
                for (size_t i = 0; i < sampler->numLights(); ++i) {
                    const Light* light = sampler->lightAt(i);
                    const Real p = sampler->pmf(ctx, light);
                    if (p > Real(0)) Ld += estimateDirect(scene, rec, wo, lambda, ctx, light, p);
                }
                return Ld;
            }

            const auto sampled = sampler->sample(ctx, sampling::Random());
            if (!sampled || sampled->p <= Real(0)) return Ld;
            return estimateDirect(scene, rec, wo, lambda, ctx, sampled->light, sampled->p);
        }

        /**
         * @brief One NEE sample of `light`, which was selected with probability `pSelect`, with MIS.
         */
        Spectrum estimateDirect(const Scene& scene,
                                const SurfaceInteraction& rec,
                                const Vector3& wo,
                                const SampledWavelengths& lambda,
                                const LightSampleContext& ctx,
                                const Light* light,
                                Real pSelect) const {
            Spectrum Ld(0.0);

            // ========== Sample Chosen Light ==========
            Point2 uLight(sampling::Random(), sampling::Random());
            auto ls = light->sampleLi(ctx, uLight, lambda);

            if (!ls || ls->pdf <= 0 || isBlack(ls->Li)) {
                return Ld;
            }

            Vector3 wi = ls->wi;
            Assert(math::isNormalized(wi) && "light sample direction must be unit length");
            Spectrum Li = ls->Li;
            const Real pdfLight = ls->pdf * pSelect;   // density of this NEE sample, w.r.t. solid angle

            // ========== Visibility Test ==========
            if (!visible(scene, rec, wi, ls->tMax, ls->pLight)) {
                return Ld;
            }

            // ========== BSDF Evaluation ==========
            Spectrum f = rec.matPtr->eval(rec, wo, wi, lambda);
            if (isBlack(f)) return Ld;

            Real cosTheta = std::abs(glm::dot(rec.n, wi));

            // ========== MIS Weight ==========
            Real w = 1.0;
            if (!ls->isDelta) {
                Real pdfBsdf = rec.matPtr->pdf(rec, wo, wi, lambda);
                if (pdfBsdf > 0) {
                    w = math::powerHeuristic(pdfLight, pdfBsdf);
                }
            }

            // ========== Final Contribution ==========
            // Li already includes the geometry term from light sampling.
            Ld = f * Li * cosTheta * w / pdfLight;

            return Ld;
        }

        /**
         * @brief Check visibility between surface point and light.
         */
        static bool visible(const Scene& scene,
            const SurfaceInteraction& ref,
            const Vector3& wi,
            Real tMax,
            const Point3& pLight) {

            // Infinite light (environment): the ray only has to escape the scene.
            if (!std::isfinite(tMax)) {
                return !scene.occluded(ref.SpawnRay(wi));
            }

            // Finite light: aim the ray exactly at the sampled point from the *offset* origin and
            // stop a relative 1e-4 short of it. An absolute epsilon is not enough here: moving the
            // origin by RAY_EPSILON shifts the ray sideways, and where it meets the light's surface
            // at a shallow angle the hit point moves along the ray by much more than that. The
            // shadow ray then hits the light itself and every area-light sample is reported blocked.
            const Point3 origin = ref.OffsetRayOrigin(wi);
            const Vector3 toLight = pLight - origin;
            const Real dist = glm::length(toLight);
            if (dist <= Real(0)) return true;

            Ray shadow(origin, toLight / dist, Real(0), ref.medium);
            shadow.time = ref.time;
            shadow.tMax = dist * (Real(1) - Real(1e-4));
            return !scene.occluded(shadow);
        }

    };

} // namespace rayt 