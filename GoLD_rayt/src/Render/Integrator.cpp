/**
 * @file src/Render/Integrator.cpp
 * @brief PathIntegrator: the render loop, the path loop (Li) and the steps it is made of.
 *
 * Li walks the vertices of one camera path. Each iteration asks one member function for one step:
 *   advanceInMedium      the ray travels in a medium: collisions, null collisions, absorption
 *   addEscapedRadiance   the ray left the scene: backdrop / environment light
 *   addSurfaceEmission   the ray hit an emitter
 *   sampleScatter        next event estimation, caustic light, the BSDF sample (and subsurface exit points)
 *   updateThroughput     the sample's weight and Russian roulette
 *   nextRay              the ray for the next iteration (and the medium it travels in)
 */

#include "pch.h"

#include "Render/Integrator.hpp"

namespace rayt {

    namespace {

#if defined(RAYT_SPECTRAL_MODE)
        constexpr int kComponents = kSpectralSamples;   // wavelengths per path
#else
        constexpr int kComponents = 3;                  // RGB channels
#endif

        /// The hero of a path in a medium whose coefficients depend on the wavelength.
        int pickHero() {
#if defined(RAYT_SPECTRAL_MODE)
            return 0;   // the wavelengths are a randomly rotated set: the first one is a uniform choice among them
#else
            return std::min(2, int(sampling::Random() * Real(3)));   // one of the three channels, chosen once per path
#endif
        }

    } // namespace

    // ------------------------------------------------------------------------------------------------
    // Render loops
    // ------------------------------------------------------------------------------------------------

    void PathIntegrator::render(const Scene& scene, Film& film) {
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
        const std::uint64_t renderCall = static_cast<std::uint64_t>(m_renderCalls++);

#pragma omp parallel for schedule(dynamic, 1)
        for (int j = 0; j < height; ++j) {
            // With --seed a row always draws the same random numbers, whichever thread renders it and in whatever order.
            sampling::SeedStream(1, renderCall, static_cast<std::uint64_t>(j));
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

    void PathIntegrator::renderOnePass(const Scene& scene, Film& film, int sampleIndex) {
        int width = film.width();
        int height = film.height();

        // Parallelization (OpenMP)
#pragma omp parallel for schedule(dynamic, 1)
        for (int j = 0; j < height; ++j) {
            sampling::SeedStream(2, static_cast<std::uint64_t>(sampleIndex), static_cast<std::uint64_t>(j));   // no-op without a run seed
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

    SampledWavelengths PathIntegrator::samplePathWavelengths() {
#if defined(RAYT_SPECTRAL_MODE)
        return SampledWavelengths::sampleVisible(sampling::Random());
#else
        return SampledWavelengths();
#endif
    }

    // ------------------------------------------------------------------------------------------------
    // The path
    // ------------------------------------------------------------------------------------------------

    Spectrum PathIntegrator::Li(Ray r, const Scene& scene, SampledWavelengths& lambda, PathAOV* aov) const {
        PathState st;

        for (int depth = 0; depth < m_maxDepth; ++depth) {
            SurfaceInteraction rec;

            const bool hit = scene.hit(r, rec);

            // ========== A captured cloud of Gaussians (docs/CAPTURE.md) ==========
            // Emits what the capture saw and absorbs; a backdrop that paths fly through (outside any medium).
            if (scene.splats() && !r.medium) {
                Vector3 rgb;
                Real T;
                scene.splats()->gather(r, hit ? rec.t : constants::INFINITY_VAL, rgb, T);
                if (T < Real(1)) {
                    st.L += st.beta * RgbToIlluminantSpectrum(rgb, lambda);
                    if (aov && !aov->valid && T < Real(0.5)) {   // mostly the cloud: its color is the denoiser's albedo
                        aov->albedo = glm::clamp(rgb / std::max(Real(1) - T, Real(1e-3)), Vector3(0.0), Vector3(1.0));
                        aov->valid = true;
                    }
                    st.beta *= Spectrum(T);
                    if (isBlack(st.beta)) break;
                }
            }

            // ========== Participating medium ==========
            if (r.medium) {
                const MediumStep step = advanceInMedium(r, scene, hit, hit ? rec.t : constants::INFINITY_VAL, st, lambda);
                if (step == MediumStep::Ended) break;
                if (step == MediumStep::Scattered) {
                    --depth;   // a collision in a medium is not a surface bounce
                    continue;
                }
            }
            if (hit) rec.medium = r.medium;   // whatever the ray hits is surrounded by the ray's medium

            // ========== The ray left the scene ==========
            if (!hit) {
                addEscapedRadiance(r, scene, lambda, aov, st);
                break;
            }

            // A boundary with no optical effect (the surface of a fog): cross it, switching the medium. Not a bounce.
            if (m_mediaNee && rec.matPtr->isNullInterface() && ++st.passThrough < 64) {
                r = nextRay(rec, r.d);
                --depth;
                continue;
            }

            // Translucent layer: with probability 1 - opacity the ray continues straight through.
            {
                const Real opacity = rec.matPtr->opacity(rec);
                if (opacity < Real(1) && sampling::Random() >= opacity && ++st.passThrough < 64) {
                    r = rec.SpawnRay(r.d);
                    --depth;   // not a scattering event
                    continue;
                }
            }

            // The BSDFs below assume unit-length directions (wo = -r.d).
            Assert(math::isNormalized(r.d) && "ray direction must be unit length");

            // A micro-relief (skin pores) bends the shading normal; the denoiser features below see it too.
            rec.matPtr->perturbShading(rec, -r.d);

            // Denoiser features: the first vertex that is not a perfect mirror / smooth glass.
            if (aov && !aov->valid && !rec.matPtr->isSpecular()) {
                aov->albedo = rec.matPtr->aovAlbedo(rec);
                aov->normal = rec.n;
                aov->valid = true;
            }

            // Dispersive materials reduce the path to a single wavelength before anything is evaluated.
            rec.matPtr->prepareWavelengths(lambda);

            addSurfaceEmission(r, scene, rec, lambda, st);

            // Remember this vertex for the MIS weight of the NEXT hit (emitter or environment).
            // It must happen AFTER the emissive-hit MIS above, which needs the PREVIOUS vertex.
            st.lastRec = rec;

            const std::optional<BSDFSample> bsdfSample = sampleScatter(r, scene, rec, lambda, st);
            if (!bsdfSample) break;
            if (!updateThroughput(rec, *bsdfSample, r.medium, depth, st)) break;

            r = nextRay(rec, bsdfSample->wi);
        }

        return st.L;
    }

    // ------------------------------------------------------------------------------------------------
    // Participating media
    // ------------------------------------------------------------------------------------------------

    PathIntegrator::MediumStep PathIntegrator::advanceInMedium(Ray& r, const Scene& scene, bool hit, Real tSurface, PathState& st,
        const SampledWavelengths& lambda) const
    {
        const Medium& med = *r.medium;

        // Purely absorbing medium (clear water): the survivor fraction of the segment is a weight, not an event.
        if (med.attenuatesContinuously()) {
            if (!hit) return MediumStep::Ended;   // a ray inside a closed boundary that escapes: a leaky mesh, nothing sensible to add
            st.beta *= med.transmittance(tSurface, lambda);
            return MediumStep::Surface;
        }

        // Scattering medium: a hero-wavelength walk with null collisions (delta tracking).
        //
        // Distances are sampled with the extinction of the hero at the densest point, sigma_bar_h (the majorant). At
        // a tentative collision x the local density q = med.density(x) decides what happens: scatter with
        // probability q albedo_h, absorb with q (1 - albedo_h), and with 1 - q nothing (a null collision: go on in
        // the same direction). q = 1 everywhere is a homogeneous medium: every tentative collision is real.
        //
        // Each wavelength i gets the weight a_i = (density of the same event if i were the hero) / (density of the
        // event as sampled), for the three kinds of event of a path segment of length t:
        //   scatter:       sigma_s,i exp(-sigma_bar_i t) / (sigma_s,h exp(-sigma_bar_h t))
        //   null:          sigma_bar_i exp(-sigma_bar_i t) / (sigma_bar_h exp(-sigma_bar_h t))
        //   surface found: exp(-sigma_bar_i t) / exp(-sigma_bar_h t)
        // (q is the same for all wavelengths and cancels; the phase function is gray and sampled exactly.)
        // Normalizing a by m = mean(rho a) is the balance heuristic over the choice of the hero.
        //
        // A medium whose coefficients are the same for every wavelength has a = 1 for all events: no hero, no weights.
        // No light is sampled inside a medium (the boundary hides the lights).
        const bool constant = med.isSpectrallyConstant();
        if (!constant && st.hero < 0) st.hero = pickHero();   // chosen once per path, when the first chromatic medium is met
        const int hero = constant ? 0 : st.hero;

        const Spectrum sT = med.sigmaTSpectrum(lambda);
        const Spectrum sS = med.sigmaSSpectrum(lambda);

        // Multiplies beta and rho by the weight of one event; false if the path cannot go on.
        auto weigh = [&](Real t, int kind) {   // kind: 0 surface, 1 null, 2 scatter
            if (constant) return true;
            Real e[kComponents], eMax = -constants::INFINITY_VAL;
            for (int i = 0; i < kComponents; ++i) { e[i] = -(sT[i] - sT[hero]) * t; eMax = std::max(eMax, e[i]); }
            Spectrum a(0.0);
            Real m = 0;
            for (int i = 0; i < kComponents; ++i) {   // exp(e - eMax): cannot overflow where the surface is far away
                const Real k = kind == 2 ? sS[i] / sS[hero] : kind == 1 ? sT[i] / sT[hero] : Real(1);
                a[i] = std::exp(e[i] - eMax) * k;
                m += st.rho[i] * a[i];
            }
            m /= Real(kComponents);
            if (!(m > Real(0)) || !std::isfinite(m)) return false;
            for (int i = 0; i < kComponents; ++i) { st.beta[i] *= a[i] / m; st.rho[i] *= a[i] / m; }
            return true;
        };

        Real tCur = 0;
        for (;;) {
            const Real tNext = tCur - std::log(Real(1) - sampling::Random()) / sT[hero];
            if (tNext >= tSurface) {
                if (!hit) return MediumStep::Ended;   // a ray that escapes a closed boundary: a leaky mesh
                return weigh(tSurface - tCur, 0) ? MediumStep::Surface : MediumStep::Ended;
            }
            if (++st.volumeEvents > m_maxVolumeEvents) { ++m_volumeCapHits; return MediumStep::Ended; }

            const Real q = med.density(r.o + r.d * tNext);
            const Real u = sampling::Random();
            if (u >= q) {   // null collision
                if (!weigh(tNext - tCur, 1)) return MediumStep::Ended;
                tCur = tNext;
                continue;
            }
            if (u * sT[hero] >= q * sS[hero]) return MediumStep::Ended;   // absorbed (probability 1 - albedo of the hero)
            if (!weigh(tNext - tCur, 2)) return MediumStep::Ended;

            // A real scattering collision at tNext.
            const Point3 p = r.o + r.d * tNext;
            if (m_mediaNee) st.L += st.beta * sampleLightsInMedium(scene, p, r.d, med, lambda);
            const Vector3 dir = med.samplePhase(r.d, Point2(sampling::Random(), sampling::Random()));
            const Real cosScatter = glm::dot(r.d, dir);
            r = Ray(p, dir, Real(0), &med);   // p is inside the medium: nothing to skip (an epsilon would hide nearby surfaces)
            st.hasLastBsdf = true;
            if (m_mediaNee) {
                // Light was sampled here: an emitter that the phase-sampled ray finds next is weighed against it (the phase function is the
                // density of the direction). A reference point without a normal: gn = 0.
                st.lastSpecular = false;
                st.lastPdf = med.phase(cosScatter);
                st.lastRec = SurfaceInteraction();
                st.lastRec.p = p;
                st.lastRec.gn = Vector3(0.0);
                st.lastRec.n = Vector3(0.0);
                st.lastRec.medium = &med;
                // Light tracing accounts for the specular chains that start at a medium vertex that the camera sees straight. A vertex behind
                // glass cannot be connected to the camera: its chains stay with the camera path.
                if (!st.sawReceiver && !st.viaSpecular) st.mediumLT = true;
                st.sawReceiver = true;
                st.fromReceiver = mediaTraced() && (st.mediumLT || st.fromReceiver);
            }
            else {
                st.lastSpecular = true;    // no NEE was done here: a light hit next is not MIS weighted
                st.fromReceiver = false;   // photons are not traced through scattering media: what follows is counted here
            }
            return MediumStep::Scattered;
        }
    }

    // ------------------------------------------------------------------------------------------------
    // Light found by the path
    // ------------------------------------------------------------------------------------------------

    void PathIntegrator::addEscapedRadiance(const Ray& r, const Scene& scene, const SampledWavelengths& lambda,
        PathAOV* aov, PathState& st) const
    {
        if (!st.hasLastBsdf && scene.hasCameraBackground()) {   // a camera ray that touched nothing: the backdrop
            st.L += st.beta * RgbToIlluminantSpectrum(scene.cameraBackground(), lambda);
            return;
        }

        if (!scene.hasEnvLight()) return;

        Spectrum envL = scene.envLight()->Le(r, lambda);

        if (aov && !aov->valid) {   // seen directly (or through glass): the background color
            aov->albedo = glm::clamp(SpectrumToRGB(envL, lambda), Vector3(0.0), Vector3(1.0));
            aov->valid = true;
        }

        if (st.hasLastBsdf && !st.lastSpecular) {   // A BSDF was evaluated before and it was not a specular bounce
            // MIS with BSDF sampling
            // Density of reaching this light by NEE = P(select it) * pdfLi,
            // evaluated for the same context the NEE sample at lastRec used.
            LightSampleContext ctx{ st.lastRec.p, st.lastRec.gn };
            Real pdfEnv = scene.envLight()->pdfLi(ctx, r.d);
            if (scene.lightSampler()) pdfEnv *= scene.lightSampler()->pmf(ctx, scene.envLight());

            Real w = 1.0;
            if (pdfEnv > 0 && st.lastPdf > 0) {
                w = math::powerHeuristic(st.lastPdf, pdfEnv);
            }
            st.L += st.beta * envL * w;
        }
        else if (!(m_caustics && st.fromReceiver && scene.envLight()->photonCaustics())) {
            // Camera ray or specular path: no MIS (an L S+ D path is left to the caustic photons)
            st.L += st.beta * envL;
        }
    }

    void PathIntegrator::addSurfaceEmission(const Ray& r, const Scene& scene, const SurfaceInteraction& rec,
        const SampledWavelengths& lambda, PathState& st) const
    {
        // wo = -r.d
        const Spectrum Le = rec.matPtr->emitted(rec, -r.d, lambda);
        if (isBlack(Le)) return;

        if (st.hasLastBsdf && !st.lastSpecular && rec.areaLight) {
            LightSampleContext ctx{ st.lastRec.p, st.lastRec.gn };
            Real pdfLight = rec.areaLight->pdfLi(ctx, r.d);
            if (scene.lightSampler()) pdfLight *= scene.lightSampler()->pmf(ctx, rec.areaLight);
            Real w = 1.0;
            if (pdfLight > 0 && st.lastPdf > 0) {
                w = math::powerHeuristic(st.lastPdf, pdfLight);
            }
            st.L += st.beta * Le * w;
        }
        else if (!(m_caustics && st.fromReceiver && rec.areaLight && rec.areaLight->photonCaustics())) {
            st.L += st.beta * Le;
        }
    }

    // ------------------------------------------------------------------------------------------------
    // Scattering at a surface
    // ------------------------------------------------------------------------------------------------

    std::optional<BSDFSample> PathIntegrator::sampleScatter(const Ray& r, const Scene& scene, SurfaceInteraction& rec,
        const SampledWavelengths& lambda, PathState& st) const
    {
        std::optional<BSDFSample> bsdfSample;
        for (int stage = 0; ; ++stage) {
            // ========== Next Event Estimation (NEE) ==========
            if (!rec.matPtr->isSpecular()) {
                st.L += st.beta * sampleAllLights(scene, rec, -r.d, lambda);
                // Caustic light (L S+ D) at this vertex from the photon map (not inside a scattering medium: no photons there)
                if (m_caustics && photonsAllowedIn(r.medium, mediaTraced()) && !(m_skipDirectCaustics && !st.hasLastBsdf)) st.L += st.beta * m_caustics->radiance(rec, -r.d, lambda);
                // The same caustic of point lights, solved exactly through a flat-faced glass body (manifold NEE)
                if (m_specularNee && !r.medium) st.L += st.beta * m_specularNee->sample(scene, rec, -r.d, lambda);
            }

            // ========== BSDF Sampling for Next Path Vertex ==========
            Point2 u(sampling::Random(), sampling::Random());
            bsdfSample = rec.matPtr->sample(rec, -r.d, u, lambda);
            if (!bsdfSample || !bsdfSample->isSubsurface() || stage > 0) break;

            SurfaceInteraction exitVertex;
            Spectrum exitWeight;
            if (!SampleSubsurfaceExit(scene, rec, lambda, exitVertex, exitWeight)) { bsdfSample.reset(); break; }
            st.beta *= bsdfSample->f * exitWeight;
            rec = exitVertex;
            st.lastRec = rec;   // later emitter hits are weighted against the exit vertex
        }
        return bsdfSample;
    }

    bool PathIntegrator::updateThroughput(const SurfaceInteraction& rec, const BSDFSample& sample, const Medium* rayMedium,
        int depth, PathState& st) const
    {
        // Monte Carlo estimator: beta_new = beta_old * (f * cos_theta / pdf)
        const Spectrum f = sample.f;
        const Real pdf = sample.pdf;
        const Vector3 wi = sample.wi;   // the new direction
        Assert(math::isNormalized(wi) && "BSDF sample direction must be unit length");

        if (sample.isSpecular()) { if (!st.sawReceiver) st.viaSpecular = true; }
        else st.sawReceiver = true;

        st.lastPdf = pdf;
        st.lastSpecular = sample.isSpecular();
        st.hasLastBsdf = true;
        // A non-delta bounce starts a chain that belongs to the photons; a delta lobe of a material that also has
        // diffuse lobes (a coat) is not a photon path (photons only continue through purely specular materials).
        if (!st.lastSpecular) { if (photonsAllowedIn(rayMedium, mediaTraced())) st.fromReceiver = true; }
        else if (!rec.matPtr->isSpecular()) st.fromReceiver = false;

        // Is this a specular (delta) bounce?
        if (sample.isSpecular()) {
            // The PDF is conceptually infinite (Dirac delta), so
            // f is taken to already hold the weight f / pdf, cancelled inside sample()
            st.beta *= f;
        }
        else {
            const Real cosTheta = std::abs(glm::dot(rec.n, wi));
            if (pdf > 1e-8f)
                st.beta *= f * cosTheta / pdf;
            else
                return false;
        }

        // Stop when the throughput has become 0
        if (isBlack(st.beta)) return false;

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
            const Real rrProb = std::clamp(math::maxComponent(st.beta), Real(0.05), Real(0.95));
            if (sampling::Random() > rrProb) return false;   // terminate this path
            st.beta /= rrProb;                               // compensate for the survival probability
        }
        return true;
    }

    Ray PathIntegrator::nextRay(SurfaceInteraction& rec, const Vector3& wi) {
        // A closed boundary with a medium inside: the new direction goes to the inside or the outside.
        // gn faces the incoming ray, frontFace says that this side is the outside of the surface.
        if (rec.matPtr->isMediumBoundary()) {
            const bool wiOnRaySide = glm::dot(wi, rec.gn) > 0;
            rec.medium = (wiOnRaySide == rec.frontFace) ? rec.matPtr->mediumOutside() : rec.matPtr->mediumInside();
        }
        if (rec.medium) {
            // The new ray travels in a medium: start exactly on the surface (tMin keeps it from hitting that
            // surface again) instead of offsetting the origin, which would skip the first collisions: the
            // epsilon (10 um) is a visible fraction of the mean free path of skin-like media (~0.2 mm).
            Ray next(rec.p, wi, constants::RAY_EPSILON, rec.medium);
            next.time = rec.time;
            return next;
        }
        return rec.SpawnRay(wi);
    }

    // ------------------------------------------------------------------------------------------------
    // Next event estimation
    // ------------------------------------------------------------------------------------------------

    Spectrum PathIntegrator::sampleAllLights(const Scene& scene, const SurfaceInteraction& rec, const Vector3& wo,
        const SampledWavelengths& lambda) const
    {
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

    Spectrum PathIntegrator::estimateDirect(const Scene& scene, const SurfaceInteraction& rec, const Vector3& wo,
        const SampledWavelengths& lambda, const LightSampleContext& ctx, const Light* light, Real pSelect) const
    {
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
        Spectrum transmittance(1.0);
        if (m_mediaNee || scene.splats()) {   // null interfaces are crossed, media (and a captured cloud) attenuate
            const Point3 origin = rec.OffsetRayOrigin(wi);
            Vector3 dirS = wi;
            Real tMaxS = constants::INFINITY_VAL;
            if (std::isfinite(ls->tMax)) {
                const Vector3 toLight = ls->pLight - origin;
                const Real dist = glm::length(toLight);
                if (dist > Real(0)) { dirS = toLight / dist; tMaxS = dist * (Real(1) - Real(1e-4)); }
            }
            // A boundary between two media: the shadow ray starts in the medium on the side it leaves to (as nextRay() decides for a path).
            const Medium* startMedium = rec.medium;
            if (rec.matPtr->isMediumBoundary())
                startMedium = ((glm::dot(wi, rec.gn) > 0) == rec.frontFace) ? rec.matPtr->mediumOutside() : rec.matPtr->mediumInside();
            transmittance = shadowTransmittance(scene, origin, startMedium, dirS, tMaxS, lambda, rec.time);
            if (isBlack(transmittance)) return Ld;
        }
        else if (!visible(scene, rec, wi, ls->tMax, ls->pLight)) {
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
        Ld = f * Li * transmittance * (cosTheta * w / pdfLight);

        return Ld;
    }

    // ------------------------------------------------------------------------------------------------
    // Next event estimation in a medium
    // ------------------------------------------------------------------------------------------------

    Spectrum PathIntegrator::sampleLightsInMedium(const Scene& scene, const Point3& p, const Vector3& dirIn, const Medium& med,
        const SampledWavelengths& lambda) const
    {
        Spectrum Ld(0.0);
        const LightSampler* sampler = scene.lightSampler();
        if (!sampler || !sampler->hasLights()) return Ld;

        LightSampleContext ctx;
        ctx.p = p;
        ctx.gn = Vector3(0.0);   // a point in a medium has no orientation
        ctx.medium = &med;

        if (sampler->samplesAll()) {
            for (size_t i = 0; i < sampler->numLights(); ++i) {
                const Light* light = sampler->lightAt(i);
                const Real pm = sampler->pmf(ctx, light);
                if (pm > Real(0)) Ld += estimateDirectInMedium(scene, p, dirIn, med, lambda, ctx, light, pm);
            }
            return Ld;
        }
        const auto sampled = sampler->sample(ctx, sampling::Random());
        if (!sampled || sampled->p <= Real(0)) return Ld;
        return estimateDirectInMedium(scene, p, dirIn, med, lambda, ctx, sampled->light, sampled->p);
    }

    Spectrum PathIntegrator::estimateDirectInMedium(const Scene& scene, const Point3& p, const Vector3& dirIn, const Medium& med,
        const SampledWavelengths& lambda, const LightSampleContext& ctx, const Light* light, Real pSelect) const
    {
        auto ls = light->sampleLi(ctx, Point2(sampling::Random(), sampling::Random()), lambda);
        if (!ls || ls->pdf <= 0 || isBlack(ls->Li)) return Spectrum(0.0);

        const Vector3 wi = ls->wi;
        const Real pdfLight = ls->pdf * pSelect;

        // The light travels along -wi into the collision and leaves along -dirIn (toward the camera): the angle between those two
        // travel directions has the cosine dot(wi, dirIn).
        const Real ph = med.phase(glm::dot(wi, dirIn));
        if (ph <= Real(0)) return Spectrum(0.0);

        Vector3 dirS = wi;
        Real tMaxS = constants::INFINITY_VAL;
        if (std::isfinite(ls->tMax)) {
            const Vector3 toLight = ls->pLight - p;
            const Real dist = glm::length(toLight);
            if (dist > Real(0)) { dirS = toLight / dist; tMaxS = dist * (Real(1) - Real(1e-4)); }
        }
        const Spectrum T = shadowTransmittance(scene, p, &med, dirS, tMaxS, lambda);
        if (isBlack(T)) return Spectrum(0.0);

        Real w = 1.0;
        if (!ls->isDelta && pdfLight > 0) w = math::powerHeuristic(pdfLight, ph);   // the phase-sampled ray has the density ph

        return ls->Li * T * (ph * w / pdfLight);
    }

    bool PathIntegrator::visible(const Scene& scene, const SurfaceInteraction& ref, const Vector3& wi, Real tMax,
        const Point3& pLight)
    {
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

} // namespace rayt
