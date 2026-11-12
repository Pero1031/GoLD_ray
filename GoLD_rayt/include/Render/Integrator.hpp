/**
 * @file Render/Integrator.hpp
 * @brief The path tracer: Integrator (interface) and PathIntegrator (unidirectional path tracing with next-event
 *        estimation, MIS, participating media, subsurface scattering and photon caustics).
 *
 * The implementation is in src/Render/Integrator.cpp; PathIntegrator::Li is the loop over path vertices and calls one
 * member function per step (medium, escaped ray, emission, scattering, throughput, next ray).
 */

#pragma once

#include <atomic>
#include <memory>
#include <optional>

#include "Core/Forward.hpp"
#include "Core/Types.hpp"
#include "Core/Ray.hpp"
#include "Core/Interaction.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Lights/Light.hpp"
#include "Materials/Material.hpp"

namespace rayt {

    class Scene;
    class Camera;
    class Film;
    class Medium;
    class CausticPhotonMap;
    class SpecularNee;

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
        void setCausticMap(const CausticPhotonMap* map);   // an empty map is ignored (a map that only traces media stores no photons)

        /// Light tracing took over the caustics on surfaces the camera sees directly (no bounce before): no photon gather there. The paths of the
        /// other receivers (seen through glass, after a diffuse bounce) still use the photons.
        void setSkipDirectCaustics(bool skip) { m_skipDirectCaustics = skip; }

        /**
         * @brief Light inside participating media: next event estimation at every scattering event, and shadow rays that cross null interfaces.
         *
         * Off (default): no light is sampled in a medium (the medium sits inside glass: only paths that find a light by chance count), and a
         * shadow ray is blocked by every surface. On: a medium that is bounded by a null interface (a fog) is lit directly. Both kinds of
         * next event estimation (at medium vertices and at surfaces) weigh the shadow ray by the transmittance of the media it crosses, a
         * phase-sampled ray that finds an emitter is MIS weighted against the estimate, and a camera path passes null interfaces without
         * counting a bounce. With a photon map that traces media (CausticPhotonParams::traceMedia) the light that a path finds after a specular
         * bounce that followed a medium vertex (L S+ M) is dropped: light tracing accounts for it.
         */
        void setMediaNee(bool on) { m_mediaNee = on; }

        /// Specular manifold NEE for point lights seen through a flat-faced glass body (null = off). It adds the caustic of
        /// those lights at every non-specular vertex, computed exactly; do not use it together with a caustic photon map.
        void setSpecularNee(const SpecularNee* nee) { m_specularNee = nee; }

        /// The render loop: spp samples per pixel, in parallel over the rows.
        void render(const Scene& scene, Film& film) override;

        /**
         * @brief Compute radiance along a ray using path tracing.
         * @param lambda Wavelengths carried by this path (spectral mode); unused in RGB mode.
         *               In/out: a dispersive material drops the secondary wavelengths (pdf = 0),
         *               so pass the same object to SpectrumToRGB afterwards.
         * @param aov    If not null, receives the denoiser features of this path (albedo, normal).
         */
        Spectrum Li(Ray r, const Scene& scene, SampledWavelengths& lambda, PathAOV* aov = nullptr) const;

        /// Progressive rendering: render one sample per pixel.
        void renderOnePass(const Scene& scene, Film& film, int sampleIndex);

    private:
        std::shared_ptr<Camera> m_camera;
        int m_maxDepth;
        int m_spp;
        int m_maxVolumeEvents = 100000;
        int m_renderCalls = 0;   // render() calls so far: with --seed every call draws its own random streams (--caustic-passes)
        mutable std::atomic<long long> m_volumeCapHits{ 0 };
        const CausticPhotonMap* m_caustics = nullptr;
        bool m_skipDirectCaustics = false;
        const SpecularNee* m_specularNee = nullptr;
        bool m_mediaNee = false;

        /// The photon map traces scattering media: a path's gather and its drop rule treat them like air.
        bool mediaTraced() const;

        /// Everything that one path carries from vertex to vertex.
        struct PathState {
            Spectrum L = Spectrum(0.0);       ///< accumulated radiance
            Spectrum beta = Spectrum(1.0);    ///< throughput (already divided by the mean of rho, see below)

            // The previous scattering event, for the MIS weight of an emitter that the next ray finds.
            Real lastPdf = 0;
            bool lastSpecular = false;
            bool hasLastBsdf = false;
            SurfaceInteraction lastRec;
            // True from a non-specular scattering vertex on until the next non-specular event, as long as only perfectly
            // specular bounces follow: light found now would be an L S+ D path, which the caustic photon map accounts for.
            bool fromReceiver = false;

            // Light tracing connects a medium vertex to the camera only when the segment between them is free of specular bounces, so the
            // camera path may drop the light it finds after a specular chain only if its FIRST receiver (a medium vertex or a non-specular
            // surface) was seen straight: viaSpecular = a specular bounce came before that receiver, sawReceiver = a receiver was met,
            // mediumLT = the first receiver was a medium vertex that was seen straight.
            bool viaSpecular = false, sawReceiver = false, mediumLT = false;

            int passThrough = 0;              ///< translucent layers passed (bounded)
            int volumeEvents = 0;             ///< collisions inside participating media (not counted in depth)

            // Scattering media (see PathIntegrator::advanceInMedium): the hero (the wavelength / channel whose extinction
            // samples the distances) and the per-wavelength path-density ratios rho_i = p_i / p_hero, kept normalized to mean 1.
            int hero = -1;
            Spectrum rho = Spectrum(1.0);
        };

        /// What a medium did to the ray in one step of the path.
        enum class MediumStep {
            Surface,     ///< the ray reaches the surface (or the next vertex) without a collision
            Scattered,   ///< a collision: the ray was replaced by one that starts at the collision point
            Ended        ///< the path is over (absorbed, a leaky boundary, the collision cap)
        };

        MediumStep advanceInMedium(Ray& r, const Scene& scene, bool hit, Real tSurface, PathState& st, const SampledWavelengths& lambda) const;

        /// Next event estimation at a scattering event of `med` at `p` (the path arrived travelling along `dirIn`): the light is chosen like at a surface.
        Spectrum sampleLightsInMedium(const Scene& scene, const Point3& p, const Vector3& dirIn, const Medium& med, const SampledWavelengths& lambda) const;
        Spectrum estimateDirectInMedium(const Scene& scene, const Point3& p, const Vector3& dirIn, const Medium& med, const SampledWavelengths& lambda,
            const LightSampleContext& ctx, const Light* light, Real pSelect) const;

        /// The ray left the scene: adds the backdrop or the environment light (with the MIS weight) and records the AOV.
        void addEscapedRadiance(const Ray& r, const Scene& scene, const SampledWavelengths& lambda, PathAOV* aov, PathState& st) const;

        /// The ray hit an emitter: adds its emission (with the MIS weight).
        void addSurfaceEmission(const Ray& r, const Scene& scene, const SurfaceInteraction& rec, const SampledWavelengths& lambda, PathState& st) const;

        /// Next event estimation (+ caustic light) at the vertex, then the BSDF sample for the next direction. A subsurface
        /// material may return a "subsurface" event instead of a direction: the path then moves to an exit point (probe
        /// rays) and the same two steps run there with the exit lobe, so `rec` can change. Null: the path ends.
        std::optional<BSDFSample> sampleScatter(const Ray& r, const Scene& scene, SurfaceInteraction& rec, const SampledWavelengths& lambda, PathState& st) const;

        /// Multiplies the throughput by the sample's weight, then Russian roulette; false: the path ends.
        bool updateThroughput(const SurfaceInteraction& rec, const BSDFSample& sample, const Medium* rayMedium, int depth, PathState& st) const;

        /// The ray that leaves `rec` in direction wi; picks the medium on the other side of a boundary that has one inside.
        static Ray nextRay(SurfaceInteraction& rec, const Vector3& wi);

        /**
         * @brief Wavelengths for one camera path.
         *
         * Spectral mode: four stratified wavelengths importance-sampled by visual sensitivity.
         * RGB mode: a dummy set (never read, no random numbers consumed).
         */
        static SampledWavelengths samplePathWavelengths();

        /**
         * @brief Next event estimation: pick one light with the scene's LightSampler, sample it, MIS.
         *
         * The light is chosen with probability p (uniform / power / spatial importance, the
         * environment light included). The estimator and both MIS densities use the SAME combined
         * density of "choose this light, then sample this direction":
         *     pdf_light_strategy = p * light.pdfLi(wi)
         * The BSDF-sampling side (Li) uses sampler->pmf(ctx, light) * pdfLi for the same context.
         */
        Spectrum sampleAllLights(const Scene& scene, const SurfaceInteraction& rec, const Vector3& wo, const SampledWavelengths& lambda) const;

        /// One NEE sample of `light`, which was selected with probability `pSelect`, with MIS.
        Spectrum estimateDirect(const Scene& scene, const SurfaceInteraction& rec, const Vector3& wo, const SampledWavelengths& lambda,
            const LightSampleContext& ctx, const Light* light, Real pSelect) const;

        /// Check visibility between surface point and light.
        static bool visible(const Scene& scene, const SurfaceInteraction& ref, const Vector3& wi, Real tMax, const Point3& pLight);
    };

} // namespace rayt
