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
 *
 * The tracing, the grid and the gathering are in src/Render/CausticPhotonMap.cpp.
 */

#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "Core/Types.hpp"
#include "Core/Interaction.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Materials/Medium.hpp"

namespace rayt {

    class Scene;
    class Camera;
    class LightSampler;
    class Material;

    /// Photons are traced and gathered in air and inside purely absorbing media (their transmittance is a weight). A scattering medium is
    /// traced only when the map traces media (CausticPhotonParams::traceMedia) and the medium is gray (the same for every wavelength).
    inline bool photonsAllowedIn(const Medium* medium, bool mediaTraced = false) {
        return !medium || medium->attenuatesContinuously() || (mediaTraced && medium->isSpectrallyConstant());
    }

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
        /// Do not store a photon whose chain is exactly two transmissions through these materials (a T T chain of a body that the specular
        /// manifold NEE solves exactly: SpecularNee::materials()), so that the two estimators do not count it twice. Empty = store them.
        std::vector<const Material*> skipTransmissionPairs;

        /// Light tracing (--caustic-trace): when set, every stored photon (L S+ D) is also connected to this camera, and its contribution is
        /// splatted into an image of traceWidth x traceHeight pixels (CausticPhotonMap::splat()). The camera is treated as a pinhole.
        const Camera* lightTraceCamera = nullptr;
        int traceWidth = 0, traceHeight = 0;

        /// Volumetric caustics (needs lightTraceCamera): a photon goes on through a gray scattering medium (delta tracking), and every real
        /// scattering event after at least one specular bounce (L S+ M) is connected to the camera. Null interfaces (a fog boundary) are
        /// crossed without counting as a bounce, and the connections cross them too. A photon that scatters before its first specular bounce
        /// is dropped: light that reaches the medium directly is left to the camera paths (PathIntegrator::setMediaNee).
        bool traceMedia = false;
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
        /// Light-traced caustic image (average radiance per pixel, linear RGB), empty unless CausticPhotonParams::lightTraceCamera was set.
        const std::vector<Vector3>& splat() const { return m_splat; }
        /// True if the photons were traced through scattering media (see CausticPhotonParams::traceMedia).
        bool tracesMedia() const { return m_tracesMedia; }
        size_t memoryBytes() const { return m_photons.capacity() * sizeof(Photon) + m_start.capacity() * sizeof(uint32_t); }

        /// Traces the photons (OpenMP) and builds the grid.
        void build(const Scene& scene, const CausticPhotonParams& P);

        /**
         * @brief Reflected caustic radiance toward `wo` at a non-specular surface point (light-like Spectrum, same units as Li()).
         * @param rec  the hit (position, geometric normal, material)
         * @param lambda the camera path's wavelengths
         */
        Spectrum radiance(const SurfaceInteraction& rec, const Vector3& wo, const SampledWavelengths& lambda) const;

    private:
        std::vector<Photon> m_photons;   // sorted by grid bucket
        std::vector<uint32_t> m_start;   // bucket b holds photons [m_start[b], m_start[b + 1])
        std::vector<Vector3> m_splat;    // light-traced image (see CausticPhotonParams::lightTraceCamera)
        uint32_t m_mask = 0;
        bool m_tracesMedia = false;
        Real m_radius = Real(0.0015);

        static Real luminance(const Vector3& c) { return Real(0.2126) * c.x + Real(0.7152) * c.y + Real(0.0722) * c.z; }
        int cellIndex(Real x) const;
        static uint32_t hashCell(int x, int y, int z);
        void buildGrid();
        static void connectToCamera(const Scene& scene, const CausticPhotonParams& P, const SurfaceInteraction& rec, const Vector3& wiLight,
            const Spectrum& flux, const SampledWavelengths& lambda, std::vector<Vector3>& splat);
        static void connectMediumToCamera(const Scene& scene, const CausticPhotonParams& P, const Point3& x, const Vector3& dirIn, const Medium& med,
            const Spectrum& flux, const SampledWavelengths& lambda, std::vector<Vector3>& splat);
        static void tracePhoton(const Scene& scene, const CausticPhotonParams& P, const LightSampler& sampler, Real invN,
            std::vector<Photon>& out, std::vector<Vector3>& splat);
    };

} // namespace rayt
