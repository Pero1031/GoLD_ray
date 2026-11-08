/**
 * @file Render/Transmittance.hpp
 * @brief Visibility that knows about participating media: the transmittance of a shadow ray.
 *
 * A shadow ray normally is blocked or not (Scene::occluded). Inside a fog the light is also attenuated along the ray, and the fog
 * is bounded by a closed surface with no optical effect of its own (Material::isNullInterface): the shadow ray crosses such a
 * surface, switches the medium, and is blocked only by other surfaces.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include "Core/Types.hpp"
#include "Core/Ray.hpp"
#include "Core/Sampling.hpp"
#include "Core/Constants.hpp"
#include "Core/Interaction.hpp"
#include "Core/SpectrumUtils.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Scene/Scene.hpp"
#include "Materials/Material.hpp"
#include "Materials/Medium.hpp"

namespace rayt {

#if defined(RAYT_SPECTRAL_MODE)
    constexpr int kSpectrumComponents = kSpectralSamples;   // wavelengths per path
#else
    constexpr int kSpectrumComponents = 3;                  // RGB channels
#endif

    /// Beer-Lambert transmittance of a straight segment inside a medium, per wavelength. Exact for a homogeneous medium; a medium with a
    /// density field is estimated by ratio tracking (unbiased, with the largest extinction of the wavelengths as the majorant).
    inline Spectrum mediumTransmittance(const Medium& m, const Point3& o, const Vector3& d, Real dist, const SampledWavelengths& lambda) {
        if (m.attenuatesContinuously()) return m.transmittance(dist, lambda);
        const Spectrum sT = m.sigmaTSpectrum(lambda);
        Spectrum T(1.0);
        if (m.isHomogeneous()) {
            for (int i = 0; i < kSpectrumComponents; ++i) T[i] = std::exp(-sT[i] * dist);
            return T;
        }
        Real majorant = 0;
        for (int i = 0; i < kSpectrumComponents; ++i) majorant = std::max(majorant, sT[i]);
        if (majorant <= Real(0)) return T;
        Real t = 0;
        for (int guard = 0; guard < 100000; ++guard) {
            t -= std::log(Real(1) - sampling::Random()) / majorant;
            if (t >= dist) break;
            const Real q = m.density(o + d * t);
            for (int i = 0; i < kSpectrumComponents; ++i) T[i] *= std::max(Real(0), Real(1) - q * sT[i] / majorant);
            if (isBlack(T)) break;
        }
        return T;
    }

    /**
     * @brief Transmittance along a shadow ray from `origin` in direction `dir` (unit length) up to distance `tMax` (infinite: out of the scene).
     *
     * `medium` is the medium the ray starts in (null: air). Null interfaces are crossed (the medium changes, the ray does not); any other
     * surface blocks the ray (the result is zero). The caller stops tMax a little before a light that it aims at.
     */
    inline Spectrum shadowTransmittance(const Scene& scene, Point3 origin, const Medium* medium, const Vector3& dir, Real tMax,
        const SampledWavelengths& lambda, Real time = Real(0))
    {
        Spectrum T(1.0);
        Real remaining = tMax;
        for (int crossings = 0; crossings < 64; ++crossings) {
            Ray ray(origin, dir, constants::RAY_EPSILON, medium);
            ray.tMax = remaining;
            ray.time = time;
            SurfaceInteraction h;
            const bool hit = scene.hit(ray, h);
            if (medium) {
                if (!hit && !std::isfinite(remaining)) return Spectrum(0.0);   // a ray that never leaves the medium: a leaky boundary
                T *= mediumTransmittance(*medium, origin, dir, hit ? h.t : remaining, lambda);
            }
            if (!hit) return T;
            if (!h.matPtr->isNullInterface()) return Spectrum(0.0);
            if (isBlack(T)) return T;
            if (h.matPtr->isMediumBoundary()) medium = h.frontFace ? h.matPtr->mediumInside() : h.matPtr->mediumOutside();   // gn faces the ray: a front-face hit enters the surface
            origin = h.p;
            if (std::isfinite(remaining)) remaining -= h.t;
        }
        return Spectrum(0.0);
    }

} // namespace rayt
