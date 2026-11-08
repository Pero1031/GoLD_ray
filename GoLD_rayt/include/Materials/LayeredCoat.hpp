/**
 * @file Materials/LayeredCoat.hpp
 * @brief A clear coat over a base surface, evaluated by a stochastic random walk between the two layers.
 *
 * Lacquer on gold, car paint, a varnished painting, oily skin: a rough dielectric interface (the coat) sits on a base that is either
 * diffuse or a rough metal. Light that enters the coat is scattered by the base, goes back up, is partly reflected down again by the
 * coat from the inside, and so on. Closed forms exist only for very special cases; here the series is summed by Monte Carlo
 * (pbrt-v4's LayeredBxDF, Jakob et al. 2014, "A Comprehensive Framework for Rendering Layered Materials"), which is exact for any pair
 * of layers and any roughness, at the price of a stochastic eval() (unbiased, with some noise).
 *
 *   - eval(): the coat's own reflection plus random walks between the layers, with next event estimation through the coat at every
 *     base scatter (multiple importance sampling between the two ways to leave).
 *   - sample(): one random walk that ends when the path leaves through the coat. The returned weight f cos / pdf is the exact throughput
 *     of the walk; the pdf reported to the integrator is a deterministic approximation (only used to weigh against light sampling:
 *     unbiased for any positive function, efficient when it is close).
 *   - pdf(): that deterministic approximation.
 *
 * The coat interface is a smooth-to-rough GGX dielectric (roughness is clamped to a small minimum: a perfectly smooth coat would
 * need delta lobes). It is lossless apart from the single-scattering energy loss of the microfacet model. An optional tint makes
 * the coat absorb: Beer-Lambert along the path inside it (amber urushi lacquer, tinted varnish).
 *
 * The base is evaluated in the coat's medium: a conductor base must be built with the incident index of the coat
 * (RoughConductor's `mediumIor` argument), otherwise its Fresnel reflectance is that of a metal in air.
 */

#pragma once

#include <algorithm>
#include <memory>
#include <optional>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Core/Interaction.hpp"
#include "Core/Sampling.hpp"
#include "Geometry/Frame.hpp"
#include "Materials/Material.hpp"
#include "Color/SpectralConversion.hpp"

namespace rayt {

    class LayeredCoat : public Material {
    public:
        struct Params {
            Real ior = Real(1.5);              ///< refractive index of the coat
            Real roughness = Real(0.1);        ///< perceptual roughness of the coat (alpha = roughness^2), clamped to >= 0.02
            Vector3 tint = Vector3(1.0);       ///< transmittance of the coat along the normal (round trip is longer): 1 = clear
            int samples = 1;                   ///< random walks per eval() (more = less noise, slower)
            int maxDepth = 12;                 ///< interface events per walk
        };

        /// A diffuse base (Lambertian reflectance `albedo`).
        LayeredCoat(const Params& p, const RgbColor& albedo);

        /// A glossy base: any reflection-only material evaluated with the coat as the incident medium (see the file comment).
        LayeredCoat(const Params& p, std::shared_ptr<const Material> base);

        Spectrum eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda, TransportMode mode) const override;

        std::optional<BSDFSample> sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u,
            const SampledWavelengths& lambda, TransportMode mode) const override;

        Real pdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda) const override;

        bool isSpecular() const override { return false; }

        Vector3 aovAlbedo(const SurfaceInteraction& rec) const override;

    private:
        Params m_p;
        Real m_alpha = Real(0.01);
        RgbColor m_albedo = RgbColor(0.5);
        std::shared_ptr<const Material> m_base;   // null = diffuse base
    };

} // namespace rayt
