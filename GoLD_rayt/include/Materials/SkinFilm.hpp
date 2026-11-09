/**
 * @file Materials/SkinFilm.hpp
 * @brief A thin liquid film (lotion water, a humectant, an oil) lying on the skin relief and pooling in its furrows (docs/SKIN_MOISTURE.md).
 *
 * The skin surface has the height h(p) (negative = a furrow or a pore). The top of the film is
 *     z(p) = smoothmax(h(p), L),     thickness t(p) = z(p) - h(p) >= 0,
 * so the film fills the low places up to the level L and leaves the ridges bare. L is found from the mean thickness the user asks for.
 * The film is smooth: the specular lobe of the skin follows the normal of z (almost flat over a filled furrow), not the normal of h.
 * The reflection of the film is the film on skin seen from the air: R = R1 + (1 - R1)^2 R2 (air-film, film-stratum corneum; no interference,
 * so the film must be thicker than a few wavelengths).
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <vector>

#include "Core/Types.hpp"
#include "Core/Fresnel.hpp"
#include "Core/Math.hpp"
#include "Materials/SkinRelief.hpp"
#include "Microfacet/GGX.hpp"

namespace rayt {

    class SkinFilm {
    public:
        struct Params {
            Real thickness = Real(5e-6);    ///< mean thickness over the skin (m)
            Real ior = Real(1.33);          ///< refractive index of the liquid (water 1.33, glycerin 1.47, oil about 1.45)
            Real roughness = Real(0.02);    ///< perceptual roughness of the film surface
            Real softness = Real(1.5e-6);   ///< width (m) over which the film's edge is rounded (the smooth maximum)
            Real wetting = Real(1e-6);      ///< thickness (m) at which the skin counts as covered (the specular lobe has changed over to the film's)
            Real skinIor = Real(1.55);      ///< refractive index of the stratum corneum below
        };

        /// What the film does at one point.
        struct Hit {
            Real thickness = 0;    ///< t (m)
            Real cover = 0;        ///< 0..1: 1 - exp(-t / wetting)
            Vector3 normal;        ///< normal of the top of the film
        };

        SkinFilm(const Params& p, std::shared_ptr<const SkinRelief> relief)
            : m_p(p), m_relief(std::move(relief)), m_alpha(MicrofacetDistribution::roughnessToAlpha(p.roughness))
        {
            solveLevel();
            // Cosine-weighted average of the film's reflectance, 2 Integral F mu dmu (as CoatedDiffuse does for its dry coat).
            constexpr int n = 64;
            Real sum = 0;
            for (int i = 0; i < n; ++i) {
                const Real mu = (Real(i) + Real(0.5)) / Real(n);
                sum += fresnelReflectance(mu) * mu;
            }
            m_fAvg = Real(2) * sum / Real(n);
        }

        const Params& params() const { return m_p; }
        Real level() const { return m_level; }
        Real alpha() const { return m_alpha; }
        Real fresnelAverage() const { return m_fAvg; }

        /// Reflectance of the film on skin seen from the air at the angle with cosine `cosTheta`.
        Real fresnelReflectance(Real cosTheta) const {
            const Real c = math::saturate(cosTheta);
            const Real r1 = fresnel::fresnelDielectric(c, Real(1), m_p.ior);
            const Real sinT = math::safe_sqrt(Real(1) - c * c) / m_p.ior;
            const Real cosT = math::safe_sqrt(Real(1) - sinT * sinT);
            const Real r2 = fresnel::fresnelDielectric(cosT, m_p.ior, m_p.skinIor);
            return r1 + (Real(1) - r1) * (Real(1) - r1) * r2;
        }

        /// Thickness, cover and top normal at p. `n` is the smooth normal of the skin, `wo` points to the viewer (the normal never turns away from it).
        Hit at(const Point3& p, const Vector3& n, const Vector3& wo) const {
            Vector3 grad;
            const Real h = m_relief->surface(p, n, grad);
            const Real d = h - m_level, root = std::sqrt(d * d + m_p.softness * m_p.softness);
            Hit hit;
            hit.thickness = Real(0.5) * (m_level - h + root);
            hit.cover = Real(1) - std::exp(-hit.thickness / m_p.wetting);
            const Real w = Real(0.5) * (Real(1) + d / root);   // slope of smoothmax(h, L) with respect to h: 1 on a bare ridge, 0 deep in a puddle
            Vector3 nf = glm::normalize(n - w * grad);
            const Real c = glm::dot(nf, wo);
            if (c < Real(0.05)) nf = glm::normalize(nf + wo * (Real(0.05) - c));
            hit.normal = nf;
            return hit;
        }

        /// Mean thickness (m) over `count` random points of a plane patch; the test of the level (should equal params().thickness).
        Real meanThickness(int count, unsigned seed = 7u) const {
            std::mt19937 rng(seed);
            std::uniform_real_distribution<double> u(0.0, 1.0);
            double sum = 0;
            for (int i = 0; i < count; ++i) sum += thicknessAt(planeHeight(Point3(u(rng) * kPatch, u(rng) * kPatch, 0.0)));
            return Real(sum / count);
        }

    private:
        static constexpr double kPatch = 0.024;   // side of the patch (m) whose heights set the level: the tile of the wrinkle field, many furrow cells

        Real planeHeight(const Point3& p) const {   // the relief on a plane with the normal z
            Vector3 g;
            return m_relief->surface(p, Vector3(0, 0, 1), g);
        }

        Real thicknessAt(Real height) const {
            const Real d = height - m_level;
            return Real(0.5) * (m_level - height + std::sqrt(d * d + m_p.softness * m_p.softness));
        }

        /// The level L at which the mean thickness over the relief is the one asked for (bisection; the mean grows with L).
        void solveLevel() {
            std::mt19937 rng(12345u);
            std::uniform_real_distribution<double> u(0.0, 1.0);
            std::vector<Real> heights(300000);   // a thin film sits in the deepest pores only: the mean is set by rare points, so many are needed
            Real lo = 1, hi = -1;
            for (Real& v : heights) {
                v = planeHeight(Point3(u(rng) * kPatch, u(rng) * kPatch, 0.0));
                lo = std::min(lo, v); hi = std::max(hi, v);
            }
            const Real target = std::max(m_p.thickness, Real(1e-9));
            Real a = lo - Real(1e-3), b = hi + target + Real(10) * m_p.softness;
            for (int it = 0; it < 80; ++it) {
                const Real mid = Real(0.5) * (a + b);
                double sum = 0;
                for (Real v : heights) {
                    const Real d = v - mid;
                    sum += Real(0.5) * (mid - v + std::sqrt(d * d + m_p.softness * m_p.softness));
                }
                if (Real(sum / double(heights.size())) < target) a = mid; else b = mid;
            }
            m_level = Real(0.5) * (a + b);
        }

        Params m_p;
        std::shared_ptr<const SkinRelief> m_relief;
        Real m_alpha;
        Real m_level = 0;
        Real m_fAvg = 0;
    };

} // namespace rayt
