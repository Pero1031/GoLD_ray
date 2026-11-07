/**
 * @file Materials/DensityField.hpp
 * @brief Density fields q(x) in [0, 1] for heterogeneous participating media.
 *
 * A heterogeneous medium has the extinction sigma_t(x, lambda) = q(x) * sigma_bar(lambda), where sigma_bar is the
 * constant upper bound (the majorant) that the integrator samples distances with. q = 1 is the densest point; where
 * q < 1 the integrator draws "null collisions" that change nothing (delta tracking, see Medium.hpp / Integrator.hpp).
 * The field must never exceed 1, otherwise the walk is biased.
 *
 * NoiseDensity is a procedural cloud: fractal value noise in world coordinates, mapped to [0, 1] by a smooth step.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "Core/Types.hpp"

namespace rayt {

    class DensityField {
    public:
        virtual ~DensityField() = default;
        /// Relative density at a world-space point, in [0, 1].
        virtual Real operator()(const Point3& p) const = 0;
    };

    /// The same density everywhere (a test: a medium with q < 1 and null collisions must look like a thinner homogeneous one).
    class ConstantDensity final : public DensityField {
    public:
        explicit ConstantDensity(Real q) : m_q(std::clamp(q, Real(0), Real(1))) {}
        Real operator()(const Point3&) const override { return m_q; }
    private:
        Real m_q;
    };

    /**
     * @brief Fractal value noise, mapped by smoothstep(lo, hi, noise).
     *
     * @param frequency  features per metre of the first octave (the size of the biggest clouds is about 1 / frequency)
     * @param octaves    number of noise layers; each has twice the frequency and half the weight
     * @param lo, hi     noise values below lo give density 0, above hi density 1 (a narrow band: sharp clouds; a wide one:
     *                   soft haze). The noise is bell shaped around 0.5, so lo = 0.45, hi = 0.65 covers roughly a third of
     *                   the volume with dense parts.
     * @param seed       another seed, another cloud
     *
     * A spherical mask (setSphereMask) fades the density to 0 towards a sphere's surface, so that a cloud inside a
     * (index-matched) sphere has soft edges instead of the sphere's silhouette.
     */
    class NoiseDensity final : public DensityField {
    public:
        NoiseDensity(Real frequency, int octaves, Real lo, Real hi, std::uint32_t seed = 1)
            : m_frequency(frequency), m_octaves(std::max(1, octaves)), m_lo(lo), m_hi(std::max(hi, lo + Real(1e-6))), m_seed(seed) {}

        Real operator()(const Point3& p) const override {
            Real sum = 0, amplitude = 1, total = 0, f = m_frequency;
            for (int o = 0; o < m_octaves; ++o) {
                sum += amplitude * valueNoise(p * f, m_seed + std::uint32_t(o) * 7919u);
                total += amplitude;
                amplitude *= Real(0.5);
                f *= Real(2);
            }
            const Real n = sum / total;
            const Real t = std::clamp((n - m_lo) / (m_hi - m_lo), Real(0), Real(1));
            Real q = t * t * (Real(3) - Real(2) * t);
            if (m_maskRadius > Real(0)) {   // 1 inside start * radius, 0 at the radius (smooth)
                const Real r = glm::length(p - m_maskCenter);
                const Real s = std::clamp((r - m_maskStart * m_maskRadius) / ((Real(1) - m_maskStart) * m_maskRadius), Real(0), Real(1));
                q *= Real(1) - s * s * (Real(3) - Real(2) * s);
            }
            return q;
        }

        /// Fades the density to 0 at `radius` from `center`; it is untouched inside `startFraction` * radius.
        NoiseDensity& setSphereMask(const Point3& center, Real radius, Real startFraction) {
            m_maskCenter = center; m_maskRadius = radius; m_maskStart = std::clamp(startFraction, Real(0), Real(0.999));
            return *this;
        }

    private:
        static Real hash01(int x, int y, int z, std::uint32_t seed) {
            std::uint32_t h = seed ^ 0x9e3779b9u;
            h ^= std::uint32_t(x) * 0x85ebca6bu; h = (h << 13) | (h >> 19); h *= 5u; h += 0xe6546b64u;
            h ^= std::uint32_t(y) * 0xc2b2ae35u; h = (h << 13) | (h >> 19); h *= 5u; h += 0xe6546b64u;
            h ^= std::uint32_t(z) * 0x27d4eb2fu; h = (h << 13) | (h >> 19); h *= 5u; h += 0xe6546b64u;
            h ^= h >> 16; h *= 0x85ebca6bu; h ^= h >> 13; h *= 0xc2b2ae35u; h ^= h >> 16;
            return Real(h) / Real(4294967296.0);
        }

        /// Value noise: random values on the integer lattice, interpolated smoothly (quintic fade). Range [0, 1).
        static Real valueNoise(const Point3& p, std::uint32_t seed) {
            const Real fx = std::floor(p.x), fy = std::floor(p.y), fz = std::floor(p.z);
            const int ix = int(fx), iy = int(fy), iz = int(fz);
            auto fade = [](Real t) { return t * t * t * (t * (t * Real(6) - Real(15)) + Real(10)); };
            const Real u = fade(p.x - fx), v = fade(p.y - fy), w = fade(p.z - fz);
            auto h = [&](int dx, int dy, int dz) { return hash01(ix + dx, iy + dy, iz + dz, seed); };
            const Real x00 = h(0, 0, 0) + u * (h(1, 0, 0) - h(0, 0, 0));
            const Real x10 = h(0, 1, 0) + u * (h(1, 1, 0) - h(0, 1, 0));
            const Real x01 = h(0, 0, 1) + u * (h(1, 0, 1) - h(0, 0, 1));
            const Real x11 = h(0, 1, 1) + u * (h(1, 1, 1) - h(0, 1, 1));
            const Real y0 = x00 + v * (x10 - x00), y1 = x01 + v * (x11 - x01);
            return y0 + w * (y1 - y0);
        }

        Real m_frequency;
        int m_octaves;
        Real m_lo, m_hi;
        std::uint32_t m_seed;
        Point3 m_maskCenter = Point3(0.0);
        Real m_maskRadius = 0, m_maskStart = Real(0.5);
    };

} // namespace rayt
