#pragma once
/**
 * @file Kernel/Media.hpp
 * @brief The participating media inside glass, as the kernel renderer reads them (docs/GPU.md).
 *
 * Materials/Medium.hpp, Materials/TabulatedPhase.hpp and Materials/DensityField.hpp step by step: the phase functions (a mixture
 * of two Henyey-Greenstein lobes, or a table such as Mie scattering), the density field of a heterogeneous medium (fractal value
 * noise with a spherical mask), and the transmittance of a straight segment (exact for a homogeneous medium, ratio tracking for
 * a density field, Render/Transmittance.hpp). With T = double they give the engine's bits; with T = float they are what the GPU
 * computes. The walk itself (PathIntegrator::advanceInMedium) is in Kernel/PathTracer.hpp.
 */

#include "Kernel/KernelMath.hpp"
#include "Kernel/KernelScene.hpp"
#include "Kernel/Materials.hpp"
#include "Kernel/Rng.hpp"

namespace rayt::kernel {

    /// A phase function's table (KScene::phaseF and phaseCdf from the medium's offset), in the precision of the evaluation.
    template <class T> struct PhaseTable {
        const T* f;     ///< density in theta at the nodes (integrates to 1)
        const T* cdf;   ///< cumulative distribution at the nodes
    };

    namespace media {

        template <class T> RAYT_HD_INLINE T kclampT(T v, T lo, T hi) { return v < lo ? lo : (hi < v ? hi : v); }   // std::clamp

        /// henyeyGreenstein: 1 / (4 pi) (1 - g^2) / (1 + g^2 - 2 g cos)^(3/2).
        template <class T> RAYT_HD_INLINE T henyeyGreenstein(T cosTheta, T g) {
            const T d = T(1) + g * g - T(2) * g * cosTheta;
            return (T(1) - g * g) / (T(4) * kPi<T>() * d * ksqrt(kmax(d, T(1e-12))));
        }

        /// A direction around `dir` with the polar cosine cosTheta and the azimuth 2 pi u (the engine's frame::Frame(dir)).
        template <class T> RAYT_HD_INLINE Vec3<T> aroundDirection(const Vec3<T>& dir, T cosTheta, T u) {
            const T sinTheta = ksqrt(kmax(T(0), T(1) - cosTheta * cosTheta));
            const T phi = T(2) * kPi<T>() * u;
            const Frame3<T> frame = Frame3<T>::fromNormal(dir);
            return frame.toWorld(Vec3<T>(sinTheta * kcos(phi), sinTheta * ksin(phi), cosTheta));
        }

        /// sampleHenyeyGreenstein.
        template <class T> RAYT_HD Vec3<T> sampleHenyeyGreenstein(const Vec3<T>& dir, T g, T ux, T uy) {
            T cosTheta;
            if (kabs(g) < T(1e-3)) cosTheta = T(1) - T(2) * ux;
            else {
                const T s = (T(1) - g * g) / (T(1) - g + T(2) * g * ux);
                cosTheta = (T(1) + g * g - s * s) / (T(2) * g);
            }
            cosTheta = kclampT(cosTheta, T(-1), T(1));
            return aroundDirection(dir, cosTheta, uy);
        }

        /// TabulatedPhase::eval.
        template <class T> RAYT_HD T tableEval(const MediumParams<T>& m, const PhaseTable<T>& tab, T cosTheta) {
            const T c = kclampT(cosTheta, T(-1), T(1));
            const T theta = kacos(c);
            const T x = theta / m.tableStep;
            const int i = kmin(int(x), int(m.tableCount) - 2);
            const T t = x - T(i);
            const T f = tab.f[i] + (tab.f[i + 1] - tab.f[i]) * t;
            const T s = kmax(ksqrt(kmax(T(1) - c * c, T(0))), T(1e-9));
            return f / (T(2) * kPi<T>() * s);
        }

        /// TabulatedPhase::sampleCosine (u in [0, 1)).
        template <class T> RAYT_HD T tableSampleCosine(const MediumParams<T>& m, const PhaseTable<T>& tab, T u) {
            const int n = int(m.tableCount);
            const T uu = kclampT(u, T(0), T(1) - T(1e-12));
            // the last node whose cumulative value is <= uu (std::upper_bound - 1)
            int lo = 0, hi = n;   // first index with cdf > uu, by bisection
            while (lo < hi) {
                const int mid = (lo + hi) / 2;
                if (uu < tab.cdf[mid]) hi = mid; else lo = mid + 1;
            }
            const int i = kclampT(lo - 1, 0, n - 2);
            const T r = uu - tab.cdf[i];
            const T f0 = tab.f[i], f1 = tab.f[i + 1];
            const T k = (f1 - f0) / m.tableStep;
            const T disc = kmax(f0 * f0 + T(2) * k * r, T(0));
            const T denom = f0 + ksqrt(disc);
            T s = denom > T(0) ? T(2) * r / denom : T(0);
            s = kclampT(s, T(0), m.tableStep);
            return kcos(T(i) * m.tableStep + s);
        }

        /// NoiseDensity's hash of a lattice point to [0, 1).
        template <class T> RAYT_HD_INLINE T hash01(int x, int y, int z, uint32_t seed) {
            uint32_t h = seed ^ 0x9e3779b9u;
            h ^= uint32_t(x) * 0x85ebca6bu; h = (h << 13) | (h >> 19); h *= 5u; h += 0xe6546b64u;
            h ^= uint32_t(y) * 0xc2b2ae35u; h = (h << 13) | (h >> 19); h *= 5u; h += 0xe6546b64u;
            h ^= uint32_t(z) * 0x27d4eb2fu; h = (h << 13) | (h >> 19); h *= 5u; h += 0xe6546b64u;
            h ^= h >> 16; h *= 0x85ebca6bu; h ^= h >> 13; h *= 0xc2b2ae35u; h ^= h >> 16;
            return T(h) / T(4294967296.0);
        }

        /// NoiseDensity::valueNoise: lattice values interpolated with the quintic fade.
        template <class T> RAYT_HD T valueNoise(const Vec3<T>& p, uint32_t seed) {
            const T fx = kfloor(p.x), fy = kfloor(p.y), fz = kfloor(p.z);
            const int ix = int(fx), iy = int(fy), iz = int(fz);
            const auto fade = [](T t) { return t * t * t * (t * (t * T(6) - T(15)) + T(10)); };
            const T u = fade(p.x - fx), v = fade(p.y - fy), w = fade(p.z - fz);
            const auto h = [&](int dx, int dy, int dz) { return hash01<T>(ix + dx, iy + dy, iz + dz, seed); };
            const T x00 = h(0, 0, 0) + u * (h(1, 0, 0) - h(0, 0, 0));
            const T x10 = h(0, 1, 0) + u * (h(1, 1, 0) - h(0, 1, 0));
            const T x01 = h(0, 0, 1) + u * (h(1, 0, 1) - h(0, 0, 1));
            const T x11 = h(0, 1, 1) + u * (h(1, 1, 1) - h(0, 1, 1));
            const T y0 = x00 + v * (x10 - x00), y1 = x01 + v * (x11 - x01);
            return y0 + w * (y1 - y0);
        }

    } // namespace media

    /// PhaseFunction::eval: the density per steradian of scattering by the angle with cosine cosTheta between the old and the new propagation.
    template <class T> RAYT_HD T phaseEval(const MediumParams<T>& m, const PhaseTable<T>& tab, T cosTheta) {
        if (m.phaseKind == kPhaseTable) return media::tableEval(m, tab, cosTheta);
        return m.w * media::henyeyGreenstein(cosTheta, m.g1) + (T(1) - m.w) * media::henyeyGreenstein(cosTheta, m.g2);
    }

    /// PhaseFunction::sample: a new propagation direction for a ray travelling along `dir`, with exactly phaseEval's density.
    template <class T> RAYT_HD Vec3<T> phaseSample(const MediumParams<T>& m, const PhaseTable<T>& tab, const Vec3<T>& dir, T ux, T uy) {
        if (m.phaseKind == kPhaseTable) return media::aroundDirection(dir, media::tableSampleCosine(m, tab, ux), uy);
        // the first number picks the lobe and, rescaled, is reused for the angle
        if (ux < m.w) return media::sampleHenyeyGreenstein(dir, m.g1, m.w > T(0) ? ux / m.w : T(0), uy);
        return media::sampleHenyeyGreenstein(dir, m.g2, (ux - m.w) / (T(1) - m.w), uy);
    }

    /// The relative density q(p) in [0, 1] (ChromaticMedium::density with NoiseDensity); 1 for a homogeneous medium.
    template <class T> RAYT_HD T mediumDensity(const MediumParams<T>& m, const Vec3<T>& p) {
        if (!m.noise) return T(1);
        T sum = T(0), amplitude = T(1), total = T(0), f = m.frequency;
        for (uint32_t o = 0; o < m.octaves; ++o) {
            sum += amplitude * media::valueNoise(p * f, m.seed + o * 7919u);
            total += amplitude;
            amplitude *= T(0.5);
            f *= T(2);
        }
        const T n = sum / total;
        const T t = media::kclampT((n - m.lo) / (m.hi - m.lo), T(0), T(1));
        T q = t * t * (T(3) - T(2) * t);
        if (m.maskRadius > T(0)) {   // 1 inside start * radius, 0 at the radius (smooth)
            const T r = klength(p - Vec3<T>(m.maskCenter[0], m.maskCenter[1], m.maskCenter[2]));
            const T s = media::kclampT((r - m.maskStart * m.maskRadius) / ((T(1) - m.maskStart) * m.maskRadius), T(0), T(1));
            q *= T(1) - s * s * (T(3) - T(2) * s);
        }
        return media::kclampT(q, T(0), T(1));
    }

    /// mediumTransmittance (Render/Transmittance.hpp): Beer-Lambert along a straight segment, per channel; exact for a homogeneous medium,
    /// ratio tracking with the largest extinction as the majorant for a density field (it draws from rng).
    RAYT_HD_INLINE Vec3<float> mediumTransmittance(const KMedium& m, const Vec3<float>& o, const Vec3<float>& d, float dist, Pcg32& rng) {
        Vec3<float> T(1.0f);
        if (!m.noise) {
            T.x = kexp(-m.sigmaT[0] * dist);
            T.y = kexp(-m.sigmaT[1] * dist);
            T.z = kexp(-m.sigmaT[2] * dist);
            return T;
        }
        const float majorant = kmax(kmax(kmax(0.0f, m.sigmaT[0]), m.sigmaT[1]), m.sigmaT[2]);
        if (majorant <= 0.0f) return T;
        float t = 0.0f;
        for (int guard = 0; guard < 100000; ++guard) {
            t -= klog(1.0f - rng.nextFloat()) / majorant;
            if (t >= dist) break;
            const float q = mediumDensity(m, o + d * t);
            T.x *= kmax(0.0f, 1.0f - q * m.sigmaT[0] / majorant);
            T.y *= kmax(0.0f, 1.0f - q * m.sigmaT[1] / majorant);
            T.z *= kmax(0.0f, 1.0f - q * m.sigmaT[2] / majorant);
            // Below 1e-20 the light that is left does not count (the sun's radiance times it is 1e-15): stopping there keeps the result
            // and the random numbers drawn independent of how a device treats subnormal floats.
            if (!(kmaxComponent(T) > 1e-20f)) return Vec3<float>(0.0f);
        }
        return T;
    }

} // namespace rayt::kernel
