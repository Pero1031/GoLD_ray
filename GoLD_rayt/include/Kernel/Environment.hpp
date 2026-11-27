#pragma once
/**
 * @file Kernel/Environment.hpp
 * @brief The environment map of the kernel renderer: lookup, importance sampling and its density (docs/GPU.md).
 *
 * Mirrors the engine's EnvMap (IO/EnvMap.hpp) step by step: the equirectangular mapping (u = (atan2(z, x) + pi) / 2 pi,
 * v = 1 - acos(y) / pi, row 0 of the image at the +y pole), the bilinear lookup in float (EnvMap::sampleBilinear: glm::mix with a
 * float weight computes in float), and Distribution2D / Distribution1D, which the engine already keeps in float. T is the
 * precision of the mapping between directions and (u, v), double in the engine; with T = double these functions give the
 * engine's values, with T = float (atan2, acos, sin and cos from Kernel/KernelMath.hpp) the same bits on the CPU and the GPU.
 * Values are the map's; the caller multiplies by KScene::envScale.
 */

#include "Kernel/KernelScene.hpp"
#include "Kernel/KernelMath.hpp"

namespace rayt::kernel {

    RAYT_HD_INLINE bool hasEnvironment(const KScene& s) { return s.envWidth > 0 && s.envHeight > 0; }

    /// EnvMap::dirToUV for a unit direction.
    template <class T> RAYT_HD_INLINE void envDirToUV(const Vec3<T>& d, T& u, T& v) {
        const T theta = kacos(kmin(kmax(d.y, T(-1)), T(1)));
        const T phi = katan2(d.z, d.x);
        u = (phi + kPi<T>()) / (T(2) * kPi<T>());
        v = T(1) - (theta / kPi<T>());
    }

    /// EnvMap::sinThetaFromV.
    template <class T> RAYT_HD_INLINE T envSinThetaFromV(T v) { return ksin(kPi<T>() * (T(1) - ksaturate(v))); }

    /// EnvMap::uvToDir.
    template <class T> RAYT_HD_INLINE Vec3<T> envUVToDir(T u, T v) {
        const T phi = (ksaturate(u) * T(2) * kPi<T>()) - kPi<T>();
        const T theta = kPi<T>() * (T(1) - ksaturate(v));
        T sinTheta, cosTheta, sinPhi, cosPhi;
        sinTheta = ksin(theta); cosTheta = kcos(theta);
        cosPhi = kcos(phi); sinPhi = ksin(phi);
        return Vec3<T>(sinTheta * cosPhi, cosTheta, sinTheta * sinPhi);
    }

    /// EnvMap::texel: wraps in x (longitude), clamps in y.
    RAYT_HD_INLINE Vec3<float> envTexel(const KScene& s, int x, int y) {
        const int w = int(s.envWidth), h = int(s.envHeight);
        x = (x % w + w) % w;
        y = kmin(kmax(y, 0), h - 1);
        const float* p = s.envPixels + 3 * (size_t(y) * size_t(w) + size_t(x));
        return Vec3<float>(p[0], p[1], p[2]);
    }

    /// EnvMap::sampleBilinear (all in float, as glm::mix with a float weight).
    RAYT_HD_INLINE Vec3<float> envBilinear(const KScene& s, float u, float v) {
        u = u - floorf(u);
        v = kmin(kmax(v, 0.0f), 1.0f);
        const float x = u * float(int(s.envWidth)) - 0.5f;
        const float y = (1.0f - v) * float(int(s.envHeight)) - 0.5f;
        const int x0 = int(floorf(x)), y0 = int(floorf(y));
        const float tx = x - floorf(x), ty = y - floorf(y);
        const Vec3<float> c00 = envTexel(s, x0, y0), c10 = envTexel(s, x0 + 1, y0), c01 = envTexel(s, x0, y0 + 1), c11 = envTexel(s, x0 + 1, y0 + 1);
        const Vec3<float> a = c00 * (1.0f - tx) + c10 * tx;   // glm::mix(c00, c10, tx)
        const Vec3<float> b = c01 * (1.0f - tx) + c11 * tx;
        return a * (1.0f - ty) + b * ty;
    }

    /// Distribution1D::sampleContinuous over func[0, n), cdf[0, n] (std::upper_bound written out).
    RAYT_HD_INLINE float envSample1D(const float* func, const float* cdf, float funcInt, int n, float u, float& pdf, int& off) {
        u = kmin(u, 0.99999994f);
        int lo = 0, hi = n + 1;   // the first index whose cdf is > u, in [0, n + 1]
        while (lo < hi) {
            const int mid = (lo + hi) / 2;
            if (cdf[mid] > u) hi = mid; else lo = mid + 1;
        }
        const int offset = kmin(kmax(lo - 1, 0), n - 1);
        off = offset;
        pdf = funcInt > 0.0f ? func[offset] / funcInt : 1.0f;
        float du = u - cdf[offset];
        const float denom = cdf[offset + 1] - cdf[offset];
        if (denom > 0.0f) du /= denom;
        else du = 0.0f;
        return (float(offset) + du) / float(n);
    }

    /// Distribution2D::sampleContinuous: (x, y) in [0, 1)^2 (y = image rows from the top) and its density.
    RAYT_HD_INLINE void envSample2D(const KScene& s, float u1, float u2, float& x, float& y, float& pdf) {
        const int w = int(s.envWidth), h = int(s.envHeight);
        float pdfV, pdfU;
        int vOff, uOff;
        y = envSample1D(s.envMarginalFunc, s.envMarginalCdf, s.envMarginalInt, h, u2, pdfV, vOff);
        x = envSample1D(s.envCondFunc + size_t(vOff) * size_t(w), s.envCondCdf + size_t(vOff) * size_t(w + 1), s.envCondInt[vOff], w, u1, pdfU, uOff);
        pdf = pdfU * pdfV;
    }

    /// Distribution2D::pdf.
    RAYT_HD_INLINE float envPdf2D(const KScene& s, float x, float y) {
        const int w = int(s.envWidth), h = int(s.envHeight);
        const int v = kmin(kmax(int(y * float(h)), 0), h - 1);
        const int u = kmin(kmax(int(x * float(w)), 0), w - 1);
        if (s.envMarginalInt == 0.0f) return 1.0f;
        if (s.envCondInt[v] == 0.0f) return 0.0f;
        const float pv = s.envMarginalFunc[v] / s.envMarginalInt;
        const float puv = s.envCondFunc[size_t(v) * size_t(w) + size_t(u)] / s.envCondInt[v];
        return pv * puv;
    }

    /// EnvMap::eval: the map's radiance in direction dir.
    template <class T> RAYT_HD_INLINE Vec3<T> envEval(const KScene& s, const Vec3<T>& dir) {
        T u, v;
        envDirToUV(knormalizeGlm(dir), u, v);
        const Vec3<float> c = envBilinear(s, float(u), float(v));
        return Vec3<T>(T(c.x), T(c.y), T(c.z));
    }

    /// EnvMap::sample: a direction drawn by the map's luminance x sin(theta), its radiance and its density over solid angle (0: none).
    template <class T> RAYT_HD_INLINE void envSample(const KScene& s, float u1, float u2, Vec3<T>& wi, Vec3<T>& le, T& pdfW) {
        pdfW = T(0);
        float x, y, pdfUVf;
        envSample2D(s, u1, u2, x, y, pdfUVf);
        const T pdfUV = T(pdfUVf);
        if (pdfUV <= T(0)) { le = Vec3<T>(0); return; }
        const T uSph = T(x), vSph = T(1) - T(y);
        wi = envUVToDir(uSph, vSph);
        const Vec3<float> c = envBilinear(s, float(uSph), float(vSph));
        le = Vec3<T>(T(c.x), T(c.y), T(c.z));
        const T sinTheta = envSinThetaFromV(vSph);
        if (sinTheta <= T(0)) { le = Vec3<T>(0); return; }
        pdfW = pdfUV / (T(2) * kPi<T>() * kPi<T>() * sinTheta);
    }

    /// EnvMap::pdf: the density of envSample over solid angle in direction wi.
    template <class T> RAYT_HD_INLINE T envPdf(const KScene& s, const Vec3<T>& wi) {
        T uSph, vSph;
        envDirToUV(knormalizeGlm(wi), uSph, vSph);
        const T vImg = T(1) - vSph;
        const T pdfUV = T(envPdf2D(s, float(uSph), float(vImg)));
        if (pdfUV <= T(0)) return T(0);
        const T sinTheta = envSinThetaFromV(vSph);
        if (sinTheta <= T(0)) return T(0);
        return pdfUV / (T(2) * kPi<T>() * kPi<T>() * sinTheta);
    }

} // namespace rayt::kernel
