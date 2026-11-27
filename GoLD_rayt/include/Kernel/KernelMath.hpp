#pragma once
/**
 * @file Kernel/KernelMath.hpp
 * @brief Scalar math and vector types for kernel code (CPU and GPU).
 *
 * Float: every operation gives the same bits on the CPU and on the GPU. +, -, *, /, sqrtf and fmaf are correctly rounded
 * on both (the GPU module is built with IEEE division and square root and without implicit fused multiply-adds, as
 * Dr.Jit emits them: an FMA only where kernel code writes kfma). sin and cos are our own float polynomials (ksincos), not
 * the C library's or CUDA's sinf / cosf, which differ in the last bits and, in CUDA, fall back to double for huge angles.
 *
 * Double (CPU only, Kernel/HostDevice.hpp): the C library's functions, so kernel code that mirrors an engine function
 * gives the engine's bits. min / max keep the argument order of std::min / std::max for the same reason.
 */

#include <math.h>
#include <stdint.h>
#include <glm/glm.hpp>

#include "Kernel/HostDevice.hpp"

namespace rayt::kernel {

    template <class T> using Vec2 = glm::vec<2, T, glm::defaultp>;
    template <class T> using Vec3 = glm::vec<3, T, glm::defaultp>;

    template <class T> RAYT_HD_INLINE constexpr T kPi() { return T(3.14159265358979323846); }
    template <class T> RAYT_HD_INLINE constexpr T kTwoPi() { return T(2) * kPi<T>(); }

    RAYT_HD_INLINE float ksqrt(float x) { return sqrtf(x); }
    RAYT_HD_INLINE float kfma(float a, float b, float c) { return fmaf(a, b, c); }
    RAYT_HD_INLINE float kcopysign(float x, float s) { return copysignf(x, s); }
    RAYT_HD_INLINE float kabs(float x) { return fabsf(x); }
    RAYT_HD_INLINE float kfloor(float x) { return floorf(x); }

    /**
     * @brief Sine and cosine in float, the same bits on every device.
     *
     * The CEPHES single-precision algorithm (S. L. Moshier), as Dr.Jit implements it for all of Mitsuba's backends
     * (references/src/drjit/include/drjit/math.h): x is reduced by multiples of pi/4 with pi/4 split into three floats
     * (Cody-Waite), then degree-7 / degree-8 polynomials on [-pi/4, pi/4]. Every step is a correctly rounded float
     * operation or an explicit FMA, so the CPU and the GPU agree bit for bit. Error about 1 ulp for |x| < 8192 (Dr.Jit
     * measures at most 6e-8 absolute there); beyond that x is first reduced by 2 pi in float, which loses accuracy but stays
     * finite. NaN and infinity give NaN.
     */
    RAYT_HD_INLINE void ksincos(float x, float& sinOut, float& cosOut) {
        float xa = fabsf(x);
        if (!(xa < 8192.0f)) {
            if (!(xa <= 3.402823466e38f)) { sinOut = cosOut = x - x; return; }   // NaN, +-infinity
            xa = kfma(-6.28318548f, floorf(xa * 0.159154943f), xa);
        }
        // j = nearest even multiple of pi/4 at or below 4|x| / pi rounded up to even: the octant pair of x.
        int32_t j = int32_t(xa * 1.27323954473516269f);
        j = (j + 1) & ~1;
        const float y0 = float(j);
        float y = kfma(-y0, 0.78515625f, xa);
        y = kfma(-y0, 2.4187564849853515625e-4f, y);
        y = kfma(-y0, 3.77489497744594108e-8f, y);
        const float z = y * y;
        float ps = kfma(kfma(-1.9515295891e-4f, z, 8.3321608736e-3f), z, -1.6666654611e-1f) * z;
        float pc = kfma(kfma(2.443315711809948e-5f, z, -1.388731625493765e-3f), z, 4.166664568298827e-2f) * z;
        ps = kfma(ps, y, y);                                  // sin y
        pc = kfma(pc, z, kfma(z, -0.5f, 1.0f));               // cos y
        const bool swap = (j & 2) != 0;                       // odd quadrant: sin and cos trade places
        const float s = swap ? pc : ps, c = swap ? ps : pc;
        const bool negSin = ((j & 4) != 0) != (x < 0.0f);     // sin is odd: the sign of x and the half turn
        const bool negCos = ((j + 2) & 4) != 0;
        sinOut = negSin ? -s : s;
        cosOut = negCos ? -c : c;
    }
    RAYT_HD_INLINE float ksin(float x) { float s, c; ksincos(x, s, c); return s; }
    RAYT_HD_INLINE float kcos(float x) { float s, c; ksincos(x, s, c); return c; }

    /**
     * @brief atan2(y, x) in float, the same bits on every device: Dr.Jit's minimax fit of atan(sqrt(z)) / sqrt(z) on [0, 1]
     *        (references/src/drjit/include/drjit/math.h; about 6e-7 absolute at most), evaluated with explicit FMAs. 0 when x = y = 0.
     */
    RAYT_HD_INLINE float katan2(float y, float x) {
        const float ax = fabsf(x), ay = fabsf(y);
        const float mn = ay < ax ? ay : ax, mx = ax < ay ? ay : ax;
        if (!(mx > 0.0f)) return 0.0f;
        const float s = mn / mx, z = s * s;
        float t = 0.0078613793713198150252f;
        t = kfma(t, z, -0.037006525670417265220f);
        t = kfma(t, z, 0.083863120428809689910f);
        t = kfma(t, z, -0.13486708938456973185f);
        t = kfma(t, z, 0.19881342388439013552f);
        t = kfma(t, z, -0.33326497518773606976f);
        t = kfma(t, z, 0.99999934166683966009f);
        t = t * s;
        if (ay > ax) t = 1.57079632679489661923f - t;
        if (x < 0.0f) t = 3.14159265358979323846f - t;
        return y < 0.0f ? -t : t;
    }

    /// acos(x) in float, the same bits on every device: the CEPHES single-precision algorithm as Dr.Jit implements it (about 1 ulp).
    RAYT_HD_INLINE float kacos(float x) {
        const float xa = fabsf(x);
        const bool big = xa > 0.5f;
        const float x1 = 0.5f * (1.0f - xa);
        const float x3 = big ? x1 : x * x;
        const float x4 = big ? sqrtf(x1) : xa;
        float z = 4.2163199048e-2f;
        z = kfma(z, x3, 2.4181311049e-2f);
        z = kfma(z, x3, 4.5470025998e-2f);
        z = kfma(z, x3, 7.4953002686e-2f);
        z = kfma(z, x3, 1.666675242e-1f);
        z = kfma(z, x3 * x4, x4);
        if (big) {
            const float z2 = z + z;
            return x < 0.0f ? 3.14159265358979323846f - z2 : z2;
        }
        return 1.57079632679489661923f - copysignf(z, x);
    }

#if !defined(RAYT_GPU_CODE)
    // CPU only: device code that calls one of these does not compile.
    inline double ksqrt(double x) { return sqrt(x); }
    inline double kfma(double a, double b, double c) { return fma(a, b, c); }
    inline double kcos(double x) { return cos(x); }
    inline double ksin(double x) { return sin(x); }
    inline double kcopysign(double x, double s) { return copysign(x, s); }
    inline double kabs(double x) { return fabs(x); }
    inline double kfloor(double x) { return floor(x); }
    inline double katan2(double y, double x) { return atan2(y, x); }
    inline double kacos(double x) { return acos(x); }
#else
    // Device code: a double argument picks these deleted overloads (instead of converting to float) and fails to compile.
    RAYT_HD double ksqrt(double) = delete;
    RAYT_HD double kfma(double, double, double) = delete;
    RAYT_HD double kcos(double) = delete;
    RAYT_HD double ksin(double) = delete;
    RAYT_HD double kcopysign(double, double) = delete;
    RAYT_HD double kabs(double) = delete;
    RAYT_HD double kfloor(double) = delete;
    RAYT_HD double katan2(double, double) = delete;
    RAYT_HD double kacos(double) = delete;
#endif

    /// std::max(a, b): b if a < b, else a.
    template <class T> RAYT_HD_INLINE T kmax(T a, T b) { return (a < b) ? b : a; }
    /// std::min(a, b): b if b < a, else a.
    template <class T> RAYT_HD_INLINE T kmin(T a, T b) { return (b < a) ? b : a; }
    /// glm::clamp(x, 0, 1) = min(max(x, 0), 1).
    template <class T> RAYT_HD_INLINE T ksaturate(T x) { return kmin(kmax(x, T(0)), T(1)); }
    /// sqrt(max(x, 0)), as math::safe_sqrt.
    template <class T> RAYT_HD_INLINE T ksafeSqrt(T x) { return ksqrt(kmax(x, T(0))); }
    template <class T> RAYT_HD_INLINE T ksqr(T x) { return x * x; }
    template <class T> RAYT_HD_INLINE T kinfinity() { return T(INFINITY); }
    /// False for NaN and +-infinity (x - x is 0 only for finite x; IEEE arithmetic, no fast math, keeps it).
    template <class T> RAYT_HD_INLINE bool kisFinite(T x) { return x - x == T(0); }
    /// 1 / x, with |x| raised to at least 1e-6 (sign kept), as math::safe_recip.
    template <class T> RAYT_HD_INLINE T ksafeRecip(T x) {
        const T eps = T(1e-6);
        if (kabs(x) < eps) return T(1) / (x >= T(0) ? eps : -eps);
        return T(1) / x;
    }

    // ---- vectors (glm's arithmetic operators are plain IEEE operations on both devices; these spell out the rest) ----

    /// (a.x b.x + a.y b.y) + a.z b.z, the order of glm::dot.
    template <class T> RAYT_HD_INLINE T kdot(const Vec3<T>& a, const Vec3<T>& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
    /// glm::cross.
    template <class T> RAYT_HD_INLINE Vec3<T> kcross(const Vec3<T>& a, const Vec3<T>& b) {
        return Vec3<T>(a.y * b.z - b.y * a.z, a.z * b.x - b.z * a.x, a.x * b.y - b.x * a.y);
    }
    template <class T> RAYT_HD_INLINE T klength(const Vec3<T>& a) { return ksqrt(kdot(a, a)); }
    /// a / |a|.
    template <class T> RAYT_HD_INLINE Vec3<T> knormalize(const Vec3<T>& a) { return a / klength(a); }
    /// glm::normalize: a times 1 / |a| (one rounding more than knormalize; kernel code that mirrors engine code uses this one).
    template <class T> RAYT_HD_INLINE Vec3<T> knormalizeGlm(const Vec3<T>& a) { return a * (T(1) / ksqrt(kdot(a, a))); }
    template <class T> RAYT_HD_INLINE T kmaxComponent(const Vec3<T>& a) { return kmax(a.x, kmax(a.y, a.z)); }

} // namespace rayt::kernel
