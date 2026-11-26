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

#if !defined(RAYT_GPU_CODE)
    // CPU only: device code that calls one of these does not compile.
    inline double ksqrt(double x) { return sqrt(x); }
    inline double kfma(double a, double b, double c) { return fma(a, b, c); }
    inline double kcos(double x) { return cos(x); }
    inline double ksin(double x) { return sin(x); }
    inline double kcopysign(double x, double s) { return copysign(x, s); }
    inline double kabs(double x) { return fabs(x); }
    inline double kfloor(double x) { return floor(x); }
#else
    // Device code: a double argument picks these deleted overloads (instead of converting to float) and fails to compile.
    RAYT_HD double ksqrt(double) = delete;
    RAYT_HD double kfma(double, double, double) = delete;
    RAYT_HD double kcos(double) = delete;
    RAYT_HD double ksin(double) = delete;
    RAYT_HD double kcopysign(double, double) = delete;
    RAYT_HD double kabs(double) = delete;
    RAYT_HD double kfloor(double) = delete;
#endif

    /// std::max(a, b): b if a < b, else a.
    template <class T> RAYT_HD_INLINE T kmax(T a, T b) { return (a < b) ? b : a; }
    /// std::min(a, b): b if b < a, else a.
    template <class T> RAYT_HD_INLINE T kmin(T a, T b) { return (b < a) ? b : a; }
    /// glm::clamp(x, 0, 1) = min(max(x, 0), 1).
    template <class T> RAYT_HD_INLINE T ksaturate(T x) { return kmin(kmax(x, T(0)), T(1)); }
    /// sqrt(max(x, 0)), as math::safe_sqrt.
    template <class T> RAYT_HD_INLINE T ksafeSqrt(T x) { return ksqrt(kmax(x, T(0))); }

} // namespace rayt::kernel
