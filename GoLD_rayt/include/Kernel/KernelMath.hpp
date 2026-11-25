#pragma once
/**
 * @file Kernel/KernelMath.hpp
 * @brief Scalar math and vector types for kernel code (CPU and GPU).
 *
 * The overloads call sqrtf / cosf / ... for float and sqrt / cos / ... for double. These names exist in the C library on
 * the CPU and as device functions in CUDA, so a call resolves on both sides. min / max keep the argument order of
 * std::min / std::max, so kernel code that mirrors an engine function gives bit-identical doubles on the CPU.
 */

#include <math.h>
#include <glm/glm.hpp>

#include "Kernel/HostDevice.hpp"

namespace rayt::kernel {

    template <class T> using Vec2 = glm::vec<2, T, glm::defaultp>;
    template <class T> using Vec3 = glm::vec<3, T, glm::defaultp>;

    template <class T> RAYT_HD_INLINE constexpr T kPi() { return T(3.14159265358979323846); }
    template <class T> RAYT_HD_INLINE constexpr T kTwoPi() { return T(2) * kPi<T>(); }

    RAYT_HD_INLINE float ksqrt(float x) { return sqrtf(x); }
    RAYT_HD_INLINE double ksqrt(double x) { return sqrt(x); }
    RAYT_HD_INLINE float kcos(float x) { return cosf(x); }
    RAYT_HD_INLINE double kcos(double x) { return cos(x); }
    RAYT_HD_INLINE float ksin(float x) { return sinf(x); }
    RAYT_HD_INLINE double ksin(double x) { return sin(x); }
    RAYT_HD_INLINE float kcopysign(float x, float s) { return copysignf(x, s); }
    RAYT_HD_INLINE double kcopysign(double x, double s) { return copysign(x, s); }
    RAYT_HD_INLINE float kabs(float x) { return fabsf(x); }
    RAYT_HD_INLINE double kabs(double x) { return fabs(x); }

    /// std::max(a, b): b if a < b, else a.
    template <class T> RAYT_HD_INLINE T kmax(T a, T b) { return (a < b) ? b : a; }
    /// std::min(a, b): b if b < a, else a.
    template <class T> RAYT_HD_INLINE T kmin(T a, T b) { return (b < a) ? b : a; }
    /// glm::clamp(x, 0, 1) = min(max(x, 0), 1).
    template <class T> RAYT_HD_INLINE T ksaturate(T x) { return kmin(kmax(x, T(0)), T(1)); }
    /// sqrt(max(x, 0)), as math::safe_sqrt.
    template <class T> RAYT_HD_INLINE T ksafeSqrt(T x) { return ksqrt(kmax(x, T(0))); }

} // namespace rayt::kernel
