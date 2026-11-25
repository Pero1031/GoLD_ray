#pragma once
/**
 * @file Kernel/Sampling.hpp
 * @brief Direction sampling, orthonormal bases and the dielectric Fresnel term for kernel code.
 *
 * Each function is the same formula, in the same order of operations, as its engine counterpart (named in its
 * comment), so in double precision on the CPU both give the same bits; --scene gpu checks that.
 */

#include "Kernel/KernelMath.hpp"

namespace rayt::kernel {

    /// Cosine-weighted direction around +z (pdf cos(theta) / pi); the engine's CosineSampleHemisphere(u), Core/Sampling.hpp.
    template <class T> RAYT_HD_INLINE Vec3<T> sampleCosineHemisphere(T u1, T u2) {
        const T z = ksqrt(kmax(T(0), T(1) - u2));
        const T phi = kTwoPi<T>() * u1;
        const T sinTheta = ksqrt(u2);
        return Vec3<T>(kcos(phi) * sinTheta, ksin(phi) * sinTheta, z);
    }

    /// Uniform direction on the unit sphere (pdf 1 / (4 pi)); the engine's UniformSampleSphere(u), Core/Sampling.hpp.
    template <class T> RAYT_HD_INLINE Vec3<T> sampleUniformSphere(T u1, T u2) {
        const T z = T(1) - T(2) * u1;
        const T r = ksafeSqrt(T(1) - z * z);
        const T phi = kTwoPi<T>() * u2;
        return Vec3<T>(r * kcos(phi), r * ksin(phi), z);
    }

    /// Orthonormal basis (s, t, n) around a unit normal n (Duff et al. 2017); the engine's frame::makeOrthonormalBasis, Geometry/Frame.hpp.
    template <class T> RAYT_HD_INLINE void orthonormalBasis(const Vec3<T>& n, Vec3<T>& s, Vec3<T>& t) {
        const T sign = kcopysign(T(1), n.z);
        const T a = -T(1) / (sign + n.z);
        const T b = n.x * n.y * a;
        s = Vec3<T>(T(1) + sign * n.x * n.x * a, sign * b, -sign * n.x);
        t = Vec3<T>(b, sign + n.y * n.y * a, -n.y);
    }

    /// Unpolarized Fresnel reflectance of a dielectric interface, cosThetaI clamped to [0, 1]; the engine's
    /// fresnel::fresnelDielectric(cosThetaI, etaI, etaT), Core/Fresnel.hpp.
    template <class T> RAYT_HD_INLINE T fresnelDielectric(T cosThetaI, T etaI, T etaT) {
        cosThetaI = ksaturate(cosThetaI);
        const T sinThetaI = ksafeSqrt(T(1) - cosThetaI * cosThetaI);
        const T sinThetaT = (etaI / etaT) * sinThetaI;
        if (sinThetaT >= T(1)) return T(1);
        const T cosThetaT = ksqrt(kmax(T(0), T(1) - sinThetaT * sinThetaT));
        const T rParl = ((etaT * cosThetaI) - (etaI * cosThetaT)) / ((etaT * cosThetaI) + (etaI * cosThetaT));
        const T rPerp = ((etaI * cosThetaI) - (etaT * cosThetaT)) / ((etaI * cosThetaI) + (etaT * cosThetaT));
        return (rParl * rParl + rPerp * rPerp) * T(0.5);
    }

} // namespace rayt::kernel
