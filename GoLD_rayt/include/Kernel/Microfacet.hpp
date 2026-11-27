#pragma once
/**
 * @file Kernel/Microfacet.hpp
 * @brief The GGX microfacet distribution and the Fresnel reflectance of a conductor, for the engine and the kernel (docs/GPU.md).
 *
 * One implementation: the engine's GGXDistribution (Microfacet/GGX.cpp) and fresnel::fresnelConductor (Core/Fresnel.hpp) call
 * these with T = double, the kernel renderer with T = float. The formulas and the order of operations are those the engine
 * had before (Walter et al. 2007; Heitz 2014 for the height-correlated Smith term; Heitz 2018 for sampling the visible normals),
 * so moving them here left the engine's pictures the same bits. Directions are in the local frame (normal +z).
 */

#include "Kernel/KernelMath.hpp"

namespace rayt::kernel {

    /// GGX normal distribution D(wh) for roughness (ax, ay); 0 below the surface.
    template <class T> RAYT_HD_INLINE T ggxD(const Vec3<T>& wh, T ax, T ay) {
        if (wh.z <= T(0)) return T(0);
        const T ax2 = ksqr(ax), ay2 = ksqr(ay);
        const T e = ksqr(wh.x) / ax2 + ksqr(wh.y) / ay2 + ksqr(wh.z);
        return T(1) / (kPi<T>() * ax * ay * ksqr(e));
    }

    /// Smith's Lambda(w) of GGX: (-1 + sqrt(1 + alpha^2 tan^2 theta)) / 2, alpha along the azimuth of w.
    template <class T> RAYT_HD_INLINE T ggxLambda(const Vec3<T>& w, T ax, T ay) {
        const T absCosTheta = kabs(w.z);
        if (absCosTheta <= T(0)) return kinfinity<T>();
        const T sin2Theta = kmax(T(0), T(1) - ksqr(absCosTheta));
        if (sin2Theta <= T(0)) return T(0);
        const T wxy2 = ksqr(w.x) + ksqr(w.y);
        const T ax2 = ksqr(ax), ay2 = ksqr(ay);
        const T alpha2 = (wxy2 > T(0)) ? (ksqr(w.x) * ax2 + ksqr(w.y) * ay2) / wxy2 : ax2;
        return T(0.5) * (ksafeSqrt(ksqr(absCosTheta) + alpha2 * sin2Theta) / absCosTheta - T(1));
    }

    /// Smith's masking of one direction, 1 / (1 + Lambda).
    template <class T> RAYT_HD_INLINE T ggxG1(const Vec3<T>& w, T ax, T ay) { return T(1) / (T(1) + ggxLambda(w, ax, ay)); }

    /// The height-correlated masking-shadowing term, 1 / (1 + Lambda(wo) + Lambda(wi)).
    template <class T> RAYT_HD_INLINE T ggxG(const Vec3<T>& wo, const Vec3<T>& wi, T ax, T ay) {
        return T(1) / (T(1) + ggxLambda(wo, ax, ay) + ggxLambda(wi, ax, ay));
    }

    /// A micro-normal drawn from the normals visible from wo (Heitz 2018), from two uniform numbers.
    template <class T> RAYT_HD_INLINE Vec3<T> ggxSampleVisibleNormal(const Vec3<T>& wo, T ax, T ay, T u1, T u2) {
        if (wo.z <= T(0)) return Vec3<T>(0, 0, 1);
        const Vec3<T> vh = knormalizeGlm(Vec3<T>(ax * wo.x, ay * wo.y, wo.z));
        const T lenSq = ksqr(vh.x) + ksqr(vh.y);
        const Vec3<T> t1v = lenSq > T(0) ? Vec3<T>(-vh.y, vh.x, T(0)) * ksafeRecip(ksafeSqrt(lenSq)) : Vec3<T>(1, 0, 0);
        const Vec3<T> t2v = kcross(vh, t1v);
        const T r = ksafeSqrt(u1);
        const T phi = kTwoPi<T>() * u2;
        const T t1 = r * kcos(phi);
        T t2 = r * ksin(phi);
        const T s = T(0.5) * (T(1) + vh.z);
        t2 = (T(1) - s) * ksafeSqrt(T(1) - ksqr(t1)) + s * t2;
        const Vec3<T> nh = t1 * t1v + t2 * t2v + ksafeSqrt(kmax(T(0), T(1) - ksqr(t1) - ksqr(t2))) * vh;
        return knormalizeGlm(Vec3<T>(ax * nh.x, ay * nh.y, kmax(T(0), nh.z)));
    }

    /// Density of ggxSampleVisibleNormal over the micro-normal wh (solid angle of wh): G1(wo) |wo . wh| D(wh) / wo.z.
    template <class T> RAYT_HD_INLINE T ggxPdfVisibleNormal(const Vec3<T>& wo, const Vec3<T>& wh, T ax, T ay) {
        if (wo.z <= T(0) || wh.z <= T(0)) return T(0);
        return ggxG1(wo, ax, ay) * kabs(kdot(wo, wh)) * ggxD(wh, ax, ay) * ksafeRecip(wo.z);
    }

    /// Unpolarized Fresnel reflectance of a conductor of complex index eta + i k (one wavelength or RGB channel); |cos| clamped to 1.
    template <class T> RAYT_HD_INLINE T fresnelConductor(T cosThetaI, T eta, T k) {
        cosThetaI = ksaturate(kabs(cosThetaI));
        const T cos2 = ksqr(cosThetaI), sin2 = T(1) - cos2;
        const T eta2 = eta * eta, k2 = k * k;
        const T t0 = eta2 - k2 - sin2;
        const T a2plusb2 = ksqrt(t0 * t0 + T(4) * eta2 * k2);
        const T t1 = a2plusb2 + cos2;
        const T a = ksafeSqrt(T(0.5) * (a2plusb2 + t0));
        const T t2 = T(2) * cosThetaI * a;
        const T rs = (t1 - t2) / (t1 + t2);
        const T t3 = cos2 * a2plusb2 + sin2 * sin2;
        const T t4 = t2 * sin2;
        const T rp = rs * (t3 - t4) / (t3 + t4);
        return T(0.5) * (rs + rp);
    }

} // namespace rayt::kernel
