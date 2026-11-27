#pragma once
/**
 * @file Kernel/Materials.hpp
 * @brief The materials of the kernel renderer: eval, pdf and sample by material type (docs/GPU.md).
 *
 * Each material mirrors an engine class step by step (the same checks, frames and order of operations), so that with T = double
 * it gives the engine's values, which --scene gpu checks on random directions, and with T = float it is what the GPU runs.
 * Directions are in world space: wo toward the viewer, wi toward the light, ns the shading normal and ng the geometric normal,
 * both turned toward wo (SurfaceInteraction's n and gn).
 *
 *   diffuse    Materials/Lambertian.hpp: f = albedo / pi above the geometric surface, cosine-weighted sampling around ns.
 *   conductor  Materials/RoughConductor.{hpp,cpp}, isotropic: GGX with the height-correlated Smith term, the conductor's Fresnel
 *              term per channel, the Kulla-Conty lobe; samples the visible normals, or the Kulla-Conty lobe with probability
 *              1 - E(mu_o).
 */

#include "Kernel/KernelScene.hpp"
#include "Kernel/Microfacet.hpp"
#include "Kernel/Sampling.hpp"

namespace rayt::kernel {

    /// A local frame (s, t, n) built as the engine's frame::Frame::buildFromW: n normalized as glm does, then Duff's basis.
    template <class T> struct Frame3 {
        Vec3<T> s, t, n;
        RAYT_HD_INLINE static Frame3 fromNormal(const Vec3<T>& nn) {
            Frame3 f;
            f.n = knormalizeGlm(nn);
            orthonormalBasis(f.n, f.s, f.t);
            return f;
        }
        RAYT_HD_INLINE Vec3<T> toLocal(const Vec3<T>& a) const { return Vec3<T>(kdot(a, s), kdot(a, t), kdot(a, n)); }
        RAYT_HD_INLINE Vec3<T> toWorld(const Vec3<T>& a) const { return a.x * s + a.y * t + a.z * n; }
    };

    template <class T> struct BsdfSample {
        Vec3<T> wi;
        Vec3<T> f;
        T pdf;
    };

    // ---------------------------------------------------------------- diffuse (Lambertian)

    template <class T> RAYT_HD_INLINE Vec3<T> diffuseEval(const Vec3<T>& albedo, const Vec3<T>& wi, const Vec3<T>& ng) {
        if (kdot(ng, wi) <= T(0)) return Vec3<T>(0);
        return albedo * (T(1) / kPi<T>());
    }

    template <class T> RAYT_HD_INLINE T diffusePdf(const Vec3<T>& wi, const Vec3<T>& ns, const Vec3<T>& ng) {
        if (kdot(ng, wi) <= T(0)) return T(0);
        const T c = kdot(ns, wi);
        if (c <= T(0)) return T(0);
        return c * (T(1) / kPi<T>());
    }

    template <class T> RAYT_HD_INLINE bool diffuseSample(const Vec3<T>& albedo, const Vec3<T>& ns, const Vec3<T>& ng, T u1, T u2, BsdfSample<T>& out) {
        const Frame3<T> fr = Frame3<T>::fromNormal(ns);
        const Vec3<T> local = sampleCosineHemisphere(u1, u2);
        out.wi = fr.toWorld(local);
        if (kdot(ng, out.wi) <= T(0)) return false;
        out.pdf = local.z * (T(1) / kPi<T>());
        out.f = albedo * (T(1) / kPi<T>());
        return true;
    }

    // ---------------------------------------------------------------- conductor (RoughConductor, isotropic)

    /// E(mu) of the conductor's alpha, linear between the table's nodes (RoughConductor::eAt).
    template <class T> RAYT_HD_INLINE T conductorE(const ConductorParams<T>& c, T mu) {
        const T x = ksaturate(mu) * T(kAlbedoMu - 1);
        const int i = kmin(int(x), kAlbedoMu - 2);
        const T f = x - T(i);
        return c.eRow[i] * (T(1) - f) + c.eRow[i + 1] * f;
    }

    template <class T> RAYT_HD Vec3<T> conductorEval(const ConductorParams<T>& c, const Vec3<T>& wo, const Vec3<T>& wi, const Vec3<T>& ns, const Vec3<T>& ng) {
        if (kdot(ng, wo) <= T(0) || kdot(ng, wi) <= T(0)) return Vec3<T>(0);
        const T cosO = kmax(T(0), kdot(ns, wo));
        const T cosI = kmax(T(0), kdot(ns, wi));
        if (cosI <= T(0) || cosO <= T(0)) return Vec3<T>(0);
        Vec3<T> wh = wi + wo;
        const T whLen2 = kdot(wh, wh);
        if (whLen2 <= T(0)) return Vec3<T>(0);
        wh *= ksafeRecip(ksafeSqrt(whLen2));
        const Frame3<T> fr = Frame3<T>::fromNormal(ns);
        const Vec3<T> woL = fr.toLocal(wo), wiL = fr.toLocal(wi), whL = fr.toLocal(wh);
        const T D = ggxD(whL, c.alpha, c.alpha);
        const T G = ggxG(woL, wiL, c.alpha, c.alpha);
        const T cosD = kabs(kdot(wi, wh));
        const Vec3<T> F(fresnelConductor(cosD, c.eta[0], c.k[0]), fresnelConductor(cosD, c.eta[1], c.k[1]), fresnelConductor(cosD, c.eta[2], c.k[2]));
        const T den = T(4) * cosI * cosO;
        if (den < T(1e-6)) return Vec3<T>(0);
        Vec3<T> f = (D * G * F) / den;
        if (c.multipleScattering) {
            const T scale = (T(1) - conductorE(c, cosO)) * (T(1) - conductorE(c, cosI)) / (kPi<T>() * (T(1) - c.eAvg));
            f += Vec3<T>(c.fms[0], c.fms[1], c.fms[2]) * scale;
        }
        return f;
    }

    template <class T> RAYT_HD T conductorPdf(const ConductorParams<T>& c, const Vec3<T>& wo, const Vec3<T>& wi, const Vec3<T>& ns, const Vec3<T>& ng) {
        if (kdot(ng, wi) <= T(0) || kdot(ng, wo) <= T(0)) return T(0);
        Vec3<T> wh = wo + wi;
        const T whLen2 = kdot(wh, wh);
        if (whLen2 <= T(0)) return T(0);
        wh = wh * ksafeRecip(ksafeSqrt(whLen2));
        const Frame3<T> fr = Frame3<T>::fromNormal(ns);
        const Vec3<T> woL = fr.toLocal(wo), whL = fr.toLocal(wh);
        const T pdfWh = ggxPdfVisibleNormal(woL, whL, c.alpha, c.alpha);
        const T dotWoWh = kabs(kdot(woL, whL));
        const T pdfGgx = (dotWoWh > T(0)) ? pdfWh / (T(4) * dotWoWh) : T(0);
        if (!c.multipleScattering) return pdfGgx;
        const T pMs = T(1) - conductorE(c, woL.z);
        const T pdfMs = kmax(kdot(ns, wi), T(0)) / kPi<T>();
        return (T(1) - pMs) * pdfGgx + pMs * pdfMs;
    }

    template <class T> RAYT_HD bool conductorSample(const ConductorParams<T>& c, const Vec3<T>& wo, const Vec3<T>& ns, const Vec3<T>& ng, T u1, T u2,
        BsdfSample<T>& out) {
        if (kdot(ng, wo) <= T(0)) return false;
        const Frame3<T> fr = Frame3<T>::fromNormal(ns);
        const Vec3<T> woL = fr.toLocal(wo);
        const T pMs = c.multipleScattering ? (T(1) - conductorE(c, woL.z)) : T(0);
        bool useMs = false;
        if (pMs > T(0)) {
            if (u1 < pMs) { useMs = true; u1 = u1 / pMs; }
            else u1 = (u1 - pMs) / (T(1) - pMs);
        }
        if (useMs) {
            const T r = ksafeSqrt(u1);
            const T phi = kTwoPi<T>() * u2;
            out.wi = fr.toWorld(Vec3<T>(r * kcos(phi), r * ksin(phi), ksafeSqrt(T(1) - u1)));
        }
        else {
            const Vec3<T> wh = fr.toWorld(ggxSampleVisibleNormal(woL, c.alpha, c.alpha, u1, u2));
            const Vec3<T> v = -wo;
            out.wi = v - T(2) * kdot(v, wh) * wh;   // math::reflectOutward
        }
        if (kdot(ng, out.wi) <= T(0)) return false;
        const Vec3<T> wiL = fr.toLocal(out.wi);
        if (woL.z == T(0) || wiL.z == T(0)) return false;
        out.pdf = conductorPdf(c, wo, out.wi, ns, ng);
        out.f = conductorEval(c, wo, out.wi, ns, ng);
        if (out.pdf <= T(1e-6) || !kisFinite(out.f.x) || !kisFinite(out.f.y) || !kisFinite(out.f.z)) return false;
        return true;
    }

    // ---------------------------------------------------------------- by type (float, what the kernel renderer calls)

    RAYT_HD_INLINE Vec3<float> materialAlbedo(const KMaterial& m) { return Vec3<float>(m.albedo[0], m.albedo[1], m.albedo[2]); }

    RAYT_HD_INLINE Vec3<float> bsdfEval(const KMaterial& m, const Vec3<float>& wo, const Vec3<float>& wi, const Vec3<float>& ns, const Vec3<float>& ng) {
        switch (m.type) {
        case kMatDiffuse: return diffuseEval(materialAlbedo(m), wi, ng);
        case kMatConductor: return conductorEval(m.conductor, wo, wi, ns, ng);
        default: return Vec3<float>(0.0f);
        }
    }

    RAYT_HD_INLINE float bsdfPdf(const KMaterial& m, const Vec3<float>& wo, const Vec3<float>& wi, const Vec3<float>& ns, const Vec3<float>& ng) {
        switch (m.type) {
        case kMatDiffuse: return diffusePdf(wi, ns, ng);
        case kMatConductor: return conductorPdf(m.conductor, wo, wi, ns, ng);
        default: return 0.0f;
        }
    }

    RAYT_HD_INLINE bool bsdfSample(const KMaterial& m, const Vec3<float>& wo, const Vec3<float>& ns, const Vec3<float>& ng, float u1, float u2, BsdfSample<float>& out) {
        switch (m.type) {
        case kMatDiffuse: return diffuseSample(materialAlbedo(m), ns, ng, u1, u2, out);
        case kMatConductor: return conductorSample(m.conductor, wo, ns, ng, u1, u2, out);
        default: return false;
        }
    }

} // namespace rayt::kernel
