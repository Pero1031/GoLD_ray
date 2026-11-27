#pragma once
/**
 * @file Kernel/Materials.hpp
 * @brief The materials of the kernel renderer: eval, pdf and sample by material type (docs/GPU.md).
 *
 * Each material mirrors an engine class step by step (the same checks, frames and order of operations), so that with T = double
 * it gives the engine's values, which --scene gpu checks on random directions, and with T = float it is what the GPU runs.
 * Directions are in world space: wo toward the viewer, wi toward the light, ns the shading normal and ng the geometric normal,
 * both turned toward wo (SurfaceInteraction's n and gn); front says that the ray came from the outside (the side of the face
 * normal, SurfaceInteraction::frontFace).
 *
 *   diffuse    Materials/Lambertian.hpp: f = albedo / pi above the geometric surface, cosine-weighted sampling around ns.
 *   conductor  Materials/RoughConductor.{hpp,cpp}, isotropic: GGX with the height-correlated Smith term, the conductor's Fresnel
 *              term per channel, the Kulla-Conty lobe; samples the visible normals, or the Kulla-Conty lobe with probability
 *              1 - E(mu_o).
 *   dielectric Materials/Dielectric.{hpp,cpp} without dispersion or a medium inside: smooth (a delta reflection or refraction,
 *              chosen by the Fresnel term) or rough (GGX reflection and transmission, Walter et al. 2007).
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

    /// A sampled direction. delta: a perfectly specular event (f is then the weight that multiplies the throughput, as the engine's
    /// BSDFSample of a specular lobe; pdf is the probability of the choice).
    template <class T> struct BsdfSample {
        Vec3<T> wi;
        Vec3<T> f;
        T pdf;
        bool delta;
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
        out.delta = false;
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
        out.delta = false;
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

    // ---------------------------------------------------------------- dielectric (Dielectric: no dispersion, isotropic, nothing inside)

    /// math::refractIncident(-wo, n, eta): false at total internal reflection.
    template <class T> RAYT_HD_INLINE bool refractOutward(const Vec3<T>& wo, const Vec3<T>& n, T eta, Vec3<T>& wi) {
        const Vec3<T> v = -wo;
        const T cosi = kdot(-v, n);
        const T sin2t = ksqr(eta) * (T(1) - ksqr(cosi));
        if (sin2t > T(1)) return false;
        const T cost = ksafeSqrt(T(1) - sin2t);
        wi = eta * v + (eta * cosi - cost) * n;
        return true;
    }

    template <class T> RAYT_HD_INLINE T dielectricFresnel(bool front, T cosAbs, T ior) {
        return front ? fresnelDielectric(cosAbs, T(1), ior) : fresnelDielectric(cosAbs, ior, T(1));
    }

    template <class T> RAYT_HD Vec3<T> dielectricEval(const DielectricParams<T>& c, bool front, const Vec3<T>& wo, const Vec3<T>& wi, const Vec3<T>& ns,
        const Vec3<T>& ng) {
        if (c.smooth) return Vec3<T>(0);
        const T ior = c.ior;
        const T cosO = kdot(ns, wo), cosI = kdot(ns, wi);
        const bool reflection = cosO * cosI > T(0);
        if (reflection) {
            if (kdot(ng, wo) <= T(0) || kdot(ng, wi) <= T(0)) return Vec3<T>(0);
        }
        else if (kdot(ng, wo) * kdot(ng, wi) > T(0)) return Vec3<T>(0);
        const T etap = front ? ior : (T(1) / ior);
        Vec3<T> wh = reflection ? wo + wi : wo + wi * etap;
        if (kdot(wh, wh) == T(0)) return Vec3<T>(0);
        wh = knormalizeGlm(wh);
        if (kdot(wh, ns) < T(0)) wh = -wh;
        if (kdot(wh, wi) * cosI <= T(0) || kdot(wh, wo) * cosO <= T(0)) return Vec3<T>(0);
        const Frame3<T> fr = Frame3<T>::fromNormal(ns);
        const Vec3<T> woL = fr.toLocal(wo), wiL = fr.toLocal(wi), whL = fr.toLocal(wh);
        const T D = ggxD(whL, c.alpha, c.alpha);
        const T G = ggxG(woL, wiL, c.alpha, c.alpha);
        const T F = dielectricFresnel(front, kabs(kdot(wo, wh)), ior);
        if (reflection) {
            const T den = kabs(T(4) * cosI * cosO);
            if (den < T(1e-8)) return Vec3<T>(0);
            return Vec3<T>(D * G * F / den);
        }
        const T dotWiWh = kdot(wi, wh), dotWoWh = kdot(wo, wh);
        const T sqrtDen = dotWiWh * etap + dotWoWh;
        const T den = ksqr(sqrtDen) * cosI * cosO;
        if (kabs(den) < T(1e-8)) return Vec3<T>(0);
        return Vec3<T>(D * G * (T(1) - F) * kabs(dotWiWh * dotWoWh / den));
    }

    template <class T> RAYT_HD T dielectricPdf(const DielectricParams<T>& c, bool front, const Vec3<T>& wo, const Vec3<T>& wi, const Vec3<T>& ns) {
        if (c.smooth) return T(0);
        const T ior = c.ior;
        const T cosO = kdot(ns, wo), cosI = kdot(ns, wi);
        const bool reflection = cosO * cosI > T(0);
        const T etap = front ? ior : (T(1) / ior);
        Vec3<T> wh = reflection ? wo + wi : wo + wi * etap;
        if (kdot(wh, wh) == T(0)) return T(0);
        wh = knormalizeGlm(wh);
        if (kdot(wh, ns) < T(0)) wh = -wh;
        if (kdot(wh, wi) * cosI <= T(0) || kdot(wh, wo) * cosO <= T(0)) return T(0);
        const T dotWoWh = kdot(wo, wh);
        if (dotWoWh == T(0)) return T(0);
        const T F = dielectricFresnel(front, kabs(dotWoWh), ior);
        const Frame3<T> fr = Frame3<T>::fromNormal(ns);
        Vec3<T> woL = fr.toLocal(wo);
        if (woL.z < T(0)) woL = -woL;
        const T pdfWh = ggxPdfVisibleNormal(woL, fr.toLocal(wh), c.alpha, c.alpha);
        if (reflection) return (pdfWh / (T(4) * kabs(dotWoWh))) * F;
        const T dotWiWh = kdot(wi, wh);
        const T sqrtDen = dotWiWh * etap + dotWoWh;
        if (sqrtDen == T(0)) return T(0);
        const T dwhDwi = kabs(dotWiWh) * ksqr(etap) / ksqr(sqrtDen);
        return (pdfWh * dwhDwi) * (T(1) - F);
    }

    /// u1, u2 for the micro-normal (u1 alone for the smooth choice), u3 for the choice between reflection and refraction of a rough
    /// surface (the engine draws that one inside Dielectric::sample).
    template <class T> RAYT_HD bool dielectricSample(const DielectricParams<T>& c, bool front, const Vec3<T>& wo, const Vec3<T>& ns, const Vec3<T>& ng,
        T u1, T u2, T u3, BsdfSample<T>& out) {
        const T ior = c.ior;
        const T cosO = kdot(ns, wo);
        const T eta = front ? (T(1) / ior) : ior;
        const T etap = front ? ior : (T(1) / ior);
        if (c.smooth) {
            const T F = dielectricFresnel(front, kabs(cosO), ior);
            out.delta = true;
            if (u1 < F) {
                const Vec3<T> v = -wo;
                out.wi = v - T(2) * kdot(v, ns) * ns;   // math::reflectOutward
                out.pdf = F;
                out.f = Vec3<T>(1);
            }
            else {
                if (!refractOutward(wo, ns, eta, out.wi)) return false;
                out.pdf = T(1) - F;
                out.f = Vec3<T>(1) / ksqr(etap);   // radiance is compressed into the smaller solid angle
            }
            return true;
        }
        out.delta = false;
        const Frame3<T> fr = Frame3<T>::fromNormal(ns);
        Vec3<T> woS = fr.toLocal(wo);
        if (woS.z < T(0)) woS = -woS;
        Vec3<T> whL = ggxSampleVisibleNormal(woS, c.alpha, c.alpha, u1, u2);
        Vec3<T> wh = fr.toWorld(whL);
        if (kdot(wh, wo) < T(0)) { wh = -wh; whL = -whL; }
        const T dotWoWh = kdot(wo, wh);
        if (dotWoWh == T(0)) return false;
        const T F = dielectricFresnel(front, kabs(dotWoWh), ior);
        if (u3 < F) {
            const Vec3<T> v = -wo;
            out.wi = v - T(2) * kdot(v, wh) * wh;   // math::reflectIncident(-wo, wh)
            if (kdot(ng, wo) <= T(0) || kdot(ng, out.wi) <= T(0)) return false;
            out.f = dielectricEval(c, front, wo, out.wi, ns, ng);
            const T pdfWh = ggxPdfVisibleNormal(woS, whL, c.alpha, c.alpha);
            out.pdf = (pdfWh / (T(4) * kabs(dotWoWh))) * F;
        }
        else {
            if (!refractOutward(wo, wh, eta, out.wi)) return false;
            if (kdot(ns, out.wi) * cosO > T(0)) return false;
            out.f = dielectricEval(c, front, wo, out.wi, ns, ng);
            const T dotWiWh = kdot(out.wi, wh);
            const T sqrtDen = dotWiWh * etap + dotWoWh;
            if (sqrtDen == T(0)) return false;
            const T dwhDwi = kabs(dotWiWh) * ksqr(etap) / ksqr(sqrtDen);
            const T pdfWh = ggxPdfVisibleNormal(woS, whL, c.alpha, c.alpha);
            out.pdf = pdfWh * dwhDwi * (T(1) - F);
        }
        if (out.pdf < T(1e-8) || !kisFinite(out.f.x) || !kisFinite(out.f.y) || !kisFinite(out.f.z)) return false;
        return true;
    }

    // ---------------------------------------------------------------- by type (float, what the kernel renderer calls)

    RAYT_HD_INLINE Vec3<float> materialAlbedo(const KMaterial& m) { return Vec3<float>(m.albedo[0], m.albedo[1], m.albedo[2]); }

    /// The material scatters only in delta lobes (no light sampling at it, no density to weigh against).
    RAYT_HD_INLINE bool bsdfIsDelta(const KMaterial& m) { return m.type == kMatDielectric && m.dielectric.smooth; }

    RAYT_HD_INLINE Vec3<float> bsdfEval(const KMaterial& m, bool front, const Vec3<float>& wo, const Vec3<float>& wi, const Vec3<float>& ns, const Vec3<float>& ng) {
        switch (m.type) {
        case kMatDiffuse: return diffuseEval(materialAlbedo(m), wi, ng);
        case kMatConductor: return conductorEval(m.conductor, wo, wi, ns, ng);
        case kMatDielectric: return dielectricEval(m.dielectric, front, wo, wi, ns, ng);
        default: return Vec3<float>(0.0f);
        }
    }

    RAYT_HD_INLINE float bsdfPdf(const KMaterial& m, bool front, const Vec3<float>& wo, const Vec3<float>& wi, const Vec3<float>& ns, const Vec3<float>& ng) {
        switch (m.type) {
        case kMatDiffuse: return diffusePdf(wi, ns, ng);
        case kMatConductor: return conductorPdf(m.conductor, wo, wi, ns, ng);
        case kMatDielectric: return dielectricPdf(m.dielectric, front, wo, wi, ns);
        default: return 0.0f;
        }
    }

    RAYT_HD_INLINE bool bsdfSample(const KMaterial& m, bool front, const Vec3<float>& wo, const Vec3<float>& ns, const Vec3<float>& ng, float u1, float u2,
        float u3, BsdfSample<float>& out) {
        out.delta = false;
        switch (m.type) {
        case kMatDiffuse: return diffuseSample(materialAlbedo(m), ns, ng, u1, u2, out);
        case kMatConductor: return conductorSample(m.conductor, wo, ns, ng, u1, u2, out);
        case kMatDielectric: return dielectricSample(m.dielectric, front, wo, ns, ng, u1, u2, u3, out);
        default: return false;
        }
    }

} // namespace rayt::kernel
