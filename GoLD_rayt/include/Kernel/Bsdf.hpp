#pragma once
/**
 * @file Kernel/Bsdf.hpp
 * @brief The kernel renderer's materials by type, in float: what the path tracer calls (docs/GPU.md).
 *
 * The materials themselves are in Kernel/Materials.hpp (diffuse, conductor, dielectric) and Kernel/Layered.hpp (a coat over a
 * base). A layered material evaluates and samples by random walks, so eval and sample take the path's random number stream;
 * the other materials do not draw from it.
 */

#include "Kernel/Layered.hpp"
#include "Kernel/Rng.hpp"

namespace rayt::kernel {

    /// The path's stream as Kernel/Layered.hpp's random number source (next: one number; next2: two, in this order).
    struct PcgRand {
        Pcg32& rng;
        RAYT_HD_INLINE float next() { return rng.nextFloat(); }
        RAYT_HD_INLINE Vec2<float> next2() {
            const float a = rng.nextFloat();
            const float b = rng.nextFloat();
            return Vec2<float>(a, b);
        }
    };

    RAYT_HD_INLINE Vec3<float> materialAlbedo(const KMaterial& m) { return Vec3<float>(m.albedo[0], m.albedo[1], m.albedo[2]); }

    /// The base of a layered material.
    RAYT_HD_INLINE LayerBase<float> layerBase(const KMaterial& m) {
        return LayerBase<float>{ materialAlbedo(m), m.layered.conductorBase ? &m.conductor : nullptr };
    }

    /// The material scatters only in delta lobes (no light sampling at it, no density to weigh against).
    RAYT_HD_INLINE bool bsdfIsDelta(const KMaterial& m) { return m.type == kMatDielectric && m.dielectric.smooth; }

    /// albedo: a diffuse material's reflectance at the hit (its albedo times its texture there, if it has one).
    RAYT_HD_INLINE Vec3<float> bsdfEval(const KMaterial& m, const Vec3<float>& albedo, bool front, const Vec3<float>& wo, const Vec3<float>& wi, const Vec3<float>& ns,
        const Vec3<float>& ng, Pcg32& rng) {
        switch (m.type) {
        case kMatDiffuse: return diffuseEval(albedo, wi, ng);
        case kMatConductor: return conductorEval(m.conductor, wo, wi, ns, ng);
        case kMatDielectric: return dielectricEval(m.dielectric, front, wo, wi, ns, ng);
        case kMatLayered: {
            PcgRand rnd{ rng };
            return layeredEval(m.layered, layerBase(m), wo, wi, ns, ng, rnd);
        }
        default: return Vec3<float>(0.0f);
        }
    }

    RAYT_HD_INLINE float bsdfPdf(const KMaterial& m, bool front, const Vec3<float>& wo, const Vec3<float>& wi, const Vec3<float>& ns, const Vec3<float>& ng) {
        switch (m.type) {
        case kMatDiffuse: return diffusePdf(wi, ns, ng);
        case kMatConductor: return conductorPdf(m.conductor, wo, wi, ns, ng);
        case kMatDielectric: return dielectricPdf(m.dielectric, front, wo, wi, ns);
        case kMatLayered: return layeredPdf(m.layered, layerBase(m), wo, wi, ns, ng);
        default: return 0.0f;
        }
    }

    /// u1, u2 (and u3 for a rough dielectric's choice) are the sample's numbers; a layered material draws the rest of its walk from rng.
    RAYT_HD_INLINE bool bsdfSample(const KMaterial& m, const Vec3<float>& albedo, bool front, const Vec3<float>& wo, const Vec3<float>& ns, const Vec3<float>& ng,
        float u1, float u2, float u3, Pcg32& rng, BsdfSample<float>& out) {
        out.delta = false;
        switch (m.type) {
        case kMatDiffuse: return diffuseSample(albedo, ns, ng, u1, u2, out);
        case kMatConductor: return conductorSample(m.conductor, wo, ns, ng, u1, u2, out);
        case kMatDielectric: return dielectricSample(m.dielectric, front, wo, ns, ng, u1, u2, u3, out);
        case kMatLayered: {
            PcgRand rnd{ rng };
            return layeredSample(m.layered, layerBase(m), wo, ns, ng, Vec2<float>(u1, u2), rnd, out);
        }
        default: return false;
        }
    }

} // namespace rayt::kernel
