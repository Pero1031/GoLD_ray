#pragma once
/**
 * @file Kernel/PathTracer.hpp
 * @brief The kernel path tracer: one camera path of one pixel sample, in float, on the CPU or the GPU (docs/GPU.md).
 *
 * The rules are the engine's PathIntegrator restricted to what KScene holds (the materials of Kernel/Bsdf.hpp, emitting
 * quads, an environment map), so the two converge to the same picture: a path has at most maxDepth surface vertices; each
 * non-delta vertex samples one light with a shadow ray (the environment with probability envSelect, by its luminance; else a
 * light chosen by power: a uniform point on a quad, the cone of a sphere, Kernel/Lights.hpp), weighted against the BSDF's density by the power heuristic, and then samples
 * the BSDF; an emitter or the environment found by a BSDF ray is weighted against that light sample the same way; Russian
 * roulette from the fourth vertex on with survival max(beta) clamped to [0.05, 0.95]. (The engine chooses among its lights by
 * another rule; any choice is unbiased.)
 *
 * Ray casting is the caller's: a Tracer has
 *     KHit closest(const V3& o, const V3& d, float tmin, float tmax) const;   // prim = kNoHit when nothing is hit
 *     bool occluded(const V3& o, const V3& d, float tmin, float tmax) const;
 * (Embree on the CPU, OptiX on the GPU). Only the triangle of a hit is taken from the tracer: the hit point is recomputed
 * from it with the kernel's own formula (triangleBarycentrics), and rays leave a surface from an origin moved off it in
 * integer steps of the float representation (offsetRayOrigin), so both devices start every ray from the same bits as long as
 * their tracers report the same triangle.
 */

#include <string.h>

#include "Kernel/Bsdf.hpp"
#include "Kernel/Camera.hpp"
#include "Kernel/Environment.hpp"
#include "Kernel/KernelScene.hpp"
#include "Kernel/Lights.hpp"
#include "Kernel/Media.hpp"
#include "Kernel/Rng.hpp"
#include "Kernel/Sampling.hpp"
#include "Kernel/Textures.hpp"

namespace rayt::kernel {

    using V3 = Vec3<float>;

    constexpr float kFloatMax = 3.402823466e38f;
    constexpr float kInvPi = 0.318309886183790671538f;

    RAYT_HD_INLINE V3 toV3(const KFloat3& a) { return V3(a.x, a.y, a.z); }

    /**
     * @brief The origin of a ray that leaves a surface at p on the side of the unit normal n.
     *
     * Wachter and Binder, "A Fast and Robust Method for Avoiding Self-Intersection" (Ray Tracing Gems, ch. 6), which Cycles
     * uses: p is moved along n by a number of float steps (ulps) proportional to n, so the offset grows with the magnitude of
     * the coordinates as the error of the hit point does; near the origin, where ulps are tiny, by a fixed small distance.
     */
    RAYT_HD_INLINE V3 offsetRayOrigin(const V3& p, const V3& n) {
        const float kOrigin = 1.0f / 32.0f, kFloatScale = 1.0f / 65536.0f, kIntScale = 256.0f;
        const int32_t ox = int32_t(kIntScale * n.x), oy = int32_t(kIntScale * n.y), oz = int32_t(kIntScale * n.z);
        const float px = kIntAsFloat(kFloatAsInt(p.x) + (p.x < 0.0f ? -ox : ox));
        const float py = kIntAsFloat(kFloatAsInt(p.y) + (p.y < 0.0f ? -oy : oy));
        const float pz = kIntAsFloat(kFloatAsInt(p.z) + (p.z < 0.0f ? -oz : oz));
        return V3(kabs(p.x) < kOrigin ? p.x + kFloatScale * n.x : px,
                  kabs(p.y) < kOrigin ? p.y + kFloatScale * n.y : py,
                  kabs(p.z) < kOrigin ? p.z + kFloatScale * n.z : pz);
    }

    /**
     * @brief Barycentrics (b1, b2) of the point where the ray o + t d meets the plane of triangle (v0, v1, v2) (Moller-Trumbore),
     *        clamped onto the triangle.
     *
     * The tracers (Embree, OptiX) report the triangle a ray hits, but each computes the barycentrics with its own arithmetic, so
     * the same hit would give hit points that differ in the last bits, and everything after it would differ a little. The
     * kernel recomputes them with this one formula on both devices (as Cycles re-intersects a triangle to get a hit point it
     * controls), so the same triangle gives the same point. The tracer has already decided that the ray hits the triangle;
     * clamping only keeps a point found at an edge on it.
     */
    RAYT_HD_INLINE void triangleBarycentrics(const V3& o, const V3& d, const V3& v0, const V3& v1, const V3& v2, float& b1, float& b2) {
        const V3 e1 = v1 - v0, e2 = v2 - v0;
        const V3 pv = kcross(d, e2);
        const float det = kdot(e1, pv);
        if (kabs(det) < 1e-30f) { b1 = b2 = 1.0f / 3.0f; return; }   // the ray lies in the plane: the centre (never reached in practice)
        const float inv = 1.0f / det;
        const V3 tv = o - v0;
        const V3 qv = kcross(tv, e1);
        b1 = kmax(kdot(tv, pv) * inv, 0.0f);
        b2 = kmax(kdot(d, qv) * inv, 0.0f);
        const float sum = b1 + b2;
        if (sum > 1.0f) { b1 = b1 / sum; b2 = b2 / sum; }
    }

    /// The light whose cdf interval holds u (binary search over cdf[0..n]).
    RAYT_HD_INLINE uint32_t pickLight(const float* cdf, uint32_t n, float u) {
        uint32_t lo = 0, hi = n;   // the answer i has cdf[i] <= u < cdf[i + 1]
        while (hi - lo > 1) {
            const uint32_t mid = (lo + hi) / 2;
            if (cdf[mid] <= u) lo = mid; else hi = mid;
        }
        return lo;
    }

    RAYT_HD_INLINE void setAov(KAov* aov, const V3& albedo, const V3& normal) {
        aov->albedo[0] = albedo.x; aov->albedo[1] = albedo.y; aov->albedo[2] = albedo.z;
        aov->normal[0] = normal.x; aov->normal[1] = normal.y; aov->normal[2] = normal.z;
        aov->valid = 1;
    }

    /// Where a tracer's hit lies: the point by the kernel's own formula (above), the face normal, the shading normal (facing the ray),
    /// whether the ray came from the outside, the material and the light.
    struct SurfacePoint {
        V3 p, faceN, ns;
        float err;   ///< a bound of p's rounding beyond what offsetRayOrigin covers (a sphere's: its centre and radius set it)
        uint32_t prim;
        float b1, b2;   ///< a triangle's barycentrics of its second and third vertex (for its texture coordinates)
        uint32_t material;
        int32_t light;
        bool front, smooth;
    };

    RAYT_HD_INLINE void surfacePoint(const KScene& s, const KHit& h, const V3& o, const V3& d, SurfacePoint& sp) {
        sp.smooth = false;
        sp.err = 0.0f;
        sp.prim = h.prim;
        sp.b1 = sp.b2 = 0.0f;
        if (h.prim < s.numTriangles) {
            const uint32_t i0 = s.indices[3 * h.prim], i1 = s.indices[3 * h.prim + 1], i2 = s.indices[3 * h.prim + 2];
            const V3 v0 = toV3(s.positions[i0]), v1 = toV3(s.positions[i1]), v2 = toV3(s.positions[i2]);
            float b1, b2;
            triangleBarycentrics(o, d, v0, v1, v2, b1, b2);   // not the tracer's own (h.b1, h.b2): the same bits on every device
            const float w0 = 1.0f - b1 - b2;
            sp.p = w0 * v0 + b1 * v1 + b2 * v2;
            sp.b1 = b1;
            sp.b2 = b2;
            sp.faceN = toV3(s.faceNormals[h.prim]);
            sp.front = kdot(d, sp.faceN) < 0.0f;
            if (s.triSmooth[h.prim]) {
                const V3 n = knormalize(w0 * toV3(s.normals[i0]) + b1 * toV3(s.normals[i1]) + b2 * toV3(s.normals[i2]));
                sp.ns = sp.front ? n : -n;
                sp.smooth = true;
            }
            sp.material = s.triMaterial[h.prim];
            sp.light = s.triLight[h.prim];
        }
        else {
            // A sphere: both tracers report the distance of kernel::intersectSphere; the point is put back on the surface (as pbrt
            // does) so that the offset of the next ray's origin need only cover the rounding of that projection.
            const KSphere& sph = s.spheres[h.prim - s.numTriangles];
            const V3 c = toV3(sph.center);
            const V3 pc = (o + h.t * d) - c;
            sp.p = c + pc * (kabs(sph.radius) / klength(pc));
            // The projection rounds at the scale of the centre and the radius, not of p: a floor made of a sphere of radius 100 centred
            // 100 below has points near y = 0 whose y is off by an ulp of 100 (8e-6), more than the integer offset of a small y moves it.
            // Leaving rays start this much further out (pbrt-v4 offsets by such a bound too); a few ulps of |c| + r.
            sp.err = 8.0f * (kmax(kmax(kabs(c.x), kabs(c.y)), kabs(c.z)) + kabs(sph.radius)) * 5.9604645e-8f;
            sp.faceN = (sp.p - c) / sph.radius;   // outward, inward for a negative radius (Geometry/Sphere.hpp)
            sp.front = kdot(d, sp.faceN) < 0.0f;
            sp.material = sph.material;
            sp.light = sph.light;
        }
        if (!sp.smooth) sp.ns = sp.front ? sp.faceN : -sp.faceN;
    }

    /// The origin of a ray that leaves the surface point toward the side of the unit normal n: beyond the point's rounding, then by the
    /// integer offset (offsetRayOrigin).
    RAYT_HD_INLINE V3 leaveSurface(const SurfacePoint& sp, const V3& n) { return offsetRayOrigin(sp.err > 0.0f ? sp.p + n * sp.err : sp.p, n); }

    /// The distance along the unit ray (o, d) to the hit: a sphere's is the one both tracers compute with kernel::intersectSphere; a
    /// triangle's is taken from the kernel's own hit point, since Embree and OptiX each compute it their own way, which differs in the
    /// last bits (a walk through a medium weighs its channels by it).
    RAYT_HD_INLINE float hitDistance(const KScene& s, const KHit& h, const V3& o, const V3& d, const SurfacePoint& sp) {
        return h.prim < s.numTriangles ? kdot(sp.p - o, d) : h.t;
    }

    /// A diffuse material's reflectance at the surface point: its albedo, times its texture at the point's uv if it has one (Lambertian::albedoAt;
    /// a triangle's uv interpolates its vertices' as TriangleMesh::uvAt, a sphere has none: (0, 0)).
    RAYT_HD_INLINE V3 surfaceAlbedo(const KScene& s, const KMaterial& m, const SurfacePoint& sp) {
        const V3 albedo = materialAlbedo(m);
        if (m.texture < 0) return albedo;
        float u = 0.0f, v = 0.0f;
        if (sp.prim < s.numTriangles) {
            const KFloat2 t0 = s.uvs[s.indices[3 * sp.prim]], t1 = s.uvs[s.indices[3 * sp.prim + 1]], t2 = s.uvs[s.indices[3 * sp.prim + 2]];
            const float w = 1.0f - sp.b1 - sp.b2;
            u = w * t0.u + sp.b1 * t1.u + sp.b2 * t2.u;
            v = w * t0.v + sp.b1 * t1.v + sp.b2 * t2.v;
        }
        const TextureArrays<float> arrays{ s.texels8, s.texelsF, s.srgbLut };
        return albedo * textureRgb(s.textures[m.texture], arrays, u, v);
    }

    /// A smooth dielectric with the same index on both sides (Dielectric::isNullInterface): with --media-nee paths and shadow rays cross it.
    RAYT_HD_INLINE bool isNullInterface(const KMaterial& m) { return m.type == kMatDielectric && m.dielectric.nullInterface != 0; }

    /// One light sample from a point p, as the path tracer chooses it (the environment with probability envSelect, else a light by power):
    /// direction wi, radiance le (times `scale`, the environment's), its density pdf (0: none), and for a light the point on it.
    struct LightChoice {
        V3 wi, le, p;
        float pdf, scale;
        bool area;    ///< a light at a point (a shadow ray stops short of it), not the environment
        bool delta;   ///< a point light: no other strategy can find it, its MIS weight is 1
    };

    RAYT_HD_INLINE void chooseLight(const KScene& s, Pcg32& rng, const V3& p, LightChoice& c) {
        const float uSel = rng.nextFloat();
        const float a = rng.nextFloat(), b = rng.nextFloat();
        c.pdf = 0.0f;
        c.scale = 1.0f;
        c.area = false;
        c.delta = false;
        if (uSel < s.envSelect) {
            float pdfW;
            envSample(s, a, b, c.wi, c.le, pdfW);
            if (pdfW > 0.0f) {
                c.pdf = s.envSelect * pdfW;
                c.scale = s.envScale;
            }
        }
        else {
            const uint32_t li = pickLight(s.lightCdf, s.numLights, (uSel - s.envSelect) / (1.0f - s.envSelect));
            const KLight& q = s.lights[li];
            LightSample<float> ls;
            if (sampleLight(q, p, a, b, ls) && ls.emits && ls.pdf > 0.0f) {
                c.wi = ls.wi;
                c.le = V3(q.radiance[0], q.radiance[1], q.radiance[2]);
                if (q.type == kLightPoint) {   // PointLight::sampleLi: the intensity over the squared distance
                    c.le = c.le * ls.invDist2;
                    c.delta = true;
                }
                c.p = ls.p;
                c.pdf = (1.0f - s.envSelect) * (s.lightCdf[li + 1] - s.lightCdf[li]) * ls.pdf;
                c.area = true;
            }
        }
    }

    /**
     * @brief shadowTransmittance (Render/Transmittance.hpp), for --media-nee: a shadow ray from org along the unit dir up to tMax (kFloatMax:
     *        out of the scene), starting in `medium` (-1: air), crosses index-matched boundaries (switching the medium) and is attenuated in
     *        media; any other surface blocks it.
     */
    template <class Tracer>
    RAYT_HD V3 shadowTransmittance(const KScene& s, const Tracer& tracer, Pcg32& rng, V3 org, const V3& dir, float tMax, int32_t medium) {
        V3 T(1.0f);
        float remaining = tMax;
        for (int crossings = 0; crossings < 64; ++crossings) {
            const KHit h = tracer.closest(org, dir, 0.0f, remaining);
            const bool hit = h.prim != kNoHit;
            SurfacePoint sp;
            if (hit) surfacePoint(s, h, org, dir, sp);
            const float dist = hit ? hitDistance(s, h, org, dir, sp) : remaining;
            if (medium >= 0) {
                if (!hit && !(remaining < kFloatMax)) return V3(0.0f);   // a ray that never leaves the medium: a leaky boundary
                T *= mediumTransmittance(s.media[medium], org, dir, dist, rng);
            }
            if (!hit) return T;
            const KMaterial& m = s.materials[sp.material];
            if (!isNullInterface(m)) return V3(0.0f);
            if (!(kmaxComponent(T) > 0.0f)) return T;
            if (m.medium >= 0) medium = sp.front ? m.medium : -1;   // a front-face hit enters the surface
            org = leaveSurface(sp, sp.front ? -sp.faceN : sp.faceN);   // on, beyond the surface
            if (remaining < kFloatMax) remaining -= dist;
        }
        return V3(0.0f);
    }

    /// The light that a scattering event at p in `medium` sends on along -dirIn, from one light sample (PathIntegrator::estimateDirectInMedium,
    /// --media-nee): the phase function instead of a BSDF, no cosine, the shadow ray through the medium and index-matched boundaries.
    template <class Tracer>
    RAYT_HD V3 lightInMedium(const KScene& s, const Tracer& tracer, Pcg32& rng, const V3& p, const V3& dirIn, int32_t medium) {
        if (!(s.numLights > 0 || hasEnvironment(s))) return V3(0.0f);
        LightChoice c;
        chooseLight(s, rng, p, c);
        if (!(c.pdf > 0.0f)) return V3(0.0f);
        const KMedium& med = s.media[medium];
        const PhaseTable<float> tab{ s.phaseF + med.tableOffset, s.phaseCdf + med.tableOffset };
        const float ph = phaseEval(med, tab, kdot(c.wi, dirIn));
        if (!(ph > 0.0f)) return V3(0.0f);
        V3 dirS = c.wi;
        float tMaxS = kFloatMax;
        if (c.area) {
            const V3 toLight = c.p - p;
            const float dist = klength(toLight);
            if (dist > 0.0f) { dirS = toLight / dist; tMaxS = dist * 0.9999f; }
        }
        const V3 T = shadowTransmittance(s, tracer, rng, p, dirS, tMaxS, medium);   // p lies inside the medium: no offset
        if (!(kmaxComponent(T) > 0.0f)) return V3(0.0f);
        const float w = c.delta ? 1.0f : c.pdf * c.pdf / (c.pdf * c.pdf + ph * ph);
        return c.le * c.scale * T * (ph * w / c.pdf);
    }

    enum : int { kStepSurface = 0, kStepScattered = 1, kStepEnded = 2 };

    /**
     * @brief The path's walk through the medium it travels in, up to the next surface (PathIntegrator::advanceInMedium).
     *
     * Delta tracking with the extinction of the hero channel at the densest point: a tentative collision is real with probability q (the
     * density there), then scatters with the hero's albedo and is absorbed otherwise; a null collision changes nothing. A medium whose
     * channels differ weighs every channel by the density of the same event had it been the hero, normalized by the running ratios rho
     * (the balance heuristic over the hero choice). With --media-nee a scattering event samples a light. Returns kStepSurface (the
     * surface the ray hit is reached; beta weighted), kStepScattered (o, d are the new ray from the collision) or kStepEnded.
     */
    template <class Tracer>
    RAYT_HD int mediumStep(const KScene& s, const Tracer& tracer, Pcg32& rng, int32_t medium, V3& o, V3& d, bool hit, float tSurface, V3& beta, V3& rho,
        int& hero, uint32_t& volumeEvents, V3& L, float& lastPdf, bool& lastDelta) {
        const KMedium& med = s.media[medium];
        const bool constant = med.constant != 0;
        if (!constant && hero < 0) hero = kmin(2, int(rng.nextFloat() * 3.0f));   // chosen once per path, at its first chromatic medium
        const int hh = constant ? 0 : hero;
        const float sT[3] = { med.sigmaT[0], med.sigmaT[1], med.sigmaT[2] };
        const float sS[3] = { med.sigmaS[0], med.sigmaS[1], med.sigmaS[2] };
        // multiplies beta and rho by the weight of one event (kind 0: surface, 1: null, 2: scatter); false if the path cannot go on
        const auto weigh = [&](float t, int kind) -> bool {
            if (constant) return true;
            float e[3], eMax = -INFINITY;
            for (int i = 0; i < 3; ++i) { e[i] = -(sT[i] - sT[hh]) * t; eMax = kmax(eMax, e[i]); }
            float a[3], m = 0.0f;
            for (int i = 0; i < 3; ++i) {
                const float k = kind == 2 ? sS[i] / sS[hh] : kind == 1 ? sT[i] / sT[hh] : 1.0f;
                a[i] = kexp(e[i] - eMax) * k;
                m += rho[i] * a[i];
            }
            m /= 3.0f;
            if (!(m > 0.0f) || !kisFinite(m)) return false;
            for (int i = 0; i < 3; ++i) { beta[i] *= a[i] / m; rho[i] *= a[i] / m; }
            return true;
        };
        float tCur = 0.0f;
        for (;;) {
            const float tNext = tCur - klog(1.0f - rng.nextFloat()) / sT[hh];
            if (tNext >= tSurface) {
                if (!hit) return kStepEnded;   // a ray that escapes a closed boundary: a leaky mesh
                return weigh(tSurface - tCur, 0) ? kStepSurface : kStepEnded;
            }
            if (++volumeEvents > s.volumeCap) return kStepEnded;
            const V3 x = o + d * tNext;
            const float q = mediumDensity(med, x);
            const float u = rng.nextFloat();
            if (u >= q) {   // a null collision
                if (!weigh(tNext - tCur, 1)) return kStepEnded;
                tCur = tNext;
                continue;
            }
            if (u * sT[hh] >= q * sS[hh]) return kStepEnded;   // absorbed (probability 1 - the hero's albedo)
            if (!weigh(tNext - tCur, 2)) return kStepEnded;
            // a real scattering collision
            if (s.mediaNee) L += beta * lightInMedium(s, tracer, rng, x, d, medium);
            const PhaseTable<float> tab{ s.phaseF + med.tableOffset, s.phaseCdf + med.tableOffset };
            const float ua = rng.nextFloat(), ub = rng.nextFloat();
            const V3 dir = phaseSample(med, tab, d, ua, ub);
            lastPdf = phaseEval(med, tab, kdot(d, dir));   // an emitter found next is weighed against the light sample made here
            lastDelta = !s.mediaNee;                       // without it no light was sampled: an emitter found next counts fully
            o = x;                                         // inside the medium: nothing to step over
            d = dir;
            return kStepScattered;
        }
    }

    /**
     * @brief The radiance a camera ray (o, d) brings back (the engine's PathIntegrator::Li), drawing its random numbers from rng.
     *
     * aov (optional, its valid = 0 on entry): receives the denoiser's features of the path, as PathIntegrator::Li records them. Recording
     * them draws no random number, so the radiance is the same bits with or without. The ray starts in air; it travels in the medium of a
     * dielectric (KMaterial::medium) once it crosses that dielectric inward, as PathIntegrator::nextRay decides.
     */
    template <class Tracer>
    RAYT_HD V3 traceRay(const KScene& s, const Tracer& tracer, Pcg32& rng, V3 o, V3 d, KAov* aov = nullptr) {
        V3 L(0.0f), beta(1.0f);
        float lastPdf = 0.0f;
        bool lastDelta = true;   // the camera ray: an emitter it finds counts fully
        int32_t medium = -1;     // the medium the ray travels in (-1: air)
        int hero = -1;           // a chromatic medium's hero channel, chosen when the path meets the first one
        V3 rho(1.0f);
        uint32_t volumeEvents = 0, passThrough = 0;
        for (uint32_t depth = 0; depth < s.maxDepth; ++depth) {
            const KHit h = tracer.closest(o, d, 0.0f, kFloatMax);
            const bool hit = h.prim != kNoHit;
            SurfacePoint sp;
            if (hit) surfacePoint(s, h, o, d, sp);
            if (medium >= 0) {
                const float tSurface = hit ? hitDistance(s, h, o, d, sp) : kFloatMax;
                const int step = mediumStep(s, tracer, rng, medium, o, d, hit, tSurface, beta, rho, hero, volumeEvents, L, lastPdf, lastDelta);
                if (step == kStepEnded) break;
                if (step == kStepScattered) { --depth; continue; }   // a collision in a medium is not a surface bounce
            }
            if (!hit) {
                // ---- the environment, weighted against the light sample of the previous vertex (PathIntegrator::addEscapedRadiance) ----
                if (hasEnvironment(s)) {
                    const V3 le = envEval(s, d);
                    if (aov && !aov->valid) {   // seen directly or through smooth glass: the environment's color (PathIntegrator::addEscapedRadiance)
                        const V3 c = le * s.envScale;
                        setAov(aov, V3(kmin(kmax(c.x, 0.0f), 1.0f), kmin(kmax(c.y, 0.0f), 1.0f), kmin(kmax(c.z, 0.0f), 1.0f)), V3(0.0f));
                    }
                    float w = 1.0f;
                    if (!lastDelta && s.envSelect > 0.0f) {
                        const float pdfL = s.envSelect * envPdf(s, d);
                        w = lastPdf * lastPdf / (lastPdf * lastPdf + pdfL * pdfL);
                    }
                    L += beta * le * s.envScale * w;
                }
                break;
            }
            const V3 p = sp.p, ns = sp.ns;
            const bool front = sp.front;
            const V3 ng = front ? sp.faceN : -sp.faceN;   // faces the incoming ray (SurfaceInteraction::setFaceNormal)
            const KMaterial& m = s.materials[sp.material];

            // ---- with --media-nee an index-matched boundary is crossed, switching the medium; not a bounce ----
            if (s.mediaNee && isNullInterface(m) && ++passThrough < 64) {
                if (m.medium >= 0) medium = front ? m.medium : -1;
                o = leaveSurface(sp, -ng);
                --depth;
                continue;
            }

            // the denoiser's features: the first vertex that is not specular, emitters included (recorded before their light, as the engine does)
            const V3 albedo = surfaceAlbedo(s, m, sp);
            if (aov && !aov->valid && !m.aovSpecular) {
                if (m.texture >= 0)   // Lambertian::aovAlbedo: the textured albedo, clamped to [0, 1]
                    setAov(aov, V3(kmin(kmax(albedo.x, 0.0f), 1.0f), kmin(kmax(albedo.y, 0.0f), 1.0f), kmin(kmax(albedo.z, 0.0f), 1.0f)), ns);
                else setAov(aov, V3(m.aovAlbedo[0], m.aovAlbedo[1], m.aovAlbedo[2]), ns);
            }

            // ---- an emitter: its light, weighted against the light sample of the previous vertex; emitters do not reflect ----
            if (m.type == kMatEmission) {
                if (front || m.twoSided) {
                    float w = 1.0f;
                    const int32_t li = sp.light;
                    if (li >= 0 && !lastDelta) {
                        const float pdfL = (1.0f - s.envSelect) * (s.lightCdf[li + 1] - s.lightCdf[li]) * lightPdf(s.lights[li], o, d, p);
                        w = lastPdf * lastPdf / (lastPdf * lastPdf + pdfL * pdfL);
                    }
                    L += beta * V3(m.emission[0], m.emission[1], m.emission[2]) * w;
                }
                break;
            }

            // ---- a scattering surface: one light sample (PathIntegrator::estimateDirect; none at a delta material, whose eval is 0),
            //      then the BSDF sample ----
            const V3 wo = -d;
            if ((s.numLights > 0 || hasEnvironment(s)) && !bsdfIsDelta(m)) {
                // The light is chosen first and the BSDF is evaluated in one place: a layered material's eval is a long random walk,
                // and on the GPU every call site is a copy of it in the one program that runs a whole path.
                LightChoice c;
                chooseLight(s, rng, p, c);
                if (c.pdf > 0.0f) {
                    const V3 wi = c.wi;
                    const V3 f = bsdfEval(m, albedo, front, wo, wi, ns, ng, rng);
                    if (kmaxComponent(f) > 0.0f) {
                        const V3 org = leaveSurface(sp, kdot(wi, ng) > 0.0f ? ng : -ng);   // on the side the light is on (rough glass transmits)
                        V3 dirS = wi;
                        float tMaxS = kFloatMax;
                        if (c.area) {
                            const V3 toL2 = c.p - org;
                            const float distB = klength(toL2);
                            dirS = toL2 / distB;
                            tMaxS = distB * 0.9999f;
                        }
                        V3 T(1.0f);
                        bool blocked;
                        if (s.mediaNee) {
                            // index-matched boundaries are crossed and media attenuate; the shadow ray starts in the medium on the side it
                            // leaves to (as the path would)
                            int32_t startMedium = medium;
                            if (m.medium >= 0) startMedium = ((kdot(wi, ng) > 0.0f) == front) ? -1 : m.medium;
                            T = shadowTransmittance(s, tracer, rng, org, dirS, tMaxS, startMedium);
                            blocked = !(kmaxComponent(T) > 0.0f);
                        }
                        else blocked = tracer.occluded(org, dirS, 0.0f, tMaxS);
                        if (!blocked) {
                            const float pdfL = c.pdf;
                            const float pdfB = bsdfPdf(m, front, wo, wi, ns, ng);
                            const float w = (pdfB > 0.0f && !c.delta) ? pdfL * pdfL / (pdfL * pdfL + pdfB * pdfB) : 1.0f;
                            const float cosT = kabs(kdot(ns, wi));
                            if (s.mediaNee) L += beta * f * c.le * c.scale * T * (cosT * w / pdfL);
                            else L += beta * f * c.le * c.scale * (cosT * w / pdfL);   // scale: the environment's (times 1 is exact for a light)
                        }
                    }
                }
            }

            const float u1 = rng.nextFloat(), u2 = rng.nextFloat(), u3 = rng.nextFloat();
            BsdfSample<float> bs;
            if (!bsdfSample(m, albedo, front, wo, ns, ng, u1, u2, u3, rng, bs)) break;   // no sample (below the geometric surface, total reflection, ...)
            if (bs.delta) beta *= bs.f;   // f is already the weight of the specular event (PathIntegrator::updateThroughput)
            else {
                if (!(bs.pdf > 1e-8f)) break;
                beta *= bs.f * kabs(kdot(ns, bs.wi)) / bs.pdf;
            }
            if (kmaxComponent(beta) <= 0.0f) break;
            if (depth >= 3) {
                const float q = kmin(kmax(kmaxComponent(beta), 0.05f), 0.95f);
                if (rng.nextFloat() > q) break;
                beta /= q;
            }
            lastPdf = bs.pdf;
            lastDelta = bs.delta;
            const V3 wi = bs.wi;
            // a boundary with a medium inside: the new direction goes to the inside or the outside (PathIntegrator::nextRay)
            if (m.medium >= 0) medium = ((kdot(wi, ng) > 0.0f) == front) ? -1 : m.medium;
            o = leaveSurface(sp, kdot(wi, ng) > 0.0f ? ng : -ng);   // a refracted ray leaves from the other side
            d = wi;
        }
        return L;
    }

    /// The camera ray of pixel (px, py) (py = 0 is the bottom row) through a uniform point of the pixel, drawn from rng, and its weight
    /// (PathIntegrator::render's camW): the pinhole draws nothing more; a thin or a real lens draws a point of the lens (Kernel/Camera.hpp).
    RAYT_HD_INLINE void cameraRay(const KScene& s, Pcg32& rng, uint32_t px, uint32_t py, V3& o, V3& d, float& weight) {
        const float u = (float(px) + rng.nextFloat()) / float(s.width);
        const float v = (float(py) + rng.nextFloat()) / float(s.height);
        if (s.camera.type == kCamPinhole) {
            o = lens::v3(s.camera.origin);
            d = knormalize(lens::v3(s.camera.lowerLeft) + u * lens::v3(s.camera.horizontal) + v * lens::v3(s.camera.vertical) - o);
            weight = 1.0f;
            return;
        }
        const float a = rng.nextFloat(), b = rng.nextFloat();
        lensCameraRay(s.camera, LensArrays<float>{ s.lensSurfaces, s.lensPupil, s.lensBladeNormals }, u, v, a, b, o, d, weight);
    }

    /// The radiance of one camera path through pixel (px, py) (py = 0 is the bottom row), sample index `sample`, times the camera's
    /// weight of the sample; aov as for traceRay. A sample the lens blocks is traced only for the features (its central ray).
    template <class Tracer>
    RAYT_HD V3 tracePath(const KScene& s, const Tracer& tracer, uint32_t px, uint32_t py, uint32_t sample, KAov* aov = nullptr) {
        Pcg32 rng = Pcg32::forSample(s.seed, py * s.width + px, sample);
        V3 o, d;
        float weight;
        cameraRay(s, rng, px, py, o, d, weight);
        if (weight <= 0.0f && !aov) return V3(0.0f);   // nothing to trace: the sample adds no light (PathIntegrator::render)
        return traceRay(s, tracer, rng, o, d, aov) * weight;
    }

} // namespace rayt::kernel
