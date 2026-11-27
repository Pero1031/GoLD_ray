#pragma once
/**
 * @file Kernel/PathTracer.hpp
 * @brief The kernel path tracer: one camera path of one pixel sample, in float, on the CPU or the GPU (docs/GPU.md).
 *
 * The rules are the engine's PathIntegrator restricted to what KScene holds (the materials of Kernel/Materials.hpp, emitting
 * quads, a uniform environment), so the two converge to the same picture: a path has at most maxDepth surface vertices; each
 * non-emitting vertex samples one light (chosen by power, a uniform point on the quad) with a shadow ray, weighted against the
 * BSDF's density by the power heuristic, and then samples the BSDF; an emitter found by a BSDF ray is weighted against that
 * light sample the same way; Russian roulette from the fourth vertex on with survival max(beta) clamped to [0.05, 0.95]. The
 * environment is only found by BSDF rays (the engine also samples it; both are unbiased).
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

#include "Kernel/KernelScene.hpp"
#include "Kernel/Materials.hpp"
#include "Kernel/Rng.hpp"
#include "Kernel/Sampling.hpp"

namespace rayt::kernel {

    using V3 = Vec3<float>;

    constexpr float kFloatMax = 3.402823466e38f;
    constexpr float kInvPi = 0.318309886183790671538f;

    RAYT_HD_INLINE V3 toV3(const KFloat3& a) { return V3(a.x, a.y, a.z); }

    RAYT_HD_INLINE int32_t kFloatAsInt(float f) {
#if defined(RAYT_GPU_CODE)
        return __float_as_int(f);
#else
        int32_t i;
        memcpy(&i, &f, sizeof(i));
        return i;
#endif
    }
    RAYT_HD_INLINE float kIntAsFloat(int32_t i) {
#if defined(RAYT_GPU_CODE)
        return __int_as_float(i);
#else
        float f;
        memcpy(&f, &i, sizeof(f));
        return f;
#endif
    }

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

    /// The radiance of one camera path through pixel (px, py) (py = 0 is the bottom row), sample index `sample`.
    template <class Tracer>
    RAYT_HD V3 tracePath(const KScene& s, const Tracer& tracer, uint32_t px, uint32_t py, uint32_t sample) {
        Pcg32 rng = Pcg32::forSample(s.seed, py * s.width + px, sample);
        const float u = (float(px) + rng.nextFloat()) / float(s.width);
        const float v = (float(py) + rng.nextFloat()) / float(s.height);
        V3 o = toV3(s.camera.origin);
        V3 d = knormalize(toV3(s.camera.lowerLeft) + u * toV3(s.camera.horizontal) + v * toV3(s.camera.vertical) - o);

        V3 L(0.0f), beta(1.0f);
        float lastPdf = 0.0f;
        bool lastDelta = true;   // the camera ray: an emitter it finds counts fully
        for (uint32_t depth = 0; depth < s.maxDepth; ++depth) {
            const KHit h = tracer.closest(o, d, 0.0f, kFloatMax);
            if (h.prim == kNoHit) {
                L += beta * V3(s.envRadiance[0], s.envRadiance[1], s.envRadiance[2]);
                break;
            }
            const uint32_t i0 = s.indices[3 * h.prim], i1 = s.indices[3 * h.prim + 1], i2 = s.indices[3 * h.prim + 2];
            const V3 v0 = toV3(s.positions[i0]), v1 = toV3(s.positions[i1]), v2 = toV3(s.positions[i2]);
            float b1, b2;
            triangleBarycentrics(o, d, v0, v1, v2, b1, b2);   // not the tracer's own (h.b1, h.b2): the same bits on every device
            const float w0 = 1.0f - b1 - b2;
            const V3 p = w0 * v0 + b1 * v1 + b2 * v2;
            const V3 faceN = toV3(s.faceNormals[h.prim]);
            const bool front = kdot(d, faceN) < 0.0f;
            const V3 ng = front ? faceN : -faceN;   // faces the incoming ray (SurfaceInteraction::setFaceNormal)
            V3 ns = ng;
            if (s.triSmooth[h.prim]) {
                const V3 n = knormalize(w0 * toV3(s.normals[i0]) + b1 * toV3(s.normals[i1]) + b2 * toV3(s.normals[i2]));
                ns = front ? n : -n;
            }
            const KMaterial& m = s.materials[s.triMaterial[h.prim]];

            // ---- an emitter: its light, weighted against the light sample of the previous vertex; emitters do not reflect ----
            if (m.type == kMatEmission) {
                if (front || m.twoSided) {
                    float w = 1.0f;
                    const int32_t li = s.triLight[h.prim];
                    if (li >= 0 && !lastDelta) {
                        const KQuadLight& q = s.lights[li];
                        const V3 dp = p - o;
                        const float cosL = kabs(kdot(d, toV3(q.normal)));
                        const float pdfL = (s.lightCdf[li + 1] - s.lightCdf[li]) * kdot(dp, dp) / (cosL * q.area);
                        w = lastPdf * lastPdf / (lastPdf * lastPdf + pdfL * pdfL);
                    }
                    L += beta * V3(m.emission[0], m.emission[1], m.emission[2]) * w;
                }
                break;
            }

            // ---- a scattering surface: one light sample (PathIntegrator::estimateDirect; none at a delta material, whose eval is 0),
            //      then the BSDF sample ----
            const V3 wo = -d;
            if (s.numLights > 0 && !bsdfIsDelta(m)) {
                const uint32_t li = pickLight(s.lightCdf, s.numLights, rng.nextFloat());
                const KQuadLight& q = s.lights[li];
                const float a = rng.nextFloat(), b = rng.nextFloat();
                const V3 pl = toV3(q.corner) + a * toV3(q.edge1) + b * toV3(q.edge2);
                const V3 toL = pl - p;
                const float dist2 = kdot(toL, toL);
                const float dist = ksqrt(dist2);
                const V3 wi = toL / dist;
                float cosL = -kdot(wi, toV3(q.normal));
                if (q.twoSided) cosL = kabs(cosL);
                if (cosL > 0.0f) {
                    const V3 f = bsdfEval(m, front, wo, wi, ns, ng);
                    if (kmaxComponent(f) > 0.0f) {
                        const V3 org = offsetRayOrigin(p, kdot(wi, ng) > 0.0f ? ng : -ng);   // on the side the light is on (rough glass transmits)
                        const V3 toL2 = pl - org;
                        const float distB = klength(toL2);
                        if (!tracer.occluded(org, toL2 / distB, 0.0f, distB * 0.9999f)) {
                            const float pdfL = (s.lightCdf[li + 1] - s.lightCdf[li]) * dist2 / (cosL * q.area);
                            const float pdfB = bsdfPdf(m, front, wo, wi, ns, ng);
                            const float w = pdfB > 0.0f ? pdfL * pdfL / (pdfL * pdfL + pdfB * pdfB) : 1.0f;
                            const float cosT = kabs(kdot(ns, wi));
                            L += beta * f * V3(q.radiance[0], q.radiance[1], q.radiance[2]) * (cosT * w / pdfL);
                        }
                    }
                }
            }

            const float u1 = rng.nextFloat(), u2 = rng.nextFloat(), u3 = rng.nextFloat();
            BsdfSample<float> bs;
            if (!bsdfSample(m, front, wo, ns, ng, u1, u2, u3, bs)) break;   // no sample (below the geometric surface, total reflection, ...)
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
            o = offsetRayOrigin(p, kdot(wi, ng) > 0.0f ? ng : -ng);   // a refracted ray leaves from the other side
            d = wi;
        }
        return L;
    }

} // namespace rayt::kernel
