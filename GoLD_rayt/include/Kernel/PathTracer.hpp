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
#include "Kernel/Environment.hpp"
#include "Kernel/KernelScene.hpp"
#include "Kernel/Lights.hpp"
#include "Kernel/Rng.hpp"
#include "Kernel/Sampling.hpp"

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

    /**
     * @brief The radiance a camera ray (o, d) brings back (the engine's PathIntegrator::Li), drawing its random numbers from rng.
     *
     * aov (optional, its valid = 0 on entry): receives the denoiser's features of the path, as PathIntegrator::Li records them. Recording
     * them draws no random number, so the radiance is the same bits with or without.
     */
    template <class Tracer>
    RAYT_HD V3 traceRay(const KScene& s, const Tracer& tracer, Pcg32& rng, V3 o, V3 d, KAov* aov = nullptr) {
        V3 L(0.0f), beta(1.0f);
        float lastPdf = 0.0f;
        bool lastDelta = true;   // the camera ray: an emitter it finds counts fully
        for (uint32_t depth = 0; depth < s.maxDepth; ++depth) {
            const KHit h = tracer.closest(o, d, 0.0f, kFloatMax);
            if (h.prim == kNoHit) {
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
            V3 p, faceN, ns;
            uint32_t material;
            int32_t lightIndex;
            bool front, smooth = false;
            if (h.prim < s.numTriangles) {
                const uint32_t i0 = s.indices[3 * h.prim], i1 = s.indices[3 * h.prim + 1], i2 = s.indices[3 * h.prim + 2];
                const V3 v0 = toV3(s.positions[i0]), v1 = toV3(s.positions[i1]), v2 = toV3(s.positions[i2]);
                float b1, b2;
                triangleBarycentrics(o, d, v0, v1, v2, b1, b2);   // not the tracer's own (h.b1, h.b2): the same bits on every device
                const float w0 = 1.0f - b1 - b2;
                p = w0 * v0 + b1 * v1 + b2 * v2;
                faceN = toV3(s.faceNormals[h.prim]);
                front = kdot(d, faceN) < 0.0f;
                if (s.triSmooth[h.prim]) {
                    const V3 n = knormalize(w0 * toV3(s.normals[i0]) + b1 * toV3(s.normals[i1]) + b2 * toV3(s.normals[i2]));
                    ns = front ? n : -n;
                    smooth = true;
                }
                material = s.triMaterial[h.prim];
                lightIndex = s.triLight[h.prim];
            }
            else {
                // A sphere: both tracers report the distance of kernel::intersectSphere; the point is put back on the surface (as pbrt
                // does) so that the offset of the next ray's origin need only cover the rounding of that projection.
                const KSphere& sp = s.spheres[h.prim - s.numTriangles];
                const V3 c = toV3(sp.center);
                const V3 pc = (o + h.t * d) - c;
                p = c + pc * (kabs(sp.radius) / klength(pc));
                faceN = (p - c) / sp.radius;   // outward, inward for a negative radius (Geometry/Sphere.hpp)
                front = kdot(d, faceN) < 0.0f;
                material = sp.material;
                lightIndex = sp.light;
            }
            const V3 ng = front ? faceN : -faceN;   // faces the incoming ray (SurfaceInteraction::setFaceNormal)
            if (!smooth) ns = ng;
            const KMaterial& m = s.materials[material];
            // the denoiser's features: the first vertex that is not specular, emitters included (recorded before their light, as the engine does)
            if (aov && !aov->valid && !m.aovSpecular) setAov(aov, V3(m.aovAlbedo[0], m.aovAlbedo[1], m.aovAlbedo[2]), ns);

            // ---- an emitter: its light, weighted against the light sample of the previous vertex; emitters do not reflect ----
            if (m.type == kMatEmission) {
                if (front || m.twoSided) {
                    float w = 1.0f;
                    const int32_t li = lightIndex;
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
                // The light: the environment with probability envSelect, else a light of the scene (a quad or a sphere) chosen by power.
                // The light is chosen first and the BSDF is evaluated in one place: a layered material's eval is a long random walk,
                // and on the GPU every call site is a copy of it in the one program that runs a whole path.
                const float uSel = rng.nextFloat();
                const float a = rng.nextFloat(), b = rng.nextFloat();
                V3 wi, le, lightP;
                float pdfL = 0.0f, scale = 1.0f;
                bool area = false;
                if (uSel < s.envSelect) {
                    float pdfW;
                    envSample(s, a, b, wi, le, pdfW);
                    if (pdfW > 0.0f) {
                        pdfL = s.envSelect * pdfW;
                        scale = s.envScale;
                    }
                }
                else {
                    const uint32_t li = pickLight(s.lightCdf, s.numLights, (uSel - s.envSelect) / (1.0f - s.envSelect));
                    const KLight& q = s.lights[li];
                    LightSample<float> ls;
                    if (sampleLight(q, p, a, b, ls) && ls.emits && ls.pdf > 0.0f) {
                        wi = ls.wi;
                        le = V3(q.radiance[0], q.radiance[1], q.radiance[2]);
                        lightP = ls.p;
                        pdfL = (1.0f - s.envSelect) * (s.lightCdf[li + 1] - s.lightCdf[li]) * ls.pdf;
                        area = true;
                    }
                }
                if (pdfL > 0.0f) {
                    const V3 f = bsdfEval(m, front, wo, wi, ns, ng, rng);
                    if (kmaxComponent(f) > 0.0f) {
                        const V3 org = offsetRayOrigin(p, kdot(wi, ng) > 0.0f ? ng : -ng);   // on the side the light is on (rough glass transmits)
                        bool blocked;
                        if (area) {
                            const V3 toL2 = lightP - org;
                            const float distB = klength(toL2);
                            blocked = tracer.occluded(org, toL2 / distB, 0.0f, distB * 0.9999f);
                        }
                        else blocked = tracer.occluded(org, wi, 0.0f, kFloatMax);
                        if (!blocked) {
                            const float pdfB = bsdfPdf(m, front, wo, wi, ns, ng);
                            const float w = pdfB > 0.0f ? pdfL * pdfL / (pdfL * pdfL + pdfB * pdfB) : 1.0f;
                            const float cosT = kabs(kdot(ns, wi));
                            L += beta * f * le * scale * (cosT * w / pdfL);   // scale: the environment's (times 1 is exact for a light)
                        }
                    }
                }
            }

            const float u1 = rng.nextFloat(), u2 = rng.nextFloat(), u3 = rng.nextFloat();
            BsdfSample<float> bs;
            if (!bsdfSample(m, front, wo, ns, ng, u1, u2, u3, rng, bs)) break;   // no sample (below the geometric surface, total reflection, ...)
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

    /// The camera ray of pixel (px, py) (py = 0 is the bottom row): the pinhole through a uniform point of the pixel, drawn from rng.
    RAYT_HD_INLINE void cameraRay(const KScene& s, Pcg32& rng, uint32_t px, uint32_t py, V3& o, V3& d) {
        const float u = (float(px) + rng.nextFloat()) / float(s.width);
        const float v = (float(py) + rng.nextFloat()) / float(s.height);
        o = toV3(s.camera.origin);
        d = knormalize(toV3(s.camera.lowerLeft) + u * toV3(s.camera.horizontal) + v * toV3(s.camera.vertical) - o);
    }

    /// The radiance of one camera path through pixel (px, py) (py = 0 is the bottom row), sample index `sample`; aov as for traceRay.
    template <class Tracer>
    RAYT_HD V3 tracePath(const KScene& s, const Tracer& tracer, uint32_t px, uint32_t py, uint32_t sample, KAov* aov = nullptr) {
        Pcg32 rng = Pcg32::forSample(s.seed, py * s.width + px, sample);
        V3 o, d;
        cameraRay(s, rng, px, py, o, d);
        return traceRay(s, tracer, rng, o, d, aov);
    }

} // namespace rayt::kernel
