#pragma once
/**
 * @file Kernel/Lights.hpp
 * @brief Spheres and the lights of the kernel renderer (docs/GPU.md).
 *
 * intersectSphere is the one ray-sphere test of the kernel: Embree calls it as user geometry and OptiX as the intersection
 * program of a custom primitive (as Mitsuba binds one routine to both), so both tracers report the same distance for a sphere.
 * The light functions mirror the engine: a quad light is sampled at a uniform point (as before), a sphere light in the cone it
 * subtends from outside and at a uniform point of its surface from inside (Geometry/Sphere.hpp sampleFrom / pdfSurface and the
 * default Shape::sampleFrom), with the densities over solid angle at the reference point.
 */

#include "Kernel/KernelScene.hpp"
#include "Kernel/KernelMath.hpp"

namespace rayt::kernel {

    /**
     * @brief The ray o + t d against a sphere: the nearest t in (tmin, tmax), or false.
     *
     * Haines et al., "Precision Improvements for Ray/Sphere Intersection" (Ray Tracing Gems, ch. 7): the discriminant from the
     * distance of the centre to the line, and the smaller root as c / q, so neither root loses its digits to cancellation.
     */
    RAYT_HD_INLINE bool intersectSphere(const Vec3<float>& o, const Vec3<float>& d, const Vec3<float>& c, float radius, float tmin, float tmax, float& t) {
        const Vec3<float> f = o - c;
        const float a = kdot(d, d);
        if (!(a > 0.0f)) return false;
        const float b = -kdot(f, d);
        const Vec3<float> l = f + (b / a) * d;           // from the centre to the point of the line closest to it
        const float r2 = radius * radius;
        const float disc = a * (r2 - kdot(l, l));
        if (disc < 0.0f) return false;
        const float q = b + kcopysign(ksqrt(disc), b);
        if (q == 0.0f) return false;
        const float cc = kdot(f, f) - r2;
        float t0 = cc / q, t1 = q / a;
        if (t1 < t0) { const float x = t0; t0 = t1; t1 = x; }
        if (t0 > tmin && t0 < tmax) { t = t0; return true; }
        if (t1 > tmin && t1 < tmax) { t = t1; return true; }
        return false;
    }

    /**
     * @brief The box both tracers put around a sphere: the exact box widened by a millionth of its size, so that the float
     *        rounding of a box test never culls a ray that intersectSphere would accept (the test itself decides every hit).
     */
    RAYT_HD_INLINE void sphereBounds(const KSphere& s, float lo[3], float hi[3]) {
        const float r = kabs(s.radius);
        const float c[3] = { s.center.x, s.center.y, s.center.z };
        const float pad = 1e-6f * (kmax(kmax(kabs(c[0]), kabs(c[1])), kabs(c[2])) + r) + 1e-30f;
        for (int i = 0; i < 3; ++i) {
            lo[i] = c[i] - r - pad;
            hi[i] = c[i] + r + pad;
        }
    }

    /// A point of a light seen from a reference point: where it is, the direction to it, and the density over solid angle (0: none).
    template <class T> struct LightSample {
        Vec3<T> p;
        Vec3<T> wi;
        T pdf;
        bool emits;   ///< the point shows its emitting side to the reference point
        T invDist2;   ///< a point light: 1 / its squared distance, the falloff of its intensity (other lights: unused)
    };

    /// PointLight::sampleLi: the direction toward the light, a delta (pdf 1 by the engine's convention), and the falloff 1 / d^2.
    template <class T> RAYT_HD_INLINE bool samplePointLight(const Vec3<T>& pos, const Vec3<T>& ref, LightSample<T>& out) {
        const Vec3<T> d = pos - ref;
        const T dist2 = kdot(d, d);
        if (dist2 <= T(0)) return false;
        const T invDist = ksafeRecip(ksafeSqrt(dist2));
        out.wi = d * invDist;
        out.p = pos;
        out.pdf = T(1);
        out.emits = true;
        out.invDist2 = ksafeRecip(dist2);
        return true;
    }

    /// Sphere::sampleFrom (with the default Shape::sampleFrom inside the sphere) and AreaLight::sampleLi's direction.
    template <class T> RAYT_HD_INLINE bool sampleSphereLight(const Vec3<T>& c, T radius, bool twoSided, const Vec3<T>& ref, T u1, T u2, LightSample<T>& out) {
        const Vec3<T> toCenter = c - ref;
        const T dc2 = kdot(toCenter, toCenter);
        const T r2 = radius * radius;
        if (dc2 <= r2 * (T(1) + T(1e-6))) {
            // inside: a uniform point of the surface, its area density turned into solid angle at ref
            const T z = T(1) - T(2) * u1;
            const T rr = ksafeSqrt(kmax(T(0), T(1) - z * z));
            const T phi = T(2) * kPi<T>() * u2;
            const Vec3<T> dir(rr * kcos(phi), rr * ksin(phi), z);
            out.p = c + dir * radius;
            const T area = T(4) * kPi<T>() * radius * radius;
            const T pdfA = ksafeRecip(area);
            const Vec3<T> dv = out.p - ref;
            const T dist2 = kdot(dv, dv);
            if (dist2 <= T(0)) return false;
            const T cosL = kdot(dir, -(dv / ksqrt(dist2)));
            if (kabs(cosL) <= T(1e-9)) return false;
            out.emits = cosL > T(0) || twoSided;
            out.pdf = pdfA * dist2 / kabs(cosL);
            out.wi = dv * ksafeRecip(ksafeSqrt(dist2));   // AreaLight::sampleLi
            return true;
        }
        const T dc = ksqrt(dc2);
        const Vec3<T> wc = toCenter / dc;
        const T sin2Max = r2 / dc2;
        const T cosMax = ksafeSqrt(T(1) - sin2Max);
        const T oneMinusCosMax = sin2Max / (T(1) + cosMax);   // stable for small spheres
        const T cosT = T(1) - u1 * oneMinusCosMax;
        const T sin2T = u1 * oneMinusCosMax * (T(1) + cosT);
        const T sinMax = ksqrt(sin2Max);
        const T cosA = sin2T / sinMax + cosT * ksafeSqrt(T(1) - sin2T / sin2Max);   // at the centre, between wc and the point (law of sines)
        const T sinA = ksafeSqrt(T(1) - cosA * cosA);
        const T phi = T(2) * kPi<T>() * u2;
        const T sign = kcopysign(T(1), wc.z);
        const T a = -T(1) / (sign + wc.z);
        const T b = wc.x * wc.y * a;
        const Vec3<T> bx(T(1) + sign * wc.x * wc.x * a, sign * b, -sign * wc.x);
        const Vec3<T> by(b, sign + wc.y * wc.y * a, -wc.y);
        const Vec3<T> n = -(bx * (sinA * kcos(phi)) + by * (sinA * ksin(phi)) + wc * cosA);
        out.p = c + n * radius;
        const Vec3<T> dv = out.p - ref;
        const T dist2 = kdot(dv, dv);
        if (dist2 <= T(0)) return false;
        out.wi = dv * ksafeRecip(ksafeSqrt(dist2));   // AreaLight::sampleLi
        out.emits = true;                              // the cone's points face ref
        out.pdf = T(1) / (T(2) * kPi<T>() * oneMinusCosMax);
        return true;
    }

    /// Sphere::pdfSurface: the density of sampleSphereLight in direction wi from ref, where wi meets the sphere at p.
    template <class T> RAYT_HD_INLINE T sphereLightPdf(const Vec3<T>& c, T radius, const Vec3<T>& ref, const Vec3<T>& wi, const Vec3<T>& p) {
        const Vec3<T> toCenter = c - ref;
        const T dc2 = kdot(toCenter, toCenter);
        const T r2 = radius * radius;
        if (dc2 > r2 * (T(1) + T(1e-6))) {
            const T sin2Max = r2 / dc2;
            const T cosMax = ksafeSqrt(T(1) - sin2Max);
            const T cosT = kdot(wi, toCenter) / ksqrt(dc2);
            if (cosT < cosMax) return T(0);
            return (T(1) + cosMax) / (T(2) * kPi<T>() * sin2Max);
        }
        const Vec3<T> dv = p - ref;
        const T dist2 = kdot(dv, dv);
        const T cosTheta = kabs(kdot((p - c) / radius, -wi));
        if (cosTheta <= T(0)) return T(0);
        return dist2 / (cosTheta * (T(4) * kPi<T>() * radius * radius));
    }

    /// A light sample of light q from ref (float, what the renderer calls).
    RAYT_HD_INLINE bool sampleLight(const KLight& q, const Vec3<float>& ref, float a, float b, LightSample<float>& out) {
        if (q.type == kLightPoint) return samplePointLight(Vec3<float>(q.center.x, q.center.y, q.center.z), ref, out);
        if (q.type == kLightSphere)
            return sampleSphereLight(Vec3<float>(q.center.x, q.center.y, q.center.z), q.radius, q.twoSided != 0, ref, a, b, out);
        const Vec3<float> corner(q.corner.x, q.corner.y, q.corner.z), e1(q.edge1.x, q.edge1.y, q.edge1.z), e2(q.edge2.x, q.edge2.y, q.edge2.z);
        const Vec3<float> n(q.normal.x, q.normal.y, q.normal.z);
        out.p = corner + a * e1 + b * e2;
        const Vec3<float> toL = out.p - ref;
        const float dist2 = kdot(toL, toL);
        out.wi = toL / ksqrt(dist2);
        float cosL = -kdot(out.wi, n);
        if (q.twoSided) cosL = kabs(cosL);
        out.emits = cosL > 0.0f;
        out.pdf = out.emits ? dist2 / (cosL * q.area) : 0.0f;
        return out.emits;
    }

    /// The density of sampleLight for light q in direction d from ref, where d meets the light at p (float).
    RAYT_HD_INLINE float lightPdf(const KLight& q, const Vec3<float>& ref, const Vec3<float>& d, const Vec3<float>& p) {
        if (q.type == kLightPoint) return 0.0f;   // a delta: no path finds it
        if (q.type == kLightSphere) return sphereLightPdf(Vec3<float>(q.center.x, q.center.y, q.center.z), q.radius, ref, d, p);
        const Vec3<float> dp = p - ref;
        const float cosL = kabs(kdot(d, Vec3<float>(q.normal.x, q.normal.y, q.normal.z)));
        return kdot(dp, dp) / (cosL * q.area);
    }

} // namespace rayt::kernel
