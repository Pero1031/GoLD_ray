#pragma once
/**
 * @file Kernel/Camera.hpp
 * @brief The kernel renderer's camera rays: a thin lens and a real lens traced surface by surface (docs/GPU.md).
 *
 * Render/Camera.hpp (generateRay, lensRay, lensToWorld, sampleBlades) and Render/LensSystem.cpp (generateRay, centralRay,
 * traceFromFilm and what it calls) step by step, for what the kernel models: a thin lens with a round or polygonal aperture, and a
 * real lens of spherical and flat surfaces and a stop (no aspheric surfaces, no dispersion: the RGB build traces n_d). With T = double
 * the rays are the engine's bits; with T = float they are what the GPU computes. (The pinhole stays in Kernel/PathTracer.hpp.)
 */

#include "Kernel/KernelMath.hpp"
#include "Kernel/KernelScene.hpp"

namespace rayt::kernel {

    /// A real lens's arrays (KScene::lensSurfaces, lensPupil, lensBladeNormals), in the precision of the evaluation.
    template <class T> struct LensArrays {
        const LensSurface<T>* surfaces;
        const T* pupil;
        const T* bladeNormals;
    };

    namespace lens {

        template <class T> RAYT_HD_INLINE Vec3<T> v3(const T* a) { return Vec3<T>(a[0], a[1], a[2]); }

        /// Refraction of the unit direction d at a surface whose unit normal n faces against d; eta = n_from / n_to (LensSystem.cpp's Refract).
        template <class T> RAYT_HD bool refract(Vec3<T>& d, const Vec3<T>& n, T eta) {
            const T cosi = -kdot(n, d);
            const T sin2t = eta * eta * (T(1) - cosi * cosi);
            if (sin2t >= T(1)) return false;
            d = knormalizeGlm(eta * d + (eta * cosi - ksqrt(T(1) - sin2t)) * n);
            return true;
        }

        /// LensSystem::insideStop for the lens's own stop (its precomputed edge normals when it has blades).
        template <class T> RAYT_HD bool insideStop(const CameraParams<T>& c, const LensArrays<T>& L, const LensSurface<T>& s, T x, T y) {
            const T r = s.apertureRadius * c.stopScale;
            if (c.stopBlades < 3) return x * x + y * y <= r * r;
            const T apothem = r * c.apothemPerRadius;
            for (uint32_t k = 0; k < c.stopBlades; ++k)
                if (x * L.bladeNormals[2 * k] + y * L.bladeNormals[2 * k + 1] > apothem) return false;
            return true;
        }

        /// LensSystem::hitSurface for a sphere or a flat surface, with the clear aperture.
        template <class T> RAYT_HD bool hitSurface(const LensSurface<T>& s, Vec3<T>& o, const Vec3<T>& d, bool backward, Vec3<T>& normal) {
            const T tMin = backward ? -kinfinity<T>() : T(1e-12);
            T t = T(0);
            if (s.kind == kLensFlat) {
                if (kabs(d.z) < T(1e-12)) return false;
                t = (s.z - o.z) / d.z;
                if (t <= tMin) return false;
                const Vec3<T> p = o + d * t;
                const T rho = ksqrt(p.x * p.x + p.y * p.y);
                if (rho > s.apertureRadius) return false;
                // the sag of a flat surface is 0 with slope 0 (Surface::sag), so the normal is the axis
                const T slope = T(0);
                normal = rho > T(0) ? knormalizeGlm(Vec3<T>(slope * p.x / rho, slope * p.y / rho, T(1))) : Vec3<T>(0, 0, 1);
                o = p;
                return true;
            }
            // a sphere through the vertex, its centre on the axis; the intersection on the vertex's side of the sphere
            const Vec3<T> centre(T(0), T(0), s.z - s.radius);
            const Vec3<T> oc = o - centre;
            const T b = kdot(oc, d), c = kdot(oc, oc) - s.radius * s.radius;
            const T disc = b * b - c;
            if (disc < T(0)) return false;
            const T root = ksqrt(disc);
            t = ((d.z > T(0)) == (s.radius > T(0))) ? -b + root : -b - root;
            if (t <= tMin) return false;
            const Vec3<T> p = o + d * t;
            const T rho = ksqrt(p.x * p.x + p.y * p.y);
            if (rho > s.apertureRadius) return false;
            normal = knormalizeGlm(p - centre);
            o = p;
            return true;
        }

        /// LensSystem::crossSurface: the stop (a plane, no refraction) or a refracting surface.
        template <class T> RAYT_HD bool crossSurface(const CameraParams<T>& c, const LensArrays<T>& L, const LensSurface<T>& s, Vec3<T>& o, Vec3<T>& d,
            T nFrom, T nTo, bool backward) {
            if (s.kind == kLensStop) {
                if (kabs(d.z) < T(1e-12)) return false;
                const T t = (s.z - o.z) / d.z;
                if (!backward && t <= T(0)) return false;
                const Vec3<T> p = o + d * t;
                if (!insideStop(c, L, s, p.x, p.y)) return false;
                o = p;
                return true;
            }
            Vec3<T> normal;
            if (!hitSurface(s, o, d, backward, normal)) return false;
            if (kdot(normal, d) > T(0)) normal = -normal;
            return refract(d, normal, nFrom / nTo);
        }

        /// LensSystem::traceFromFilm: through every surface from the back to the front (in place; false if blocked).
        template <class T> RAYT_HD bool traceFromFilm(const CameraParams<T>& c, const LensArrays<T>& L, Vec3<T>& o, Vec3<T>& d) {
            const int n = int(c.numSurfaces);
            for (int i = n - 1; i >= 0; --i) {
                // a negative spacing puts surface i behind surface i + 1, which the ray has just left
                const bool backward = i + 1 < n && L.surfaces[i].z < L.surfaces[i + 1].z;
                const T nTo = i > 0 ? L.surfaces[i - 1].nBehind : T(1);
                if (!crossSurface(c, L, L.surfaces[i], o, d, L.surfaces[i].nBehind, nTo, backward)) return false;
            }
            return true;
        }

        /// LensSystem::generateRay: from the film point (x, y) through the point (ux, uy) of the exit pupil's bounds (lens frame).
        template <class T> RAYT_HD bool generateRay(const CameraParams<T>& c, const LensArrays<T>& L, T x, T y, T ux, T uy, Vec3<T>& o, Vec3<T>& d, T& weight) {
            weight = T(0);
            const T r = khypot(x, y);
            if (c.numBands == 0 || r >= c.filmRadius) return false;
            const uint32_t band = kmin(uint32_t(r / c.filmRadius * T(c.numBands)), c.numBands - 1);
            const T* b = L.pupil + 5 * band;   // minX, minY, maxX, maxY, empty
            if (b[4] != T(0)) return false;
            // a point of the box, turned to the film point's angle (the box is for a film point on +x)
            const T px = b[0] + (b[2] - b[0]) * ux, py = b[1] + (b[3] - b[1]) * uy;
            const T cs = r > T(0) ? x / r : T(1), sn = r > T(0) ? y / r : T(0);
            const T z = L.surfaces[c.numSurfaces - 1].z - c.filmZ;
            o = Vec3<T>(x, y, c.filmZ);
            d = knormalizeGlm(Vec3<T>(px * cs - py * sn - x, px * sn + py * cs - y, z));
            const T cos4 = d.z * d.z * d.z * d.z;
            if (!traceFromFilm(c, L, o, d)) return false;
            weight = cos4 * (b[2] - b[0]) * (b[3] - b[1]) / (z * z) / c.centreIrradiance;
            return true;
        }

        /// LensSystem::centralRay: the middle of the pupil's box first, then a 7 x 7 grid over it, nearest to the middle first.
        template <class T> RAYT_HD bool centralRay(const CameraParams<T>& c, const LensArrays<T>& L, T x, T y, Vec3<T>& o, Vec3<T>& d) {
            const T r = khypot(x, y);
            if (c.numBands == 0 || r >= c.filmRadius) return false;
            for (int ring = 0; ring <= 3; ++ring)
                for (int a = -ring; a <= ring; ++a)
                    for (int bb = -ring; bb <= ring; ++bb) {
                        if (kmax(a < 0 ? -a : a, bb < 0 ? -bb : bb) != ring) continue;
                        T w;
                        if (generateRay(c, L, x, y, T(0.5) + T(a) / T(7), T(0.5) + T(bb) / T(7), o, d, w)) return true;
                    }
            return false;
        }

        /// Camera::sampleBlades: a uniform point of the polygonal aperture (one of its triangles, then a point in it).
        template <class T> RAYT_HD Vec3<T> sampleBlades(const CameraParams<T>& c, T ux, T uy) {
            const T scaled = ux * T(c.blades);
            const int k = kmin(int(scaled), int(c.blades) - 1);
            const T a = scaled - T(k);
            const T r = ksqrt(uy);
            const T step = kTwoPi<T>() / T(c.blades);
            const T a0 = c.bladeRotation + step * T(k);
            const T a1 = a0 + step;
            const T x = c.bladeRadius * r * ((T(1) - a) * kcos(a0) + a * kcos(a1));
            const T y = c.bladeRadius * r * ((T(1) - a) * ksin(a0) + a * ksin(a1));
            return Vec3<T>(x, y, T(0));
        }

    } // namespace lens

    /**
     * @brief The camera ray of film point (s, t) and lens sample (lu, lv), for a thin lens or a real lens (Camera::generateRay at the
     *        reference wavelength): origin o, unit direction d, and the weight of the sample (the real lens: its share of the light, 0 when
     *        the lens blocks it; the ray is then the pixel's central one, for the denoiser's features).
     */
    template <class T> RAYT_HD void lensCameraRay(const CameraParams<T>& c, const LensArrays<T>& L, T s, T t, T lu, T lv, Vec3<T>& o, Vec3<T>& d, T& weight) {
        const Vec3<T> U = lens::v3(c.u), V = lens::v3(c.v), W = lens::v3(c.w);
        if (c.type == kCamLens) {
            // Camera::lensRay: the film is turned 180 degrees (the lens forms an upside-down picture)
            const T x = -(ksaturate(s) - T(0.5)) * c.filmW, y = -(ksaturate(t) - T(0.5)) * c.filmH;
            Vec3<T> lo(0), ld(0, 0, 1);
            T w = T(0);
            const bool passed = lens::generateRay(c, L, x, y, lu, lv, lo, ld, w);
            if (!passed) lens::centralRay(c, L, x, y, lo, ld);
            weight = passed ? w : T(0);
            // Camera::lensToWorld
            const Vec3<T> axis = -W;
            o = lens::v3(c.origin) + U * lo.x + V * lo.y + axis * (lo.z - c.frontPrincipalZ);
            d = knormalizeGlm(U * ld.x + V * ld.y + axis * ld.z);
            return;
        }
        // the thin lens: the chief ray from the lens centre to the point of the plane of focus, moved to a point of the aperture
        const Vec3<T> target = lens::v3(c.lowerLeft) + ksaturate(s) * lens::v3(c.horizontal) + ksaturate(t) * lens::v3(c.vertical);
        Vec3<T> rd;
        if (c.blades >= 3) rd = lens::sampleBlades(c, lu, lv);
        else {
            const T r = ksqrt(lu), theta = kTwoPi<T>() * lv;   // sampling::UniformSampleDisk
            rd = c.lensRadius * Vec3<T>(r * kcos(theta), r * ksin(theta), T(0));
        }
        const Vec3<T> offset = U * rd.x + V * rd.y;
        o = lens::v3(c.origin) + offset;
        d = knormalizeGlm(target - o);
        weight = T(1);
    }

} // namespace rayt::kernel
