#pragma once
/**
 * @file Kernel/Layered.hpp
 * @brief A clear coat over a diffuse or a rough metal base, as the engine's LayeredCoat evaluates it (docs/GPU.md).
 *
 * Materials/LayeredCoat.{hpp,cpp} step by step: the coat is a rough GGX dielectric interface (its own GGX in the coat's local
 * frame, with pbrt-v4's sampling of the visible normals), the base a Lambertian or the kernel's isotropic RoughConductor (built
 * with the coat as its incident medium), and the light between them is summed by random walks (pbrt-v4's LayeredBxDF): eval()
 * is the coat's reflection plus `samples` walks with next-event estimation at both layers, sample() one walk that leaves through
 * the coat (its weight is exact; its pdf the deterministic stand-in pdf() returns). With T = double it gives the engine's bits
 * when it draws the engine's random numbers; with T = float it is what the GPU runs.
 *
 * Random numbers come from `Rand`, which has T next() and Vec2<T> next2() (the engine's sampling::Random and Random2D). The draws
 * happen in the engine's order (LayeredCoat.cpp names each one before the call that uses it).
 */

#include "Kernel/Materials.hpp"

namespace rayt::kernel {

    /// The base under the coat: a Lambertian of reflectance `albedo`, or a rough conductor (`conductor` not null).
    template <class T> struct LayerBase {
        Vec3<T> albedo;
        const ConductorParams<T>* conductor;
    };

    namespace layered {

        template <class T> RAYT_HD_INLINE T maxOf(const Vec3<T>& s) {
            T m = T(0);
            m = kmax(m, s.x);
            m = kmax(m, s.y);
            m = kmax(m, s.z);
            return m;
        }

        template <class T> RAYT_HD_INLINE T powerHeuristic(T f, T g) {
            if (kisInf(f)) return T(1);
            const T f2 = f * f, g2 = g * g;
            return (f2 + g2) > T(0) ? f2 / (f2 + g2) : T(0);
        }

        template <class T> RAYT_HD_INLINE T kclamp(T v, T lo, T hi) { return v < lo ? lo : (hi < v ? hi : v); }   // std::clamp

        /// Unpolarized Fresnel reflectance; cosThetaI < 0: the light comes from the inside (eta inverted).
        template <class T> RAYT_HD T frDielectric(T cosThetaI, T eta) {
            cosThetaI = kclamp(cosThetaI, T(-1), T(1));
            if (cosThetaI < T(0)) { eta = T(1) / eta; cosThetaI = -cosThetaI; }
            const T sin2I = T(1) - cosThetaI * cosThetaI;
            const T sin2T = sin2I / (eta * eta);
            if (sin2T >= T(1)) return T(1);
            const T cosT = ksqrt(T(1) - sin2T);
            const T rParl = (eta * cosThetaI - cosT) / (eta * cosThetaI + cosT);
            const T rPerp = (cosThetaI - eta * cosT) / (cosThetaI + eta * cosT);
            return T(0.5) * (rParl * rParl + rPerp * rPerp);
        }

        /// Refraction of wi about n (either side): false on total internal reflection.
        template <class T> RAYT_HD bool refract(const Vec3<T>& wi, Vec3<T> n, T eta, T& etap, Vec3<T>& wt) {
            T cosI = kdot(n, wi);
            if (cosI < T(0)) { eta = T(1) / eta; cosI = -cosI; n = -n; }
            const T sin2I = kmax(T(0), T(1) - cosI * cosI);
            const T sin2T = sin2I / (eta * eta);
            if (sin2T >= T(1)) return false;
            const T cosT = ksqrt(T(1) - sin2T);
            wt = -wi / eta + (cosI / eta - cosT) * n;
            etap = eta;
            return true;
        }

        /// LayeredCoat's isotropic GGX in the local frame.
        template <class T> struct Ggx {
            T a;
            RAYT_HD T D(const Vec3<T>& wm) const {
                const T cos2 = wm.z * wm.z;
                if (cos2 <= T(0)) return T(0);
                const T tan2 = (T(1) - cos2) / cos2;
                if (kisInf(tan2)) return T(0);
                const T cos4 = cos2 * cos2;
                const T e = tan2 / (a * a);
                return T(1) / (kPi<T>() * a * a * cos4 * (T(1) + e) * (T(1) + e));
            }
            RAYT_HD T lambda(const Vec3<T>& w) const {
                const T cos2 = w.z * w.z;
                if (cos2 <= T(0)) return T(0);
                const T tan2 = (T(1) - cos2) / cos2;
                if (kisInf(tan2)) return T(0);
                return (ksqrt(T(1) + a * a * tan2) - T(1)) / T(2);
            }
            RAYT_HD T G1(const Vec3<T>& w) const { return T(1) / (T(1) + lambda(w)); }
            RAYT_HD T G(const Vec3<T>& wo, const Vec3<T>& wi) const { return T(1) / (T(1) + lambda(wo) + lambda(wi)); }
            RAYT_HD T pdfWm(const Vec3<T>& w, const Vec3<T>& wm) const { return G1(w) / kabs(w.z) * D(wm) * kabs(kdot(w, wm)); }
            RAYT_HD Vec3<T> sampleWm(const Vec3<T>& w, const Vec2<T>& u) const {
                Vec3<T> wh = knormalizeGlm(Vec3<T>(a * w.x, a * w.y, w.z));
                if (wh.z < T(0)) wh = -wh;
                const Vec3<T> t1 = (wh.z < T(0.99999)) ? knormalizeGlm(kcross(Vec3<T>(0, 0, 1), wh)) : Vec3<T>(1, 0, 0);
                const Vec3<T> t2 = kcross(wh, t1);
                const T r = ksqrt(u.x), phi = T(2) * kPi<T>() * u.y;
                T px = r * kcos(phi), py = r * ksin(phi);
                const T h = ksqrt(T(1) - px * px);
                py = (T(1) - (T(1) + wh.z) / T(2)) * h + (T(1) + wh.z) / T(2) * py;
                const T pz = ksqrt(kmax(T(0), T(1) - px * px - py * py));
                const Vec3<T> nh = px * t1 + py * t2 + pz * wh;
                return knormalizeGlm(Vec3<T>(a * nh.x, a * nh.y, kmax(T(1e-6), nh.z)));
            }
        };

        template <class T> struct TopSample {
            T f;
            Vec3<T> wi;
            T pdf;
            bool reflection;
        };

        /// The coat: a rough dielectric interface in the local frame (air above, z > 0; the coat below).
        template <class T> struct Top {
            T eta;
            Ggx<T> g;

            RAYT_HD T f(const Vec3<T>& wo, const Vec3<T>& wi, bool radiance) const {
                const T cosO = wo.z, cosI = wi.z;
                const bool reflect = cosI * cosO > T(0);
                T etap = T(1);
                if (!reflect) etap = cosO > T(0) ? eta : T(1) / eta;
                Vec3<T> wm = wi * etap + wo;
                if (cosI == T(0) || cosO == T(0) || kdot(wm, wm) == T(0)) return T(0);
                wm = knormalizeGlm(wm);
                if (wm.z < T(0)) wm = -wm;
                if (kdot(wm, wi) * cosI < T(0) || kdot(wm, wo) * cosO < T(0)) return T(0);
                const T F = frDielectric(kdot(wo, wm), eta);
                if (reflect) return g.D(wm) * g.G(wo, wi) * F / kabs(T(4) * cosI * cosO);
                const T s = kdot(wi, wm) + kdot(wo, wm) / etap;
                const T denom = s * s * cosI * cosO;
                T ft = g.D(wm) * (T(1) - F) * g.G(wo, wi) * kabs(kdot(wi, wm) * kdot(wo, wm) / denom);
                if (radiance) ft /= etap * etap;
                return ft;
            }

            RAYT_HD T pdf(const Vec3<T>& wo, const Vec3<T>& wi, bool allowR, bool allowT) const {
                const T cosO = wo.z, cosI = wi.z;
                const bool reflect = cosI * cosO > T(0);
                T etap = T(1);
                if (!reflect) etap = cosO > T(0) ? eta : T(1) / eta;
                Vec3<T> wm = wi * etap + wo;
                if (cosI == T(0) || cosO == T(0) || kdot(wm, wm) == T(0)) return T(0);
                wm = knormalizeGlm(wm);
                if (wm.z < T(0)) wm = -wm;
                if (kdot(wm, wi) * cosI < T(0) || kdot(wm, wo) * cosO < T(0)) return T(0);
                const T R = frDielectric(kdot(wo, wm), eta), Tr = T(1) - R;
                const T pr = allowR ? R : T(0), pt = allowT ? Tr : T(0);
                if (pr == T(0) && pt == T(0)) return T(0);
                if (reflect) return g.pdfWm(wo, wm) / (T(4) * kabs(kdot(wo, wm))) * pr / (pr + pt);
                const T s = kdot(wi, wm) + kdot(wo, wm) / etap;
                const T dwmDwi = kabs(kdot(wi, wm)) / (s * s);
                return g.pdfWm(wo, wm) * dwmDwi * pt / (pr + pt);
            }

            RAYT_HD bool sample(const Vec3<T>& wo, T uc, const Vec2<T>& u, bool radiance, bool allowR, bool allowT, TopSample<T>& s) const {
                if (wo.z == T(0)) return false;
                const Vec3<T> wm = g.sampleWm(wo, u);
                const T R = frDielectric(kdot(wo, wm), eta), Tr = T(1) - R;
                const T pr = allowR ? R : T(0), pt = allowT ? Tr : T(0);
                if (pr == T(0) && pt == T(0)) return false;
                if (uc < pr / (pr + pt)) {
                    s.wi = T(2) * kdot(wo, wm) * wm - wo;   // reflect wo about wm
                    if (wo.z * s.wi.z <= T(0)) return false;
                    s.pdf = g.pdfWm(wo, wm) / (T(4) * kabs(kdot(wo, wm))) * pr / (pr + pt);
                    s.f = g.D(wm) * g.G(wo, s.wi) * R / (T(4) * s.wi.z * wo.z);
                    s.reflection = true;
                    s.f = kabs(s.f);
                    return true;
                }
                T etap = T(1);
                Vec3<T> wi;
                if (!refract(wo, wm, eta, etap, wi)) return false;
                if (wo.z * wi.z > T(0) || wi.z == T(0)) return false;
                const T sdot = kdot(wi, wm) + kdot(wo, wm) / etap;
                const T denom = sdot * sdot;
                const T dwmDwi = kabs(kdot(wi, wm)) / denom;
                s.wi = wi;
                s.pdf = g.pdfWm(wo, wm) * dwmDwi * pt / (pr + pt);
                s.f = Tr * g.D(wm) * g.G(wo, wi) * kabs(kdot(wi, wm) * kdot(wo, wm) / (wi.z * wo.z * denom));
                if (radiance) s.f /= etap * etap;
                s.reflection = false;
                return true;
            }
        };

        /// The base in the local frame (z > 0 faces the coat). The engine evaluates a material base with n = gn = +z.
        template <class T> RAYT_HD Vec3<T> baseF(const LayerBase<T>& b, const Vec3<T>& wo, const Vec3<T>& wi) {
            if (wo.z <= T(0) || wi.z <= T(0)) return Vec3<T>(0);
            if (!b.conductor) return b.albedo * (T(1) / kPi<T>());
            const Vec3<T> z(0, 0, 1);
            return conductorEval(*b.conductor, wo, wi, z, z);
        }
        template <class T> RAYT_HD T basePdf(const LayerBase<T>& b, const Vec3<T>& wo, const Vec3<T>& wi) {
            if (wo.z <= T(0) || wi.z <= T(0)) return T(0);
            if (!b.conductor) return wi.z / kPi<T>();
            const Vec3<T> z(0, 0, 1);
            return conductorPdf(*b.conductor, wo, wi, z, z);
        }
        template <class T> RAYT_HD bool baseSample(const LayerBase<T>& b, const Vec3<T>& wo, const Vec2<T>& u, Vec3<T>& f, Vec3<T>& wi, T& pdf) {
            if (wo.z <= T(0)) return false;
            if (!b.conductor) {
                const T r = ksqrt(u.x), phi = T(2) * kPi<T>() * u.y;
                wi = Vec3<T>(r * kcos(phi), r * ksin(phi), ksqrt(kmax(T(0), T(1) - u.x)));
                pdf = wi.z / kPi<T>();
                f = b.albedo * (T(1) / kPi<T>());
                return pdf > T(0);
            }
            const Vec3<T> z(0, 0, 1);
            BsdfSample<T> bs{};
            if (!conductorSample(*b.conductor, wo, z, z, u.x, u.y, bs) || bs.pdf <= T(0) || bs.wi.z <= T(0)) return false;
            f = bs.f; wi = bs.wi; pdf = bs.pdf;
            return true;
        }

        /// Beer-Lambert across the coat in direction w (sigma: optical depth along the normal). A clear coat skips the exponentials
        /// (exp(-0) is 1, so the product is the same bits).
        template <class T> RAYT_HD_INLINE Vec3<T> transmittance(const LayeredParams<T>& p, const Vec3<T>& w) {
            if (!p.absorbing) return Vec3<T>(1);
            const T inv = T(1) / kmax(kabs(w.z), T(1e-3));
            return Vec3<T>(kexp(-p.sigma[0] * inv), kexp(-p.sigma[1] * inv), kexp(-p.sigma[2] * inv));
        }

        /// The layers' frame: z along the shading normal (frame::Frame::buildFromW; the engine uses the surface tangent where there is
        /// one, which only turns the frame about the normal: the coat and the bases are isotropic).
        template <class T> RAYT_HD_INLINE Frame3<T> frame(const Vec3<T>& ns) { return Frame3<T>::fromNormal(ns); }

    } // namespace layered

    /// LayeredCoat::pdf: a deterministic stand-in for the density of the walk (the coat's reflection lobe mixed with a cosine lobe,
    /// and a coat-like lobe for a glossy base).
    template <class T> RAYT_HD T layeredPdf(const LayeredParams<T>& p, const LayerBase<T>& base, const Vec3<T>& wo, const Vec3<T>& wi, const Vec3<T>& ns,
        const Vec3<T>& ng) {
        if (kdot(ng, wo) <= T(0) || kdot(ng, wi) <= T(0)) return T(0);
        const Frame3<T> fr = layered::frame(ns);
        const Vec3<T> o = fr.toLocal(wo), i = fr.toLocal(wi);
        if (o.z <= T(0) || i.z <= T(0)) return T(0);
        const layered::Top<T> top{ p.ior, layered::Ggx<T>{ p.alpha } };
        const T pTop = top.pdf(o, i, true, false);
        const T pCos = i.z / kPi<T>();
        const T rr = layered::kclamp(layered::frDielectric(o.z, p.ior), T(0.05), T(0.95));
        const T pBase = base.conductor ? T(0.5) * (pTop + pCos) : pCos;
        return rr * pTop + (T(1) - rr) * pBase;
    }

    /// LayeredCoat::eval (radiance transport).
    template <class T, class Rand> RAYT_HD Vec3<T> layeredEval(const LayeredParams<T>& p, const LayerBase<T>& base, const Vec3<T>& wo, const Vec3<T>& wi,
        const Vec3<T>& ns, const Vec3<T>& ng, Rand& rnd) {
        using namespace layered;
        if (kdot(ng, wo) <= T(0) || kdot(ng, wi) <= T(0)) return Vec3<T>(0);
        const Frame3<T> fr = frame(ns);
        const Vec3<T> o = fr.toLocal(wo), i = fr.toLocal(wi);
        if (o.z <= T(0) || i.z <= T(0)) return Vec3<T>(0);
        const Top<T> top{ p.ior, Ggx<T>{ p.alpha } };
        const int N = kmax(int(p.samples), 1);

        Vec3<T> f = Vec3<T>(top.f(o, i, true) * T(N));   // the coat's own reflection
        for (int s = 0; s < N; ++s) {
            // In from outside along wo, and the "virtual light" that leaves along wi (reversed transport).
            const Vec2<T> uo = rnd.next2();
            const T co = rnd.next();
            TopSample<T> wos{};
            if (!top.sample(o, co, uo, true, false, true, wos) || wos.f == T(0) || wos.pdf == T(0) || wos.wi.z == T(0)) continue;
            const Vec2<T> ui = rnd.next2();
            const T ci = rnd.next();
            TopSample<T> wis{};
            if (!top.sample(i, ci, ui, false, false, true, wis) || wis.f == T(0) || wis.pdf == T(0) || wis.wi.z == T(0)) continue;

            Vec3<T> beta = Vec3<T>(wos.f * kabs(wos.wi.z) / wos.pdf);
            Vec3<T> w = wos.wi;   // travelling down
            for (int depth = 0; depth < int(p.maxDepth); ++depth) {
                const T peak = maxOf(beta);
                if (depth > 3 && peak < T(0.25)) {   // Russian roulette
                    const T q = kmax(T(0), T(1) - peak);
                    if (rnd.next() < q) break;
                    beta = beta * (T(1) / (T(1) - q));
                }
                beta = beta * transmittance(p, w);   // down to the base

                // next-event estimation: from the base straight out through the coat along the pre-sampled direction
                const Vec3<T> fb = baseF(base, -w, -wis.wi);
                if (maxOf(fb) > T(0)) {
                    const T wt = powerHeuristic(wis.pdf, basePdf(base, -w, -wis.wi));
                    f = f + beta * fb * transmittance(p, wis.wi) * (kabs(wis.wi.z) * wt * wis.f / wis.pdf);
                }

                // scatter at the base
                const Vec2<T> ub = rnd.next2();
                Vec3<T> bf, bwi;
                T bpdf;
                if (!baseSample(base, -w, ub, bf, bwi, bpdf) || maxOf(bf) <= T(0)) break;
                beta = beta * bf * (kabs(bwi.z) / bpdf);
                w = bwi;   // travelling up

                // next-event estimation: leave through the coat towards wi
                const T fExit = top.f(-w, i, true);
                if (fExit > T(0)) {
                    const T exitPdf = top.pdf(i, -w, false, true);
                    const T wt = powerHeuristic(bpdf, exitPdf);
                    f = f + beta * transmittance(p, w) * (fExit * wt);
                }

                // up to the coat and back down by internal reflection
                beta = beta * transmittance(p, w);
                const Vec2<T> ur = rnd.next2();
                const T cr = rnd.next();
                TopSample<T> rs{};
                if (!top.sample(-w, cr, ur, true, true, false, rs) || rs.f == T(0) || rs.pdf == T(0) || rs.wi.z == T(0)) break;
                beta = beta * (rs.f * kabs(rs.wi.z) / rs.pdf);
                w = rs.wi;
            }
        }
        return f * (T(1) / T(N));
    }

    /// LayeredCoat::sample (radiance transport): u for the coat's first micro-normal, the rest from rnd. f is the walk's exact weight
    /// times pdf / cos, so that f cos / pdf is that weight.
    template <class T, class Rand> RAYT_HD bool layeredSample(const LayeredParams<T>& p, const LayerBase<T>& base, const Vec3<T>& wo, const Vec3<T>& ns,
        const Vec3<T>& ng, const Vec2<T>& u, Rand& rnd, BsdfSample<T>& out) {
        using namespace layered;
        out.delta = false;
        if (kdot(ng, wo) <= T(0)) return false;
        const Frame3<T> fr = frame(ns);
        const Vec3<T> o = fr.toLocal(wo);
        if (o.z <= T(0)) return false;
        const Top<T> top{ p.ior, Ggx<T>{ p.alpha } };

        Vec3<T> exitDir;
        Vec3<T> beta;
        const T c0 = rnd.next();
        TopSample<T> first{};
        if (!top.sample(o, c0, u, true, true, true, first) || first.f == T(0) || first.pdf == T(0)) return false;
        beta = Vec3<T>(first.f * kabs(first.wi.z) / first.pdf);
        if (first.reflection) exitDir = first.wi;   // reflected by the coat itself
        else {
            Vec3<T> w = first.wi;   // travelling down
            bool left = false;
            for (int depth = 0; depth < int(p.maxDepth) && !left; ++depth) {
                const T peak = maxOf(beta);
                if (depth > 3 && peak < T(0.25)) {
                    const T q = kmax(T(0), T(1) - peak);
                    if (rnd.next() < q) return false;
                    beta = beta * (T(1) / (T(1) - q));
                }
                if (w.z == T(0)) return false;
                beta = beta * transmittance(p, w);   // down to the base
                const Vec2<T> ub = rnd.next2();
                Vec3<T> bf, bwi;
                T bpdf;
                if (!baseSample(base, -w, ub, bf, bwi, bpdf) || maxOf(bf) <= T(0)) return false;
                beta = beta * bf * (kabs(bwi.z) / bpdf);
                w = bwi;   // travelling up
                beta = beta * transmittance(p, w);   // up to the coat
                const Vec2<T> ut = rnd.next2();
                const T ct = rnd.next();
                TopSample<T> ts{};
                if (!top.sample(-w, ct, ut, true, true, true, ts) || ts.f == T(0) || ts.pdf == T(0)) return false;
                beta = beta * (ts.f * kabs(ts.wi.z) / ts.pdf);
                if (!ts.reflection) { exitDir = ts.wi; left = true; }
                else w = ts.wi;
            }
            if (!left) return false;
        }
        if (exitDir.z <= T(0)) return false;
        out.wi = fr.toWorld(exitDir);
        if (kdot(ng, out.wi) <= T(0)) return false;
        out.pdf = layeredPdf(p, base, wo, out.wi, ns, ng);
        if (out.pdf <= T(1e-8) || !kisFinite(beta.x) || !kisFinite(beta.y) || !kisFinite(beta.z)) return false;
        out.f = beta * (out.pdf / exitDir.z);   // f cos / pdf = beta: the exact weight of the walk
        return true;
    }

} // namespace rayt::kernel
