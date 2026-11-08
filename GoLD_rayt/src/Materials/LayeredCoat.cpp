// LayeredCoat: a rough dielectric coat over a diffuse or glossy base, by a stochastic random walk (see Materials/LayeredCoat.hpp).
#include "pch.h"

#include "Materials/LayeredCoat.hpp"
#include "Microfacet/Distribution.hpp"

namespace rayt {

    namespace {

#if defined(RAYT_SPECTRAL_MODE)
        constexpr int kN = kSpectralSamples;
#else
        constexpr int kN = 3;
#endif

        Real maxOf(const Spectrum& s) {
            Real m = Real(0);
            for (int i = 0; i < kN; ++i) m = std::max(m, s[i]);
            return m;
        }

        /// Power heuristic with exponent 2 (one sample of each strategy).
        Real powerHeuristic(Real f, Real g) {
            if (std::isinf(f)) return Real(1);
            const Real f2 = f * f, g2 = g * g;
            return (f2 + g2) > Real(0) ? f2 / (f2 + g2) : Real(0);
        }

        /// Unpolarized Fresnel reflectance of a dielectric interface; cosThetaI < 0 means the light comes from the inside (eta is inverted).
        Real frDielectric(Real cosThetaI, Real eta) {
            cosThetaI = std::clamp(cosThetaI, Real(-1), Real(1));
            if (cosThetaI < Real(0)) { eta = Real(1) / eta; cosThetaI = -cosThetaI; }
            const Real sin2I = Real(1) - cosThetaI * cosThetaI;
            const Real sin2T = sin2I / (eta * eta);
            if (sin2T >= Real(1)) return Real(1);
            const Real cosT = std::sqrt(Real(1) - sin2T);
            const Real rParl = (eta * cosThetaI - cosT) / (eta * cosThetaI + cosT);
            const Real rPerp = (cosThetaI - eta * cosT) / (cosThetaI + eta * cosT);
            return Real(0.5) * (rParl * rParl + rPerp * rPerp);
        }

        /// Refraction of wi about n (n on the side of wi or not): returns false on total internal reflection. etap = relative index used.
        bool refract(const Vector3& wi, Vector3 n, Real eta, Real& etap, Vector3& wt) {
            Real cosI = glm::dot(n, wi);
            if (cosI < Real(0)) { eta = Real(1) / eta; cosI = -cosI; n = -n; }
            const Real sin2I = std::max(Real(0), Real(1) - cosI * cosI);
            const Real sin2T = sin2I / (eta * eta);
            if (sin2T >= Real(1)) return false;
            const Real cosT = std::sqrt(Real(1) - sin2T);
            wt = -wi / eta + (cosI / eta - cosT) * n;
            etap = eta;
            return true;
        }

        /// Isotropic GGX (Trowbridge-Reitz) in the local frame, with visible-normal sampling (Heitz 2018).
        struct Ggx {
            Real a;
            Real D(const Vector3& wm) const {
                const Real cos2 = wm.z * wm.z;
                if (cos2 <= Real(0)) return Real(0);
                const Real tan2 = (Real(1) - cos2) / cos2;
                if (std::isinf(tan2)) return Real(0);
                const Real cos4 = cos2 * cos2;
                const Real e = tan2 / (a * a);
                return Real(1) / (constants::PI * a * a * cos4 * (Real(1) + e) * (Real(1) + e));
            }
            Real lambda(const Vector3& w) const {
                const Real cos2 = w.z * w.z;
                if (cos2 <= Real(0)) return Real(0);
                const Real tan2 = (Real(1) - cos2) / cos2;
                if (std::isinf(tan2)) return Real(0);
                return (std::sqrt(Real(1) + a * a * tan2) - Real(1)) / Real(2);
            }
            Real G1(const Vector3& w) const { return Real(1) / (Real(1) + lambda(w)); }
            Real G(const Vector3& wo, const Vector3& wi) const { return Real(1) / (Real(1) + lambda(wo) + lambda(wi)); }
            /// Density of the visible normals as seen from w.
            Real pdfWm(const Vector3& w, const Vector3& wm) const {
                return G1(w) / std::abs(w.z) * D(wm) * std::abs(glm::dot(w, wm));
            }
            Vector3 sampleWm(const Vector3& w, const Point2& u) const {
                Vector3 wh = glm::normalize(Vector3(a * w.x, a * w.y, w.z));
                if (wh.z < Real(0)) wh = -wh;
                const Vector3 t1 = (wh.z < Real(0.99999)) ? glm::normalize(glm::cross(Vector3(0, 0, 1), wh)) : Vector3(1, 0, 0);
                const Vector3 t2 = glm::cross(wh, t1);
                const Real r = std::sqrt(u.x), phi = Real(2) * constants::PI * u.y;
                Real px = r * std::cos(phi), py = r * std::sin(phi);
                const Real h = std::sqrt(Real(1) - px * px);
                py = (Real(1) - (Real(1) + wh.z) / Real(2)) * h + (Real(1) + wh.z) / Real(2) * py;
                const Real pz = std::sqrt(std::max(Real(0), Real(1) - px * px - py * py));
                const Vector3 nh = px * t1 + py * t2 + pz * wh;
                return glm::normalize(Vector3(a * nh.x, a * nh.y, std::max(Real(1e-6), nh.z)));
            }
        };

        struct TopSample {
            Real f = 0;
            Vector3 wi{};
            Real pdf = 0;
            bool reflection = false;
        };

        /// The coat: a rough dielectric interface in the local frame (air above, z > 0; the coat medium below).
        struct Top {
            Real eta;
            Ggx g;

            Real f(const Vector3& wo, const Vector3& wi, bool radiance) const {
                const Real cosO = wo.z, cosI = wi.z;
                const bool reflect = cosI * cosO > Real(0);
                Real etap = Real(1);
                if (!reflect) etap = cosO > Real(0) ? eta : Real(1) / eta;
                Vector3 wm = wi * etap + wo;
                if (cosI == Real(0) || cosO == Real(0) || glm::dot(wm, wm) == Real(0)) return Real(0);
                wm = glm::normalize(wm);
                if (wm.z < Real(0)) wm = -wm;
                if (glm::dot(wm, wi) * cosI < Real(0) || glm::dot(wm, wo) * cosO < Real(0)) return Real(0);
                const Real F = frDielectric(glm::dot(wo, wm), eta);
                if (reflect) return g.D(wm) * g.G(wo, wi) * F / std::abs(Real(4) * cosI * cosO);
                const Real s = glm::dot(wi, wm) + glm::dot(wo, wm) / etap;
                const Real denom = s * s * cosI * cosO;
                Real ft = g.D(wm) * (Real(1) - F) * g.G(wo, wi) * std::abs(glm::dot(wi, wm) * glm::dot(wo, wm) / denom);
                if (radiance) ft /= etap * etap;
                return ft;
            }

            Real pdf(const Vector3& wo, const Vector3& wi, bool allowR, bool allowT) const {
                const Real cosO = wo.z, cosI = wi.z;
                const bool reflect = cosI * cosO > Real(0);
                Real etap = Real(1);
                if (!reflect) etap = cosO > Real(0) ? eta : Real(1) / eta;
                Vector3 wm = wi * etap + wo;
                if (cosI == Real(0) || cosO == Real(0) || glm::dot(wm, wm) == Real(0)) return Real(0);
                wm = glm::normalize(wm);
                if (wm.z < Real(0)) wm = -wm;
                if (glm::dot(wm, wi) * cosI < Real(0) || glm::dot(wm, wo) * cosO < Real(0)) return Real(0);
                const Real R = frDielectric(glm::dot(wo, wm), eta), T = Real(1) - R;
                const Real pr = allowR ? R : Real(0), pt = allowT ? T : Real(0);
                if (pr == Real(0) && pt == Real(0)) return Real(0);
                if (reflect) return g.pdfWm(wo, wm) / (Real(4) * std::abs(glm::dot(wo, wm))) * pr / (pr + pt);
                const Real s = glm::dot(wi, wm) + glm::dot(wo, wm) / etap;
                const Real dwmDwi = std::abs(glm::dot(wi, wm)) / (s * s);
                return g.pdfWm(wo, wm) * dwmDwi * pt / (pr + pt);
            }

            std::optional<TopSample> sample(const Vector3& wo, Real uc, const Point2& u, bool radiance, bool allowR, bool allowT) const {
                if (wo.z == Real(0)) return std::nullopt;
                const Vector3 wm = g.sampleWm(wo, u);
                const Real R = frDielectric(glm::dot(wo, wm), eta), T = Real(1) - R;
                const Real pr = allowR ? R : Real(0), pt = allowT ? T : Real(0);
                if (pr == Real(0) && pt == Real(0)) return std::nullopt;
                TopSample s;
                if (uc < pr / (pr + pt)) {
                    s.wi = Real(2) * glm::dot(wo, wm) * wm - wo;                     // reflect wo about wm
                    if (wo.z * s.wi.z <= Real(0)) return std::nullopt;
                    s.pdf = g.pdfWm(wo, wm) / (Real(4) * std::abs(glm::dot(wo, wm))) * pr / (pr + pt);
                    s.f = g.D(wm) * g.G(wo, s.wi) * R / (Real(4) * s.wi.z * wo.z);
                    s.reflection = true;
                    s.f = std::abs(s.f);
                    return s;
                }
                Real etap = Real(1);
                Vector3 wi;
                if (!refract(wo, wm, eta, etap, wi)) return std::nullopt;
                if (wo.z * wi.z > Real(0) || wi.z == Real(0)) return std::nullopt;
                const Real sdot = glm::dot(wi, wm) + glm::dot(wo, wm) / etap;
                const Real denom = sdot * sdot;
                const Real dwmDwi = std::abs(glm::dot(wi, wm)) / denom;
                s.wi = wi;
                s.pdf = g.pdfWm(wo, wm) * dwmDwi * pt / (pr + pt);
                s.f = T * g.D(wm) * g.G(wo, wi) * std::abs(glm::dot(wi, wm) * glm::dot(wo, wm) / (wi.z * wo.z * denom));
                if (radiance) s.f /= etap * etap;
                s.reflection = false;
                return s;
            }
        };

        struct BaseSample {
            Spectrum f;
            Vector3 wi;
            Real pdf = 0;
        };

        /// The base under the coat, in the local frame (z > 0 faces the coat): Lambertian or any reflection-only material.
        struct Base {
            const Material* mat;
            Spectrum R;
            const SampledWavelengths& lambda;
            SurfaceInteraction rec;

            Base(const Material* m, const Spectrum& r, const SampledWavelengths& l) : mat(m), R(r), lambda(l) {
                rec.n = rec.gn = Vector3(0, 0, 1);
                rec.dpdu = Vector3(1, 0, 0);   // the local x axis is the tangent of the base's anisotropy
            }
            Spectrum f(const Vector3& wo, const Vector3& wi) const {
                if (wo.z <= Real(0) || wi.z <= Real(0)) return Spectrum(0.0);
                if (!mat) return R * (Real(1) / constants::PI);
                return mat->eval(rec, wo, wi, lambda, TransportMode::Radiance);
            }
            Real pdf(const Vector3& wo, const Vector3& wi) const {
                if (wo.z <= Real(0) || wi.z <= Real(0)) return Real(0);
                if (!mat) return wi.z / constants::PI;
                return mat->pdf(rec, wo, wi, lambda);
            }
            std::optional<BaseSample> sample(const Vector3& wo, const Point2& u) const {
                if (wo.z <= Real(0)) return std::nullopt;
                BaseSample s;
                if (!mat) {
                    const Real r = std::sqrt(u.x), phi = Real(2) * constants::PI * u.y;
                    s.wi = Vector3(r * std::cos(phi), r * std::sin(phi), std::sqrt(std::max(Real(0), Real(1) - u.x)));
                    s.pdf = s.wi.z / constants::PI;
                    s.f = R * (Real(1) / constants::PI);
                    if (s.pdf <= Real(0)) return std::nullopt;
                    return s;
                }
                auto bs = mat->sample(rec, wo, u, lambda, TransportMode::Radiance);
                if (!bs || bs->pdf <= Real(0) || bs->wi.z <= Real(0)) return std::nullopt;
                s.f = bs->f; s.wi = bs->wi; s.pdf = bs->pdf;
                return s;
            }
        };

        /// The local frame of the layers: z along the shading normal, x along the surface tangent (the base's anisotropy axis), if there is one.
        frame::Frame surfaceFrame(const SurfaceInteraction& rec) {
            frame::Frame f;
            if (glm::dot(rec.dpdu, rec.dpdu) > Real(0)) f.buildFromNormalAndTangent(rec.n, rec.dpdu);
            else f.buildFromW(rec.n);
            return f;
        }

        /// Transmittance of the coat along a straight path across it in direction w (Beer-Lambert; sigma = optical depth along the normal).
        Spectrum transmittance(const Spectrum& sigma, const Vector3& w) {
            Spectrum t;
            const Real inv = Real(1) / std::max(std::abs(w.z), Real(1e-3));
            for (int i = 0; i < kN; ++i) t[i] = std::exp(-sigma[i] * inv);
            return t;
        }

        Spectrum coatDepth(const Vector3& tint, const SampledWavelengths& lambda) {
            const Spectrum T = RgbColor(tint).reflectance(lambda);
            Spectrum s;
            for (int i = 0; i < kN; ++i) s[i] = -std::log(std::clamp(T[i], Real(1e-4), Real(1)));
            return s;
        }

    } // namespace

    LayeredCoat::LayeredCoat(const Params& p, const RgbColor& albedo) : m_p(p), m_albedo(albedo) {
        m_alpha = MicrofacetDistribution::roughnessToAlpha(std::max(p.roughness, Real(0.02)));
    }

    LayeredCoat::LayeredCoat(const Params& p, std::shared_ptr<const Material> base) : m_p(p), m_base(std::move(base)) {
        m_alpha = MicrofacetDistribution::roughnessToAlpha(std::max(p.roughness, Real(0.02)));
    }

    Vector3 LayeredCoat::aovAlbedo(const SurfaceInteraction& rec) const {
        if (m_base) return m_base->aovAlbedo(rec);
        return glm::clamp(m_albedo.rgb(), Vector3(0.0), Vector3(1.0));
    }

    Real LayeredCoat::pdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const SampledWavelengths&) const {
        if (glm::dot(rec.gn, wo) <= Real(0) || glm::dot(rec.gn, wi) <= Real(0)) return Real(0);
        const frame::Frame frame = surfaceFrame(rec);
        const Vector3 o = frame.worldToLocal(wo), i = frame.worldToLocal(wi);
        if (o.z <= Real(0) || i.z <= Real(0)) return Real(0);
        const Top top{ m_p.ior, Ggx{ m_alpha } };
        // A deterministic stand-in for the density of the random walk (used for weighing against light sampling only): the coat's own
        // reflection lobe, mixed with what comes out of the base (a cosine lobe, plus a coat-like lobe if the base is glossy).
        const Real pTop = top.pdf(o, i, true, false);
        const Real pCos = i.z / constants::PI;
        const Real rr = std::clamp(frDielectric(o.z, m_p.ior), Real(0.05), Real(0.95));
        const Real pBase = m_base ? Real(0.5) * (pTop + pCos) : pCos;
        return rr * pTop + (Real(1) - rr) * pBase;
    }

    Spectrum LayeredCoat::eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
        const SampledWavelengths& lambda, TransportMode mode) const
    {
        if (glm::dot(rec.gn, wo) <= Real(0) || glm::dot(rec.gn, wi) <= Real(0)) return Spectrum(0.0);
        const frame::Frame frame = surfaceFrame(rec);
        const Vector3 o = frame.worldToLocal(wo), i = frame.worldToLocal(wi);
        if (o.z <= Real(0) || i.z <= Real(0)) return Spectrum(0.0);

        const bool radiance = (mode == TransportMode::Radiance);
        const Top top{ m_p.ior, Ggx{ m_alpha } };
        const Base base(m_base.get(), m_albedo.reflectance(lambda), lambda);
        const Spectrum sigma = coatDepth(m_p.tint, lambda);
        const int N = std::max(m_p.samples, 1);

        Spectrum f = Spectrum(top.f(o, i, radiance) * Real(N));    // the coat's own reflection

        for (int s = 0; s < N; ++s) {
            // In from outside along wo, and the "virtual light" that leaves along wi (reversed transport).
            const auto wos = top.sample(o, sampling::Random(), sampling::Random2D(), radiance, false, true);
            if (!wos || wos->f == Real(0) || wos->pdf == Real(0) || wos->wi.z == Real(0)) continue;
            const auto wis = top.sample(i, sampling::Random(), sampling::Random2D(), !radiance, false, true);
            if (!wis || wis->f == Real(0) || wis->pdf == Real(0) || wis->wi.z == Real(0)) continue;

            Spectrum beta = Spectrum(wos->f * std::abs(wos->wi.z) / wos->pdf);
            Vector3 w = wos->wi;                                   // travelling down
            for (int depth = 0; depth < m_p.maxDepth; ++depth) {
                const Real peak = maxOf(beta);
                if (depth > 3 && peak < Real(0.25)) {                // Russian roulette
                    const Real q = std::max(Real(0), Real(1) - peak);
                    if (sampling::Random() < q) break;
                    beta = beta * (Real(1) / (Real(1) - q));
                }
                beta = beta * transmittance(sigma, w);               // down to the base

                // Next event estimation: from the base straight out through the coat along the pre-sampled wis direction.
                const Spectrum fb = base.f(-w, -wis->wi);
                if (maxOf(fb) > Real(0)) {
                    const Real wt = powerHeuristic(wis->pdf, base.pdf(-w, -wis->wi));
                    f = f + beta * fb * transmittance(sigma, wis->wi) * (std::abs(wis->wi.z) * wt * wis->f / wis->pdf);
                }

                // Scatter at the base.
                const auto bs = base.sample(-w, sampling::Random2D());
                if (!bs || maxOf(bs->f) <= Real(0)) break;
                beta = beta * bs->f * (std::abs(bs->wi.z) / bs->pdf);
                w = bs->wi;                                          // travelling up

                // Next event estimation: leave through the coat towards wi along the sampled direction.
                const Real fExit = top.f(-w, i, radiance);
                if (fExit > Real(0)) {
                    // The other strategy would have sampled this very path from the outside: wi in, travelling along -w. Its density is
                    // the one of sampling the inside direction -w from wi (not of sampling wi from -w: the coat is not symmetric).
                    const Real exitPdf = top.pdf(i, -w, false, true);
                    const Real wt = powerHeuristic(bs->pdf, exitPdf);
                    f = f + beta * transmittance(sigma, w) * (fExit * wt);
                }

                // Up to the coat and back down by internal reflection.
                beta = beta * transmittance(sigma, w);
                const auto rs = top.sample(-w, sampling::Random(), sampling::Random2D(), radiance, true, false);
                if (!rs || rs->f == Real(0) || rs->pdf == Real(0) || rs->wi.z == Real(0)) break;
                beta = beta * (rs->f * std::abs(rs->wi.z) / rs->pdf);
                w = rs->wi;
            }
        }
        return f * (Real(1) / Real(N));
    }

    std::optional<BSDFSample> LayeredCoat::sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u,
        const SampledWavelengths& lambda, TransportMode mode) const
    {
        if (glm::dot(rec.gn, wo) <= Real(0)) return std::nullopt;
        const frame::Frame frame = surfaceFrame(rec);
        const Vector3 o = frame.worldToLocal(wo);
        if (o.z <= Real(0)) return std::nullopt;

        const bool radiance = (mode == TransportMode::Radiance);
        const Top top{ m_p.ior, Ggx{ m_alpha } };
        const Base base(m_base.get(), m_albedo.reflectance(lambda), lambda);
        const Spectrum sigma = coatDepth(m_p.tint, lambda);

        Vector3 exitDir;
        Spectrum beta;
        auto first = top.sample(o, sampling::Random(), u, radiance, true, true);
        if (!first || first->f == Real(0) || first->pdf == Real(0)) return std::nullopt;
        beta = Spectrum(first->f * std::abs(first->wi.z) / first->pdf);
        if (first->reflection) {
            exitDir = first->wi;                                    // reflected by the coat itself
        }
        else {
            Vector3 w = first->wi;                                  // travelling down
            bool left = false;
            for (int depth = 0; depth < m_p.maxDepth && !left; ++depth) {
                const Real peak = maxOf(beta);
                if (depth > 3 && peak < Real(0.25)) {
                    const Real q = std::max(Real(0), Real(1) - peak);
                    if (sampling::Random() < q) return std::nullopt;
                    beta = beta * (Real(1) / (Real(1) - q));
                }
                if (w.z == Real(0)) return std::nullopt;
                beta = beta * transmittance(sigma, w);              // down to the base
                const auto bs = base.sample(-w, sampling::Random2D());
                if (!bs || maxOf(bs->f) <= Real(0)) return std::nullopt;
                beta = beta * bs->f * (std::abs(bs->wi.z) / bs->pdf);
                w = bs->wi;                                         // travelling up
                beta = beta * transmittance(sigma, w);              // up to the coat
                const auto ts = top.sample(-w, sampling::Random(), sampling::Random2D(), radiance, true, true);
                if (!ts || ts->f == Real(0) || ts->pdf == Real(0)) return std::nullopt;
                beta = beta * (ts->f * std::abs(ts->wi.z) / ts->pdf);
                if (!ts->reflection) { exitDir = ts->wi; left = true; }
                else w = ts->wi;
            }
            if (!left) return std::nullopt;
        }
        if (exitDir.z <= Real(0)) return std::nullopt;

        BSDFSample out;
        out.wi = frame.localToWorld(exitDir);
        if (glm::dot(rec.gn, out.wi) <= Real(0)) return std::nullopt;
        out.pdf = pdf(rec, wo, out.wi, lambda);
        if (out.pdf <= Real(1e-8) || math::hasNonFinite(beta)) return std::nullopt;
        out.f = beta * (out.pdf / exitDir.z);                        // f cos / pdf = beta: the exact weight of the walk
        out.flags = BxDFFlags::Reflection | BxDFFlags::Glossy;
        return out;
    }

} // namespace rayt
