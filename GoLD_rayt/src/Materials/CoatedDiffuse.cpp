// Member functions of CoatedDiffuse that are defined here (not in Materials/CoatedDiffuse.hpp) so that changing them recompiles one file.
//
// Two normals: the specular lobe follows specNormal(rec) (the top of a liquid film where there is one, else the shading normal rec.n);
// the diffuse lobe and the subsurface event follow rec.n (the skin below the film). Without a film the two are the same vector.
#include "pch.h"

#include "Materials/CoatedDiffuse.hpp"
#include "Core/Sampling.hpp"

namespace rayt {

    Spectrum CoatedDiffuse::eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
        const SampledWavelengths& lambda, TransportMode /*mode*/) const
    {
        if (glm::dot(rec.gn, wo) <= 0 || glm::dot(rec.gn, wi) <= 0) return Spectrum(0.0);
        const Vector3& ns = specNormal(rec);
        const bool split = glm::dot(rec.ns, rec.ns) > Real(0);
        const Real cosO = glm::dot(rec.n, wo), cosI = glm::dot(rec.n, wi);
        const Real cosOs = split ? glm::dot(ns, wo) : cosO, cosIs = split ? glm::dot(ns, wi) : cosI;

        const Spectrum R = albedoAt(rec).reflectance(lambda);
        Spectrum f(0.0);
        if (cosOs > 0 && cosIs > 0) {
            f = specular(rec, wo, wi, R, cosOs, cosIs);
            // The integrator weighs by |cos(rec.n, wi)|; a lobe about the film's normal wants cos(ns, wi), so change one for the other.
            if (split) f = f * (cosIs / std::max(std::abs(cosI), Real(1e-3)));
        }
        if (!m_profile && cosO > 0 && cosI > 0)   // with subsurface scattering the diffuse lobe is the subsurface event
            f = f + diffuse(R, std::max(cosOs, Real(0)), std::max(cosIs, Real(0)), rec);
        return f;
    }

    std::optional<BSDFSample> CoatedDiffuse::sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u,
        const SampledWavelengths& lambda, TransportMode mode) const
    {
        if (glm::dot(rec.gn, wo) <= 0) return std::nullopt;

        const Vector3& ns = specNormal(rec);
        frame::Frame frame(ns);
        const Vector3 wo_local = frame.worldToLocal(wo);
        if (wo_local.z <= 0) return std::nullopt;

        const Real pS = lobeProbability(rec, wo_local.z);
        BSDFSample s;
        bool specularLobe = false;
        if (u.x < pS) {
            const GGXDistribution dist(alphaAt(rec), alphaAt(rec));
            const Vector3 wh = frame.localToWorld(dist.sample_wh(wo_local, Point2(u.x / pS, u.y)));
            s.wi = math::reflectOutward(wo, wh);
            specularLobe = true;
        }
        else if (m_profile) {
            // The light passes the coat (probability 1 - F) and enters: the integrator picks the exit.
            // Choosing it with probability 1 - pS gives the weight (1 - F) / (1 - pS), exactly 1 unless pS was clamped.
            s.f = Spectrum((Real(1) - coatF(rec, wo_local.z)) / (Real(1) - pS));
            s.pdf = Real(1);
            s.wi = rec.n;
            s.flags = BxDFFlags::Subsurface;
            return s;
        }
        else {
            const Real ux = (u.x - pS) / (Real(1) - pS);
            const Real r = math::safe_sqrt(ux), phi = Real(2) * constants::PI * u.y;
            const frame::Frame skin(rec.n);
            s.wi = skin.localToWorld(Vector3(r * std::cos(phi), r * std::sin(phi), math::safe_sqrt(Real(1) - ux)));
        }
        if (glm::dot(rec.gn, s.wi) <= 0 || glm::dot(specularLobe ? ns : rec.n, s.wi) <= 0) return std::nullopt;

        s.pdf = pdf(rec, wo, s.wi, lambda);
        s.f = eval(rec, wo, s.wi, lambda, mode);
        if (s.pdf <= Real(1e-6) || math::hasNonFinite(s.f)) return std::nullopt;
        s.flags = BxDFFlags::Reflection | BxDFFlags::Glossy;
        return s;
    }

    Real CoatedDiffuse::pdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
        const SampledWavelengths& /*lambda*/) const
    {
        if (glm::dot(rec.gn, wo) <= 0 || glm::dot(rec.gn, wi) <= 0) return Real(0);
        const Vector3& ns = specNormal(rec);
        const Real cosO = glm::dot(rec.n, wo), cosI = glm::dot(rec.n, wi);
        const Real cosOs = glm::dot(ns, wo), cosIs = glm::dot(ns, wi);
        if (cosOs <= 0) return Real(0);

        Vector3 wh = wo + wi;
        const Real len2 = glm::dot(wh, wh);
        if (len2 <= Real(0)) return Real(0);
        wh /= std::sqrt(len2);

        frame::Frame frame(ns);
        const GGXDistribution dist(alphaAt(rec), alphaAt(rec));
        const Vector3 wo_local = frame.worldToLocal(wo), wh_local = frame.worldToLocal(wh);
        const Real dotOH = std::abs(glm::dot(wo_local, wh_local));
        const Real pdfSpec = (cosIs > 0 && dotOH > 0) ? dist.pdf(wo_local, wh_local) / (Real(4) * dotOH) : Real(0);

        const Real pS = lobeProbability(rec, cosOs);
        if (m_profile) return pS * pdfSpec;   // the rest of the probability is the subsurface event (no direction)
        const Real pdfDiffuse = (cosI > 0 && cosO > 0) ? cosI / constants::PI : Real(0);
        return pS * pdfSpec + (Real(1) - pS) * pdfDiffuse;
    }

} // namespace rayt
