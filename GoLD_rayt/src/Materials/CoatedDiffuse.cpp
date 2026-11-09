// Member functions of CoatedDiffuse that are defined here (not in Materials/CoatedDiffuse.hpp) so that changing them recompiles one file.
#include "pch.h"

#include "Materials/CoatedDiffuse.hpp"
#include "Core/Sampling.hpp"

namespace rayt {

    Spectrum CoatedDiffuse::eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
        const SampledWavelengths& lambda, TransportMode /*mode*/) const
    {
        if (glm::dot(rec.gn, wo) <= 0 || glm::dot(rec.gn, wi) <= 0) return Spectrum(0.0);
        const Real cosO = glm::dot(rec.n, wo), cosI = glm::dot(rec.n, wi);
        if (cosO <= 0 || cosI <= 0) return Spectrum(0.0);

        const Spectrum R = albedoAt(rec).reflectance(lambda);
        if (m_profile) return specular(rec, wo, wi, R, cosO, cosI);   // the diffuse lobe is the subsurface event
        return diffuse(R, cosO, cosI, weightAt(rec)) + specular(rec, wo, wi, R, cosO, cosI);
    }

    std::optional<BSDFSample> CoatedDiffuse::sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u,
        const SampledWavelengths& lambda, TransportMode mode) const
    {
        if (glm::dot(rec.gn, wo) <= 0) return std::nullopt;

        frame::Frame frame(rec.n);
        const Vector3 wo_local = frame.worldToLocal(wo);
        if (wo_local.z <= 0) return std::nullopt;

        const Real pS = lobeProbability(rec, wo_local.z);
        BSDFSample s;
        if (u.x < pS) {
            const GGXDistribution dist(alphaAt(rec), alphaAt(rec));
            const Vector3 wh = frame.localToWorld(dist.sample_wh(wo_local, Point2(u.x / pS, u.y)));
            s.wi = math::reflectOutward(wo, wh);
        }
        else if (m_profile) {
            // The light passes the coat (probability 1 - w F) and enters: the integrator picks the exit.
            // Choosing it with probability 1 - pS gives the weight (1 - w F) / (1 - pS), exactly 1 unless pS was clamped.
            s.f = Spectrum((Real(1) - coatFresnel(wo_local.z, weightAt(rec))) / (Real(1) - pS));
            s.pdf = Real(1);
            s.wi = rec.n;
            s.flags = BxDFFlags::Subsurface;
            return s;
        }
        else {
            const Real ux = (u.x - pS) / (Real(1) - pS);
            const Real r = math::safe_sqrt(ux), phi = Real(2) * constants::PI * u.y;
            s.wi = frame.localToWorld(Vector3(r * std::cos(phi), r * std::sin(phi), math::safe_sqrt(Real(1) - ux)));
        }
        if (glm::dot(rec.gn, s.wi) <= 0 || glm::dot(rec.n, s.wi) <= 0) return std::nullopt;

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
        const Real cosO = glm::dot(rec.n, wo), cosI = glm::dot(rec.n, wi);
        if (cosO <= 0 || cosI <= 0) return Real(0);

        Vector3 wh = wo + wi;
        const Real len2 = glm::dot(wh, wh);
        if (len2 <= Real(0)) return Real(0);
        wh /= std::sqrt(len2);

        frame::Frame frame(rec.n);
        const GGXDistribution dist(alphaAt(rec), alphaAt(rec));
        const Vector3 wo_local = frame.worldToLocal(wo), wh_local = frame.worldToLocal(wh);
        const Real dotOH = std::abs(glm::dot(wo_local, wh_local));
        const Real pdfSpec = dotOH > 0 ? dist.pdf(wo_local, wh_local) / (Real(4) * dotOH) : Real(0);

        const Real pS = lobeProbability(rec, cosO);
        if (m_profile) return pS * pdfSpec;   // the rest of the probability is the subsurface event (no direction)
        return pS * pdfSpec + (Real(1) - pS) * cosI / constants::PI;
    }

} // namespace rayt
