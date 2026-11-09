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

        // A layer of foundation on top: the mirror-like lobe is seen through the light that crosses the layer without a collision (both ways),
        // the diffuse lobes through all the light that crosses it, with the bounces between skin and layer summed.
        // The layer has its own normal (it levels the relief): what it lets through depends on the directions against that normal, and where
        // a direction is below the layer's plane nothing reaches the skin.
        const Real cosOl = m_cosmetic ? glm::dot(layerNormal(rec), wo) : cosO, cosIl = m_cosmetic ? glm::dot(layerNormal(rec), wi) : cosI;
        const bool layered = m_cosmetic && cosOl > 0 && cosIl > 0;
        Real straight = m_cosmetic ? Real(0) : Real(1);
        Spectrum through(m_cosmetic ? 0.0 : 1.0);
        if (layered) {
            straight = m_cosmetic->straight(cosOl) * m_cosmetic->straight(cosIl);
            through = RgbColor(m_cosmetic->through(cosOl)).reflectance(lambda) * RgbColor(m_cosmetic->through(cosIl)).reflectance(lambda);
            const Spectrum back = RgbColor(m_cosmetic->backDiffuse()).reflectance(lambda);
            for (int i = 0; i < kSssChannels; ++i) through[i] /= std::max(Real(1) - std::clamp(back[i] * R[i], Real(0), Real(1)), Real(0.05));
        }

        Spectrum f(0.0);
        if (cosOs > 0 && cosIs > 0) {
            f = specular(rec, wo, wi, R, cosOs, cosIs);
            // The integrator weighs by |cos(rec.n, wi)|; a lobe about the film's normal wants cos(ns, wi), so change one for the other.
            if (split) f = f * (cosIs / std::max(std::abs(cosI), Real(1e-3)));
            f = f * straight;
        }
        if (cosO > 0 && cosI > 0) {
            Spectrum below(0.0);
            if (!m_profile)   // with subsurface scattering the diffuse lobe is the subsurface event
                below = diffuse(R, std::max(cosOs, Real(0)), std::max(cosIs, Real(0)), rec);
            if (rec.surfScatter > Real(0))   // the white scatter of flaking skin: local and diffuse, in both cases
                below = below + Spectrum(scatterLobe(rec, std::max(cosOs, Real(0)), std::max(cosIs, Real(0))));
            f = f + below * through;
        }
        if (layered) f = f + cosmeticReflection(rec, wo, wi, lambda);
        return f;
    }

    Spectrum CoatedDiffuse::cosmeticReflection(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const SampledWavelengths& lambda) const {
        const frame::Frame skin = layerFrame(rec);
        const Vector3 o = skin.worldToLocal(wo), i = skin.worldToLocal(wi);
        if (o.z <= Real(1e-4) || i.z <= Real(1e-4)) return Spectrum(0.0);
        // The walk follows the light: it arrives travelling along -wi and is seen leaving along wo; its samples average to f cos(theta_o).
        // The integrator weighs by |cos(rec.n, wi)|; the layer's lobe wants the cosine with its own normal, so change one for the other.
        const Real swap = i.z / std::max(std::abs(glm::dot(rec.n, wi)), Real(1e-3));
        return m_cosmetic->layer().reflectionSampleT<Spectrum>(-i, o, m_cosmeticDiffuser.reflectance(lambda), m_cosmeticPlatelet.reflectance(lambda)) * (swap / o.z);
    }

    Real CoatedDiffuse::cosmeticPdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi) const {
        const frame::Frame skin = layerFrame(rec);
        const Vector3 o = skin.worldToLocal(wo), i = skin.worldToLocal(wi);
        if (o.z <= 0 || i.z <= 0) return Real(0);
        const CosmeticLayer& layer = m_cosmetic->layer();
        const Real q = std::min(layer.plateletShare(-o), Real(0.9));
        return q * layer.plateletPhase(-o, i) + (Real(1) - q) * i.z / constants::PI;
    }

    std::optional<BSDFSample> CoatedDiffuse::sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& uIn,
        const SampledWavelengths& lambda, TransportMode mode) const
    {
        if (glm::dot(rec.gn, wo) <= 0) return std::nullopt;

        const Vector3& ns = specNormal(rec);
        frame::Frame frame(ns);
        const Vector3 wo_local = frame.worldToLocal(wo);
        if (wo_local.z <= 0) return std::nullopt;

        BSDFSample s;

        // A layer of foundation: its own reflection is one more lobe, chosen with a probability that follows how much the layer reflects.
        Point2 u = uIn;
        Real pL = Real(0);
        const Real cosOn = m_cosmetic ? glm::dot(layerNormal(rec), wo) : Real(0);   // against the layer's own normal
        if (m_cosmetic && cosOn > 0) {
            pL = m_cosmetic->lobeProbability(cosOn);
            if (pL > 0 && u.x < pL) {
                const Real ux = u.x / pL;
                const frame::Frame skin = layerFrame(rec);
                const Vector3 o = skin.worldToLocal(wo);
                const CosmeticLayer& layer = m_cosmetic->layer();
                const Real q = std::min(layer.plateletShare(-o), Real(0.9));
                Vector3 il;
                if (ux < q) il = layer.samplePlatelet(-o, ux / q, u.y);   // the mirror reflection off a platelet, seen from the viewer's side
                else {
                    const Real v = (ux - q) / (Real(1) - q), r = math::safe_sqrt(v), phi = Real(2) * constants::PI * u.y;
                    il = Vector3(r * std::cos(phi), r * std::sin(phi), math::safe_sqrt(Real(1) - v));
                }
                if (il.z <= 0) return std::nullopt;
                s.wi = skin.localToWorld(il);
                if (glm::dot(rec.gn, s.wi) <= 0) return std::nullopt;
                s.pdf = pdf(rec, wo, s.wi, lambda);
                s.f = eval(rec, wo, s.wi, lambda, mode);
                if (s.pdf <= Real(1e-6) || math::hasNonFinite(s.f)) return std::nullopt;
                s.flags = BxDFFlags::Reflection | BxDFFlags::Glossy;
                return s;
            }
            if (pL > 0) u.x = (u.x - pL) / (Real(1) - pL);
        }

        const Real pS = lobeProbability(rec, wo_local.z);
        bool specularLobe = false;
        if (u.x < pS) {
            const GGXDistribution dist(alphaAt(rec), alphaAt(rec));
            const Vector3 wh = frame.localToWorld(dist.sample_wh(wo_local, Point2(u.x / pS, u.y)));
            s.wi = math::reflectOutward(wo, wh);
            specularLobe = true;
        }
        else if (m_profile) {
            // The light passes the coat (probability 1 - F) and enters. The white scatter of flaking skin takes the share rho of it, as a cosine lobe about n;
            // the rest (1 - rho) is the subsurface event: the integrator picks the exit.
            // Choosing the subsurface event with probability (1 - pS)(1 - rho) gives the weight (1 - F)(1 - rho) / ((1 - pS)(1 - rho)) = (1 - F) / (1 - pS).
            const Real rho = scatterAt(rec);
            const Real ux = (u.x - pS) / (Real(1) - pS);
            if (ux >= rho) {
                s.f = Spectrum((Real(1) - coatF(rec, wo_local.z)) / (Real(1) - pS));
                // under a layer of foundation: the light that comes out here crosses the layer towards the viewer (and this event was chosen with 1 - pL)
                if (m_cosmetic && cosOn > 0) s.f = s.f * RgbColor(m_cosmetic->through(cosOn)).reflectance(lambda) * (Real(1) / (Real(1) - pL));
                s.pdf = Real(1);
                s.wi = rec.n;
                s.flags = BxDFFlags::Subsurface;
                return s;
            }
            const Real r = math::safe_sqrt(ux / rho), phi = Real(2) * constants::PI * u.y;
            const frame::Frame skin(rec.n);
            s.wi = skin.localToWorld(Vector3(r * std::cos(phi), r * std::sin(phi), math::safe_sqrt(Real(1) - ux / rho)));
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
        if (!m_cosmetic) return skinPdf(rec, wo, wi);
        const Real cosOn = glm::dot(layerNormal(rec), wo);
        const Real pL = cosOn > 0 ? m_cosmetic->lobeProbability(cosOn) : Real(0);
        return pL * (pL > 0 ? cosmeticPdf(rec, wo, wi) : Real(0)) + (Real(1) - pL) * skinPdf(rec, wo, wi);
    }

    Real CoatedDiffuse::skinPdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi) const {
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
        if (m_profile) {   // the rest of the probability is the subsurface event (no direction), except the white scatter's share: a cosine lobe about n
            const Real pdfScatter = (cosI > 0 && cosO > 0) ? scatterAt(rec) * cosI / constants::PI : Real(0);
            return pS * pdfSpec + (Real(1) - pS) * pdfScatter;
        }
        const Real pdfDiffuse = (cosI > 0 && cosO > 0) ? cosI / constants::PI : Real(0);
        return pS * pdfSpec + (Real(1) - pS) * pdfDiffuse;
    }

} // namespace rayt
