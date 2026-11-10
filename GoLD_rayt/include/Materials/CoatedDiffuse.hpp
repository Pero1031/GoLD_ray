/**
 * @file Materials/CoatedDiffuse.hpp
 * @brief Diffuse base under a glossy GGX layer, with a metallic blend (a compact "plastic / principled" material).
 *
 * Lobes (all share one GGX roughness):
 *  - diffuse   : (1 - m) R / pi  *  (1 - w F(mu_o)) (1 - w F(mu_i)) / (1 - w F_avg)
 *                R = albedo (constant x texture); the factors take away what the coat reflects,
 *                so the diffuse lobe never returns more than it receives (Kulla-Conty style layering).
 *  - specular  : F(theta_d) D G / (4 mu_o mu_i)  with
 *                F = (1 - m) w F_dielectric(theta_d, IOR 1.5) + m * Schlick(F0 = R, theta_d).
 *
 * w (specularWeight, 0..1) scales the dielectric coat (0 = matte), m (metallic, 0..1) blends toward a
 * metal whose reflectance color is the albedo. Roughness is perceptual: alpha = roughness^2.
 *
 * Sampling picks the specular lobe (GGX VNDF) or the diffuse lobe (cosine) with probabilities that
 * follow their weights; pdf() returns the same mixture, so MIS stays consistent.
 * Like single-scattering GGX it loses some energy at high roughness (see RoughConductor for the
 * compensated metal).
 */

#pragma once

#include <algorithm>
#include <memory>
#include <optional>

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Constants.hpp"
#include "Core/Interaction.hpp"
#include "Core/Math.hpp"
#include "Core/Fresnel.hpp"
#include "Core/Texture.hpp"
#include "Geometry/Frame.hpp"
#include "Materials/Material.hpp"
#include "Materials/Subsurface.hpp"
#include "Materials/SubsurfaceWalk.hpp"
#include "Materials/SkinRelief.hpp"
#include "Materials/SkinPigment.hpp"
#include "Materials/SkinFilm.hpp"
#include "Materials/SkinScatter.hpp"
#include "Materials/CosmeticLayer.hpp"
#include "Materials/DensityField.hpp"
#include "Microfacet/GGX.hpp"
#include "Color/SpectralConversion.hpp"

namespace rayt {

    class CoatedDiffuse : public Material {
    public:
        CoatedDiffuse(const RgbColor& albedo, Real roughness, Real specularWeight = Real(1), Real metallic = Real(0),
            std::shared_ptr<const Texture2D> texture = nullptr)
            : m_albedo(albedo), m_tex(std::move(texture))
            , m_w(std::clamp(specularWeight, Real(0), Real(1))), m_m(std::clamp(metallic, Real(0), Real(1)))
            , m_alpha(MicrofacetDistribution::roughnessToAlpha(roughness))
        {
            // Cosine-averaged dielectric Fresnel (IOR 1.5): 2 Integral F(mu) mu dmu.
            constexpr int n = 64;
            Real sum = 0;
            for (int i = 0; i < n; ++i) {
                const Real mu = (Real(i) + Real(0.5)) / Real(n);
                sum += fresnel::fresnelDielectric(mu, Real(1), Real(1.5)) * mu;
            }
            m_fAvg = Real(2) * sum / Real(n);
        }

        /**
         * @brief Turns the diffuse base into a subsurface scatterer (skin).
         *
         * The coat (the specular lobe, weight w) stays a local reflection. The light that passes the
         * coat enters the surface and leaves at another point of any surface sharing `profile`; the
         * diffuse color becomes the albedo of the radial profile. Set before rendering; metallic is ignored.
         */
        void setSubsurface(std::shared_ptr<const SubsurfaceProfile> profile) {
            m_profile = std::move(profile);
            if (m_profile) m_exit = std::make_unique<SubsurfaceExit>(m_w, m_profile->eta());
            else m_exit.reset();
            m_walk.reset();
        }

        /**
         * @brief Follows the subsurface light by a random walk through the inside of the object instead of the BSSRDF's jump along
         *        the surface (Materials/SubsurfaceWalk.hpp). Call after setSubsurface(); needs a closed surface with outward normals.
         *
         * The material becomes the boundary of a medium made from its albedo and the profile's scatter distance. The light that passes
         * the coat enters with a cosine distribution; a walk that reaches the surface from inside leaves through the BSSRDF's exit lobe.
         * The layers on the skin (pigment, flakes, emulsion, foundation) act as they do with the BSSRDF: factors where the light goes
         * in and where it comes out. A thick part looks as with the BSSRDF; a thin part lets through what its thickness allows.
         * @param g mean cosine of the phase function (0 = isotropic; skin is about 0.8). It does not change the color or the spread
         *          of a thick part.
         */
        void setSubsurfaceWalk(Real g = Real(0)) {
            if (m_profile && !m_profile->layered()) m_walk = std::make_shared<SubsurfaceWalkMedium>(m_albedo, m_profile, g);
        }
        const SubsurfaceWalkMedium* subsurfaceWalk() const { return m_walk.get(); }
        const Medium* mediumInside() const override { return m_walk.get(); }

        /// A micro-relief that bends the shading normal (skin furrows and pores, Materials/SkinRelief.hpp). Set before rendering.
        void setRelief(std::shared_ptr<const SkinRelief> relief) { m_relief = std::move(relief); }

        /// A liquid film that pools in the furrows of the relief (Materials/SkinFilm.hpp; needs the relief it was built on). The specular lobe then
        /// follows the smooth top of the film, with the film's Fresnel reflectance and roughness, where the film covers the skin. Set before rendering.
        void setFilm(std::shared_ptr<const SkinFilm> film) { m_film = std::move(film); }

        /// A white local scatter of flaking skin (Materials/SkinScatter.hpp): a third lobe, a Lambertian reflection under the coat that takes its share of the light
        /// that enters (what is left goes to the diffuse lobe or the subsurface event). A liquid film (setFilm) cuts it by index matching. Set before rendering.
        void setScatter(std::shared_ptr<const SkinScatter> scatter) { m_scatter = std::move(scatter); }

        /// The oil droplets of an emulsion in the liquid film (needs the film): a white layer on top of the skin whose transport optical depth is `tauMean` across the film of mean
        /// thickness and grows and shrinks with the local thickness of the film (it pools in the furrows). Reflects R = 0.75 tau / (1 + 0.75 tau) of the light that arrives (two-stream, no
        /// absorption); it is combined with the flakes below it into one white scatter (rec.surfScatter). Set before rendering.
        void setEmulsion(Real tauMean) { m_emulTau = std::max(tauMean, Real(0)); }

        /// A layer of foundation make-up on top of everything (Materials/CosmeticLayer.hpp, docs/COSMETICS.md). The reflection of the layer itself is
        /// computed by its random walk at every evaluation; the skin is seen through the layer's tables: the mirror-like lobe through the light that
        /// crosses the layer without a collision, the diffuse lobes and the subsurface light through all the light that crosses it, with the light
        /// that bounces between skin and layer summed (1 / (1 - Rback A)).
        /// `leveling` (0..1): how far the layer levels the relief of the skin (setRelief). 1: its surface follows the smooth normal, the furrows and pores
        /// do not show in its reflection (a paste fills them); 0: it follows the bent normal of the skin. The skin's own lobes keep the relief. Set before rendering.
        void setCosmetic(std::shared_ptr<const CosmeticCoat> coat, Real leveling = Real(1)) {
            m_cosmeticLeveling = std::clamp(leveling, Real(0), Real(1));
            m_cosmetic = std::move(coat);
            if (m_cosmetic) {
                m_cosmeticDiffuser = RgbColor(m_cosmetic->layer().params().diffuserAlbedo);
                m_cosmeticPlatelet = RgbColor(m_cosmetic->layer().params().plateletAlbedo);
            }
        }

        /// An uneven skin color (freckles, uneven pigment, redness: Materials/SkinPigment.hpp), and how much the pores of the relief darken the skin
        /// (0..1; needs the relief). With subsurface scattering the factor is applied where the light goes in and where it comes out (square root each).
        void setPigment(std::shared_ptr<const SkinPigment> pigment, Real poreDarkness = Real(0)) {
            m_pigment = std::move(pigment);
            m_poreDark = std::clamp(poreDarkness, Real(0), Real(1));
        }

        Spectrum subsurfaceExitTint(const SurfaceInteraction& pi, const SampledWavelengths& lambda) const override {
            Spectrum tint(1.0);
            if (hasPigment()) tint = RgbColor(glm::sqrt(pigmentFactor(pi.p, pi.gn))).reflectance(lambda);
            if (m_scatter || m_emulTau > Real(0)) {   // the light that comes out passes the white scatter (flakes, emulsion) too: (1 - rho), and what it sends back is scattered once more (1 / (1 - rho A))
                const Real rho = std::clamp(scatterTotal(pi.p, pi.gn), Real(0), Real(0.95));
                if (rho > Real(0)) {
                    const Spectrum A = m_albedo.reflectance(lambda);
                    Spectrum f(1.0);
                    for (int i = 0; i < kSssChannels; ++i) f[i] = (Real(1) - rho) / std::max(Real(1) - rho * std::clamp(A[i], Real(0), Real(1)), Real(0.05));
                    tint = tint * f;
                }
            }
            if (m_cosmetic) {   // the light that comes out crosses the foundation (its diffuse transmittance); what the layer sends back down is scattered by the skin once more
                const Spectrum T = RgbColor(m_cosmetic->throughDiffuse()).reflectance(lambda);
                const Spectrum Rb = RgbColor(m_cosmetic->backDiffuse()).reflectance(lambda);
                const Spectrum A = m_albedo.reflectance(lambda);
                Spectrum f(1.0);
                for (int i = 0; i < kSssChannels; ++i) f[i] = T[i] / std::max(Real(1) - std::clamp(Rb[i] * A[i], Real(0), Real(1)), Real(0.05));
                tint = tint * f;
            }
            return tint;
        }

        /// Oily and dry patches (sebum): the weight of the coat and its roughness change slowly over the skin, by up to +-`amount` (0..1) of their values,
        /// the oily patches having the stronger and sharper highlight. `frequency`: patches per metre.
        void setGlossVariation(Real amount, Real frequency, std::uint32_t seed = 1) {
            m_glossAmount = std::clamp(amount, Real(0), Real(0.95));
            m_gloss = m_glossAmount > Real(0) ? std::make_shared<NoiseDensity>(frequency, 2, Real(0.3), Real(0.7), seed) : nullptr;
        }

        void perturbShading(SurfaceInteraction& rec, const Vector3& wo) const override {
            if (m_walk && rec.medium == m_walk.get()) {
                // A random walk reached the surface from inside: turn the record into the outside view of the exit point (normals pointing out,
                // as the BSSRDF's exit vertex has them), so that the lights are chosen and the shadow rays start on the outer side.
                // eval / sample / pdf recognise it by wo lying below gn. The exit lobe follows the smooth normal, as the BSSRDF's does.
                rec.gn = -rec.gn;
                rec.n = -rec.n;
                rec.frontFace = true;
                return;
            }
            if (m_gloss) {
                const Real g = Real(2) * (*m_gloss)(rec.p) - Real(1);   // -1 (dry) .. 1 (oily)
                rec.glossScale = Real(1) + m_glossAmount * g;
                rec.roughScale = Real(1) - Real(0.6) * m_glossAmount * g;
            }
            const Vector3 smoothN = rec.n;
            if (m_relief) rec.n = m_relief->perturb(rec.p, smoothN, wo);
            Real filmThickness = 0;
            if (m_film && m_relief) {   // the specular lobe follows the top of the film where there is one, the bent skin normal elsewhere
                const SkinFilm::Hit f = m_film->at(rec.p, smoothN, wo);
                filmThickness = f.thickness;
                rec.filmCover = f.cover;
                Vector3 ns = glm::normalize(rec.n * (Real(1) - f.cover) + f.normal * f.cover);
                const Real c = glm::dot(ns, wo);
                if (c < Real(0.05)) ns = glm::normalize(ns + wo * (Real(0.05) - c));
                rec.ns = ns;
            }
            if (m_scatter || m_emulTau > Real(0)) rec.surfScatter = scatterTotal(rec.p, smoothN, filmThickness);
            if (m_cosmetic && m_relief && m_cosmeticLeveling > Real(0)) {   // the make-up levels the relief: its normal goes from the bent one towards the smooth one
                Vector3 nl = glm::normalize(rec.n * (Real(1) - m_cosmeticLeveling) + smoothN * m_cosmeticLeveling);
                const Real c = glm::dot(nl, wo);
                if (c < Real(0.05)) nl = glm::normalize(nl + wo * (Real(0.05) - c));
                rec.nl = nl;
            }
            if (hasPigment()) {   // computed once per hit; with SSS the light passes the pigment twice, so it is applied half going in and half coming out
                const Vector3 f = pigmentFactor(rec.p, rec.gn);
                rec.albedoTint = m_profile ? glm::sqrt(f) : f;
            }
        }

        const SubsurfaceProfile* subsurface() const override { return m_profile.get(); }
        const Material* subsurfaceExit() const override { return m_exit.get(); }

        void subsurfaceParams(const SurfaceInteraction& rec, const SampledWavelengths& lambda,
            SubsurfaceParams& out) const override
        {
            const Spectrum A = albedoAt(rec).reflectance(lambda);
            const Spectrum ell = m_profile->scatterDistance(lambda);
            for (int i = 0; i < kSssChannels; ++i) {
                out.albedo[i] = std::clamp(A[i], Real(0), Real(1));
                out.d[i] = std::max(ell[i], Real(1e-6)) / sss::shape(out.albedo[i]);
            }
        }

        Spectrum eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda, TransportMode /*mode*/) const override;

        std::optional<BSDFSample> sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u,
            const SampledWavelengths& lambda, TransportMode mode) const override;

        Real pdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& /*lambda*/) const override;

        bool isSpecular() const override { return false; }

        /// The base color (also the reflectance color of the metallic blend).
        Vector3 aovAlbedo(const SurfaceInteraction& rec) const override {
            return glm::clamp(albedoAt(rec).rgb(), Vector3(0.0), Vector3(1.0));
        }

    private:
        std::shared_ptr<const SkinRelief> m_relief;
        std::shared_ptr<const SkinFilm> m_film;
        std::shared_ptr<const SkinScatter> m_scatter;
        Real m_emulTau = Real(0);
        std::shared_ptr<const CosmeticCoat> m_cosmetic;
        Real m_cosmeticLeveling = Real(1);
        /// The normal of the make-up layer at this hit: the levelled one where perturbShading set it, else the shading normal.
        const Vector3& layerNormal(const SurfaceInteraction& rec) const { return glm::dot(rec.nl, rec.nl) > Real(0) ? rec.nl : rec.n; }
        /// The frame of the make-up layer: its normal, and the tangent of the shape where it gives one. With an arbitrary tangent the direction a
        /// platelet tilt leans to would jump along lines on the surface (a visible seam in the sheen).
        frame::Frame layerFrame(const SurfaceInteraction& rec) const {
            frame::Frame f;
            if (glm::dot(rec.dpdu, rec.dpdu) > Real(0)) f.buildFromNormalAndTangent(layerNormal(rec), rec.dpdu);
            else f.buildFromW(layerNormal(rec));
            return f;
        }
        RgbColor m_cosmeticDiffuser{ Vector3(1.0) }, m_cosmeticPlatelet{ Vector3(1.0) };   // the albedos of the layer's two particles

        /// One sample of the reflection of the foundation layer alone (a BSDF value; noisy, its mean is the layer's BSDF).
        Spectrum cosmeticReflection(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const SampledWavelengths& lambda) const;
        /// The density with which sample() draws wi from the layer's lobe: the mirror reflection off the platelets a first collision would meet, and a cosine lobe.
        Real cosmeticPdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi) const;
        /// The pdf of the lobes of the skin itself (what pdf() returned before there was a layer).
        Real skinPdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi) const;

        /// Random-walk mode, seen from outside along wo: the share of the light that enters the medium (`through`: past the foundation, the coat
        /// and the white scatter) and the probability with which sample() picks that event. False where nothing enters.
        bool walkEntry(const SurfaceInteraction& rec, const Vector3& wo, const SampledWavelengths& lambda, Spectrum& through, Real& probability) const;
        /// True if this hit is the end of a random walk inside the object (perturbShading turned the normals outward, so wo is below gn).
        bool walkExit(const SurfaceInteraction& rec, const Vector3& wo) const { return m_walk && glm::dot(rec.gn, wo) < Real(0); }

        /// The white scatter at p: the emulsion layer (on top) and the flakes (below) as two non-absorbing layers, R = Re + (1 - Re)^2 Rf / (1 - Re Rf).
        /// `localThickness` < 0: not known (the exit of a subsurface path): look it up.
        Real scatterTotal(const Point3& p, const Vector3& n, Real localThickness = Real(-1)) const {
            Real rf = 0;
            if (m_scatter) rf = m_scatter->reflectance(p, m_film ? m_film->params().ior : Real(1), m_film ? m_film->params().thickness : Real(0));
            Real re = 0;
            if (m_emulTau > Real(0) && m_film) {
                if (localThickness < Real(0)) localThickness = m_film->at(p, n, n).thickness;
                const Real tau = m_emulTau * localThickness / std::max(m_film->params().thickness, Real(1e-12));
                re = Real(0.75) * tau / (Real(1) + Real(0.75) * tau);
            }
            return re + (Real(1) - re) * (Real(1) - re) * rf / std::max(Real(1) - re * rf, Real(0.05));
        }
        std::shared_ptr<const SkinPigment> m_pigment;
        std::shared_ptr<const NoiseDensity> m_gloss;
        Real m_glossAmount = Real(0);
        Real m_poreDark = Real(0);

        bool hasPigment() const { return m_pigment || (m_poreDark > Real(0) && m_relief); }

        /// The reflectance multiplier of the pigment and of the pores at p.
        Vector3 pigmentFactor(const Point3& p, const Vector3& n) const {
            Vector3 f(1.0);
            if (m_pigment) f *= m_pigment->factor(p);
            if (m_poreDark > Real(0) && m_relief) f *= Vector3(1.0) - m_poreDark * m_relief->poreMask(p, n) * Vector3(0.30, 0.42, 0.45);
            return f;
        }

        RgbColor albedoAt(const SurfaceInteraction& rec) const {
            if (!m_tex && rec.albedoTint == Vector3(1.0)) return m_albedo;
            Vector3 c = m_albedo.rgb() * rec.albedoTint;
            if (m_tex) c *= m_tex->rgb(rec.uv.x, rec.uv.y);
            return RgbColor(c);
        }

        static Spectrum schlick(const Spectrum& f0, Real cosTheta) {
            const Real k = std::pow(Real(1) - std::clamp(cosTheta, Real(0), Real(1)), Real(5));
            return f0 + (Spectrum(1.0) - f0) * k;
        }

        /// The weight of the coat and the roughness at this hit (a skin with oily and dry patches varies them; otherwise the constants).
        Real weightAt(const SurfaceInteraction& rec) const { return std::clamp(m_w * rec.glossScale, Real(0), Real(1)); }
        Real alphaAt(const SurfaceInteraction& rec) const {
            const Real a = std::clamp(m_alpha * rec.roughScale, Real(constants::ALPHA_EPSILON), Real(1));
            if (!m_film || rec.filmCover <= Real(0)) return a;
            return (Real(1) - rec.filmCover) * a + rec.filmCover * std::max(m_film->alpha(), Real(constants::ALPHA_EPSILON));   // the film's surface is smooth
        }

        Real coatFresnel(Real cosTheta, Real w) const { return w * fresnel::fresnelDielectric(cosTheta, Real(1), Real(1.5)); }

        /// Reflectance of the top at this hit: the dry coat, or, where the film covers the skin, the film's.
        Real coatF(const SurfaceInteraction& rec, Real cosTheta) const {
            const Real dry = coatFresnel(cosTheta, weightAt(rec));
            if (!m_film || rec.filmCover <= Real(0)) return dry;
            return (Real(1) - rec.filmCover) * dry + rec.filmCover * m_film->fresnelReflectance(cosTheta);
        }
        Real fAvgAt(const SurfaceInteraction& rec) const {
            const Real dry = weightAt(rec) * m_fAvg;
            if (!m_film || rec.filmCover <= Real(0)) return dry;
            return (Real(1) - rec.filmCover) * dry + rec.filmCover * m_film->fresnelAverage();
        }

        /// The normal of the specular lobe: the top of the film where one was set at this hit, else the shading normal.
        const Vector3& specNormal(const SurfaceInteraction& rec) const { return glm::dot(rec.ns, rec.ns) > Real(0) ? rec.ns : rec.n; }

        /// The share of the entering light that the white scatter of flaking skin takes (0 where there is none).
        Real scatterAt(const SurfaceInteraction& rec) const { return std::clamp(rec.surfScatter, Real(0), Real(0.95)); }

        /// cosO, cosI: cosines with the normal of the specular lobe (the Fresnel reflectance of the top is a function of them).
        /// The base takes what the white scatter leaves on the way in, (1 - rho), and again on the way out, and what the scatter sends back down
        /// bounces between it and the base: (1 - rho) / (1 - rho R). Together with the scatter's own reflection rho this is the albedo
        /// rho + (1 - rho)^2 R / (1 - rho R) of a layer that absorbs nothing over the base, the same as the subsurface path's (subsurfaceExitTint).
        Spectrum diffuse(const Spectrum& R, Real cosO, Real cosI, const SurfaceInteraction& rec) const {
            const Real rho = scatterAt(rec);
            const Real k = (Real(1) - m_m) * (Real(1) - rho) * (Real(1) - coatF(rec, cosO)) * (Real(1) - coatF(rec, cosI))
                / (Real(1) - fAvgAt(rec)) / constants::PI;
            Spectrum f = R * k;
            if (rho > Real(0))
                for (int i = 0; i < kSssChannels; ++i) f[i] *= (Real(1) - rho) / std::max(Real(1) - rho * std::clamp(R[i], Real(0), Real(1)), Real(0.05));
            return f;
        }

        /// The white scatter: a Lambertian reflection of reflectance rho under the coat, normalised like the diffuse lobe.
        Real scatterLobe(const SurfaceInteraction& rec, Real cosO, Real cosI) const {
            return (Real(1) - m_m) * scatterAt(rec) * (Real(1) - coatF(rec, cosO)) * (Real(1) - coatF(rec, cosI))
                / (Real(1) - fAvgAt(rec)) / constants::PI;
        }

        Spectrum specular(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const Spectrum& R,
            Real cosO, Real cosI) const
        {
            Vector3 wh = wo + wi;
            const Real len2 = glm::dot(wh, wh);
            if (len2 <= Real(0)) return Spectrum(0.0);
            wh /= std::sqrt(len2);

            frame::Frame frame(specNormal(rec));
            const GGXDistribution dist(alphaAt(rec), alphaAt(rec));
            const Vector3 wo_l = frame.worldToLocal(wo), wi_l = frame.worldToLocal(wi), wh_l = frame.worldToLocal(wh);
            const Real D = dist.D(wh_l), G = dist.G(wo_l, wi_l);
            const Real cosD = std::abs(glm::dot(wi, wh));
            const Spectrum F = Spectrum(coatF(rec, cosD) * (Real(1) - m_m)) + schlick(R, cosD) * m_m;
            return F * (D * G / (Real(4) * cosO * cosI));
        }

        /// Probability of sampling the specular lobe: its weight against the diffuse lobe's.
        Real lobeProbability(const SurfaceInteraction& rec, Real cosO) const {
            const Real F = coatF(rec, cosO);
            if (m_profile) return std::clamp(F, Real(0.05), Real(0.95));   // reflected by the coat, else enters
            const Real sW = (Real(1) - m_m) * F + m_m * Real(0.6);
            const Vector3 c = albedoAt(rec).rgb();   // wavelength independent: the same value in eval(), sample() and pdf()
            const Real dW = (Real(1) - m_m) * std::max({ c.x, c.y, c.z }) * (Real(1) - F);
            return std::clamp(sW / std::max(sW + dW, Real(1e-6)), Real(0.05), Real(0.95));
        }

        RgbColor m_albedo;
        std::shared_ptr<const Texture2D> m_tex;
        Real m_w, m_m, m_alpha;
        Real m_fAvg = 0;
        std::shared_ptr<const SubsurfaceProfile> m_profile;   // null = plain coated diffuse
        std::unique_ptr<SubsurfaceExit> m_exit;                // exit vertex material (owns the coat weight / IOR)
        std::shared_ptr<const SubsurfaceWalkMedium> m_walk;   // set: the subsurface light walks through this medium (setSubsurfaceWalk)
    };

} // namespace rayt
