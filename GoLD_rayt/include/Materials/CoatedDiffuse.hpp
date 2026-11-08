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
#include "Materials/SkinRelief.hpp"
#include "Materials/SkinPigment.hpp"
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
        }

        /// A micro-relief that bends the shading normal (skin furrows and pores, Materials/SkinRelief.hpp). Set before rendering.
        void setRelief(std::shared_ptr<const SkinRelief> relief) { m_relief = std::move(relief); }

        /// An uneven skin color (freckles, uneven pigment, redness: Materials/SkinPigment.hpp), and how much the pores of the relief darken the skin
        /// (0..1; needs the relief). With subsurface scattering the factor is applied where the light goes in and where it comes out (square root each).
        void setPigment(std::shared_ptr<const SkinPigment> pigment, Real poreDarkness = Real(0)) {
            m_pigment = std::move(pigment);
            m_poreDark = std::clamp(poreDarkness, Real(0), Real(1));
        }

        Spectrum subsurfaceExitTint(const SurfaceInteraction& pi, const SampledWavelengths& lambda) const override {
            if (!hasPigment()) return Spectrum(1.0);
            return RgbColor(glm::sqrt(pigmentFactor(pi.p))).reflectance(lambda);
        }

        void perturbShading(SurfaceInteraction& rec, const Vector3& wo) const override {
            if (m_relief) rec.n = m_relief->perturb(rec.p, rec.n, wo);
            if (hasPigment()) {   // computed once per hit; with SSS the light passes the pigment twice, so it is applied half going in and half coming out
                const Vector3 f = pigmentFactor(rec.p);
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
        std::shared_ptr<const SkinPigment> m_pigment;
        Real m_poreDark = Real(0);

        bool hasPigment() const { return m_pigment || (m_poreDark > Real(0) && m_relief); }

        /// The reflectance multiplier of the pigment and of the pores at p.
        Vector3 pigmentFactor(const Point3& p) const {
            Vector3 f(1.0);
            if (m_pigment) f *= m_pigment->factor(p);
            if (m_poreDark > Real(0) && m_relief) f *= Vector3(1.0) - m_poreDark * m_relief->poreMask(p) * Vector3(0.30, 0.42, 0.45);
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

        Real coatFresnel(Real cosTheta) const { return m_w * fresnel::fresnelDielectric(cosTheta, Real(1), Real(1.5)); }

        Spectrum diffuse(const Spectrum& R, Real cosO, Real cosI) const {
            const Real k = (Real(1) - m_m) * (Real(1) - coatFresnel(cosO)) * (Real(1) - coatFresnel(cosI))
                / (Real(1) - m_w * m_fAvg) / constants::PI;
            return R * k;
        }

        Spectrum specular(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const Spectrum& R,
            Real cosO, Real cosI) const
        {
            Vector3 wh = wo + wi;
            const Real len2 = glm::dot(wh, wh);
            if (len2 <= Real(0)) return Spectrum(0.0);
            wh /= std::sqrt(len2);

            frame::Frame frame(rec.n);
            const GGXDistribution dist(m_alpha, m_alpha);
            const Vector3 wo_l = frame.worldToLocal(wo), wi_l = frame.worldToLocal(wi), wh_l = frame.worldToLocal(wh);
            const Real D = dist.D(wh_l), G = dist.G(wo_l, wi_l);
            const Real cosD = std::abs(glm::dot(wi, wh));
            const Spectrum F = Spectrum(coatFresnel(cosD) * (Real(1) - m_m)) + schlick(R, cosD) * m_m;
            return F * (D * G / (Real(4) * cosO * cosI));
        }

        /// Probability of sampling the specular lobe: its weight against the diffuse lobe's.
        Real lobeProbability(const SurfaceInteraction& rec, Real cosO) const {
            if (m_profile) return std::clamp(coatFresnel(cosO), Real(0.05), Real(0.95));   // reflected by the coat, else enters
            const Real sW = (Real(1) - m_m) * coatFresnel(cosO) + m_m * Real(0.6);
            const Vector3 c = albedoAt(rec).rgb();   // wavelength independent: the same value in eval(), sample() and pdf()
            const Real dW = (Real(1) - m_m) * std::max({ c.x, c.y, c.z }) * (Real(1) - coatFresnel(cosO));
            return std::clamp(sW / std::max(sW + dW, Real(1e-6)), Real(0.05), Real(0.95));
        }

        RgbColor m_albedo;
        std::shared_ptr<const Texture2D> m_tex;
        Real m_w, m_m, m_alpha;
        Real m_fAvg = 0;
        std::shared_ptr<const SubsurfaceProfile> m_profile;   // null = plain coated diffuse
        std::unique_ptr<SubsurfaceExit> m_exit;                // exit vertex material (owns the coat weight / IOR)
    };

} // namespace rayt
