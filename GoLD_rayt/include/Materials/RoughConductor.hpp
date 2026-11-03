/**
 * @file Materials/RoughConductor.hpp
 * @brief Rough conductor material based on Microfacet theory (Cook-Torrance).
 * Implements GGX distribution with VNDF sampling and complex Fresnel terms.
 */

#pragma once

#include <algorithm>
#include <complex>
#include <memory>
#include <optional>
#include <vector>

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Constants.hpp"
#include "Core/Interaction.hpp"
#include "Core/Math.hpp"
#include "Core/Fresnel.hpp"

#include "Geometry/Frame.hpp"
#include "Materials/Material.hpp"
#include "Microfacet/GGX.hpp"
#include "Microfacet/GGXAlbedo.hpp"
#include "Color/ComplexIorTable.hpp"
#include "Color/IorRgbApprox.hpp"

namespace rayt {

    /**
     * @class RoughConductor
     * @brief Models a rough conductive surface using the Cook-Torrance BRDF.
     * * The BRDF is defined as: fr = (D * G * F) / (4 * (n.wi) * (n.wo)).
     * Features:
     * - Trowbridge-Reitz (GGX) Normal Distribution Function.
     * - Smith Shadowing-Masking function.
     * - Fresnel equations for conductors using complex refractive indices.
     * - Visible Normal Distribution Function (VNDF) importance sampling.
     */
    class RoughConductor : public Material {
#if defined(RAYT_SPECTRAL_MODE)
        // Complex IOR resampled on a 1 nm grid for fast per-wavelength lookup.
        static constexpr int kLutFirstNm = 360;
        static constexpr int kLutLastNm = 830;
        std::vector<Real> etaLut;   // Index of Refraction (Real) per nm
        std::vector<Real> kLut;     // Extinction Coefficient per nm
#else
        Spectrum eta;   // Index of Refraction (Real), RGB
        Spectrum k;     // Extinction Coefficient, RGB
#endif
        Real alpha_x;   // Roughness X
        Real alpha_y;   // Roughness Y (same as X if isotropic)
        Vector3 m_f0Rgb = Vector3(1.0);   // normal-incidence reflectance in RGB (denoiser albedo AOV)

        // ---- Multiple-scattering compensation (Kulla & Conty 2017) ----
        // Single-scattering GGX loses energy that grows with roughness. A diffuse-like lobe
        //     f_ms = F_ms (1 - E(mu_o)) (1 - E(mu_i)) / (pi (1 - E_avg))
        //     F_ms = F_avg^2 E_avg / (1 - F_avg (1 - E_avg))
        // adds it back (E = directional albedo of the single-scattering lobe, E_avg its cosine
        // average, F_avg the cosine-averaged Fresnel reflectance, per wavelength / RGB channel).
        // Anisotropic roughness uses the geometric-mean alpha for E (an approximation).
        bool m_msActive = false;                              // true if compensation is on and needed
        bool m_msEnabled = true;
        Real m_eAvg = 0;
        std::array<Real, GGXAlbedoTable::kMu> m_eRow{};       // E(mu) at the table's mu nodes for this alpha
#if defined(RAYT_SPECTRAL_MODE)
        std::vector<Real> m_fmsLut;                           // F_ms per nm (same grid as etaLut)
#else
        Spectrum m_fms = Spectrum(0.0);
#endif

        /// Cosine-weighted average Fresnel reflectance: 2 Integral F(mu) mu dmu.
        static Real averageFresnel(Real eta, Real k) {
            constexpr int n = 64;
            Real sum = 0;
            for (int i = 0; i < n; ++i) {
                const Real mu = (Real(i) + Real(0.5)) / Real(n);
                sum += fresnel::fresnelConductor(mu, eta, k) * mu;
            }
            return Real(2) * sum / Real(n);
        }

        static Real msFresnel(Real fAvg, Real eAvg) {
            return fAvg * fAvg * eAvg / (Real(1) - fAvg * (Real(1) - eAvg));
        }

        /// Interpolates E(mu) from the per-material row.
        Real eAt(Real mu) const {
            const Real x = std::clamp(mu, Real(0), Real(1)) * Real(GGXAlbedoTable::kMu - 1);
            const int i = std::min(static_cast<int>(x), GGXAlbedoTable::kMu - 2);
            const Real f = x - Real(i);
            return m_eRow[static_cast<size_t>(i)] * (Real(1) - f) + m_eRow[static_cast<size_t>(i) + 1] * f;
        }

        void initMultipleScattering() {
            const Real alphaEff = std::sqrt(alpha_x * alpha_y);
            GGXAlbedoTable::get().row(alphaEff, m_eRow, m_eAvg);
            m_msActive = m_msEnabled && (Real(1) - m_eAvg) > Real(1e-3);
            if (!m_msActive) return;
#if defined(RAYT_SPECTRAL_MODE)
            m_fmsLut.resize(etaLut.size());
            for (size_t i = 0; i < etaLut.size(); ++i)
                m_fmsLut[i] = msFresnel(averageFresnel(etaLut[i], kLut[i]), m_eAvg);
#else
            m_fms = Spectrum(msFresnel(averageFresnel(eta.x, k.x), m_eAvg),
                msFresnel(averageFresnel(eta.y, k.y), m_eAvg),
                msFresnel(averageFresnel(eta.z, k.z), m_eAvg));
#endif
        }

        /// The compensation lobe (BSDF value, cosine not included).
        Spectrum msBrdf(Real cosO, Real cosI, const SampledWavelengths& lambda) const {
            const Real scale = (Real(1) - eAt(cosO)) * (Real(1) - eAt(cosI))
                / (constants::PI * (Real(1) - m_eAvg));
#if defined(RAYT_SPECTRAL_MODE)
            Spectrum s;
            for (int i = 0; i < kSpectralSamples; ++i) s[i] = lutAt(m_fmsLut, lambda[i]) * scale;
            return s;
#else
            (void)lambda;
            return m_fms * scale;
#endif
        }

        // Convert perceptual roughness to linear alpha using an anisotropic mapping.
        void initRoughness(Real roughness, Real anisotropy) {
            Real aspect = math::safe_sqrt(1.0 - anisotropy * 0.9);
            alpha_x = MicrofacetDistribution::roughnessToAlpha(roughness / aspect);
            alpha_y = MicrofacetDistribution::roughnessToAlpha(roughness * aspect);
        }

#if defined(RAYT_SPECTRAL_MODE)
        static Real lutAt(const std::vector<Real>& lut, Real lambda_nm) {
            const Real x = std::clamp(lambda_nm, Real(kLutFirstNm), Real(kLutLastNm)) - Real(kLutFirstNm);
            const int i = std::min(static_cast<int>(x), static_cast<int>(lut.size()) - 2);
            const Real f = x - Real(i);
            return lut[static_cast<size_t>(i)] * (Real(1) - f) + lut[static_cast<size_t>(i) + 1] * f;
        }
#endif

        // Fresnel reflectance for each wavelength carried by the path.
        Spectrum fresnelAt(Real cosThetaD, const SampledWavelengths& lambda) const {
#if defined(RAYT_SPECTRAL_MODE)
            Spectrum F;
            for (int i = 0; i < kSpectralSamples; ++i) {
                F[i] = fresnel::fresnelConductor(cosThetaD, lutAt(etaLut, lambda[i]), lutAt(kLut, lambda[i]));
            }
            return F;
#else
            (void)lambda;
            return fresnel::fresnelConductor(cosThetaD, eta, k);
#endif
        }

    public:
#if !defined(RAYT_SPECTRAL_MODE)
        /**
         * @brief Constructs a rough conductor from RGB optical constants (RGB mode only).
         * @param eta Real part of IOR.
         * @param k Imaginary part of IOR.
         * @param roughness Roughness value [0, 1].
         * @param anisotropy Anisotropy factor [-1, 1] (0 for isotropic).
         */
        RoughConductor(const Spectrum& eta, const Spectrum& k, Real roughness, Real anisotropy = 0.0)
            : eta(eta), k(k) {
            m_f0Rgb = fresnel::fresnelConductor(Real(1), eta, k);
            initRoughness(roughness, anisotropy);
            initMultipleScattering();
        }
#endif

        /**
         * @brief Constructs a rough conductor from a tabulated complex IOR (eta(l), k(l)).
         *
         * Spectral mode: the Fresnel term is evaluated per sampled wavelength from the table.
         * RGB mode: the table is reduced to RGB eta/k by CMF integration (integrateIorToRgb_CMF).
         * @param table Complex IOR table (wavelengths in nm).
         * @param roughness Roughness value [0, 1].
         * @param anisotropy Anisotropy factor [-1, 1] (0 for isotropic).
         */
        RoughConductor(const std::shared_ptr<const color::ComplexIorTable>& table,
            Real roughness, Real anisotropy = 0.0) {
            {
                // RGB reduction of the table: only used for the albedo AOV in spectral mode.
                Vector3 etaRgb(1.0), kRgb(0.0);
                color::integrateIorToRgb_CMF(*table, etaRgb, kRgb);
                m_f0Rgb = fresnel::fresnelConductor(Real(1), etaRgb, kRgb);
            }
#if defined(RAYT_SPECTRAL_MODE)
            const int n = kLutLastNm - kLutFirstNm + 1;
            etaLut.resize(static_cast<size_t>(n));
            kLut.resize(static_cast<size_t>(n));
            for (int i = 0; i < n; ++i) {
                const auto [e, kk] = table->evaluateNK(Real(kLutFirstNm + i));
                etaLut[static_cast<size_t>(i)] = e;
                kLut[static_cast<size_t>(i)] = kk;
            }
#else
            Vector3 etaRgb(1.0), kRgb(0.0);
            color::integrateIorToRgb_CMF(*table, etaRgb, kRgb);
            eta = etaRgb;
            k = kRgb;
#endif
            initRoughness(roughness, anisotropy);
            initMultipleScattering();
        }

        /**
         * @brief Changes the roughness after construction (the GUI roughness slider).
         *
         * Recomputes the multiple-scattering table row. Not thread safe: call it between render
         * passes, never while rays are being traced.
         * @param roughness Roughness value [0, 1].
         * @param anisotropy Anisotropy factor [-1, 1] (0 for isotropic).
         */
        void setRoughness(Real roughness, Real anisotropy = 0.0) {
            initRoughness(roughness, anisotropy);
            initMultipleScattering();
        }

        /// Turns the multiple-scattering compensation on / off (on by default). Off = plain single-scattering GGX.
        void setMultipleScattering(bool enabled) {
            m_msEnabled = enabled;
            initMultipleScattering();
        }

        /**
         * @brief Evaluates the BRDF for a given pair of directions.
         * @param rec Surface interaction details (position, normal, etc.).
         * @param wo Outgoing direction (towards camera).
         * @param wi Incident direction (towards light).
         * @param lambda Wavelengths carried by the current path.
         * @param mode Transport mode (Radiance or Importance).
         * @return The evaluated BRDF spectrum.
         */
        Spectrum eval(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda,
            TransportMode mode) const override {

            // Ensure both directions are in the same hemisphere as the geometric normal.
            if (glm::dot(rec.gn, wo) <= 0 || glm::dot(rec.gn, wi) <= 0) return Spectrum(0.0);

            Real cosThetaO = std::max(Real(0), glm::dot(rec.n, wo));
            Real cosThetaI = std::max(Real(0), glm::dot(rec.n, wi));

            if (cosThetaI <= 0 || cosThetaO <= 0) return Spectrum(0.0);

            // 1. Calculate the half-vector and ensure it's valid.
            Vector3 wh = wi + wo;
            Real whLen2 = glm::dot(wh, wh);
            if (whLen2 <= Real(0)) return Spectrum(0.0);
            wh *= math::safe_recip(math::safe_sqrt(whLen2));

            // 2. Transform world vectors to the local shading frame.
            // The Distribution class assumes operations in Tangent Space.
            // We need to transform World Space vectors to Local Space for D() and G().
            frame::Frame frame(rec.n);
            Vector3 wo_local = frame.worldToLocal(wo);
            Vector3 wi_local = frame.worldToLocal(wi);
            Vector3 wh_local = frame.worldToLocal(wh);

            // 3. Create Distribution (On-the-fly is fine as it's lightweight)
            GGXDistribution dist(alpha_x, alpha_y);

            // 4. Evaluate Terms
            Real D = dist.D(wh_local);
            Real G = dist.G(wo_local, wi_local);

            Real cosThetaD = std::abs(glm::dot(wi, wh));
            Spectrum F = fresnelAt(cosThetaD, lambda);

            // 5. Cook-Torrance Formula
            Real denominator = 4.0 * cosThetaI * cosThetaO;
            if (denominator < 1e-6) return Spectrum(0.0);

            Spectrum f = (D * G * F) / denominator;
            if (m_msActive) f += msBrdf(cosThetaO, cosThetaI, lambda);
            return f;
        }

        /**
         * @brief Importance samples the BRDF using the Visible Normal Distribution Function (VNDF).
         * @param rec Surface interaction details.
         * @param wo Outgoing direction.
         * @param u Uniform random sample in [0, 1]^2.
         * @param lambda Wavelengths carried by the current path.
         * @return A BSDFSample containing the sampled direction, weight, and PDF, or nullopt on failure.
         */
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec,
            const Vector3& wo,
            const Point2& u,
            const SampledWavelengths& lambda,
            TransportMode mode) const override {

            // Reflection only
            if (glm::dot(rec.gn, wo) <= 0)
                return std::nullopt;

            BSDFSample bsdfSample;

            // 1. Setup Frame & Distribution
            rayt::frame::Frame frame(rec.n);
            Vector3 wo_local = frame.worldToLocal(wo);
            GGXDistribution dist(alpha_x, alpha_y);

            // Choose the lobe: GGX single scattering, or the multiple-scattering lobe (cosine-weighted)
            // with probability 1 - E(mu_o), the energy the GGX lobe cannot account for.
            const Real pMs = m_msActive ? (Real(1) - eAt(wo_local.z)) : Real(0);
            Point2 uLobe = u;
            bool useMs = false;
            if (pMs > Real(0)) {
                if (u.x < pMs) { useMs = true; uLobe.x = u.x / pMs; }
                else { uLobe.x = (u.x - pMs) / (Real(1) - pMs); }
            }

            if (useMs) {
                const Real r = math::safe_sqrt(uLobe.x);
                const Real phi = Real(2) * constants::PI * uLobe.y;
                const Vector3 wi_ms(r * std::cos(phi), r * std::sin(phi), math::safe_sqrt(Real(1) - uLobe.x));
                bsdfSample.wi = frame.localToWorld(wi_ms);
            }
            else {
                // 2. Sample a visible microfacet normal (wh) and transform to world space.
                Vector3 wh_local = dist.sample_wh(wo_local, uLobe);
                Vector3 wh = frame.localToWorld(wh_local);

                // 3. Reflect wo about wh to determine the incident direction wi.
                bsdfSample.wi = math::reflectOutward(wo, wh);
            }

            if (glm::dot(rec.gn, bsdfSample.wi) <= 0)
                return std::nullopt;

            Vector3 wi_local = frame.worldToLocal(bsdfSample.wi);

            // 4. Validate vectors against the local shading frame.
            if (wo_local.z == 0 || wi_local.z == 0) return std::nullopt;

            // 5. The PDF of the sampled direction: the mixture of both lobes (see pdf()).
            bsdfSample.pdf = pdf(rec, wo, bsdfSample.wi, lambda);

            // 6. Evaluate the sample weight (throughput).
            bsdfSample.f = eval(rec, wo, bsdfSample.wi, lambda, mode);

            if (bsdfSample.pdf <= 1e-6 || math::hasNonFinite(bsdfSample.f)) return std::nullopt;

            // Flags 
            bsdfSample.flags =
                BxDFFlags::Reflection |
                BxDFFlags::Glossy;


            return bsdfSample;
        }

        /**
         * @brief Computes the probability density function (PDF) for a given pair of directions.
         * @param rec Surface interaction details.
         * @param wo Outgoing direction.
         * @param wi Incident direction.
         * @return The PDF value with respect to solid angle.
         */
        Real pdf(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda) const override {

            if (glm::dot(rec.gn, wi) <= 0 || glm::dot(rec.gn, wo) <= 0) return 0.0;

            Vector3 wh = wo + wi;
            Real whLen2 = glm::dot(wh, wh);
            if (whLen2 <= Real(0.0)) return Real(0.0);
            wh = wh * math::safe_recip(math::safe_sqrt(whLen2));

            rayt::frame::Frame frame(rec.n);
            Vector3 wo_local = frame.worldToLocal(wo);
            Vector3 wh_local = frame.worldToLocal(wh);

            GGXDistribution dist(alpha_x, alpha_y);
            Real pdf_wh = dist.pdf(wo_local, wh_local);

            // Transform PDF from half-vector to solid angle
            Real dot_wo_wh = std::abs(glm::dot(wo_local, wh_local));
            const Real pdfGgx = (dot_wo_wh > Real(0.0)) ? pdf_wh / (Real(4) * dot_wo_wh) : Real(0.0);

            // Mixture with the multiple-scattering lobe (cosine-weighted), chosen with probability 1 - E(mu_o).
            if (!m_msActive) return pdfGgx;
            const Real pMs = Real(1) - eAt(wo_local.z);
            const Real pdfMs = std::max(glm::dot(rec.n, wi), Real(0)) / constants::PI;
            return (Real(1) - pMs) * pdfGgx + pMs * pdfMs;
        }

        /// @brief Returns false as this is a glossy material, not a perfect specular mirror.
        bool isSpecular() const override { return false; } // It is Glossy, not delta-Specular

        /// Normal-incidence reflectance (the color of the metal).
        Vector3 aovAlbedo(const SurfaceInteraction& /*rec*/) const override {
            return glm::clamp(m_f0Rgb, Vector3(0.0), Vector3(1.0));
        }

    };

} // namespace rayt