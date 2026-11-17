/**
 * @file Materials/MeasuredHpBRDF.hpp
 * @brief Measured hyperspectral-polarimetric BRDF (hpBRDF, Moon et al., SIGGRAPH Asia 2025), or any PCA-compressed
 *        BRDF table (.hpbpca v2, made by brdfkit: MERL tables are RGB and unpolarized). Header only; the decoder is Materials/hpbrdf_pca.hpp.
 *
 * What the table holds: an isotropic pBRDF as a 4x4 Mueller matrix over the Rusinkiewicz angles (phi_d, theta_d, theta_h)
 * and 68 wavelengths (414-950 nm), compressed to a 40-float latent vector per angular bin (PCA of log m00 and of the
 * normalized Mueller elements). The decoder interpolates the latent vector trilinearly and decodes it per wavelength.
 *
 * Spectral build : eval() decodes m00 at each wavelength of the path.
 * RGB build      : an RGB table (CMF-weighted m00 over the 68 bands, equal-energy white) is baked at load; eval() interpolates it.
 *
 * Polarization (docs/POLARIZATION.md): scatterMueller() returns the measured Mueller matrix divided by m00, rotated from the
 * table's reference frames (x in the plane of reflection, Mitsuba 3 convention) into the renderer's (s, p, k) bases with
 * pol::basisChange. The sign conventions were checked against pol::setInterface (Fresnel of Al) at the specular configuration:
 * m01 > 0, m22, m33 < 0, m23 > 0, m32 < 0 on both sides, so no sign flip is applied.
 *
 * Sampling: a GGX VNDF lobe (alpha = alphaSample, a user parameter) mixed with a cosine lobe (weight 0.1), as in Mitsuba's
 * measured_polarized plugin. Any alpha gives a correct image; only the noise changes.
 *
 * Limits: wavelengths outside 414-950 nm are clamped to the nearest band. About 25 % of the angular bins (grazing angles) were
 * not measured; they hold the nearest measured bin. m00 near the specular peak can be 10-18 % below the raw table (the latent
 * interpolation is a log-space interpolation).
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "Color/ColorSpace.hpp"
#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Core/Interaction.hpp"
#include "Core/Math.hpp"
#include "Core/Polarization.hpp"
#include "Geometry/Frame.hpp"
#include "Materials/Material.hpp"
#include "Microfacet/GGX.hpp"

#include "Materials/hpbrdf_pca.hpp"

namespace rayt {

    class MeasuredHpBRDF : public Material {
    public:
        static constexpr int kMaxCoef = 64;              // latent size the stack buffers allow (the tables use 40)
        static constexpr Real kCosineWeight = Real(0.1); // share of the cosine lobe in the sampling mixture

        MeasuredHpBRDF() = default;

        /**
         * @brief Loads a .hpbpca table.
         * @param path        the table (relative to the working directory, like every path of a scene file)
         * @param alphaSample GGX alpha of the sampling lobe; set it near the apparent roughness of the material
         * @param err         receives the reason when loading fails
         */
        static std::shared_ptr<MeasuredHpBRDF> load(const std::string& path, Real alphaSample = Real(0.1), std::string* err = nullptr) {
            auto m = std::make_shared<MeasuredHpBRDF>();
            if (!m->m_table.load(path, err)) return nullptr;
            if (m->m_table.k() > uint32_t(kMaxCoef)) { if (err) *err = "latent vector longer than kMaxCoef"; return nullptr; }
            m->m_alpha = std::max(alphaSample, Real(1e-3));
            m->m_path = path;
            m->init();
            return m;
        }

        const hpbrdf::PCATable& table() const { return m_table; }
        const std::string& path() const { return m_path; }

        // ------------------------------------------------------------------ Material
        Spectrum eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda, TransportMode mode = TransportMode::Radiance) const override {
            if (glm::dot(rec.gn, wo) <= 0 || glm::dot(rec.gn, wi) <= 0) return Spectrum(0.0);
            const frame::Frame f(rec.n);
            const Vector3 woL = f.worldToLocal(wo), wiL = f.worldToLocal(wi);
            if (woL.z <= 0 || wiL.z <= 0) return Spectrum(0.0);
            // light arrives from wi and leaves toward wo (radiance); the importance walk sees the reverse.
            const bool radiance = mode == TransportMode::Radiance;
            const hpbrdf::Geometry g = hpbrdf::geometry(toVec(radiance ? wiL : woL), toVec(radiance ? woL : wiL));

            Spectrum out(0.0);
#if defined(RAYT_SPECTRAL_MODE)
            float c[kMaxCoef];
            m_table.coeffs(g.phi_d, g.theta_d, g.theta_h, c);
            for (int i = 0; i < kSpectralSamples; ++i) out[i] = std::max(Real(0), Real(m_table.m00(c, channelOf(float(lambda[i])))));
#else
            (void)lambda;
            const std::array<float, 3> rgb = rgbAt(g.phi_d, g.theta_d, g.theta_h);
            for (int i = 0; i < 3; ++i) out[i] = std::max(Real(0), Real(rgb[size_t(i)]));
#endif
            return out;
        }

        std::optional<BSDFSample> sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u,
            const SampledWavelengths& lambda, TransportMode mode = TransportMode::Radiance) const override {
            if (glm::dot(rec.gn, wo) <= 0) return std::nullopt;
            const frame::Frame f(rec.n);
            const Vector3 woL = f.worldToLocal(wo);
            if (woL.z <= 0) return std::nullopt;

            BSDFSample s;
            if (u.x < kCosineWeight) {                               // cosine lobe
                const Real ux = u.x / kCosineWeight;
                const Real r = math::safe_sqrt(ux), phi = Real(2) * constants::PI * u.y;
                s.wi = f.localToWorld(Vector3(r * std::cos(phi), r * std::sin(phi), math::safe_sqrt(Real(1) - ux)));
            }
            else {                                                   // GGX visible normals
                const Point2 ug((u.x - kCosineWeight) / (Real(1) - kCosineWeight), u.y);
                GGXDistribution dist(m_alpha, m_alpha);
                const Vector3 wh = f.localToWorld(dist.sample_wh(woL, ug));
                s.wi = math::reflectOutward(wo, wh);
            }
            if (glm::dot(rec.gn, s.wi) <= 0) return std::nullopt;
            s.pdf = pdf(rec, wo, s.wi, lambda);
            if (!(s.pdf > Real(1e-7))) return std::nullopt;
            s.f = eval(rec, wo, s.wi, lambda, mode);
            if (math::hasNonFinite(s.f)) return std::nullopt;
            s.flags = BxDFFlags::Reflection | BxDFFlags::Glossy;
            return s;
        }

        Real pdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& /*lambda*/) const override {
            if (glm::dot(rec.gn, wo) <= 0 || glm::dot(rec.gn, wi) <= 0) return Real(0);
            const frame::Frame f(rec.n);
            const Vector3 woL = f.worldToLocal(wo), wiL = f.worldToLocal(wi);
            if (woL.z <= 0 || wiL.z <= 0) return Real(0);
            Vector3 whL = woL + wiL;
            const Real len2 = glm::dot(whL, whL);
            if (len2 <= Real(0)) return Real(0);
            whL *= math::safe_recip(math::safe_sqrt(len2));
            GGXDistribution dist(m_alpha, m_alpha);
            const Real dotWoWh = std::abs(glm::dot(woL, whL));
            const Real pdfGgx = dotWoWh > Real(0) ? dist.pdf(woL, whL) / (Real(4) * dotWoWh) : Real(0);
            const Real pdfCos = wiL.z / constants::PI;
            return kCosineWeight * pdfCos + (Real(1) - kCosineWeight) * pdfGgx;
        }

        /// The measured Mueller matrix (divided by m00) of the scattering wi -> wo, in the renderer's s bases (docs/POLARIZATION.md).
        bool scatterMueller(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda, pol::ScatterEvent& out) const override {
            if (!m_table.polarized) return false;                       // an unpolarized table (MERL): depolarizing bounce
            if (glm::dot(rec.gn, wo) <= 0 || glm::dot(rec.gn, wi) <= 0) return false;
            const frame::Frame f(rec.n);
            const Vector3 woL = f.worldToLocal(wo), wiL = f.worldToLocal(wi);
            if (woL.z <= 0 || wiL.z <= 0) return false;
            Vector3 wh = wo + wi;
            const Real len2 = glm::dot(wh, wh);
            if (len2 <= Real(0)) return false;
            wh *= math::safe_recip(math::safe_sqrt(len2));

            // table geometry; its Stokes reference vectors x_in (beam along -wi) and x_out (beam along wo) lie in the plane of reflection
            const hpbrdf::Geometry g = hpbrdf::geometry(toVec(wiL), toVec(woL));
            const Vector3 xIn = f.localToWorld(toVector3(g.x_in)), xOut = f.localToWorld(toVector3(g.x_out));
            out.s = pol::planeNormalS(wo, wh);

            float c[kMaxCoef];
            m_table.coeffs(g.phi_d, g.theta_d, g.theta_h, c);
            pol::Mueller Mt;   // the table's matrix / m00, per component
            for (int i = 0; i < pol::kN; ++i) {
#if defined(RAYT_SPECTRAL_MODE)
                const float lam = channelOf(float(lambda[i]));
#else
                const float lam = m_table.rgb ? float(i) : kRgbLambda[size_t(i)];
#endif
                const hpbrdf::Mueller M = m_table.mueller(c, lam);
                const float m00 = M.m[0][0];
                if (!(m00 > 0.f)) { Mt.e[0][0][i] = Real(1); continue; }   // nothing measured: a depolarizer
                for (int r = 0; r < 4; ++r)
                    for (int cc = 0; cc < 4; ++cc) Mt.e[r][cc][i] = Real(M.m[r][cc] / m00);
            }
            (void)lambda;
            // Stokes_out(s) = R(x_out -> s | wo) . Mt . R(s -> x_in | -wi) . Stokes_in(s)
            out.m = pol::basisChange(xOut, out.s, wo) * Mt * pol::basisChange(out.s, xIn, -wi);
            return true;
        }

        bool isSpecular() const override { return false; }

        /// Directional albedo at normal incidence, in RGB (denoiser AOV); computed once at load.
        Vector3 aovAlbedo(const SurfaceInteraction& /*rec*/) const override { return m_albedoRgb; }

    private:
        hpbrdf::PCATable m_table;
        std::string m_path;
        Real m_alpha = Real(0.1);
        Vector3 m_albedoRgb = Vector3(1.0);
        std::vector<float> m_rgbW;                 // [band][3]: CMF-based weights of the 68 bands for R, G, B (each channel sums to 1)
        static constexpr std::array<float, 3> kRgbLambda = { 610.f, 550.f, 465.f };   // wavelengths standing for R, G, B in the RGB build's Mueller matrix
#if !defined(RAYT_SPECTRAL_MODE)
        std::vector<float> m_rgb;                  // baked [phi][thd][thh][3]
#endif

        /// channel coordinate for the decoder: the wavelength itself, or the nearest RGB channel of an RGB table (spectral build fallback)
        float channelOf(float lambda_nm) const {
            if (!m_table.rgb) return lambda_nm;
            return lambda_nm < 490.f ? 2.f : (lambda_nm < 580.f ? 1.f : 0.f);
        }
        static hpbrdf::Vec3 toVec(const Vector3& v) { return hpbrdf::Vec3(float(v.x), float(v.y), float(v.z)); }
        static Vector3 toVector3(const hpbrdf::Vec3& v) { return Vector3(Real(v.x), Real(v.y), Real(v.z)); }

        // CIE 1931 2-degree color matching functions (the 1 nm table of the film)
        static void cmf(double lam, double& x, double& y, double& z) {
            const color::CIE_XYZ c = color::evaluateCIE_XYZ_1931(Real(lam));
            x = double(c.X); y = double(c.Y); z = double(c.Z);
        }

        /// m00 over the table's bands -> linear sRGB with the band weights
        std::array<float, 3> rgbOfSpectrum(const float* c) const {
            std::array<float, 3> rgb = { 0.f, 0.f, 0.f };
            if (m_table.rgb) { for (int ch = 0; ch < 3; ++ch) rgb[size_t(ch)] = m_table.m00_band(c, uint32_t(ch)); return rgb; }
            for (uint32_t w = 0; w < m_table.n_ch; ++w) {
                const float v = m_table.m00_band(c, w);
                for (int ch = 0; ch < 3; ++ch) rgb[size_t(ch)] += m_rgbW[size_t(w) * 3 + size_t(ch)] * v;
            }
            return rgb;
        }

        void init() {
            // band weights: XYZ -> linear sRGB (D65), normalized per channel so that a flat spectrum is white (wavelength tables only)
            if (!m_table.rgb) {
            const uint32_t n = m_table.n_ch;
            m_rgbW.assign(size_t(n) * 3, 0.f);
            double sum[3] = { 0, 0, 0 };
            for (uint32_t w = 0; w < n; ++w) {
                double x, y, z; cmf(double(m_table.channels[w]), x, y, z);
                const double r = 3.2406 * x - 1.5372 * y - 0.4986 * z;
                const double g = -0.9689 * x + 1.8758 * y + 0.0415 * z;
                const double b = 0.0557 * x - 0.2040 * y + 1.0570 * z;
                m_rgbW[size_t(w) * 3 + 0] = float(r); m_rgbW[size_t(w) * 3 + 1] = float(g); m_rgbW[size_t(w) * 3 + 2] = float(b);
                sum[0] += r; sum[1] += g; sum[2] += b;
            }
            for (uint32_t w = 0; w < n; ++w)
                for (int ch = 0; ch < 3; ++ch) m_rgbW[size_t(w) * 3 + size_t(ch)] = float(m_rgbW[size_t(w) * 3 + size_t(ch)] / sum[ch]);

            }
#if !defined(RAYT_SPECTRAL_MODE)
            // bake the RGB table: one exp per band per bin (about 3 M bins; a few seconds, parallel when OpenMP is on)
            const long long nb = (long long) m_table.n_phi * m_table.n_thd * m_table.n_thh;
            m_rgb.assign(size_t(nb) * 3, 0.f);
            const uint32_t k = m_table.k();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (long long b = 0; b < nb; ++b) {
                const std::array<float, 3> rgb = rgbOfSpectrum(&m_table.coef[size_t(b) * k]);
                for (int ch = 0; ch < 3; ++ch) m_rgb[size_t(b) * 3 + size_t(ch)] = rgb[size_t(ch)];
            }
#endif
            // directional albedo at normal incidence: integral of m00 cos over the hemisphere (coarse quadrature), in RGB
            {
                double acc[3] = { 0, 0, 0 };
                const int nt = 24, np = 16;
                float c[kMaxCoef];
                for (int it = 0; it < nt; ++it) {
                    const double th = (it + 0.5) / nt * (constants::PI / 2), ct = std::cos(th), st = std::sin(th);
                    const double dw = st * (constants::PI / 2 / nt) * (2 * constants::PI / np);
                    for (int ip = 0; ip < np; ++ip) {
                        const double ph = (ip + 0.5) / np * 2 * constants::PI;
                        const hpbrdf::Geometry g = hpbrdf::geometry(hpbrdf::Vec3(float(st * std::cos(ph)), float(st * std::sin(ph)), float(ct)),
                                                                    hpbrdf::Vec3(0.f, 0.f, 1.f));
                        m_table.coeffs(g.phi_d, g.theta_d, g.theta_h, c);
                        const std::array<float, 3> rgb = rgbOfSpectrum(c);
                        for (int ch = 0; ch < 3; ++ch) acc[ch] += double(rgb[size_t(ch)]) * ct * dw;
                    }
                }
                m_albedoRgb = glm::clamp(Vector3(Real(acc[0]), Real(acc[1]), Real(acc[2])), Vector3(0.0), Vector3(1.0));
            }
        }

#if !defined(RAYT_SPECTRAL_MODE)
        std::array<float, 3> rgbAt(float phi_d, float theta_d, float theta_h) const {
            uint32_t ip, id, ih; float tp, td, th;
            hpbrdf::PCATable::locate(m_table.phi_d, phi_d, ip, tp);
            hpbrdf::PCATable::locate(m_table.theta_d, theta_d, id, td);
            hpbrdf::PCATable::locate(m_table.theta_h, theta_h, ih, th);
            std::array<float, 3> out = { 0.f, 0.f, 0.f };
            for (int dp = 0; dp < 2; ++dp) for (int dd = 0; dd < 2; ++dd) for (int dh = 0; dh < 2; ++dh) {
                const float w = (dp ? tp : 1 - tp) * (dd ? td : 1 - td) * (dh ? th : 1 - th);
                if (w == 0.f) continue;
                const size_t b = ((size_t(ip + dp) * m_table.n_thd + (id + dd)) * m_table.n_thh + (ih + dh)) * 3;
                for (int ch = 0; ch < 3; ++ch) out[size_t(ch)] += w * m_rgb[b + size_t(ch)];
            }
            return out;
        }
#endif
    };

} // namespace rayt
