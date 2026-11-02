/**
 * @file Materials/Subsurface.hpp
 * @brief Subsurface scattering (SSS): the radial profile, the profile shared by a group of
 *        surfaces, and the pseudo material of the point where the light leaves.
 *
 * Model (BSSRDF, as in PBRT v3 / Christensen & Burley 2015):
 *     S(po, wo, pi, wi) = (1 - F(wo)) * Sp(|po - pi|) * Sw(wi)
 *  - Sp(r) is the radial profile: how much of the light that enters at po leaves at distance r
 *    along the surface. Burley's normalized diffusion, per color channel c:
 *        Sp_c(r) = A_c (exp(-r / d_c) + exp(-r / 3 d_c)) / (8 pi d_c r)
 *    which integrates to A_c over the surface (A = the diffuse color). The length scale is
 *    d_c = l_c / s(A_c) (l = scatter distance, s the shape function below).
 *  - Sw is a cosine-like lobe at the exit point with the Fresnel transmission, normalized so that
 *    it integrates to 1 (the SubsurfaceExit material below).
 *  - The entry factor (1 - F(wo)) is the probability of passing the coat (CoatedDiffuse picks the
 *    specular lobe with the complementary probability).
 *
 * Exit points are found by probe rays (Render/SubsurfaceProbe.hpp): only surfaces whose material
 * returns the SAME SubsurfaceProfile pointer take part, so one profile object is shared by
 * every skin material (face, neck, body) and clothes are never exit points.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <memory>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Core/Math.hpp"
#include "Core/Fresnel.hpp"
#include "Core/Sampling.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Geometry/Frame.hpp"
#include "Materials/Material.hpp"
#include "Color/SpectralConversion.hpp"

namespace rayt {

    /// Number of color channels (RGB) or wavelengths (spectral) the profile is evaluated for.
#if defined(RAYT_SPECTRAL_MODE)
    constexpr int kSssChannels = kSpectralSamples;
#else
    constexpr int kSssChannels = 3;
#endif

    namespace sss {

        /// Burley's fit: the profile length scale is d = scatterDistance / shape(albedo).
        inline Real shape(Real albedo) {
            const Real a = std::clamp(albedo, Real(0), Real(1));
            const Real e = std::abs(a - Real(0.8));
            return Real(1.85) - a + Real(7) * e * e * e;
        }

        /// Profile value per unit area (integrates to 1 over the plane); the area density of exit points.
        inline Real profile(Real r, Real d) {
            r = std::max(r, Real(1e-7));
            return (std::exp(-r / d) + std::exp(-r / (Real(3) * d))) / (Real(8) * constants::PI * d * r);
        }

        /// Fraction of the profile inside radius r.
        inline Real cdf(Real r, Real d) {
            const Real x = r / d;
            return Real(1) - Real(0.25) * std::exp(-x) - Real(0.75) * std::exp(-x / Real(3));
        }

        /// Inverse of cdf(): radius for the uniform number u in [0, 1).
        inline Real sampleRadius(Real u, Real d) {
            // Solve 1 - 0.25 e^-x - 0.75 e^(-x/3) = u for x = r / d (Newton, kept inside a bracket).
            Real lo = 0, hi = Real(80);
            Real x = std::clamp(Real(3) * std::log(Real(0.75) / std::max(Real(1) - u, Real(1e-12))), Real(0.01), hi);
            for (int i = 0; i < 24; ++i) {
                const Real f = Real(1) - Real(0.25) * std::exp(-x) - Real(0.75) * std::exp(-x / Real(3)) - u;
                (f > 0 ? hi : lo) = x;
                const Real fp = Real(0.25) * std::exp(-x) + Real(0.25) * std::exp(-x / Real(3));
                Real xn = x - f / fp;
                if (std::abs(xn - x) < Real(1e-12) * (Real(1) + x)) { x = xn; break; }   // converged (before the bracket test: a tiny step can land on the bracket's edge)
                if (!(xn > lo && xn < hi)) xn = Real(0.5) * (lo + hi);
                x = xn;
            }
            return x * d;
        }

    } // namespace sss

    /**
     * @brief Profile shared by a group of skin-like surfaces: the scatter distance of each color
     *        channel (meters) and the index of refraction of the medium.
     */
    class SubsurfaceProfile {
    public:
        /// @param scatterDistance mean free path per RGB channel in meters (skin: about 3.7 / 1.4 / 0.7 mm)
        SubsurfaceProfile(const Vector3& scatterDistance, Real eta = Real(1.4))
            : m_distance(scatterDistance), m_eta(eta) {}

        Real eta() const { return m_eta; }
        const Vector3& scatterDistanceRgb() const { return m_distance; }

        /// Scatter distance per channel (RGB mode) or per sampled wavelength (the RGB triple is upsampled).
        Spectrum scatterDistance(const SampledWavelengths& lambda) const {
#if defined(RAYT_SPECTRAL_MODE)
            const Real m = std::max({ m_distance.x, m_distance.y, m_distance.z, Real(1e-9) });
            return RgbColor(m_distance / m).reflectance(lambda) * m;
#else
            (void)lambda;
            return Spectrum(m_distance);
#endif
        }

    private:
        Vector3 m_distance;
        Real m_eta;
    };

    /**
     * @brief The exit vertex of a subsurface path: a cosine lobe with the Fresnel transmission.
     *
     *     Sw(wi) = (1 - w F(cos_i)) / (c pi),   c = 1 - 2 w Integral F(mu) mu dmu
     *
     * so that the integral of Sw cos over the hemisphere is exactly 1 (w = weight of the coat).
     * It does not depend on the outgoing direction, which is why the path can simply continue
     * from the exit point as if it had arrived there.
     */
    class SubsurfaceExit : public Material {
    public:
        SubsurfaceExit(Real coatWeight, Real eta) : m_w(coatWeight), m_eta(eta) {
            constexpr int n = 128;
            Real sum = 0;
            for (int i = 0; i < n; ++i) {
                const Real mu = (Real(i) + Real(0.5)) / Real(n);
                sum += fresnel::fresnelDielectric(mu, Real(1), m_eta) * mu;
            }
            m_c = Real(1) - Real(2) * m_w * sum / Real(n);
        }

        Spectrum eval(const SurfaceInteraction& rec, const Vector3& /*wo*/, const Vector3& wi,
            const SampledWavelengths& /*lambda*/, TransportMode /*mode*/) const override
        {
            if (glm::dot(rec.gn, wi) <= 0) return Spectrum(0.0);
            const Real cosI = glm::dot(rec.n, wi);
            if (cosI <= 0) return Spectrum(0.0);
            return Spectrum(weight(cosI) / constants::PI);
        }

        Real pdf(const SurfaceInteraction& rec, const Vector3& /*wo*/, const Vector3& wi,
            const SampledWavelengths& /*lambda*/) const override
        {
            if (glm::dot(rec.gn, wi) <= 0) return Real(0);
            return std::max(glm::dot(rec.n, wi), Real(0)) / constants::PI;
        }

        std::optional<BSDFSample> sample(const SurfaceInteraction& rec, const Vector3& /*wo*/, const Point2& u,
            const SampledWavelengths& lambda, TransportMode mode) const override
        {
            frame::Frame frame(rec.n);
            BSDFSample s;
            s.wi = frame.localToWorld(sampling::CosineSampleHemisphere(u));
            if (glm::dot(rec.gn, s.wi) <= 0) return std::nullopt;
            s.pdf = pdf(rec, Vector3(0.0), s.wi, lambda);
            s.f = eval(rec, Vector3(0.0), s.wi, lambda, mode);
            if (s.pdf <= Real(1e-8)) return std::nullopt;
            s.flags = BxDFFlags::Diffuse | BxDFFlags::Reflection;
            return s;
        }

    private:
        /// (1 - w F(cos_i)) / c
        Real weight(Real cosI) const {
            return (Real(1) - m_w * fresnel::fresnelDielectric(cosI, Real(1), m_eta)) / m_c;
        }

        Real m_w, m_eta, m_c = 1;
    };

} // namespace rayt
