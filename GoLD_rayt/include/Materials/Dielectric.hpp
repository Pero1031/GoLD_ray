/**
 * @file Materials/Dielectric.hpp
 * @brief Dielectric material (Glass, Water, Diamond).
 * * Supports both smooth and rough surfaces using Microfacet theory (Walter et al. 2007).
 * * Ported and adapted from PBRT v4 DielectricBxDF.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <string>

#include "Core/Types.hpp"
#include "Color/Sellmeier.hpp"
#include "Core/Fresnel.hpp"
#include "Core/Math.hpp"
#include "Core/Constants.hpp"
#include "Core/Forward.hpp"

#include "Geometry/Frame.hpp"
#include "Materials/Material.hpp"
#include "Materials/Medium.hpp"
#include "Microfacet/GGX.hpp"

namespace rayt {

    /**
     * @brief Dielectric (glass-like) material, smooth or rough, optionally dispersive.
     *
     * Dispersion (spectral mode only): the index of refraction depends on the wavelength via the
     * Cauchy relation  n(l) = A + B / l^2  (l in micrometres), built from the refractive index at
     * the d line (587.6 nm) and the Abbe number V_d = (n_d - 1) / (n_F - n_C).
     * A path that meets a dispersive surface keeps only its first ("hero") wavelength, because
     * different wavelengths take different directions (see prepareWavelengths()).
     * In RGB mode, or when no Abbe number is given, the IOR is a constant.
     */
    class Dielectric : public Material {
        Real iorRef;       // IOR at the d line (587.6 nm); the constant IOR when not dispersive
        bool m_hasSellmeier = false;
        Sellmeier m_sell;
        Real cauchyA = 0;  // Cauchy coefficient A (dimensionless)
        Real cauchyB = 0;  // Cauchy coefficient B [um^2]; 0 means "not dispersive"
        Real alpha_x;      // Roughness X
        Real alpha_y;      // Roughness Y
        std::shared_ptr<const Medium> m_medium;   // medium inside the surface (null = a plain glass boundary)
        bool m_hasBoundary = false;               // setBoundaryMedia(): inside and outside media given explicitly
        const Medium* m_inside = nullptr;
        const Medium* m_outside = nullptr;

        Real m_exteriorIor = Real(1);             // index of the medium on the outside (normal side) of the surface: 1 = air, 1.333 = water

        // The RELATIVE index of refraction (inside / outside) to use for this path: the hero wavelength's IOR in spectral mode, divided by the
        // exterior index. Every formula of the BSDF (Snell, Fresnel, the radiance scaling of refraction) depends on the ratio only.
        Real iorAt(const SampledWavelengths& lambda) const {
#if defined(RAYT_SPECTRAL_MODE)
            if (m_hasSellmeier) return m_sell.n(lambda[0] * Real(1e-3)) / m_exteriorIor;
            if (cauchyB != Real(0)) {
                const Real um = lambda[0] * Real(1e-3);   // nm -> um
                return (cauchyA + cauchyB / (um * um)) / m_exteriorIor;
            }
#else
            (void)lambda;
#endif
            return iorRef / m_exteriorIor;
        }

    public:
        /**
         * @param ior         IOR at the d line (587.6 nm).
         * @param roughness   Perceptual roughness (0 = perfectly smooth).
         * @param anisotropy  Anisotropy factor.
         * @param abbeNumber  Abbe number V_d (typical: crown glass ~60, flint ~30, dense flint ~25,
         *                    diamond ~55). Smaller = stronger dispersion. 0 = no dispersion.
         */
        /// Dispersion from Sellmeier coefficients (instead of an Abbe number). The IOR of the RGB build, and of non-dispersive uses, becomes n(587.6 nm).
        void setSellmeier(const Sellmeier& s) { m_sell = s; m_hasSellmeier = true; iorRef = s.n(Real(0.5876)); }
        const Sellmeier* sellmeier() const { return m_hasSellmeier ? &m_sell : nullptr; }

        /**
         * @brief The index of refraction on the OUTSIDE (the side the normal points to), default 1 (air).
         *
         * Nested dielectrics the Mitsuba way: a surface between two media carries both indices, and the mesh is split where the
         * outside changes. A drop of water on glass is two surfaces: the cap (water 1.333 inside, air outside) and the contact disc
         * (water inside, glass 1.5 outside; the glass has no surface of its own under the drop). The BSDF uses the ratio inside /
         * outside; dispersion applies to the inside only. The medium pointers (setMedium / setBoundaryMedia) are independent of this.
         * The specular manifold NEE (Render/SpecularNee.hpp) still assumes air outside its bodies.
         */
        void setExteriorIor(Real n) { m_exteriorIor = std::max(n, Real(1e-3)); }
        Real exteriorIor() const { return m_exteriorIor; }

        /// Fills the inside of this closed surface with a participating medium (volumetric random walk).
        /// The mesh must be closed with outward normals; the outside is vacuum.
        void setMedium(std::shared_ptr<const Medium> medium) { m_medium = std::move(medium); }
        const Medium* mediumInside() const override { return m_hasBoundary ? m_inside : m_medium.get(); }

        /// For a closed surface that stands inside another medium (a prism in a fog): the medium that fills it (null = vacuum) and the one around it.
        /// The media must outlive the material. Replaces setMedium().
        void setBoundaryMedia(const Medium* inside, const Medium* outside) { m_hasBoundary = true; m_inside = inside; m_outside = outside; }
        bool isMediumBoundary() const override { return m_hasBoundary || m_medium != nullptr; }
        const Medium* mediumOutside() const override { return m_hasBoundary ? m_outside : nullptr; }

        Dielectric(Real ior, Real roughness, Real anisotropy = 0.0, Real abbeNumber = 0.0) : iorRef(ior) {
            Real aspect = std::sqrt(1.0 - anisotropy * 0.9);
            alpha_x = MicrofacetDistribution::roughnessToAlpha(roughness / aspect);
            alpha_y = MicrofacetDistribution::roughnessToAlpha(roughness * aspect);

            if (abbeNumber > Real(0) && ior > Real(1)) {
                // Fraunhofer lines in micrometres: d (587.6 nm), F (486.1 nm), C (656.3 nm).
                constexpr Real lambdaD = Real(0.5876), lambdaF = Real(0.4861), lambdaC = Real(0.6563);
                // n_F - n_C = (n_d - 1) / V_d = B * (1/lambda_F^2 - 1/lambda_C^2)
                cauchyB = (ior - Real(1)) /
                    (abbeNumber * (Real(1) / (lambdaF * lambdaF) - Real(1) / (lambdaC * lambdaC)));
                cauchyA = ior - cauchyB / (lambdaD * lambdaD);
            }
        }

        /// Index of refraction at a wavelength in nm (the constant IOR when not dispersive or in RGB mode).
        Real iorAtNm(Real nm) const {
#if defined(RAYT_SPECTRAL_MODE)
            if (m_hasSellmeier) return m_sell.n(nm * Real(1e-3));
            if (cauchyB != Real(0)) {
                const Real um = nm * Real(1e-3);
                return cauchyA + cauchyB / (um * um);
            }
#else
            (void)nm;
#endif
            return iorRef;
        }

        /// True if the IOR depends on the wavelength (spectral mode with an Abbe number).
        bool isDispersive() const {
#if defined(RAYT_SPECTRAL_MODE)
            return m_hasSellmeier || cauchyB != Real(0);
#else
            return false;
#endif
        }

        /// Dispersive surfaces keep a single wavelength per path (no-op otherwise / in RGB mode).
        void prepareWavelengths(SampledWavelengths& lambda) const override {
#if defined(RAYT_SPECTRAL_MODE)
            if (isDispersive()) lambda.terminateSecondary();
#else
            (void)lambda;
#endif
        }

        // ---------------------------------------------------------------------
        // Evaluation (BSDF)
        // ---------------------------------------------------------------------
        Spectrum eval(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda,
            TransportMode mode) const override;

        // ---------------------------------------------------------------------
        // Sampling
        // ---------------------------------------------------------------------
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec,
            const Vector3& wo,
            const Point2& u,
            const SampledWavelengths& lambda,
            TransportMode mode) const override;

        // ---------------------------------------------------------------------
        // PDF
        // ---------------------------------------------------------------------
        Real pdf(const SurfaceInteraction& rec,
            const Vector3& wo, const Vector3& wi,
            const SampledWavelengths& lambda) const override;

        bool isSpecular() const override { return isSmooth() && !m_hasFar; }
        /// A smooth surface with the same index on both sides (1 against air, or 1.333 against water) bends and reflects nothing.
        bool isNullInterface() const override { return iorRef == m_exteriorIor && isSmooth() && !m_hasFar && !isDispersive(); }

        /// The roughness grows with the distance from the y axis, from this material's own roughness up to `farRoughness` (perceptual, as in the
        /// constructor) between the radii `start` and `end` (metres, in the x-z plane). For a sea: the mesh cannot resolve the small waves far away,
        /// and the slope they would add is given back as roughness (the sun glitter keeps its width toward the horizon).
        void setDistanceRoughness(Real farRoughness, Real start, Real end) {
            m_farAlpha = MicrofacetDistribution::roughnessToAlpha(farRoughness);
            m_farStart = start;
            m_farEnd = std::max(end, start + Real(1e-3));
            m_hasFar = true;
        }

    private:
        bool isSmooth() const { return alpha_x < 0.001 && alpha_y < 0.001; }

        bool m_hasFar = false;
        Real m_farAlpha = Real(0), m_farStart = Real(0), m_farEnd = Real(1);

        /// The roughness at this hit (the material's, or the distance-dependent one).
        void alphasAt(const SurfaceInteraction& rec, Real& ax, Real& ay) const {
            ax = alpha_x; ay = alpha_y;
            if (!m_hasFar) return;
            const Real r = std::hypot(rec.p.x, rec.p.z);
            Real t = std::clamp((r - m_farStart) / (m_farEnd - m_farStart), Real(0), Real(1));
            t = t * t * (Real(3) - Real(2) * t);
            ax += (std::max(ax, m_farAlpha) - ax) * t;
            ay += (std::max(ay, m_farAlpha) - ay) * t;
        }
        bool smoothAt(const SurfaceInteraction& rec) const {
            Real ax, ay;
            alphasAt(rec, ax, ay);
            return ax < 0.001 && ay < 0.001;
        }
    };

} // namespace rayt
