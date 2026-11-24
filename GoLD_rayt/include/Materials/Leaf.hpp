/**
 * @file Materials/Leaf.hpp
 * @brief A thin leaf whose diffuse reflection and transmission come from the PROSPECT model (docs/LEAF.md).
 *
 * After Kruk, Noizet and Lucas, "A Biophysical Approach to Foliage Rendering: Time-Dependent BSDF for Aging Leaves" (EGSR 2026),
 * with one change: their BSDF adds a GGX lobe for the cuticle to PROSPECT's reflectance, which already holds the reflection of the
 * first face (1 - tav(alpha, n)), so that light was counted twice. Here the diffuse lobes carry only the light that entered the leaf.
 *
 * On each side (the geometric normal decides which side the light is on):
 *  - specular      F(theta_d, n) D G / (4 mu_o mu_i), GGX, the plates' index n(lambda); the upper and the lower face have their own roughness
 *  - reflection    Rd (1 - E(mu_o)) (1 - E(mu_i)) / (pi (1 - E_avg)^2)       Rd = R - (1 - tav(90, n)): PROSPECT's reflectance less the first face
 *  - transmission  T  (1 - E_o(mu_o)) (1 - E_i(mu_i)) / (pi (1 - E_o,avg) (1 - E_i,avg))   (towards the other side, whose face may be rougher)
 * with R, T from PROSPECT for isotropic light (alpha = 90 degrees) and E(mu) the directional albedo of a face's GGX lobe (tabulated at the
 * index of 550 nm), E_avg its cosine average. The light the GGX lobe does not reflect goes into the leaf, as in the paper's (1 - Es(wi)); the
 * second factor and the normalisation make the lobes reciprocal and make them reflect exactly Rd and transmit exactly T under light that is
 * uniform over all directions. The leaf as a whole then reflects Rd + E_avg instead of PROSPECT's Rd + (1 - tav(90, n)): the GGX lobe loses
 * a little light at high roughness. Since R + T <= 1 and E_avg <= 1 - tav(90, n), Rd + T <= 1 - E_avg: no direction gets more than it gives.
 *
 * PROSPECT itself is symmetric: both faces have the same diffuse reflectance. Real leaves are lighter underneath (the spongy mesophyll);
 * not modelled.
 */
#pragma once
#include <memory>
#include <optional>
#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Interaction.hpp"
#include "Materials/Material.hpp"
#include "Materials/Prospect.hpp"
#include "Core/Texture.hpp"
#include <array>
#include <vector>

namespace rayt {

    class Leaf : public Material {
    public:
        /// `roughnessTop` / `roughnessBottom`: perceptual GGX roughness of the upper (front-face) and the lower (back-face) side, alpha = roughness^2.
        /// `pigmentMap` (optional): the contents per texel, an RGBA image read as numbers (tools/data/leaf_scan.py pigments): chlorophyll / 80,
        /// carotenoids / 20, anthocyanins / 40 (ug/cm^2), brown pigments (0..1); N, water and dry matter stay those of `pigments`.
        /// `scale` multiplies chlorophyll, carotenoids, anthocyanins and brown pigments (of the map or of `pigments`).
        Leaf(std::shared_ptr<const leaf::ProspectTable> table, const leaf::Pigments& pigments, Real roughnessTop = Real(0.55),
            Real roughnessBottom = Real(0.75), std::shared_ptr<const Texture2D> pigmentMap = nullptr,
            const std::array<Real, 4>& scale = { 1, 1, 1, 1 });

        Spectrum eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const SampledWavelengths& lambda,
            TransportMode mode = TransportMode::Radiance) const override;
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u, const SampledWavelengths& lambda,
            TransportMode mode = TransportMode::Radiance) const override;
        Real pdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const SampledWavelengths& lambda) const override;
        Vector3 aovAlbedo(const SurfaceInteraction& rec) const override;
        /// eval() without the GGX lobe: the diffuse reflection or transmission alone (for the checks).
        Spectrum evalDiffuse(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const SampledWavelengths& lambda) const;

        const leaf::Pigments& pigments() const { return m_p; }
        /// The contents at a hit (the map's, or the constant ones), with the scale applied.
        leaf::Pigments pigmentsAt(const SurfaceInteraction& rec) const;
        /// Linear sRGB of the leaf's reflectance R and transmittance T under D65 (PROSPECT, isotropic light).
        Vector3 rgbReflectance() const { return m_rgbR; }
        Vector3 rgbTransmittance() const { return m_rgbT; }

        /// What the lobes need at the path's wavelengths, per channel: Rd, T, t12 = tav(90, n) and the plates' index n.
        struct Optics {
            Spectrum Rd, T, t12, n;
        };
        Optics opticsAt(const SurfaceInteraction& rec, const SampledWavelengths& lambda) const;

    private:
        Real alphaFor(const SurfaceInteraction& rec) const { return rec.frontFace ? m_alphaTop : m_alphaBottom; }
        /// The directional albedo E(mu) of a face's GGX lobe (with the Fresnel term at the index of 550 nm) and its cosine average.
        struct CuticleTable {
            static constexpr int kSize = 65;
            Real e[kSize] = {};
            Real avg = 0;
            Real at(Real mu) const;
        };
        static CuticleTable cuticleTable(Real alpha, Real n);
        Real lobeProbability(const SurfaceInteraction& rec, const Optics& o, Real cosO, Real& pReflect) const;
        Spectrum specular(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const Optics& o) const;
        Spectrum diffuseLobe(const SurfaceInteraction& rec, const Optics& o, Real cosO, Real cosI, bool reflect) const;

        std::shared_ptr<const leaf::ProspectTable> m_table;
        leaf::Pigments m_p;
        Real m_alphaTop, m_alphaBottom;
        Vector3 m_rgbR, m_rgbT;              // R and T under D65 (for the AOV and the log)
        Vector3 m_rgbRd, m_rgbTd;            // the RGB build's lobes: Rd and T as linear sRGB
        Real m_n550 = Real(1.47);
        CuticleTable m_esTop, m_esBottom;
        std::shared_ptr<const Texture2D> m_map;
        std::array<Real, 4> m_scale{ 1, 1, 1, 1 };
        // The RGB build with a map: Rd and T as linear sRGB on a grid over the four map channels (interpolated per hit).
        static constexpr int kGrid[4] = { 11, 7, 11, 6 };
        std::vector<Vector3> m_gridRd, m_gridT;
        void buildRgbGrid();
        void rgbFromGrid(const std::array<Real, 4>& texel, Vector3& rd, Vector3& t) const;
    };

} // namespace rayt
