/**
 * @file Materials/ThinFilm.hpp
 * @brief Thin-film interference: soap bubbles, oil on water, coated or tarnished metal.
 *
 * A smooth surface covered by one thin transparent layer (index nF, thickness d of the order of the wavelength). Light reflected at the top
 * of the layer and at its bottom interferes, and the result depends on  d, the angle and - strongly - the wavelength: the colors of a soap
 * bubble. An RGB renderer can not show this (the reflectance oscillates many times inside one color channel); the spectral build
 * evaluates it at the wavelengths of the path.
 *
 *   ambient (nI) | film (nF, d) | substrate (nT, possibly complex)
 *
 * The reflectance is Airy's formula for each polarization, averaged (unpolarized light):
 *     r = (r01 + r12 e^{2 i beta}) / (1 + r01 r12 e^{2 i beta}),     beta = 2 pi nF d cos(thetaF) / lambda,     R = (|r_s|^2 + |r_p|^2) / 2
 * with the Fresnel amplitudes r01 (ambient to film) and r12 (film to substrate) written with complex indices, so that a metal substrate
 * (nT = eta + i k) needs no special case. See docs/THIN_FILM.md for the derivation.
 *
 * Three substrates:  Free        a free film in air on both sides (a soap bubble): reflects R and transmits T = 1 - R straight through.
 *                    Dielectric  a transparent substrate of index n2 (oil on water): reflects R, transmits T = 1 - R, refracted.
 *                    Conductor   a metal under the film (anodized titanium, lacquered gold): reflects R, nothing is transmitted.
 * The film itself does not absorb (T = 1 - R) and has no dispersion; the substrate of a conductor uses its tabulated (eta, k).
 *
 * The thickness may vary over the surface: a constant, a linear gradient along an axis (a bubble drained by gravity), and a noise
 * modulation (swirls). In RGB mode the reflectance is integrated against the D65 color of the film sensor, i.e. the color the film
 * has under white light; everything that depends on the spectrum of the light is then lost (that is the point of the comparison).
 */

#pragma once

#include <complex>
#include <memory>
#include <optional>
#include <vector>

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Interaction.hpp"
#include "Color/ComplexIorTable.hpp"
#include "Materials/DensityField.hpp"
#include "Materials/Material.hpp"

namespace rayt {

    namespace thinfilm {

        using Complex = std::complex<Real>;

        /**
         * @brief Unpolarized reflectance of ambient | film | substrate.
         * @param cosI      cosine of the angle of incidence in the ambient medium, in [0, 1]
         * @param nI        index of the ambient medium (real)
         * @param nF        index of the film (real)
         * @param dNm       thickness of the film in nm (0 gives the plain interface ambient | substrate)
         * @param nT        index of the substrate: eta + i k (k = 0 for a dielectric)
         * @param lambdaNm  wavelength in nm
         */
        Real reflectance(Real cosI, Real nI, Real nF, Real dNm, Complex nT, Real lambdaNm);

    } // namespace thinfilm

    /// Thickness of the film over the surface, in nm.
    ///
    /// A soap film drains: it is thin at the top (down to a black film) and thickens downwards, in horizontal bands that swirl where the liquid
    /// circulates. The model: a gradient along an axis with a power-law profile (`exponent` > 1 keeps the top thin for longer), the position bent
    /// by a smooth displacement field (`warp`, in metres: the bands wobble and swirl), and a fractal noise (`noiseAmplitude`, stretched by
    /// `noiseStretch`: a stretch of (1, 3, 1) makes features thin along y, i.e. horizontal streaks).
    struct FilmThickness {
        Real base = Real(400);                  ///< constant thickness, or the thickness at the start of the gradient
        bool gradient = false;                  ///< a gradient along `axis` (world space) ...
        Vector3 axis = Vector3(0, 1, 0);
        Real s0 = 0, s1 = 1;                    ///< ... from `base` at dot(p, axis) = s0 to `top` at s1 (clamped outside)
        Real top = Real(400);
        Real exponent = Real(1);                ///< profile: d = base + (top - base) * t^exponent
        Real warp = 0;                          ///< displacement of the position, up to +- this many metres (0 = none)
        std::shared_ptr<const DensityField> warpX, warpY, warpZ;
        Real noiseAmplitude = 0;                ///< +- nm added by a fractal noise (swirls); 0 = none
        Vector3 noiseStretch = Vector3(1.0);    ///< the noise is evaluated at the position multiplied by this
        std::shared_ptr<const DensityField> noise;

        Real at(const Point3& p) const {
            Point3 q = p;
            if (warp != Real(0) && warpX && warpY && warpZ)
                q += warp * Vector3(Real(2) * (*warpX)(p) - Real(1), Real(2) * (*warpY)(p) - Real(1), Real(2) * (*warpZ)(p) - Real(1));
            Real d = base;
            if (gradient) {
                const Real t = std::clamp((glm::dot(q, axis) - s0) / (s1 - s0), Real(0), Real(1));
                d = base + (top - base) * std::pow(t, exponent);
            }
            if (noise && noiseAmplitude != Real(0)) d += noiseAmplitude * (Real(2) * (*noise)(Point3(q.x * noiseStretch.x, q.y * noiseStretch.y, q.z * noiseStretch.z)) - Real(1));
            return std::max(d, Real(0));
        }
    };

    class ThinFilm : public Material {
    public:
        enum class Substrate { Free, Dielectric, Conductor };

        /// A free film (the substrate is the ambient medium on the other side): a soap bubble.
        ThinFilm(Real filmIor, const FilmThickness& thickness);

        /// A transparent substrate of index n2 under the film (oil on water: filmIor 1.45, n2 1.33).
        ThinFilm(Real filmIor, const FilmThickness& thickness, Real substrateIor);

        /// A metal under the film (its tabulated eta(lambda), k(lambda)).
        ThinFilm(Real filmIor, const FilmThickness& thickness, const std::shared_ptr<const color::ComplexIorTable>& metal);

        Spectrum eval(const SurfaceInteraction&, const Vector3&, const Vector3&, const SampledWavelengths&, TransportMode) const override {
            return Spectrum(0.0);   // delta lobes only
        }
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u,
            const SampledWavelengths& lambda, TransportMode mode) const override;
        Real pdf(const SurfaceInteraction&, const Vector3&, const Vector3&, const SampledWavelengths&) const override { return Real(0); }
        bool isSpecular() const override { return true; }
        Vector3 aovAlbedo(const SurfaceInteraction&) const override { return Vector3(0.9); }

        /// Reflectance and transmittance of the film at a point, per sampled wavelength (RGB mode: per channel). Used by sample() and the checks.
        void reflectTransmit(const Point3& p, Real cosI, bool entering, const SampledWavelengths& lambda, Spectrum& R, Spectrum& T) const;

        /// Where the film is: opacity q(p) in [0, 1] (stochastic transparency: a ray passes through where there is no film). An irregular patch of
        /// oil on clean water is a noise field with a narrow band. Without it the film covers the whole surface.
        void setCoverage(std::shared_ptr<const DensityField> coverage) { m_coverage = std::move(coverage); }
        Real opacity(const SurfaceInteraction& rec) const override { return m_coverage ? (*m_coverage)(rec.p) : Real(1); }

        Substrate substrate() const { return m_substrate; }

    private:
        Real m_nF;
        FilmThickness m_thickness;
        std::shared_ptr<const DensityField> m_coverage;
        Substrate m_substrate = Substrate::Free;
        Real m_n2 = Real(1);                       ///< substrate index (Dielectric); 1 for a free film
        std::vector<Real> m_eta, m_k;              ///< metal: eta, k on a 1 nm grid from kFirstNm
        static constexpr int kFirstNm = 360, kLastNm = 830;

        Real reflectanceAt(Real cosI, bool entering, Real d, Real lambdaNm) const;

#if !defined(RAYT_SPECTRAL_MODE)
        // RGB mode: the film color under white light (R integrated against the D65 sensor weights, 81 wavelengths) is tabulated once over
        // (cos theta, thickness), one table per side (entering / leaving), and interpolated per hit: integrating at every hit is far too slow.
        static constexpr int kTableCos = 64;
        static constexpr Real kTableStep = Real(4);          // nm
        static constexpr Real kTableMax = Real(6000);        // nm; thicker films are integrated directly
        int m_tableD = 0;                                    // number of thickness nodes
        std::vector<Vector3> m_table[2];                     // [entering ? 0 : 1][iCos * m_tableD + iD]
        void buildRgbTable();
        Vector3 rgbReflectance(Real cosI, bool entering, Real d) const;
#endif
    };

} // namespace rayt
