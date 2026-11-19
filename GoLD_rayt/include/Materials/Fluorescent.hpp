/**
 * @file Materials/Fluorescent.hpp
 * @brief A diffuse surface that fluoresces, described by a measured bispectral matrix (docs/FLUORESCENCE.md).
 *
 * A bispectral (Donaldson) matrix D(lambda_o, lambda_i) gives the radiance factor at the emission wavelength lambda_o for light of
 * the excitation wavelength lambda_i: the diagonal is the ordinary reflectance R, the rest is fluorescence. Each column stands for an
 * excitation band as wide as the column step (10 nm in the RIT database), so under a light of spectral irradiance E a diffuse surface
 * sends
 *     L(lambda_o) = (1 / pi) [ R(lambda_o) E(lambda_o) + Sum_c D(lambda_o, c) Ebar(c) ],
 * where Ebar(c) is E averaged over the band of column c (ASTM E2152 writes the same sum with E sampled at the column wavelengths).
 * Rows between the measured emission wavelengths are interpolated linearly.
 *
 * Path tracing from the camera: the path arrives with the wavelength lambda_o that the camera records. reemit() chooses reflection
 * with probability R / (R + F) and fluorescence with F / (R + F), where F(lambda_o) = Sum_c D(lambda_o, c); for fluorescence it picks
 * a column with probability D(lambda_o, c) / F and a uniform wavelength in its band. The path goes on at that wavelength, and the
 * vertex is a Lambertian surface of albedo R + F, which can exceed 1 at lambda_o: a fluorescent color can be brighter than a white.
 * The integrator follows one wavelength from the first such surface on (the secondary wavelengths are dropped).
 *
 * The RGB build cannot shift wavelengths: there the surface is a Lambertian of the color it shows under D65, fluorescence that the
 * ultraviolet of D65 excites included (its total radiance factor under D65).
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Core/Types.hpp"
#include "Materials/Material.hpp"

namespace rayt {

    /// A bispectral radiance-factor matrix: rows = emission wavelengths, columns = excitation wavelengths (both evenly spaced, nm).
    struct BispectralMatrix {
        std::string title;
        Real emFirst = Real(380), emStep = Real(10);
        int nEm = 0;
        Real exFirst = Real(300), exStep = Real(10);
        int nEx = 0;
        std::vector<Real> d;   ///< nEm x nEx, row-major

        Real at(int r, int c) const { return d[static_cast<size_t>(r) * static_cast<size_t>(nEx) + static_cast<size_t>(c)]; }
        Real em(int r) const { return emFirst + emStep * Real(r); }
        Real ex(int c) const { return exFirst + exStep * Real(c); }
        /// Column of the excitation wavelength equal to row r's emission wavelength (the reflectance), or -1.
        int diagonalColumn(int r) const;
    };

    /// What loading removed from a measured matrix (see LoadBfcMatrix).
    struct BispectralCleaning {
        int negative = 0;          ///< entries below 0 (measurement noise) set to 0
        double negativeSum = 0;
        int antiStokes = 0;        ///< entries whose excitation is 20 nm or more above the emission set to 0
        double antiStokesSum = 0;  ///< their sum (positive part)
        double noise = 0;          ///< RMS of those entries (they should be 0): the matrix's measurement noise
        int belowNoise = 0;        ///< fluorescence entries below 3 x noise set to 0
        double belowNoiseSum = 0;
    };

    /**
     * @brief Reads a Labsphere BFC-450 matrix file (".BFC", as in the RIT bispectral database of Gonzalez & Fairchild 2000).
     *
     * The file is text: comment lines (the fourth one, after ';', is the sample's name), then a header row "r:c:" with the excitation
     * wavelengths, then one row per emission wavelength. Entries whose excitation lies 20 nm or more above the emission become 0 (light
     * does not fluoresce to shorter wavelengths; next to the diagonal, +10 nm, the instrument's 10 nm bands overlap, so that column is
     * kept). Their RMS is the matrix's noise: fluorescence entries below 3 x noise become 0 too, so that a surface that does not
     * fluoresce does not glow under a light that is mostly ultraviolet, and negative entries of the diagonal become 0.
     */
    bool LoadBfcMatrix(const std::string& path, BispectralMatrix& out, BispectralCleaning* cleaning, std::string* err);

    class Fluorescent : public Material {
    public:
        /// @param fluorescence false: the same surface without its fluorescence (the diagonal alone), for comparison.
        explicit Fluorescent(std::shared_ptr<const BispectralMatrix> m, bool fluorescence = true);

        Spectrum eval(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const SampledWavelengths& lambda,
                      TransportMode mode = TransportMode::Radiance) const override;
        std::optional<BSDFSample> sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u, const SampledWavelengths& lambda,
                                         TransportMode mode = TransportMode::Radiance) const override;
        Real pdf(const SurfaceInteraction& rec, const Vector3& wo, const Vector3& wi, const SampledWavelengths& lambda) const override;
        Vector3 aovAlbedo(const SurfaceInteraction& rec) const override;

        bool shiftsWavelength() const override;
        void reemit(SurfaceInteraction& rec, SampledWavelengths& lambda, Real u0, Real u1) const override;

        /// R(lambda_o): the reflectance (the diagonal), interpolated between rows and held at the ends.
        Real reflectance(Real lambdaO) const;
        /// F(lambda_o): all fluorescence that leaves at lambda_o per unit of excitation (0 outside the rows, or when switched off).
        Real fluorescence(Real lambdaO) const;
        /// An excitation wavelength for light leaving at lambda_o, with density excitationDensity(lambda_o, .) / F(lambda_o).
        Real sampleExcitation(Real lambdaO, Real u) const;
        /// b(lambda_o, lambda_i): the fluorescence density per nm of excitation (D / column width), 0 on the diagonal column.
        Real excitationDensity(Real lambdaO, Real lambdaI) const;
        /// The linear sRGB color under D65 (with its ultraviolet), used by the RGB build and as the denoiser's albedo.
        const Vector3& rgbUnderD65() const { return m_rgbD65; }
        const BispectralMatrix& matrix() const { return *m_m; }

    private:
        std::shared_ptr<const BispectralMatrix> m_m;
        bool m_fluo;
        std::vector<Real> m_R, m_F;   ///< per row: reflectance and the sum of the fluorescence columns
        Vector3 m_rgbD65 = Vector3(0.5);

        /// Row r and weight t of lambda_o (lambda_o = (1 - t) em(r) + t em(r + 1)); false outside the rows.
        bool rowOf(Real lambdaO, int& r, Real& t) const;
        /// Fluorescence weight of column c for lambda_o between rows r and r + 1.
        Real columnWeight(int r, Real t, int c) const;
    };

} // namespace rayt
