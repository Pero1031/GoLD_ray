/**
 * @file Materials/Prospect.hpp
 * @brief The PROSPECT leaf model: the reflectance and transmittance of a leaf from its pigments and its structure.
 *
 * PROSPECT (Jacquemoud and Baret 1990) treats a leaf as a pile of N absorbing plates with flat dielectric faces (Allen et al. 1969,
 * Stokes 1862). One plate transmits the fraction tau(k) = (1 - k) e^-k + k^2 E1(k) of the isotropic light that enters it, where
 * the absorption k(lambda) = sum_i K_i(lambda) C_i / N comes from the specific absorption K_i and the content C_i of each
 * constituent. PROSPECT-D (Feret et al. 2017) has chlorophylls, carotenoids, anthocyanins, brown pigments, water and dry matter;
 * the coefficient table here is the one of PROSPECT-PRO (Feret et al. 2021), whose pigment columns are those of PROSPECT-D.
 *
 * The table (dataSpec_PRO_v2.txt of the authors' R package, 400-2500 nm, 1 nm, GPL-3) is not part of this repository:
 * tools/data/leaf_optics.py fetch downloads it to assets/Leaf/prospect/ (git-ignored). Below 400 nm the 400 nm row is used.
 *
 * Units (as in PROSPECT): chlorophyll a+b, carotenoids and anthocyanins in ug/cm^2, brown pigments in arbitrary units (0..1),
 * equivalent water thickness (water) and leaf mass per area (dry matter) in g/cm^2 (= cm for water).
 */
#pragma once
#include <memory>
#include <string>
#include <vector>
#include "Core/Types.hpp"

namespace rayt::leaf {

    /// What a leaf is made of. The defaults are a mature green leaf (the PROSPECT-D defaults).
    struct Pigments {
        Real N = Real(1.5);           ///< structure: the number of plates (1-1.5 young or thin, 1.5-2.5 mature, more: senescent)
        Real chl = Real(40);          ///< chlorophyll a + b (ug/cm^2)
        Real car = Real(8);           ///< carotenoids (ug/cm^2)
        Real ant = Real(0);           ///< anthocyanins (ug/cm^2)
        Real brown = Real(0);         ///< brown pigments (arbitrary units)
        Real ewt = Real(0.01);        ///< equivalent water thickness (g/cm^2)
        Real lma = Real(0.009);       ///< leaf mass per area, dry matter (g/cm^2)
    };

    /// The specific absorption coefficients and the refractive index of the plates, per nm.
    class ProspectTable {
    public:
        /// Reads a PROSPECT coefficient table (tab-separated, header "lambda nrefrac sac_chl sac_car sac_ant sac_brown sac_ewt sac_lma ...").
        /// Returns null and sets `err` if the file cannot be read.
        static std::shared_ptr<const ProspectTable> load(const std::string& path, std::string* err = nullptr);
        /// The table at its usual place, assets/Leaf/prospect/dataSpec_PRO_v2.txt (relative to the working folder), loaded once.
        static std::shared_ptr<const ProspectTable> standard(std::string* err = nullptr);
        static const char* standardPath() { return "assets/Leaf/prospect/dataSpec_PRO_v2.txt"; }

        struct Row {
            Real n = Real(1.5);                      ///< refractive index of the plates
            Real chl = 0, car = 0, ant = 0, brown = 0, ewt = 0, lma = 0;   ///< specific absorption coefficients
            Real tav90 = Real(0.9);                  ///< transmissivity of a face for isotropic light (alpha = 90 degrees)
            Real tav40 = Real(0.95);                 ///< ... for light within 40 degrees of the normal (the measurements' geometry)
        };
        /// The coefficients at a wavelength (nm), linear between the rows, clamped to the table's range.
        Row at(Real lambdaNm) const;
        Real lambdaMin() const { return m_lambda0; }
        Real lambdaMax() const { return m_lambda0 + Real(m_rows.size() - 1); }

    private:
        Real m_lambda0 = 400;
        std::vector<Row> m_rows;   // one per nm from m_lambda0
    };

    /// The result of PROSPECT at one wavelength.
    struct ReflTrans {
        Real R = 0;        ///< reflectance of the leaf (directional-hemispherical for light within alpha of the normal)
        Real T = 0;        ///< transmittance
        Real rFace = 0;    ///< the part of R reflected by the first face, 1 - tav(alpha, n) (a mirror-like reflection)
        Real t12 = 0;      ///< transmissivity of a face for isotropic light, tav(90, n)
        Real n = 0;        ///< refractive index of the plates
    };

    /// PROSPECT at one wavelength for light within `alphaDeg` of the normal: 90 = isotropic light (what a rendering needs),
    /// 40 = the geometry of the integrating-sphere measurements PROSPECT was calibrated on (its default).
    ReflTrans prospect(const ProspectTable& table, const Pigments& p, Real lambdaNm, Real alphaDeg = Real(90));

    /// Transmissivity of a flat dielectric face of index n for isotropic light within alphaDeg of the normal (Stern 1964, Allen 1973).
    Real tav(Real alphaDeg, Real n);

    /// The exponential integral E1(x) = Integral_x^inf e^-t / t dt for x > 0 (series below 1, continued fraction above).
    Real expintE1(Real x);

} // namespace rayt::leaf
