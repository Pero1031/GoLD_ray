/**
 * @file Color/Gamut.hpp
 * @brief Output gamuts (sRGB, Display P3, Rec.2020) and the conversion from the film's linear sRGB.
 *
 * The film always holds *linear sRGB* (BT.709 primaries, D65 white). The spectral integrator converts
 * each path's spectrum to XYZ and then to these primaries, and that conversion does not clip: a color
 * outside the sRGB triangle (a candle-lit gold, a saturated LED) comes out with a negative channel.
 * Because sRGB, Display P3 and Rec.2020 all share the D65 white point, going to another set of primaries
 * is one exact 3x3 matrix on those linear values - no chromatic adaptation, and no information lost as
 * long as nothing was clipped before. So the output stage can pick the gamut late, and the renderer does
 * not change. (Clipping happens only after the conversion, when the values are written to an 8-bit file.)
 *
 * The three matrices are the usual sRGB -> target products M_target^-1 * M_sRGB, with the primaries
 *   sRGB    R (0.640, 0.330)  G (0.300, 0.600)  B (0.150, 0.060)
 *   P3      R (0.680, 0.320)  G (0.265, 0.690)  B (0.150, 0.060)
 *   Rec2020 R (0.708, 0.292)  G (0.170, 0.797)  B (0.131, 0.046)     white D65 (0.3127, 0.3290) for all.
 * tools/test/gamut_check.py recomputes them from these chromaticities and compares.
 */

#pragma once

#include <string>

#include "Core/Types.hpp"
#include "Color/ChromaticAdaptation.hpp"

namespace rayt::color {

    enum class Gamut { SRGB, DisplayP3, Rec2020 };

    /// "srgb" | "p3" | "rec2020" (also "display-p3"); returns false for anything else.
    inline bool parseGamut(const std::string& s, Gamut& out) {
        if (s == "srgb") { out = Gamut::SRGB; return true; }
        if (s == "p3" || s == "display-p3" || s == "displayp3") { out = Gamut::DisplayP3; return true; }
        if (s == "rec2020" || s == "bt2020") { out = Gamut::Rec2020; return true; }
        return false;
    }

    inline const char* gamutName(Gamut g) {
        switch (g) {
        case Gamut::SRGB: return "sRGB";
        case Gamut::DisplayP3: return "Display P3";
        case Gamut::Rec2020: return "Rec.2020";
        }
        return "?";
    }

    /// The matrix linear sRGB -> linear RGB of `g` (rows as in the header comment; identity for sRGB).
    inline const Matrix3& matrixFromLinearSRGB(Gamut g) {
        static const Matrix3 p3 = matrixFromRows(
            Real(0.8224621), Real(0.1775380), Real(0.0000000),
            Real(0.0331941), Real(0.9668058), Real(0.0000000),
            Real(0.0170827), Real(0.0723974), Real(0.9105199));
        static const Matrix3 rec2020 = matrixFromRows(
            Real(0.6274040), Real(0.3292820), Real(0.0433136),
            Real(0.0690970), Real(0.9195400), Real(0.0113612),
            Real(0.0163916), Real(0.0880132), Real(0.8955950));
        static const Matrix3 identity(Real(1));
        switch (g) {
        case Gamut::DisplayP3: return p3;
        case Gamut::Rec2020: return rec2020;
        case Gamut::SRGB:
        default: return identity;
        }
    }

    /// Linear sRGB (may be negative or above 1) -> linear RGB in the primaries of `g`. Not clamped.
    inline Vector3 fromLinearSRGB(Gamut g, const Vector3& c) {
        if (g == Gamut::SRGB) return c;
        return matrixFromLinearSRGB(g) * c;
    }

    /// The inverse: linear RGB of `g` -> linear sRGB (negative where the color is outside the sRGB triangle).
    inline Vector3 toLinearSRGB(Gamut g, const Vector3& c) {
        static const Matrix3 p3Inv = glm::inverse(matrixFromLinearSRGB(Gamut::DisplayP3));
        static const Matrix3 rec2020Inv = glm::inverse(matrixFromLinearSRGB(Gamut::Rec2020));
        switch (g) {
        case Gamut::DisplayP3: return p3Inv * c;
        case Gamut::Rec2020: return rec2020Inv * c;
        case Gamut::SRGB:
        default: return c;
        }
    }

    /// H.273 colour primaries code (for the PNG cICP chunk): 1 = BT.709/sRGB, 12 = Display P3 (D65), 9 = BT.2020.
    inline int gamutPrimariesCode(Gamut g) {
        switch (g) {
        case Gamut::DisplayP3: return 12;
        case Gamut::Rec2020: return 9;
        case Gamut::SRGB:
        default: return 1;
        }
    }

} // namespace rayt::color
