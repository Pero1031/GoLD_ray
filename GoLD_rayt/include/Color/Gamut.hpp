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
 * tools/test/gamut_matrices.py recomputes them from these chromaticities and compares.
 */

#pragma once

#include <string>

#include "Core/Types.hpp"

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

    /// Linear sRGB (may be negative or above 1) -> linear RGB in the primaries of `g`. Not clamped.
    inline Vector3 fromLinearSRGB(Gamut g, const Vector3& c) {
        switch (g) {
        case Gamut::DisplayP3:
            return Vector3(
                Real(0.8224621) * c.x + Real(0.1775380) * c.y,
                Real(0.0331941) * c.x + Real(0.9668058) * c.y,
                Real(0.0170827) * c.x + Real(0.0723974) * c.y + Real(0.9105199) * c.z);
        case Gamut::Rec2020:
            return Vector3(
                Real(0.6274040) * c.x + Real(0.3292820) * c.y + Real(0.0433136) * c.z,
                Real(0.0690970) * c.x + Real(0.9195400) * c.y + Real(0.0113612) * c.z,
                Real(0.0163916) * c.x + Real(0.0880132) * c.y + Real(0.8955950) * c.z);
        case Gamut::SRGB:
        default:
            return c;
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
