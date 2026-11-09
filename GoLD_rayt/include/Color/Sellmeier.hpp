/**
 * @file Color/Sellmeier.hpp
 * @brief Sellmeier dispersion of optical glasses, with a few presets.
 */

#pragma once

#include <cmath>
#include <string>

#include "Core/Types.hpp"

namespace rayt {

    /// Sellmeier dispersion  n^2(l) = 1 + sum_i B_i l^2 / (l^2 - C_i)  (l in micrometres, C_i in um^2). Two- and three-term forms (B3 = 0 for two).
    struct Sellmeier {
        Real b[3] = { 0, 0, 0 };
        Real c[3] = { 0, 0, 0 };
        Real n(Real um) const {
            const Real l2 = um * um;
            Real s = 1;
            for (int i = 0; i < 3; ++i) if (b[i] != Real(0)) s += b[i] * l2 / (l2 - c[i]);
            return std::sqrt(s);
        }
    };

    /// Sellmeier coefficients of some glasses ("diamond", "sapphire", "bk7", "silica"). Returns false for an unknown name.
    /// Diamond: the two-term fit of Peter, Z. Phys. 15 (1923) (n_d 2.417); sapphire: Malitson 1962, ordinary ray (1.768); BK7: Schott (1.5168);
    /// fused silica: Malitson 1965 (1.4585). Written from the usual published values: check n_d and the Abbe number (--scene dispersion).
    inline bool sellmeierPreset(const std::string& name, Sellmeier& out) {
        if (name == "diamond") { out = Sellmeier{ { Real(0.3306), Real(4.3356), 0 }, { Real(0.1750 * 0.1750), Real(0.1060 * 0.1060), 0 } }; return true; }
        if (name == "sapphire") { out = Sellmeier{ { Real(1.4313493), Real(0.65054713), Real(5.3414021) }, { Real(0.0726631 * 0.0726631), Real(0.1193242 * 0.1193242), Real(18.028251 * 18.028251) } }; return true; }
        if (name == "bk7") { out = Sellmeier{ { Real(1.03961212), Real(0.231792344), Real(1.01046945) }, { Real(0.00600069867), Real(0.0200179144), Real(103.560653) } }; return true; }
        if (name == "silica") { out = Sellmeier{ { Real(0.6961663), Real(0.4079426), Real(0.8974794) }, { Real(0.0684043 * 0.0684043), Real(0.1162414 * 0.1162414), Real(9.896161 * 9.896161) } }; return true; }
        return false;
    }

} // namespace rayt
