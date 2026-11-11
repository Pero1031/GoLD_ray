/**
 * @file Materials/KubelkaMunkLayer.hpp
 * @brief A coloured layer of paint-like material in the Kubelka-Munk two-flux model (docs/COSMETICS.md, 4.2).
 *
 * The layer is given by its intrinsic colour R_inf (the reflectance of an opaque film of it) and its opacity X = S D (the scattering
 * coefficient times the thickness). With a = (R_inf + 1 / R_inf) / 2 and b = sqrt(a^2 - 1) (per channel):
 *     R = sinh(bX) / (a sinh(bX) + b cosh(bX)),    T = b / (a sinh(bX) + b cosh(bX))
 * is the reflectance of the layer alone (nothing behind it) and its transmittance, both for diffuse light.
 * On a surface of reflectance A the stack reflects  R + T^2 A / (1 - R A)  (Kubelka 1954; Doi et al. 2005; Yoshida and Okiyama 2021).
 * Without absorption (R_inf = 1) the limit is R = X / (1 + X), T = 1 / (1 + X).
 */

#pragma once

#include <algorithm>
#include <cmath>

#include "Core/Types.hpp"

namespace rayt {

    struct KubelkaMunkLayer {
        Vector3 R{ 0.0 };   ///< reflectance of the layer alone
        Vector3 T{ 1.0 };   ///< transmittance

        KubelkaMunkLayer() = default;

        /// `intrinsic`: the colour of an opaque film (each channel in (0, 1]); `opacity` X >= 0 (0 = no layer).
        KubelkaMunkLayer(const Vector3& intrinsic, Real opacity) {
            const Real x = std::max(opacity, Real(0));
            for (int c = 0; c < 3; ++c) {
                const Real rInf = std::clamp(intrinsic[c], Real(0.01), Real(1));
                const Real a = Real(0.5) * (rInf + Real(1) / rInf);
                const Real b2 = a * a - Real(1);
                if (b2 < Real(1e-9)) { R[c] = x / (Real(1) + x); T[c] = Real(1) / (Real(1) + x); continue; }
                const Real b = std::sqrt(b2);
                const Real bx = b * x;
                const Real sh = std::sinh(bx), ch = std::cosh(bx);
                const Real den = a * sh + b * ch;
                R[c] = sh / den;
                T[c] = b / den;
            }
        }

        /// Reflectance of this layer lying on a surface of reflectance A (per channel).
        Vector3 over(const Vector3& A) const {
            Vector3 out;
            for (int c = 0; c < 3; ++c) out[c] = R[c] + T[c] * T[c] * A[c] / std::max(Real(1) - R[c] * A[c], Real(1e-4));
            return out;
        }
    };

} // namespace rayt
