/**
 * @file Color/ChromaticAdaptation.hpp
 * @brief CAT16 chromatic adaptation (von Kries in the CAM16 cone-like space).
 *
 * Reference: C. Li, Z. Li, Z. Wang, Y. Xu, M. R. Luo, G. Cui, M. Melgosa, M. H. Brill, M. Pointer,
 * "Comprehensive color solutions: CAM16, CAT16, and CAM16-UCS", Color Research & Application 42(6), 2017.
 *
 * What it does: the XYZ of a color is taken to a cone-like space (L, M, S) by the matrix M16, each of the three channels
 * is scaled by (destination white / source white), and the result is taken back to XYZ. With D = 1 the source white lands
 * exactly on the destination white; D < 1 adapts only partly (the channel scale becomes D * ratio + 1 - D).
 *
 * Why not divide RGB by the light's RGB: the sRGB primaries are not cone-like, and the white of a very reddish light
 * (a candle, 1850 K) has a negative blue channel in sRGB, so the division is meaningless there. The white of any real light
 * has three positive L, M, S values, so this ratio is always defined. The whole transform is linear, so out-of-gamut
 * (negative) colors pass through it unharmed.
 */

#pragma once

#include "Core/Types.hpp"

namespace rayt::color {

    using Matrix3 = glm::mat<3, 3, Real, glm::defaultp>;

    /// A 3x3 matrix from its rows (glm stores columns).
    inline Matrix3 matrixFromRows(Real a, Real b, Real c, Real d, Real e, Real f, Real g, Real h, Real i) {
        Matrix3 m(Real(0));
        m[0][0] = a; m[1][0] = b; m[2][0] = c;
        m[0][1] = d; m[1][1] = e; m[2][1] = f;
        m[0][2] = g; m[1][2] = h; m[2][2] = i;
        return m;
    }

    /// CAM16's matrix: XYZ -> (L, M, S) (Li et al. 2017).
    inline const Matrix3& matrixM16() {
        static const Matrix3 m = matrixFromRows(
            Real(0.401288), Real(0.650173), Real(-0.051461),
            Real(-0.250268), Real(1.204414), Real(0.045854),
            Real(-0.002079), Real(0.048952), Real(0.953127));
        return m;
    }

    /**
     * @brief XYZ -> XYZ matrix that adapts from white `srcWhite` to white `dstWhite` (both XYZ, any common scale).
     * @param degree D in [0,1]: 1 = complete adaptation.
     * @pre L, M, S of the source white are positive (true for any real white).
     */
    inline Matrix3 cat16Matrix(const Vector3& srcWhite, const Vector3& dstWhite, Real degree = Real(1)) {
        const Matrix3& M = matrixM16();
        const Vector3 ls = M * srcWhite;
        const Vector3 ld = M * dstWhite;
        Matrix3 G(Real(0));
        for (int c = 0; c < 3; ++c) G[c][c] = degree * (ld[c] / ls[c]) + (Real(1) - degree);
        return glm::inverse(M) * G * M;
    }

    /// True when all three cone-like responses of an XYZ white are positive (CAT16 can adapt from it).
    inline bool cat16WhiteValid(const Vector3& white) {
        const Vector3 lms = matrixM16() * white;
        return lms.x > Real(0) && lms.y > Real(0) && lms.z > Real(0);
    }

} // namespace rayt::color
