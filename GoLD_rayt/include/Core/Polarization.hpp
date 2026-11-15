/**
 * @file Core/Polarization.hpp
 * @brief Stokes vectors, Mueller matrices and the complex Fresnel amplitudes (docs/POLARIZATION.md).
 *
 * Conventions
 * - A ray travelling along k carries its Stokes vector (S0, S1, S2, S3) in a basis (x, y, k) with y = k x x (right handed):
 *   E = Ex x + Ey y, S0 = |Ex|^2 + |Ey|^2, S1 = |Ex|^2 - |Ey|^2, S2 = 2 Re(Ex Ey*), S3 = -2 Im(Ex Ey*)   (time factor exp(-i w t)).
 *   S1 > 0: polarized along x, S2 > 0: along the diagonal x + y.
 * - At an interface the basis is (s, p, k) with s perpendicular to the plane of incidence and p = k x s, on both sides of the
 *   interface. The Fresnel amplitudes below are those of Born & Wolf in that basis (r_p = -r_s at normal incidence).
 * - Light flows from the source to the camera; the integrator walks the other way, so the chain matrices are multiplied on the right.
 *
 * Only what the path tracer needs: ideal smooth interfaces (a diagonal Jones matrix in the (s, p) basis) and basis rotations.
 */

#pragma once

#include <cmath>
#include <complex>

#include "Core/Types.hpp"
#include "Core/Math.hpp"

namespace rayt::pol {

    using Complex = std::complex<Real>;

#if defined(RAYT_SPECTRAL_MODE)
    constexpr int kN = kSpectralSamples;   // components of a Spectrum (wavelengths of a path)
#else
    constexpr int kN = 3;                  // RGB channels
#endif

    /// 4 x 4 Mueller matrix; every entry is a Spectrum (one value per wavelength / channel).
    struct Mueller {
        Spectrum e[4][4];

        Mueller() { for (int r = 0; r < 4; ++r) for (int c = 0; c < 4; ++c) e[r][c] = Spectrum(0.0); }

        static Mueller identity() {
            Mueller m;
            for (int i = 0; i < 4; ++i) m.e[i][i] = Spectrum(1.0);
            return m;
        }
    };

    inline Mueller operator*(const Mueller& a, const Mueller& b) {
        Mueller r;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                Spectrum s(0.0);
                for (int k = 0; k < 4; ++k) s += a.e[i][k] * b.e[k][j];
                r.e[i][j] = s;
            }
        return r;
    }

    /// The Stokes basis rotation: the same polarization state, expressed in a basis whose x axis is rotated by psi (cos2psi, sin2psi given).
    inline Mueller rotation(Real cos2psi, Real sin2psi) {
        Mueller m = Mueller::identity();
        m.e[1][1] = Spectrum(cos2psi);  m.e[1][2] = Spectrum(sin2psi);
        m.e[2][1] = Spectrum(-sin2psi); m.e[2][2] = Spectrum(cos2psi);
        return m;
    }

    /**
     * @brief Rotation from the basis A = (a, k x a) to the basis B = (b, k x b), both perpendicular to k.
     *
     * The angle psi of b in the basis A: cos psi = b . a, sin psi = b . (k x a).
     */
    inline Mueller basisChange(const Vector3& a, const Vector3& b, const Vector3& k) {
        const Real c = glm::dot(b, a);
        const Real s = glm::dot(b, glm::cross(k, a));
        return rotation(c * c - s * s, Real(2) * c * s);
    }

    /// Amplitude coefficients of an interface for the s and p components, for light arriving on the side of relative index 1.
    struct Amplitudes {
        Complex rs, rp;   ///< reflection
        Complex ts, tp;   ///< transmission (field amplitude; the energy factor is the integrator's business)
    };

    /**
     * @brief Fresnel amplitude coefficients.
     * @param cosThetaI cosine of the angle of incidence (on the side of the light).
     * @param eta       index of the second medium relative to the first: n2 / n1, complex for a conductor (eta = n + i k, k >= 0).
     *
     * Total internal reflection gives |rs| = |rp| = 1 with a phase difference (the Fresnel rhomb).
     */
    inline Amplitudes fresnelAmplitudes(Real cosThetaI, Complex eta) {
        const Real c = math::saturate(std::abs(cosThetaI));
        const Real sin2 = Real(1) - c * c;
        const Complex cosT = [&] {
            Complex ct = std::sqrt(Complex(1) - Complex(sin2) / (eta * eta));
            const Complex ec = eta * ct;
            // The wave must propagate into the second medium or decay in it: Im(eta cos) > 0, or Re > 0 when it is real.
            if (std::imag(ec) < Real(0) || (std::imag(ec) == Real(0) && std::real(ec) < Real(0))) ct = -ct;
            return ct;
        }();
        const Complex ec = eta * cosT;
        Amplitudes a;
        a.rs = (c - ec) / (c + ec);
        a.rp = (eta * c - cosT) / (eta * c + cosT);
        a.ts = Real(2) * c / (c + ec);
        a.tp = Real(2) * c / (eta * c + cosT);
        return a;
    }

    /**
     * @brief Mueller matrix of an interface that multiplies the s and p amplitudes by `as` and `ap` (reflection or transmission),
     *        divided by its first element: M00 = 1.
     *
     * The first element (the unpolarized gain, 1/2 (|as|^2 + |ap|^2), times an energy factor for transmission) is the scalar the
     * integrator already applies; this matrix is the polarization-dependent part on top of it. Writes the entries of component i.
     */
    inline void setInterface(Mueller& m, int i, Complex as, Complex ap) {
        const Real a2 = std::norm(as), p2 = std::norm(ap);
        const Real n = Real(0.5) * (a2 + p2);
        if (!(n > Real(0))) {   // nothing passes: any matrix will do, the scalar throughput is 0
            for (int r = 0; r < 4; ++r) m.e[r][r][i] = Real(1);
            return;
        }
        const Real inv = Real(1) / n;
        const Complex c = as * std::conj(ap) * inv;
        m.e[0][0][i] = Real(1);
        m.e[0][1][i] = m.e[1][0][i] = Real(0.5) * (a2 - p2) * inv;
        m.e[1][1][i] = Real(1);
        m.e[2][2][i] = m.e[3][3][i] = std::real(c);
        m.e[2][3][i] = std::imag(c);
        m.e[3][2][i] = -std::imag(c);
    }

    /// The same amplitudes for every component (a dielectric whose index does not depend on the wavelength of the path).
    inline Mueller interfaceMueller(Complex as, Complex ap) {   // not "interface": a macro in the Windows headers
        Mueller m;
        for (int i = 0; i < kN; ++i) setInterface(m, i, as, ap);
        return m;
    }

    /**
     * @brief What a material tells the integrator about one specular event it has just sampled.
     *
     * `m` is the polarization-dependent part of the event (M00 = 1), `s` the unit vector perpendicular to the plane of incidence
     * (the s axis of both bases: the one of the outgoing ray, along wo, and the one of the incoming ray, along -wi).
     */
    struct SpecularEvent {
        Mueller m;
        Vector3 s = Vector3(1, 0, 0);
    };

    /// The s axis for an event: perpendicular to the plane of wo and the normal (any perpendicular when they are parallel).
    inline Vector3 planeNormalS(const Vector3& wo, const Vector3& n) {
        Vector3 s = glm::cross(wo, n);
        const Real len = glm::length(s);
        if (len > Real(1e-6)) return s / len;
        // Normal incidence: the plane is undefined; the same axis serves both bases, so any perpendicular does.
        const Vector3 a = std::abs(wo.x) < Real(0.9) ? Vector3(1, 0, 0) : Vector3(0, 1, 0);
        return glm::normalize(glm::cross(wo, a));
    }

} // namespace rayt::pol
