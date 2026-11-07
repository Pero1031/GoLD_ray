// Thin-film interference (see Materials/ThinFilm.hpp and docs/THIN_FILM.md).
#include "pch.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>

#include "Materials/ThinFilm.hpp"
#include "Color/SpectralConversion.hpp"
#include "Color/Illuminants.hpp"
#include "Core/Math.hpp"

namespace rayt {

    namespace thinfilm {

        namespace {
            /// sqrt(1 - s^2) for the cosine of the refraction angle, on the branch where the wave does not grow into the layer:
            /// Im(n c) >= 0, and Re(n c) >= 0 when that is zero. (The principal root alone fails for an evanescent wave:
            /// the sign of a zero imaginary part decides the branch.)
            Complex cosine(Complex s, Complex n) {
                Complex c = std::sqrt(Complex(1) - s * s);
                const Complex nc = n * c;
                const Real tol = Real(1e-12) * (std::abs(nc) + Real(1e-30));
                if (nc.imag() < -tol || (std::abs(nc.imag()) <= tol && nc.real() < 0)) c = -c;
                return c;
            }
        }

        Real reflectance(Real cosI, Real nI, Real nF, Real dNm, Complex nT, Real lambdaNm) {
            cosI = std::clamp(cosI, Real(0), Real(1));
            const Real sinI = std::sqrt(std::max(Real(0), Real(1) - cosI * cosI));
            const Complex n0(nI, 0), n1(nF, 0);
            const Complex sin1 = n0 * sinI / n1;     // Snell: n sin(theta) is the same in every layer
            const Complex sin2 = n0 * sinI / nT;
            const Complex cos0(cosI, 0);
            const Complex cos1 = cosine(sin1, n1);
            const Complex cos2 = cosine(sin2, nT);

            // Fresnel amplitudes, ambient -> film (01) and film -> substrate (12), s and p polarization (the same convention at both)
            const Complex r01s = (n0 * cos0 - n1 * cos1) / (n0 * cos0 + n1 * cos1);
            const Complex r01p = (n1 * cos0 - n0 * cos1) / (n1 * cos0 + n0 * cos1);
            const Complex r12s = (n1 * cos1 - nT * cos2) / (n1 * cos1 + nT * cos2);
            const Complex r12p = (nT * cos1 - n1 * cos2) / (nT * cos1 + n1 * cos2);

            // phase thickness of one pass through the film; the round trip is 2 beta
            const Real twoPi = Real(2) * Real(3.14159265358979323846);
            const Complex beta = Complex(twoPi * dNm / lambdaNm, 0) * n1 * cos1;
            const Complex e2 = std::exp(Complex(0, 2) * beta);

            const Complex rs = (r01s + r12s * e2) / (Complex(1) + r01s * r12s * e2);
            const Complex rp = (r01p + r12p * e2) / (Complex(1) + r01p * r12p * e2);
            return std::clamp(Real(0.5) * (std::norm(rs) + std::norm(rp)), Real(0), Real(1));
        }

    } // namespace thinfilm

    namespace {
#if defined(RAYT_SPECTRAL_MODE)
        constexpr int kLanes = kSpectralSamples;
#else
        constexpr int kLanes = 3;

        /// RGB mode: weights w_k (3 channels) on a 5 nm grid with  sum_k w_k = (1, 1, 1): the film sensor under D65, so that a reflectance of 1
        /// gives (1, 1, 1) and a thin film gives the color it has under white light.
        struct RgbWeights {
            static constexpr int kN = 81;   // 380 .. 780 nm
            std::array<Real, kN> nm;
            std::array<Vector3, kN> w;
            RgbWeights() {
                Vector3 sum(0.0);
                for (int k = 0; k < kN; ++k) {
                    nm[static_cast<size_t>(k)] = Real(380 + 5 * k);
                    w[static_cast<size_t>(k)] = color::rgbSensorWeights(nm[static_cast<size_t>(k)]) * color::evaluateD65(nm[static_cast<size_t>(k)]);
                    sum += w[static_cast<size_t>(k)];
                }
                for (auto& v : w) v = Vector3(v.x / sum.x, v.y / sum.y, v.z / sum.z);
            }
        };
        const RgbWeights& rgbWeights() { static const RgbWeights t; return t; }
#endif
    }

    ThinFilm::ThinFilm(Real filmIor, const FilmThickness& thickness) : m_nF(filmIor), m_thickness(thickness) {
#if !defined(RAYT_SPECTRAL_MODE)
        buildRgbTable();
#endif
    }

    ThinFilm::ThinFilm(Real filmIor, const FilmThickness& thickness, Real substrateIor)
        : m_nF(filmIor), m_thickness(thickness), m_substrate(Substrate::Dielectric), m_n2(substrateIor) {
#if !defined(RAYT_SPECTRAL_MODE)
        buildRgbTable();
#endif
    }

    ThinFilm::ThinFilm(Real filmIor, const FilmThickness& thickness, const std::shared_ptr<const color::ComplexIorTable>& metal)
        : m_nF(filmIor), m_thickness(thickness), m_substrate(Substrate::Conductor) {
        const int n = kLastNm - kFirstNm + 1;
        m_eta.resize(static_cast<size_t>(n));
        m_k.resize(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            const auto [e, k] = metal->evaluateNK(Real(kFirstNm + i));
            m_eta[static_cast<size_t>(i)] = e;
            m_k[static_cast<size_t>(i)] = k;
        }
#if !defined(RAYT_SPECTRAL_MODE)
        buildRgbTable();
#endif
    }

    Real ThinFilm::reflectanceAt(Real cosI, bool entering, Real d, Real lambdaNm) const {
        Real nI = Real(1);
        thinfilm::Complex nT(1, 0);
        switch (m_substrate) {
        case Substrate::Free: break;
        case Substrate::Dielectric:
            if (entering) nT = thinfilm::Complex(m_n2, 0);
            else { nI = m_n2; nT = thinfilm::Complex(1, 0); }
            break;
        case Substrate::Conductor: {
            const Real x = std::clamp(lambdaNm, Real(kFirstNm), Real(kLastNm)) - Real(kFirstNm);
            const int i = std::min(static_cast<int>(x), static_cast<int>(m_eta.size()) - 2);
            const Real f = x - Real(i);
            const Real eta = m_eta[static_cast<size_t>(i)] * (Real(1) - f) + m_eta[static_cast<size_t>(i) + 1] * f;
            const Real k = m_k[static_cast<size_t>(i)] * (Real(1) - f) + m_k[static_cast<size_t>(i) + 1] * f;
            nT = thinfilm::Complex(eta, k);
            break;
        }
        }
        return thinfilm::reflectance(cosI, nI, m_nF, d, nT, lambdaNm);
    }

#if !defined(RAYT_SPECTRAL_MODE)
    namespace {
        Vector3 integrateRgb(const ThinFilm& f, const std::function<Real(Real)>& R) {
            (void)f;
            const RgbWeights& t = rgbWeights();
            Vector3 acc(0.0);
            for (int k = 0; k < RgbWeights::kN; ++k) acc += t.w[static_cast<size_t>(k)] * R(t.nm[static_cast<size_t>(k)]);
            return acc;
        }
    }

    void ThinFilm::buildRgbTable() {
        Real dMax = m_thickness.base;
        if (m_thickness.gradient) dMax = std::max(dMax, m_thickness.top);
        if (m_thickness.noise) dMax += m_thickness.noiseAmplitude;
        dMax = std::min(std::max(dMax, Real(0)) + kTableStep, kTableMax);
        m_tableD = static_cast<int>(std::ceil(dMax / kTableStep)) + 2;
        const int sides = (m_substrate == Substrate::Dielectric) ? 2 : 1;
        for (int side = 0; side < sides; ++side) {
            m_table[side].assign(static_cast<size_t>(kTableCos) * static_cast<size_t>(m_tableD), Vector3(0.0));
            for (int ic = 0; ic < kTableCos; ++ic) {
                const Real cosI = Real(ic) / Real(kTableCos - 1);
                for (int id = 0; id < m_tableD; ++id) {
                    const Real d = Real(id) * kTableStep;
                    m_table[side][static_cast<size_t>(ic) * static_cast<size_t>(m_tableD) + static_cast<size_t>(id)] =
                        integrateRgb(*this, [&](Real nm) { return reflectanceAt(cosI, side == 0, d, nm); });
                }
            }
        }
    }

    Vector3 ThinFilm::rgbReflectance(Real cosI, bool entering, Real d) const {
        if (d >= Real(m_tableD - 1) * kTableStep) {   // beyond the table: integrate directly (rare, slow)
            return integrateRgb(*this, [&](Real nm) { return reflectanceAt(cosI, entering, d, nm); });
        }
        const std::vector<Vector3>& tab = m_table[(m_substrate == Substrate::Dielectric && !entering) ? 1 : 0];
        const Real x = std::clamp(cosI, Real(0), Real(1)) * Real(kTableCos - 1);
        const int i = std::min(static_cast<int>(x), kTableCos - 2);
        const Real fx = x - Real(i);
        const Real y = std::max(d, Real(0)) / kTableStep;
        const int j = std::min(static_cast<int>(y), m_tableD - 2);
        const Real fy = y - Real(j);
        auto at = [&](int a, int b) { return tab[static_cast<size_t>(a) * static_cast<size_t>(m_tableD) + static_cast<size_t>(b)]; };
        return (at(i, j) * (Real(1) - fy) + at(i, j + 1) * fy) * (Real(1) - fx) + (at(i + 1, j) * (Real(1) - fy) + at(i + 1, j + 1) * fy) * fx;
    }
#endif

    void ThinFilm::reflectTransmit(const Point3& p, Real cosI, bool entering, const SampledWavelengths& lambda, Spectrum& R, Spectrum& T) const {
        const Real d = m_thickness.at(p);
        R = Spectrum(0.0);
#if defined(RAYT_SPECTRAL_MODE)
        for (int i = 0; i < kLanes; ++i) R[i] = reflectanceAt(cosI, entering, d, lambda[i]);
#else
        (void)lambda;
        const Vector3 acc = rgbReflectance(cosI, entering, d);
        R = Spectrum(acc.x, acc.y, acc.z);
#endif
        T = Spectrum(0.0);
        if (m_substrate != Substrate::Conductor) {
            for (int i = 0; i < kLanes; ++i) T[i] = std::max(Real(0), Real(1) - R[i]);   // the film does not absorb
        }
    }

    std::optional<BSDFSample> ThinFilm::sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u,
        const SampledWavelengths& lambda, TransportMode mode) const {
        const Real cosI = std::abs(glm::dot(rec.n, wo));
        const bool entering = rec.frontFace;

        Spectrum R, T;
        reflectTransmit(rec.p, cosI, entering, lambda, R, T);

        // One random choice for all wavelengths of the path: reflect with the mean reflectance, and give every wavelength its own weight
        // R_i / pR (reflection) or T_i / (1 - pR) (transmission). The expectation is R_i and T_i for each wavelength.
        Real pR = 0;
        for (int i = 0; i < kLanes; ++i) pR += R[i];
        pR /= Real(kLanes);
        if (m_substrate == Substrate::Conductor) pR = Real(1);
        pR = std::clamp(pR, Real(0), Real(1));

        BSDFSample s;
        if (u.x < pR) {
            s.wi = math::reflectOutward(wo, rec.n);
            s.f = R / pR;
            s.pdf = pR;
            s.flags = BxDFFlags::Specular | BxDFFlags::Reflection;
            return s;
        }

        // transmission: straight through a free film, refracted by the substrate otherwise
        const Real n2 = (m_substrate == Substrate::Dielectric) ? m_n2 : Real(1);
        const Real eta = entering ? Real(1) / n2 : n2;     // eta_i / eta_t
        const Real etap = entering ? n2 : Real(1) / n2;    // eta_t / eta_i
        Vector3 wi;
        if (!math::refractOutward(wo, rec.n, eta, wi)) return std::nullopt;   // total internal reflection: R is 1 there, so this is rare
        s.wi = wi;
        s.pdf = Real(1) - pR;
        s.f = T / (Real(1) - pR);
        if (mode == TransportMode::Radiance) s.f /= etap * etap;               // radiance changes with the index (as in Dielectric)
        s.flags = BxDFFlags::Specular | BxDFFlags::Transmission;
        return s;
    }

} // namespace rayt
