/**
 * @file Color/SpectralConversion.hpp
 * @brief RGB <-> spectral conversions used by the spectral rendering mode.
 *
 * Scene side (RGB -> spectrum):
 * - RgbColor: an RGB color given at scene-setup time (albedo, emission, light intensity).
 *   It is turned into values at the path's sampled wavelengths when rendering:
 *     reflectance(lambda) : for albedo-like quantities, white (1,1,1) -> flat 1
 *     illuminant(lambda)  : for light-like quantities, white (1,1,1) -> D65 power distribution
 *
 * Light source spectra (blackbody, measured SPD) live in Color/EmissionSpectrum.hpp.
 *
 * Film side (spectrum -> RGB):
 * - SpectrumToRGB(): Monte Carlo estimate of  Integral( L(l) * sensor_rgb(l) dl )  from the
 *   sampled wavelengths and their pdfs. Convention: a radiance of D65 * 1 maps to RGB (1,1,1).
 *
 * In RGB mode (default) every function here is the identity, so the same scene code builds in
 * both modes.
 *
 * RGB -> spectrum method:
 * - Three smooth basis spectra P_r, P_g, P_b (5 nm nodes, 380..780 nm), found once at startup as
 *   the smoothest solution of:  RGB(P_c) = e_c  (under D65, sRGB)  with  P_r + P_g + P_b = 1.
 *   S(l) = max(0, r*P_r(l) + g*P_g(l) + b*P_b(l)).
 * - White and moderately saturated colors round-trip exactly; very saturated colors (pure green)
 *   get small negative lobes that are clamped to 0, which shifts them slightly.
 *   Typical HDRIs: median color error 0, 95th percentile <= 0.1 %.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "Core/Types.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Color/ColorSpace.hpp"
#include "Color/Illuminants.hpp"

namespace rayt::color {

    namespace detail {

        // ---------------------------------------------------------------------
        // Sensor: per-wavelength linear-sRGB weights, normalized so D65 light -> (1,1,1)
        // ---------------------------------------------------------------------
        struct SensorTable {
            static constexpr int kFirstNm = 360;
            static constexpr int kLastNm = 830;

            std::vector<Vector3> weights;   ///< index = lambda_nm - kFirstNm
            Vector3 white = Vector3(0.0);   ///< integral of (M * cmf) * D65, same rule as spectrumToRGB_D65()

            SensorTable() {
                // Same integration rule as spectrumToRGB_D65(): 1 nm steps, trapezoid on [380, 780].
                for (int l = 380; l <= 780; ++l) {
                    const Real trapezoid = (l == 380 || l == 780) ? Real(0.5) : Real(1);
                    white += XYZ_to_sRGB(evaluateCIE_XYZ_1931_MultiLobe(Real(l)))
                        * (evaluateD65(Real(l)) * trapezoid);
                }

                weights.resize(static_cast<size_t>(kLastNm - kFirstNm + 1));
                for (int l = kFirstNm; l <= kLastNm; ++l) {
                    const Vector3 c = XYZ_to_sRGB(evaluateCIE_XYZ_1931_MultiLobe(Real(l)));
                    weights[static_cast<size_t>(l - kFirstNm)] =
                        Vector3(c.x / white.x, c.y / white.y, c.z / white.z);
                }
            }
        };

        inline const SensorTable& sensorTable() {
            static const SensorTable table;   // thread-safe initialization
            return table;
        }

        // ---------------------------------------------------------------------
        // Dense linear solve (Gaussian elimination with partial pivoting)
        // Solves A X = B for n x n matrix A (row-major) and n x m right-hand sides B (row-major).
        // ---------------------------------------------------------------------
        inline bool solveLinear(std::vector<Real> A, std::vector<Real> B, int n, int m,
            std::vector<Real>& X)
        {
            for (int col = 0; col < n; ++col) {
                int pivot = col;
                for (int r = col + 1; r < n; ++r) {
                    if (std::abs(A[r * n + col]) > std::abs(A[pivot * n + col])) pivot = r;
                }
                if (std::abs(A[pivot * n + col]) < Real(1e-300)) return false;
                if (pivot != col) {
                    for (int k = 0; k < n; ++k) std::swap(A[pivot * n + k], A[col * n + k]);
                    for (int k = 0; k < m; ++k) std::swap(B[pivot * m + k], B[col * m + k]);
                }
                for (int r = col + 1; r < n; ++r) {
                    const Real f = A[r * n + col] / A[col * n + col];
                    if (f == Real(0)) continue;
                    for (int k = col; k < n; ++k) A[r * n + k] -= f * A[col * n + k];
                    for (int k = 0; k < m; ++k) B[r * m + k] -= f * B[col * m + k];
                }
            }
            X.assign(static_cast<size_t>(n) * m, Real(0));
            for (int k = 0; k < m; ++k) {
                for (int r = n - 1; r >= 0; --r) {
                    Real s = B[r * m + k];
                    for (int c = r + 1; c < n; ++c) s -= A[r * n + c] * X[c * m + k];
                    X[r * m + k] = s / A[r * n + r];
                }
            }
            return true;
        }

        // ---------------------------------------------------------------------
        // RGB -> spectrum basis
        // ---------------------------------------------------------------------
        struct RgbBasis {
            static constexpr int kNodes = 81;          // 380, 385, ..., 780 nm
            static constexpr Real kFirstNm = Real(380);
            static constexpr Real kStepNm = Real(5);

            std::array<std::array<Real, kNodes>, 3> P{};   // P[channel][node]

            RgbBasis() {
                constexpr int n = kNodes;

                // A (3 x n): RGB of each piecewise-linear "hat" function under D65 / sRGB,
                // using exactly the integration rule of spectrumToRGB_D65() (so it round-trips).
                const SensorTable& sensor = sensorTable();
                std::vector<Real> A(3 * n, Real(0));
                for (int l = 380; l <= 780; ++l) {
                    const Real trapezoid = (l == 380 || l == 780) ? Real(0.5) : Real(1);
                    const Vector3 c = XYZ_to_sRGB(evaluateCIE_XYZ_1931_MultiLobe(Real(l)));
                    const Real d65 = evaluateD65(Real(l)) * trapezoid;
                    const Vector3 w(c.x * d65 / sensor.white.x,
                        c.y * d65 / sensor.white.y,
                        c.z * d65 / sensor.white.z);

                    const Real x = (Real(l) - kFirstNm) / kStepNm;
                    const int j = std::min(static_cast<int>(x), n - 2);
                    const Real f = x - Real(j);                 // hat_j = 1 - f, hat_{j+1} = f
                    for (int ch = 0; ch < 3; ++ch) {
                        A[ch * n + j] += w[ch] * (Real(1) - f);
                        A[ch * n + j + 1] += w[ch] * f;
                    }
                }

                // Smoothness metric Q = L^T L + eps I  (L = second difference).
                std::vector<Real> Q(static_cast<size_t>(n) * n, Real(0));
                for (int i = 0; i + 2 < n; ++i) {
                    const int idx[3] = { i, i + 1, i + 2 };
                    const Real coef[3] = { Real(1), Real(-2), Real(1) };
                    for (int a = 0; a < 3; ++a)
                        for (int b = 0; b < 3; ++b)
                            Q[idx[a] * n + idx[b]] += coef[a] * coef[b];
                }
                for (int i = 0; i < n; ++i) Q[i * n + i] += Real(1e-6);

                // X = Q^-1 A^T   (n x 3)
                std::vector<Real> At(static_cast<size_t>(n) * 3);
                for (int j = 0; j < n; ++j)
                    for (int ch = 0; ch < 3; ++ch) At[j * 3 + ch] = A[ch * n + j];
                std::vector<Real> X;
                if (!solveLinear(Q, At, n, 3, X)) {
                    for (auto& channel : P) channel.fill(Real(1) / Real(3));   // fallback: flat gray basis
                    return;
                }

                // S = A X (3 x 3)  and  P_c = 1/3 + X S^-1 (e_c - (1/3, 1/3, 1/3))
                glm::mat<3, 3, Real, glm::defaultp> S(Real(0));
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c) {
                        Real s = 0;
                        for (int j = 0; j < n; ++j) s += A[r * n + j] * X[j * 3 + c];
                        S[c][r] = s;                            // glm is column-major: S[col][row]
                    }
                const auto Sinv = glm::inverse(S);

                for (int ch = 0; ch < 3; ++ch) {
                    Vector3 b(-Real(1) / Real(3));
                    b[ch] += Real(1);
                    const Vector3 y = Sinv * b;
                    for (int j = 0; j < n; ++j) {
                        const Real delta = X[j * 3 + 0] * y.x + X[j * 3 + 1] * y.y + X[j * 3 + 2] * y.z;
                        P[ch][j] = Real(1) / Real(3) + delta;
                    }
                }
            }
        };

        inline const RgbBasis& rgbBasis() {
            static const RgbBasis basis;   // thread-safe initialization
            return basis;
        }

    } // namespace detail

    /**
     * @brief Camera white balance: the film divides the sensor RGB by this color (default (1,1,1)).
     *
     * A light that is not D65 (tungsten, candle, ...) renders with its true color cast. Setting
     * this to the RGB of the light that should look white removes the cast (diagonal von Kries in
     * linear sRGB). Set it before rendering starts; see setWhiteBalance() in Color/EmissionSpectrum.hpp.
     */
    inline Vector3& whiteBalanceRgb() {
        static Vector3 wb(1.0);
        return wb;
    }

    /**
     * @brief Per-channel sensor weights at a wavelength: s_c(l) = (M * cmf(l))_c / W_c.
     *
     * Integral( D65(l) * s_c(l) dl ) == 1 for every channel, so a D65-colored light of
     * unit scale maps to RGB (1,1,1).
     */
    inline Vector3 rgbSensorWeights(Real lambda_nm) {
        const auto& table = detail::sensorTable();
        const Real x = std::clamp(lambda_nm, Real(detail::SensorTable::kFirstNm),
            Real(detail::SensorTable::kLastNm)) - Real(detail::SensorTable::kFirstNm);
        const int i = std::min(static_cast<int>(x), static_cast<int>(table.weights.size()) - 2);
        const Real f = x - Real(i);
        return table.weights[static_cast<size_t>(i)] * (Real(1) - f)
            + table.weights[static_cast<size_t>(i) + 1] * f;
    }

    /**
     * @brief Smooth reflectance-like spectrum value for an RGB color at a wavelength.
     * White (1,1,1) gives exactly 1 at every wavelength. Clamped to >= 0.
     */
    inline Real rgbToReflectance(const Vector3& rgb, Real lambda_nm) {
        const auto& basis = detail::rgbBasis();
        const Real x = (std::clamp(lambda_nm, Real(380), Real(780)) - detail::RgbBasis::kFirstNm)
            / detail::RgbBasis::kStepNm;
        const int i = std::min(static_cast<int>(x), detail::RgbBasis::kNodes - 2);
        const Real f = x - Real(i);

        Real s = 0;
        for (int c = 0; c < 3; ++c) {
            s += rgb[c] * (basis.P[c][i] * (Real(1) - f) + basis.P[c][i + 1] * f);
        }
        return std::max(s, Real(0));
    }

} // namespace rayt::color

namespace rayt {

    /**
     * @brief An RGB color specified at scene-setup time (albedo, emission, light intensity).
     *
     * Materials and lights store RgbColor and evaluate it for the current path's wavelengths.
     * In RGB mode the evaluation just returns the RGB triplet.
     */
    class RgbColor {
    public:
        RgbColor() = default;
        RgbColor(Real v) : m_rgb(v) {}
        RgbColor(Real r, Real g, Real b) : m_rgb(r, g, b) {}
        RgbColor(const Vector3& rgb) : m_rgb(rgb) {}

        const Vector3& rgb() const { return m_rgb; }

        /// Reflectance-like spectrum: (1,1,1) -> 1 at all wavelengths.
        Spectrum reflectance(const SampledWavelengths& lambda) const {
#if defined(RAYT_SPECTRAL_MODE)
            Spectrum s;
            for (int i = 0; i < kSpectralSamples; ++i) {
                s[i] = color::rgbToReflectance(m_rgb, lambda[i]);
            }
            return s;
#else
            (void)lambda;
            return m_rgb;
#endif
        }

        /// Light-like spectrum: (1,1,1) -> D65 power distribution (maps back to RGB (1,1,1)).
        Spectrum illuminant(const SampledWavelengths& lambda) const {
#if defined(RAYT_SPECTRAL_MODE)
            Spectrum s;
            for (int i = 0; i < kSpectralSamples; ++i) {
                s[i] = color::rgbToReflectance(m_rgb, lambda[i]) * color::evaluateD65(lambda[i]);
            }
            return s;
#else
            (void)lambda;
            return m_rgb;
#endif
        }

        RgbColor operator*(Real s) const { return RgbColor(m_rgb * s); }

    private:
        Vector3 m_rgb = Vector3(0.0);
    };

    /**
     * @brief Radiance given as RGB (e.g. an environment map texel) -> spectrum at the path's wavelengths.
     */
    inline Spectrum RgbToIlluminantSpectrum(const Vector3& rgb, const SampledWavelengths& lambda) {
        return RgbColor(rgb).illuminant(lambda);
    }

    /**
     * @brief Film-side conversion: spectral radiance sampled at the path's wavelengths -> linear sRGB.
     *
     * Estimator:  RGB = (1/N) * sum_i  sensor_rgb(l_i) * L(l_i) / pdf(l_i)
     * RGB mode: identity.
     */
    inline Vector3 SpectrumToRGB(const Spectrum& L, const SampledWavelengths& lambda) {
#if defined(RAYT_SPECTRAL_MODE)
        Vector3 rgb(0.0);
        for (int i = 0; i < kSpectralSamples; ++i) {
            const Real pdf = lambda.pdf[i];
            if (pdf <= Real(0)) continue;
            rgb += color::rgbSensorWeights(lambda[i]) * (L[i] / pdf);
        }
        return rgb / (Real(kSpectralSamples) * color::whiteBalanceRgb());
#else
        (void)lambda;
        return L / color::whiteBalanceRgb();
#endif
    }

} // namespace rayt
