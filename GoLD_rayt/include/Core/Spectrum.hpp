/**
 * @file Core/Spectrum.hpp
 * @brief Unified spectrum representation supporting both RGB and spectral modes.
 *
 * This header provides a compile-time switchable spectrum implementation:
 * - RGB mode: Fast 3-channel rendering (production, default)
 * - Spectral mode: Physically accurate wavelength sampling (research/quality)
 *
 * Design inspired by PBRT v4's SampledSpectrum and Mitsuba 3's spectrum system.
 *
 * @note To enable spectral mode, define RAYT_SPECTRAL_MODE before including this header
 *       or add -DRAYT_SPECTRAL_MODE to compiler flags.
 */

#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <optional>
#include <glm/glm.hpp>

#include "Core/Types.hpp"
#include "Core/Forward.hpp"
#include "Core/Constants.hpp"
#include "Core/Assert.hpp"
#include "Core/SampledWavelengths.hpp"

namespace rayt {

    // ========================================================================
    // Configuration
    // ========================================================================

    enum class SpectrumMode {
        RGB,        ///< Fast RGB rendering
        Spectral    ///< Physically accurate spectral rendering
    };

    // Global spectrum mode (compile-time constant)
#if defined(RAYT_SPECTRAL_MODE)
    constexpr SpectrumMode kSpectrumMode = SpectrumMode::Spectral;
#else
    constexpr SpectrumMode kSpectrumMode = SpectrumMode::RGB;
#endif

    // ========================================================================
    // RGB Spectrum (Fast Mode - Default)
    // ========================================================================

    class RGBSpectrum {
    public:
        Vector3 v;

        // Constructors
        RGBSpectrum(Real x = 0) : v(x, x, x) {}
        RGBSpectrum(Real r, Real g, Real b) : v(r, g, b) {}
        explicit RGBSpectrum(const Vector3& rgb) : v(rgb) {}

        // Kept self-contained to avoid circular includes
        // Query
        bool isBlack() const {
            return v.x <= Real(0) && v.y <= Real(0) && v.z <= Real(0);
        }

        bool hasNaNs() const {
            return !std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z);
        }

        Real maxComponent() const {
            return std::max({ v.x, v.y, v.z });
        }

        // Component access (0=R, 1=G, 2=B). Same API as SampledSpectrumT.
        Real& operator[](int i) { return v[i]; }
        Real operator[](int i) const { return v[i]; }

        // Operators
        RGBSpectrum& operator+=(const RGBSpectrum& o) { v += o.v; return *this; }
        RGBSpectrum& operator*=(const RGBSpectrum& o) { v *= o.v; return *this; }
        RGBSpectrum& operator*=(Real s) { v *= s; return *this; }
        //RGBSpectrum& operator/=(Real s) { v /= s; return *this; }

        RGBSpectrum& operator/=(Real s) {
#ifndef NDEBUG
            Assert(std::isfinite(s));
            Assert(s != Real(0));
#endif
            if (!std::isfinite(s) || s == Real(0)) [[unlikely]] {
                v = Vector3(0);
                return *this;
                }
            v /= s;
            return *this;
        }

        friend RGBSpectrum operator+(RGBSpectrum a, const RGBSpectrum& b) {
            return a += b;
        }
        friend RGBSpectrum operator-(RGBSpectrum a, const RGBSpectrum& b) {
            a.v -= b.v; return a;
        }
        friend RGBSpectrum operator*(RGBSpectrum a, const RGBSpectrum& b) {
            return a *= b;
        }
        friend RGBSpectrum operator*(RGBSpectrum a, Real s) {
            return a *= s;
        }
        friend RGBSpectrum operator*(Real s, RGBSpectrum a) {
            return a *= s;
        }
        friend RGBSpectrum operator/(RGBSpectrum a, Real s) {
            return a /= s;
        }

        // Conversion to rendering format
        Vector3 toRGB() const {
            return v;
        }

        // Sanitize for robustness
        RGBSpectrum sanitize() const {
#ifndef NDEBUG
            Assert(!hasNaNs());
#endif
            if (hasNaNs()) [[unlikely]] return RGBSpectrum(0.0);
            return RGBSpectrum(glm::max(v, Vector3(0.0)));
        }

        // API consistency with spectral mode
        static RGBSpectrum fromWavelengths(const SampledWavelengths& /*lambda*/) {
            return RGBSpectrum(1.0);  // White
        }
    };

    // ========================================================================
    // Sampled Spectrum (Spectral Mode)
    // ========================================================================

    template<int N>
    class SampledSpectrumT {
    public:
        std::array<Real, N> c{};

        // Constructors
        SampledSpectrumT(Real x = 0) { c.fill(x); }

        // Query
        bool isBlack() const {
            for (Real x : c) if (x > Real(0)) return false;
            return true;
        }

        bool hasNaNs() const {
            for (Real x : c) {
                if (!std::isfinite(x)) return true;
            }
            return false;
        }

        Real maxComponent() const {
            return *std::max_element(c.begin(), c.end());
        }

        // Component access (one value per sampled wavelength). Same API as RGBSpectrum.
        Real& operator[](int i) { return c[i]; }
        Real operator[](int i) const { return c[i]; }

        // Operators
        SampledSpectrumT& operator+=(const SampledSpectrumT& o) {
            for (int i = 0; i < N; ++i) c[i] += o.c[i];
            return *this;
        }

        SampledSpectrumT& operator*=(const SampledSpectrumT& o) {
            for (int i = 0; i < N; ++i) c[i] *= o.c[i];
            return *this;
        }

        SampledSpectrumT& operator*=(Real s) {
            for (int i = 0; i < N; ++i) c[i] *= s;
            return *this;
        }

        /*SampledSpectrumT& operator/=(Real s) {
            Real inv = Real(1) / s;
            return (*this) *= inv;
        }*/
        SampledSpectrumT& operator/=(Real s) {
#ifndef NDEBUG
            Assert(std::isfinite(s));
            Assert(s != Real(0));
#endif
            if (!std::isfinite(s) || s == Real(0)) [[unlikely]] {
                // Policy choice:
                // A) Collapse to black (safe)
                c.fill(Real(0));
                return *this;
                // B) Return unchanged (quieter, but it hides the problem)
                // return *this;
                }
            Real inv = Real(1) / s;
            for (int i = 0; i < N; ++i) c[i] *= inv;
            return *this;
        }

        friend SampledSpectrumT operator+(SampledSpectrumT a, const SampledSpectrumT& b) {
            return a += b;
        }
        friend SampledSpectrumT operator-(SampledSpectrumT a, const SampledSpectrumT& b) {
            for (int i = 0; i < N; ++i) a.c[i] -= b.c[i];
            return a;
        }
        friend SampledSpectrumT operator*(SampledSpectrumT a, const SampledSpectrumT& b) {
            return a *= b;
        }
        friend SampledSpectrumT operator*(SampledSpectrumT a, Real s) {
            return a *= s;
        }
        friend SampledSpectrumT operator*(Real s, SampledSpectrumT a) {
            return a *= s;
        }
        friend SampledSpectrumT operator/(SampledSpectrumT a, Real s) {
            return a /= s;
        }

        // Element-wise division (a lane that divides by zero becomes 0).
        SampledSpectrumT& operator/=(const SampledSpectrumT& o) {
            for (int i = 0; i < N; ++i) c[i] = (o.c[i] != Real(0)) ? c[i] / o.c[i] : Real(0);
            return *this;
        }
        friend SampledSpectrumT operator/(SampledSpectrumT a, const SampledSpectrumT& b) {
            return a /= b;
        }

        Real average() const {
            Real sum = 0;
            for (Real x : c) sum += x;
            return sum / Real(N);
        }

        /**
         * @brief Placeholder conversion (grayscale average). Do NOT use for rendering.
         *
         * The real conversion needs the sampled wavelengths and their pdfs:
         * use SpectrumToRGB(L, lambda) from Color/SpectralConversion.hpp.
         */
        Vector3 toRGB() const {
            Real sum = 0;
            for (Real x : c) sum += x;
            Real avg = sum / Real(N);
            return Vector3(avg, avg, avg);
        }

        // Sanitize for robustness
        SampledSpectrumT sanitize() const {
#ifndef NDEBUG
            Assert(!hasNaNs());
#endif
            if (hasNaNs()) [[unlikely]] {
                return SampledSpectrumT(0);
                }
            SampledSpectrumT result;
            for (int i = 0; i < N; ++i) {
                result.c[i] = std::max(c[i], Real(0));
            }
            return result;
        }

        /**
         * @brief Create spectrum from wavelength-dependent function.
         * TODO: Implement actual spectral evaluation
         */
        static SampledSpectrumT fromWavelengths(const SampledWavelengths& lambda) {
            SampledSpectrumT s;
            for (int i = 0; i < N; ++i) {
                s.c[i] = Real(1);  // Placeholder: constant SPD
            }
            return s;
        }
    };

    // ========================================================================
    // Public Spectrum Type (Mode-Dependent)
    // ========================================================================

#if defined(RAYT_SPECTRAL_MODE)
    using Spectrum = SampledSpectrumT<kSpectralSamples>;
#endif
    // NOTE: In RGB mode (default) `Spectrum` is the glm::vec3 alias in Core/Types.hpp and the
    //       helpers (isBlack / HasInvalidValues / Sanitize) come from Core/SpectrumUtils.hpp.
    //       RGBSpectrum is compiled but not aliased yet (candidate for a later migration).

#if defined(RAYT_SPECTRAL_MODE)
    // ========================================================================
    // Utility Functions (Backwards Compatible with Existing Code)
    // ========================================================================

    /**
     * @brief Checks if the spectrum contributes no energy to the scene.
     *
     * @note Maintains compatibility with existing SpectrumUtil.hpp
     */
    inline bool isBlack(const Spectrum& s) {
        return s.isBlack();
    }

    /**
     * @brief Checks for invalid numerical states (NaN or Infinity).
     *
     * @note Maintains compatibility with existing SpectrumUtil.hpp
     */
    inline bool HasInvalidValues(const Spectrum& s) {
        return s.hasNaNs();
    }

    /**
     * @brief Ensures a spectrum is physically valid and numerically safe.
     *
     * Clamps negative values to zero and recovers from NaNs/Infs by returning black.
     *
     * @note Maintains compatibility with existing SpectrumUtil.hpp
     */
    inline Spectrum Sanitize(const Spectrum& s) {
        return s.sanitize();
    }

    // Overloads so generic code (Russian roulette, NaN checks) works for both Spectrum types.
    namespace math {
        template<int N>
        inline Real maxComponent(const SampledSpectrumT<N>& s) { return s.maxComponent(); }

        template<int N>
        inline bool hasNonFinite(const SampledSpectrumT<N>& s) { return s.hasNaNs(); }
    }
#endif // RAYT_SPECTRAL_MODE

    // NOTE: detail::lerp and SampledWavelengths live in Core/SampledWavelengths.hpp.

} // namespace rayt