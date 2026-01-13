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

    // Number of wavelength samples (used only in Spectral mode)
    constexpr int kSpectralSamples = 4;  // PBRT v4 uses 4 for Hero wavelength sampling

    // ========================================================================
    // Wavelength Sampling (Forward Declaration)
    // ========================================================================

    /**
     * @brief Represents a set of sampled wavelengths for spectral rendering.
     *
     * In spectral mode, each path traces specific wavelengths chosen either:
     * - Uniformly across the visible spectrum, or
     * - Using "hero wavelength" sampling (one primary + correlated secondaries)
     */
    class SampledWavelengths {
    public:
        // Wavelengths in nanometers [360, 830]
        std::array<Real, kSpectralSamples> lambda;

        // PDF for sampling these wavelengths (for MIS)
        std::array<Real, kSpectralSamples> pdf;

        SampledWavelengths() {
            lambda.fill(Real(550));  // Default: green (mid-spectrum)
            pdf.fill(Real(1) / (constants::LAMBDA_MAX - constants::LAMBDA_MIN));
        }

        /**
         * @brief Sample wavelengths uniformly across visible spectrum.
         * @param u Random value in [0,1)
         */
        static SampledWavelengths sampleUniform(Real u);

        /**
         * @brief Sample using "hero wavelength" strategy (PBRT v4 approach).
         */
        static SampledWavelengths sampleHero(Real u);

        /**
         * @brief Get wavelength at index i.
         */
        Real operator[](int i) const { return lambda[i]; }
    };

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

        // èzä¬éQè∆Çñhé~Ç∑ÇÈÇΩÇﬂ
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
                // Ç±Ç±ÇÕï˚êjÇ≈ëIëÅF
                // A) çïÇ…í◊Ç∑Åià¿ëSÅj
                c.fill(Real(0));
                return *this;
                // B) âΩÇ‡ÇµÇ»Ç¢Ç≈ï‘Ç∑ÅiÇÊÇËÅgê√Ç©ÅhÇæÇ™åÎñÇâªÇ∑Åj
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

        /**
         * @brief Convert sampled spectrum to RGB.
         *
         * TODO: Implement proper spectral to RGB conversion using CIE color matching
         * For now, use average as grayscale
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
#else
    using Spectrum = RGBSpectrum;
#endif

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

    // ========================================================================
    // SampledWavelengths Implementation (out-of-line to avoid circular deps)
    // ========================================================================

    // ========================================================================
    // Helper: Simple linear interpolation (avoids Math.hpp dependency)
    // ========================================================================

    namespace detail {
        template<typename T>
        inline T lerp(const T& a, const T& b, Real t) {
            return a + t * (b - a);
        }
    }

    // ========================================================================
    // SampledWavelengths Implementation
    // ========================================================================

    inline SampledWavelengths SampledWavelengths::sampleUniform(Real u) {
        SampledWavelengths swl;

        for (int i = 0; i < kSpectralSamples; ++i) {
            Real offset = (u + Real(i)) / Real(kSpectralSamples);
            offset = offset - std::floor(offset);  // Wrap to [0,1)
            swl.lambda[i] = detail::lerp(constants::LAMBDA_MIN, constants::LAMBDA_MAX, offset);
            swl.pdf[i] = Real(1) / (constants::LAMBDA_MAX - constants::LAMBDA_MIN);
        }

        return swl;
    }

    inline SampledWavelengths SampledWavelengths::sampleHero(Real u) {
        SampledWavelengths swl;

        // Hero wavelength (uniformly sampled)
        swl.lambda[0] = detail::lerp(constants::LAMBDA_MIN, constants::LAMBDA_MAX, u);
        swl.pdf[0] = Real(1) / (constants::LAMBDA_MAX - constants::LAMBDA_MIN);

        // Secondary wavelengths (evenly spaced around hero)
        Real delta = (constants::LAMBDA_MAX - constants::LAMBDA_MIN) / Real(kSpectralSamples);
        for (int i = 1; i < kSpectralSamples; ++i) {
            swl.lambda[i] = swl.lambda[0] + Real(i) * delta;

            // Wrap around spectrum (handle overflow)
            if (swl.lambda[i] > constants::LAMBDA_MAX) {
                swl.lambda[i] -= (constants::LAMBDA_MAX - constants::LAMBDA_MIN);
            }

            swl.pdf[i] = swl.pdf[0];
        }

        return swl;
    }

} // namespace rayt