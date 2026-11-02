/**
 * @file Core/SampledWavelengths.hpp
 * @brief Wavelength sampling container and utilities for spectral rendering.
 */

#pragma once

#include <array>
#include <cmath>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"

namespace rayt {

    // Number of wavelength samples (used only in Spectral mode)
    constexpr int kSpectralSamples = 4;  // PBRT v4 uses 4 for Hero wavelength sampling

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
         * @brief Sample wavelengths with a density that follows the eye's sensitivity.
         *
         * Wavelengths the eye barely sees (deep red / violet) are sampled less often, which
         * greatly reduces color noise compared to uniform sampling. The four wavelengths are
         * stratified (u, u+1/4, u+2/4, u+3/4 wrapped), each with its own pdf.
         * @param u Random value in [0,1)
         */
        static SampledWavelengths sampleVisible(Real u);

        /**
         * @brief Get wavelength at index i.
         */
        Real operator[](int i) const { return lambda[i]; }
    };

    namespace detail {
        template<typename T>
        inline T lerp(const T& a, const T& b, Real t) {
            return a + t * (b - a);
        }
    }

    /**
     * @brief Density used by sampleVisible(), normalized on [360, 830] nm (0 outside).
     * Closed form from PBRT v4: p(l) = 0.0039398042 / cosh^2(0.0072 * (l - 538)).
     */
    inline Real visibleWavelengthPdf(Real lambda_nm) {
        if (lambda_nm < Real(360) || lambda_nm > Real(830)) return Real(0);
        const Real c = std::cosh(Real(0.0072) * (lambda_nm - Real(538)));
        return Real(0.0039398042) / (c * c);
    }

    /**
     * @brief Inverse CDF of visibleWavelengthPdf(): maps u in [0,1] to a wavelength in [360, 830] nm.
     */
    inline Real sampleVisibleWavelength(Real u) {
        return Real(538) - Real(138.888889) * std::atanh(Real(0.85691062) - Real(1.82750197) * u);
    }

    inline SampledWavelengths SampledWavelengths::sampleVisible(Real u) {
        SampledWavelengths swl;

        for (int i = 0; i < kSpectralSamples; ++i) {
            Real up = u + Real(i) / Real(kSpectralSamples);
            if (up >= Real(1)) up -= Real(1);   // wrap into [0,1)
            swl.lambda[i] = sampleVisibleWavelength(up);
            swl.pdf[i] = visibleWavelengthPdf(swl.lambda[i]);
        }

        return swl;
    }

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