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