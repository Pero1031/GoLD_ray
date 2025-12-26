#pragma once

/**
 * @file Distribution1D.hpp
 * @brief 1D Probability Distribution for Importance Sampling.
 *
 * This file defines the Distribution1D class, which implements the Inverse Transform Sampling method.
 * It constructs a Cumulative Distribution Function (CDF) from a discrete data array (e.g., texture luminance)
 * and allows for continuous sampling proportional to the data values.
 * This is essential for Monte Carlo integration, specifically for environment map importance sampling.
 */

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>

namespace rayt {

    /**
     * @struct Distribution1D
     * @brief A helper class for 1D importance sampling.
     * * Handles the creation of CDF from a discrete function and provides
     * inverse transform sampling methods. The domain is assumed to be [0, 1].
     */
    struct Distribution1D {
        std::vector<float> func; ///< The piecewise constant function values (PDF * integral).
        std::vector<float> cdf;  ///< The Cumulative Distribution Function.
        float funcInt;           ///< The integral of the function over [0, 1].

        /**
         * @brief Constructs the distribution from a data array.
         * @param f Pointer to the data array.
         * @param n Number of elements.
         */
        Distribution1D(const float* f, int n) : func(f, f + n) {
            cdf.resize(n + 1);
            cdf[0] = 0;

            // Compute integral (CDF)
            // Assumes the domain is [0, 1], so the width of each bin is 1/n.
            for (int i = 0; i < n; ++i) {
                cdf[i + 1] = cdf[i] + func[i] / n; 
            }

            funcInt = cdf[n]; 

            // Normalize CDF
            if (funcInt == 0) {
                for (int i = 1; i < n + 1; ++i) cdf[i] = float(i) / float(n);
            }
            else {
                for (int i = 1; i < n + 1; ++i) cdf[i] /= funcInt;
            }
        }

        /**
         * @brief Samples the distribution using the inverse CDF method.
         * * @param u Uniform random number in [0, 1).
         * @param[out] pdf The probability density at the sampled location.
         * @param[out] off The integer index of the sampled bin.
         * @return The continuous sampled coordinate in [0, 1].
         */
        float sampleContinuous(float u, float& pdf, int& off) const { // FIX: Changed return type from int to float
            int n = static_cast<int>(func.size());

            // Handle edge case u=1.0 slightly inside to avoid out of bounds
            u = std::min(u, 0.99999994f);

            // Find the index such that cdf[offset] <= u < cdf[offset+1]
            auto ptr = std::upper_bound(cdf.begin(), cdf.end(), u);
            int offset = std::clamp(static_cast<int>(ptr - cdf.begin()) - 1, 0, n - 1);

            off = offset; // Export the index

            // Compute the PDF
            // p(x) = f(x) / Integral(f)
            if (funcInt > 0.0f) {
                pdf = func[offset] / funcInt;
            }
            else {
                pdf = 1.0f; // Uniform distribution case
            }

            // Compute the offset within the bin
            // du becomes 0..1 representing position inside the bin
            float du = u - cdf[offset];
            float denom = cdf[offset + 1] - cdf[offset];

            if (denom > 0.0f) {
                du /= denom;
            }
            else {
                // Should not happen if normalized correctly, but fail-safe
                du = 0.0f;
            }

            // Map back to [0, 1] domain
            return (offset + du) / float(n);
        }

        /**
         * @brief Returns the number of elements.
         */
        int count() const { return static_cast<int>(func.size()); }
    };

} // namespace rayt