#pragma once

/**
 * @file Distribution2D.hpp
 * @brief Generic 2D Probability Distribution for Monte Carlo Integration.
 *
 * This file defines the Distribution2D structure, which implements hierarchical
 * sampling for 2D functions. It decomposes a 2D distribution p(u, v) into
 * a marginal distribution p(v) and conditional distributions p(u|v).
 * * It is commonly used for importance sampling of environment maps,
 * but can be applied to any 2D tabulated data.
 */

#include <vector>
#include <memory>
#include <cmath>

#include "Distribution1D.hpp"
#include "Core/Types.hpp"

namespace rayt {

    /**
     * @struct Distribution2D
     * @brief A class for 2D importance sampling using a hierarchical approach.
     *
     * It decomposes the 2D distribution p(u, v) into a marginal distribution p(v)
     * and conditional distributions p(u|v).
     * Typically used for sampling directions from an environment map (IBL).
     */
    struct Distribution2D {
        std::vector<std::unique_ptr<Distribution1D>> pConditionalV; ///< Conditional distributions p(u|v) for each row.
        std::unique_ptr<Distribution1D> pMarginal;                  ///< Marginal distribution p(v) for selecting rows.

        /**
         * @brief Constructs a 2D distribution from raw floating-point data.
         *
         * @param data The flat array of 2D data (e.g., luminance values). Size must be width * height.
         * @param width The width of the 2D data.
         * @param height The height of the 2D data.
         */
        Distribution2D(const float* data, int width, int height) {
            // 1. Build conditional distributions p(u|v) for each row v.
            for (int v = 0; v < height; ++v) {
                pConditionalV.emplace_back(std::make_unique<Distribution1D>(&data[v * width], width));
            }

            // 2. Compute the marginal integrals to build p(v).
            // Each element of marginalFunc represents the integral of one row.
            std::vector<float> marginalFunc;
            for (int v = 0; v < height; ++v) {
                marginalFunc.push_back(pConditionalV[v]->funcInt);
            }
            pMarginal = std::make_unique<Distribution1D>(marginalFunc.data(), height);
        }

        /**
         * @brief Samples a continuous 2D coordinate based on the distribution.
         *
         * First samples the v-coordinate using the marginal distribution,
         * then samples the u-coordinate using the conditional distribution of the selected row.
         *
         * @param u Two uniform random numbers in [0, 1).
         * @param[out] uv The sampled 2D coordinates in [0, 1].
         * @param[out] pdf The probability density value at the sampled point.
         */
        void sampleContinuous(const Point2& u, Point2& uv, float& pdf) const {
            float pdfV, pdfU;
            int vOff, uOff;

            // 1. Sample v from p(v)
            float v = pMarginal->sampleContinuous(u.y, pdfV, vOff);

            // 2. Sample u from p(u|v)
            float x = pConditionalV[vOff]->sampleContinuous(u.x, pdfU, uOff);

            uv = Point2(x, v);

            // Joint PDF: p(u, v) = p(u|v) * p(v)
            pdf = pdfU * pdfV;   
        }

        /**
         * @brief Evaluates the probability density function (PDF) at a given coordinate.
         *
         * @param uv The 2D coordinates in [0, 1].
         * @return The probability density value.
         */
        float pdf(const Point2& uv) const {
            // Note: Assuming pMarginal and pConditionalV are valid and match the dimensions.
            int height = (int)pMarginal->func.size();
            int width = (int)pConditionalV[0]->func.size();

            // Clamp coordinates to valid range indices
            int v = std::clamp(int(uv.y * height), 0, height - 1);
            int u = std::clamp(int(uv.x * width), 0, width - 1);

            // Handle edge case where the entire image is black (integral is 0)
            if (pMarginal->funcInt == 0.0f) return 1.0f; // Return uniform PDF
            if (pConditionalV[v]->funcInt == 0.0f) return 0.0f;

            // Compute p(v) and p(u|v)
            // Since Distribution1D::funcInt is the integral (sum / N),
            // func[i] / funcInt directly gives the probability density.
            float pv = pMarginal->func[v] / pMarginal->funcInt;                 
            float puv = pConditionalV[v]->func[u] / pConditionalV[v]->funcInt;   

            return pv * puv;
        }

    };

} // namespace rayt

