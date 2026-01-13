/**
 * @file Film.hpp
 * @brief Film class for storing per-pixel accumulated radiance (sum buffer).
 *
 * Film represents the camera sensor buffer used by the renderer.
 * It stores high-dynamic-range (HDR) radiance values per pixel as an
 * accumulation (sum) over samples/passes. Any normalization (1/spp),
 * tone mapping, and display encoding (gamma/sRGB) are performed at output time
 * (e.g., UI display or file saving).
 */

#pragma once

#include <string>
#include <vector>
#include <memory>
#include <algorithm>

#include "Core/Core.hpp"  // Include core definitions (Spectrum, Real, etc.)

namespace rayt {

    /**
    * @brief Accumulation buffer for rendered radiance.
    *
    * The internal pixel buffer stores the sum of radiance contributions:
    *   sum(x,y) = Σ_k L_k(x,y)
    *
    * To obtain the unbiased estimate (average), scale by 1/spp at output:
    *   avg = sum / spp
    *
    * This design is convenient for progressive rendering where the application
    * controls the current sample count.
    */
    class Film {
    public:
        /**
         * @brief Initializes the film with a given resolution.
         * @param width  Horizontal resolution in pixels.
         * @param height Vertical resolution in pixels.
         */
        Film(int width, int height);

        /**
         * @brief Sets the spectral radiance for a specific pixel coordinate.
         * * In advanced integrators, this may be extended to an 'addSample' method
         * to facilitate multi-sample accumulation and reconstruction filtering.
         * @param x        Pixel x-coordinate.
         * @param y        Pixel y-coordinate.
         * @param radiance The physical intensity (radiance) to store.
         */
        void setPixel(int x, int y, const Spectrum& radiance);

        /// Adds radiance to the accumulation buffer at (x,y) (progressive rendering).
        void addPixel(int x, int y, const Spectrum& radiance);

        /// Returns pointer to the accumulation (sum) buffer (linear HDR).
        const Spectrum* getData() const { return m_pixels.data(); }

        // 蓄積バッファをクリアする
        void clear();

        /**
         * @brief Returns the width of the film in pixels.
         */
        int width() const { return m_width; }

        /**
         * @brief Returns the height of the film in pixels.
         */
        int height() const { return m_height; }

        /**
        * @brief Saves the film to a file.
        *
        * @note If the film stores accumulated sums (recommended), use the overload with
        *       an explicit scale (e.g., scale = 1/spp) for correct exposure.
        */
        void save(const std::string& filename) const;

        /**
        * @brief Saves the film to a file with an explicit scale applied.
        *
        * Typical usage:
        *   - Progressive rendering: scale = 1 / currentSpp
        *   - HDR analysis: scale = 1 (raw sum) or 1/spp depending on intent
         */
        void save(const std::string& filename, float scale) const;

    private:
        int m_width;
        int m_height;

        /**
         * @brief Internal buffer for raw pixel data.
         * * Using 'Spectrum' ensures that no physical light information is lost
         * during the rendering process, allowing for flexible post-processing.
         */
        std::vector<Spectrum> m_pixels;
    };

} // namespace rayt