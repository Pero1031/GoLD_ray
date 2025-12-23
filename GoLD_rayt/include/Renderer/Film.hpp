#pragma once

/**
 * @file Film.hpp
 * @brief Film class for capturing and storing rendered spectral radiance.
 * * The Film class represents the image sensor of the camera. It captures
 * the high-dynamic-range (HDR) radiance values for each pixel and manages
 * the final image output with appropriate post-processing (tone mapping, gamma).
 */

#include <string>
#include <vector>
#include <memory>
#include <algorithm>

#include "Core/Core.hpp"  // Include core definitions (Spectrum, Real, etc.)

namespace rayt {

    /**
     * @brief Models the light-sensing device in a simulated camera.
     * * It stores raw 'Spectrum' data for each pixel to preserve physical intensity,
     * which is essential for accurate light transport and HDR output formats.
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

        /**
         * @brief Returns the width of the film in pixels.
         */
        int width() const { return m_width; }

        /**
         * @brief Returns the height of the film in pixels.
         */
        int height() const { return m_height; }

        /**
         * @brief Saves the current film data to a file.
         * * Automatic format handling based on extension:
         * - ".hdr": Saves raw linear radiance (Radiance HDR format).
         * - ".png" / ".jpg": Performs tone mapping and gamma correction (LDR).
         * @param filename Path to the output file.
         */
        void save(const std::string& filename) const;

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