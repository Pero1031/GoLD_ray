/*
* @file src/Render/Fil.cpp
* @brief Film class implementation (pixel storage, accumulation, and image output).
* 
* This file implements the Film class which stores per-pixel radiance values
* in linear RGB space and provides routines to save the result as HDR (float)
* or LDR (8-bit) image formats via stb_image_write.
*/

#include "pch.h"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>
#include <vector>

#include "Core/Math.hpp"
#include "Render/Film.hpp"
#include "Color/ColorTransform.hpp"

// stb_image_write implementation is assumed to be compiled
// in a separate translation unit (e.g., ImageIO.cpp).
#include "stb_image_write.h"

namespace rayt {

    /**
     * @brief Constructs a Film with the given resolution.
     *
     * The Film stores radiance values per pixel in linear color space.
     * All pixels are initialized to black (zero radiance).
     *
     * @param width  Image width in pixels
     * @param height Image height in pixels
     */
    Film::Film(int width, int height)
        : m_width(width), m_height(height) {
        // Initialize all pixels to black (0, 0, 0).
        m_pixels.resize(static_cast<size_t>(width) * static_cast<size_t>(height), Spectrum(0.0));
    }

    /**
     * @brief Sets the radiance value of a single pixel.
     *
     * This function overwrites the existing pixel value.
     * Bounds checking is performed to avoid invalid memory access.
     *
     * @param x         Pixel x-coordinate
     * @param y         Pixel y-coordinate
     * @param radiance  Radiance value to store (linear space)
     */
    void Film::setPixel(int x, int y, const Spectrum& radiance) {

        // Boundary check to prevent segmentation faults.
        if (x < 0 || x >= m_width || y < 0 || y >= m_height) {
            return;
        }

        // Store the value directly.
        // Image coordinates assume (0,0) at the top-left corner.
        m_pixels[static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)] = radiance;
    }

    /**
     * @brief Adds radiance to a pixel (for progressive rendering).
     *
     * This function accumulates radiance values using +=,
     * which is typically used in Monte Carlo integration
     * and progressive / iterative rendering.
     *
     * @param x         Pixel x-coordinate
     * @param y         Pixel y-coordinate
     * @param radiance  Radiance contribution to add
     */
    void Film::addPixel(int x, int y, const Spectrum& radiance) {
        if (x < 0 || x >= m_width || y < 0 || y >= m_height) {
            return;
        }

        // Accumulate radiance
        m_pixels[static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)] += radiance;
    }

    /**
     * @brief Clears all accumulated pixel values.
     *
     * Resets the film to a black image.
     * Useful when restarting a progressive render.
     */
    void Film::clear() {
        std::fill(m_pixels.begin(), m_pixels.end(), Spectrum(0.0));
    }

    /**
     * @brief Converts a normalized float value in [0, 1] to an 8-bit byte.
     *
     * Values are saturated to [0, 1] for safety and then converted to
     * [0, 255] using round-to-nearest.
     *
     * @param x Normalized value (expected in [0,1])
     * @return 8-bit value in [0,255]
     */
    static inline unsigned char toByte(Real x) {
        // x is expected in [0,1] but saturate to be safe
        x = rayt::math::saturate(x);
        return static_cast<unsigned char>(Real(255.0) * x + Real(0.5)); // round-to-nearest
    }

    /**
     * @brief Saves the film to an image file (default scale = 1.0).
     * This overload forwards to save(filename, scale).
     *
     * @param filename Output file path. Extension selects the format.
     */
    void Film::save(const std::string& filename) const {
        save(filename, Real(1.0));
    }

    /**
     * @brief Saves the film to an image file with a multiplicative scale.
     *
     * The file extension determines the output format:
     * - ".hdr": saved as 32-bit float RGB (linear), without tone mapping.
     * - ".png/.bmp/.jpg": saved as 8-bit RGB after display transform.
     *
     * For LDR output, the pipeline is:
     *   exposure -> tone map (Reinhard) -> sRGB OETF -> clamp -> 8-bit
     *
     * @param filename Output file path.
     * @param scale    Multiply all stored radiance by this factor before saving.
     */
    void Film::save(const std::string& filename, float scale) const {

        /**
         * @brief Parse and normalize extension to lowercase.
         *
         * We do a simple extension parse by taking the substring after the last '.'.
         * If the filename has no '.', behavior is undefined (same as original code).
         */
        std::string ext = filename.substr(filename.find_last_of('.') + 1);
        std::transform(ext.begin(), ext.end(), ext.begin(),
            [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

        const Real s = static_cast<Real>(scale);

        if (ext == "hdr") {
            // --------------------------------------------------
            // HDR Output (always pack to float buffer safely)
            // --------------------------------------------------
            std::vector<float> hdrData(static_cast<size_t>(m_width) * static_cast<size_t>(m_height) * 3);

            const int numPixels = static_cast<int>(m_pixels.size());

            #pragma omp parallel for schedule(static)
            for (int i = 0; i < numPixels; ++i) {
                Spectrum p = m_pixels[i] * s;

                // NaN/Inf -> 0
                if (math::hasNonFinite(p.r)) p.r = Real(0);
                if (math::hasNonFinite(p.g)) p.g = Real(0);
                if (math::hasNonFinite(p.b)) p.b = Real(0);

                hdrData[i * 3 + 0] = static_cast<float>(p.r);
                hdrData[i * 3 + 1] = static_cast<float>(p.g);
                hdrData[i * 3 + 2] = static_cast<float>(p.b);
            }

            stbi_write_hdr(filename.c_str(), m_width, m_height, 3, hdrData.data());
            std::cout << "[Film] Saved HDR image: " << filename << std::endl;
            return;
        }

        // --------------------------------------------------
        // LDR Output (PNG / BMP / JPG)
        // --------------------------------------------------
        std::vector<unsigned char> outputData(static_cast<size_t>(m_width) * static_cast<size_t>(m_height) * 3);

        // TODO: UI/設定から渡す（今は固定）
        const Real exposureStops = Real(0);
        const int numPixels = static_cast<int>(m_pixels.size());

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < numPixels; ++i) {
            Spectrum pixel = m_pixels[i] * s;

            // New pipeline: exposure -> Reinhard -> sRGB encode -> clamp
            pixel = color::toDisplaySRGB_Reinhard(pixel, exposureStops);

            outputData[i * 3 + 0] = toByte(pixel.r);
            outputData[i * 3 + 1] = toByte(pixel.g);
            outputData[i * 3 + 2] = toByte(pixel.b);
        }

        if (ext == "png") {
            stbi_write_png(filename.c_str(), m_width, m_height, 3, outputData.data(), m_width * 3);
        }
        else if (ext == "bmp") {
            stbi_write_bmp(filename.c_str(), m_width, m_height, 3, outputData.data());
        }
        else if (ext == "jpg") {
            stbi_write_jpg(filename.c_str(), m_width, m_height, 3, outputData.data(), 90);
        }
        else {
            std::cerr << "[Film] Error: Unsupported file extension: " << ext << std::endl;
            return;
        }

        std::cout << "[Film] Saved LDR image: " << filename << std::endl;
    }

} // namespace rayt