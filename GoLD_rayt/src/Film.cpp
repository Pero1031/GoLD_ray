#include "pch.h"

#include "Core/Math.hpp"
#include "Renderer/Film.hpp"
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
        m_pixels.resize(width * height, Spectrum(0.0));
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
        m_pixels[y * m_width + x] = radiance;
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
        m_pixels[y * m_width + x] += radiance;
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
     * @brief Saves the film to an image file.
     *
     * The output format is determined by the file extension:
     * - HDR (.hdr): stored as linear floating-point RGB
     * - LDR (.png, .bmp, .jpg): tone-mapped and gamma-corrected
     *
     * @param filename Output file path
     */
    void Film::save(const std::string& filename) const {

        // Check file extension to determine output format.
        std::string ext = filename.substr(filename.find_last_of(".") + 1);

        // Convert extension to lowercase
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

        if (ext == "hdr") {
            // --------------------------------------------------
            // HDR Output
            // --------------------------------------------------
            // Stores raw linear radiance values (float RGB).
            // Best suited for research, analysis, and relighting.
            // Spectrum is assumed to be a contiguous vec3<float>.

            stbi_write_hdr(filename.c_str(), m_width, m_height, 3, reinterpret_cast<const float*>(m_pixels.data()));

            std::cout << "[Film] Saved HDR image: " << filename << std::endl;
        }
        else {
            // --------------------------------------------------
            // LDR Output (PNG / BMP / JPG)
            // --------------------------------------------------
            // Applies display gamma correction and quantization.

            std::vector<unsigned char> outputData(m_width * m_height * 3);

            for (int i = 0; i < m_width * m_height; ++i) {
                Spectrum pixel = m_pixels[i];

                // Convert from linear space to display gamma (gamma 2.2)
                pixel = color::toDisplayGamma22(pixel);

                // Quantize from [0,1] to [0,255]
                outputData[i * 3 + 0] = static_cast<unsigned char>(255.99 * rayt::math::saturate(pixel.r));
                outputData[i * 3 + 1] = static_cast<unsigned char>(255.99 * rayt::math::saturate(pixel.g));
                outputData[i * 3 + 2] = static_cast<unsigned char>(255.99 * rayt::math::saturate(pixel.b));
            }

            if (ext == "png") {
                stbi_write_png(filename.c_str(), m_width, m_height, 3, outputData.data(), m_width * 3);
            }
            else if (ext == "bmp") {
                stbi_write_bmp(filename.c_str(), m_width, m_height, 3, outputData.data());
            }
            else if (ext == "jpg") {
                stbi_write_jpg(filename.c_str(), m_width, m_height, 3, outputData.data(), 90); // Quality 90
            }
            else {
                std::cerr << "[Film] Error: Unsupported file extension: " << ext << std::endl;
                return;
            }

            std::cout << "[Film] Saved LDR image: " << filename << std::endl;
        }
    }

    /**
     * @brief Saves the film to an image file with a scaling factor.
     *
     * This function is typically used to apply normalization,
     * such as dividing by the number of samples per pixel (SPP).
     *
     * @param filename Output file path
     * @param scale    Scaling factor applied to all pixels
     */
    void Film::save(const std::string& filename, float scale) const {

        // Check file extension to determine output format.
        std::string ext = filename.substr(filename.find_last_of(".") + 1);

        // Convert to lowercase for comparison
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

        if (ext == "hdr") {
            // --------------------------------------------------
            // HDR Output with Scaling
            // --------------------------------------------------
            // A temporary buffer is required because scaling
            // must be applied before writing contiguous data.
            std::vector<float> hdrData(m_width * m_height * 3);

            for (size_t i = 0; i < m_pixels.size(); ++i) {
                Spectrum s = m_pixels[i] * static_cast<Real>(scale);
                hdrData[i * 3 + 0] = s.r;
                hdrData[i * 3 + 1] = s.g;
                hdrData[i * 3 + 2] = s.b;
            }

            stbi_write_hdr(filename.c_str(), m_width, m_height, 3, hdrData.data());
            std::cout << "[Film] Saved HDR image: " << filename << std::endl;
        }
        else {
            // --------------------------------------------------
            // LDR Output with Scaling
            // --------------------------------------------------
            std::vector<unsigned char> outputData(m_width * m_height * 3);

            for (int i = 0; i < m_width * m_height; ++i) {
                
                // Apply scaling (e.g., 1 / spp)
                Spectrum pixel = m_pixels[i] * static_cast<Real>(scale);

                pixel = color::toDisplayGamma22(pixel);

                outputData[i * 3 + 0] = static_cast<unsigned char>(255.99 * math::saturate(pixel.r));
                outputData[i * 3 + 1] = static_cast<unsigned char>(255.99 * math::saturate(pixel.g));
                outputData[i * 3 + 2] = static_cast<unsigned char>(255.99 * math::saturate(pixel.b));
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
    }

} // namespace rayt