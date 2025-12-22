#include "pch.h"

#include "Renderer/Film.hpp"
//#include "Core/Utils.hpp" // For linearToGamma, saturate, etc.
#include "Renderer/ColorTransform.hpp"
#include "Core/Math.hpp"

// Assuming STB_IMAGE_WRITE_IMPLEMENTATION is defined in ImageIO.cpp
#include "stb_image_write.h"

namespace rayt {

    Film::Film(int width, int height)
        : m_width(width), m_height(height) {
        // Initialize all pixels to black (0, 0, 0).
        m_pixels.resize(width * height, Spectrum(0.0));
    }

    void Film::setPixel(int x, int y, const Spectrum& radiance) {
        // Boundary check to prevent segmentation faults.
        if (x < 0 || x >= m_width || y < 0 || y >= m_height) {
            return;
        }
        // Store the value directly.
        // Note: Coordinate system usually assumes (0,0) is top-left for images.
        m_pixels[y * m_width + x] = radiance;
    }

    void Film::save(const std::string& filename) const {
        // Check file extension to determine output format.
        std::string ext = filename.substr(filename.find_last_of(".") + 1);

        // Convert to lowercase for comparison
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

        if (ext == "hdr") {
            // --- HDR Output ---
            // Save raw linear float data. Best for research and analysis.
            // stbi_write_hdr expects contiguous float array (3 floats per pixel).
            // Since Spectrum is glm::vec3 (3 floats), we can cast the pointer safely.
            stbi_write_hdr(filename.c_str(), m_width, m_height, 3, reinterpret_cast<const float*>(m_pixels.data()));

            std::cout << "[Film] Saved HDR image: " << filename << std::endl;
        }
        else {
            // --- LDR Output (PNG, BMP, JPG) ---
            // Requires Tone Mapping and Gamma Correction.

            std::vector<unsigned char> outputData(m_width * m_height * 3);

            for (int i = 0; i < m_width * m_height; ++i) {
                Spectrum pixel = m_pixels[i];

                // 1. Tone Mapping (Optional but recommended)
                // Simple Reinhard tone mapping: L_d = L / (1 + L)
                // This prevents high-intensity values (e.g., > 1.0) from being clamped to white instantly.
                // You can comment this out if you want to see raw clipping.
                pixel = pixel / (pixel + Spectrum(1.0));

                // 2. Gamma Correction (Linear -> sRGB)
                // Monitors expect sRGB data (approx. gamma 2.2).
                pixel.r = renderer::linearToGamma(pixel.r);
                pixel.g = renderer::linearToGamma(pixel.g);
                pixel.b = renderer::linearToGamma(pixel.b);

                // 3. Quantization (0.0-1.0 -> 0-255)
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

} // namespace rayt