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
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "Core/Math.hpp"
#include "Render/Film.hpp"
#include "Color/ColorSpace.hpp"
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
        m_pixels.resize(static_cast<size_t>(width) * static_cast<size_t>(height), Vector3(0.0));
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
    void Film::setPixel(int x, int y, const Vector3& radiance) {

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
    void Film::addPixel(int x, int y, const Vector3& radiance) {
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
        std::fill(m_pixels.begin(), m_pixels.end(), Vector3(0.0));
    }

    void Film::enableAOVs() {
        const size_t n = static_cast<size_t>(m_width) * static_cast<size_t>(m_height);
        m_albedo.assign(n, Vector3(0.0));
        m_normal.assign(n, Vector3(0.0));
    }

    void Film::setAOV(int x, int y, const Vector3& albedo, const Vector3& normal) {
        if (!hasAOVs() || x < 0 || x >= m_width || y < 0 || y >= m_height) {
            return;
        }
        const size_t i = static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x);
        m_albedo[i] = albedo;
        m_normal[i] = normal;
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

    /// CRC-32 of the PNG spec (polynomial 0xEDB88320), over a chunk's type and data.
    static uint32_t pngCrc32(const unsigned char* data, size_t n) {
        uint32_t crc = 0xFFFFFFFFu;
        for (size_t i = 0; i < n; ++i) {
            crc ^= data[i];
            for (int k = 0; k < 8; ++k) crc = (crc >> 1) ^ (0xEDB88320u & (0u - (crc & 1u)));
        }
        return ~crc;
    }

    /**
     * @brief Splices a cICP chunk (ITU-T H.273 code points) into an encoded PNG, right after IHDR (it must come before PLTE / IDAT).
     *
     * stb_image_write cannot add chunks. Code points: primaries 12 (Display P3) or 9 (BT.2020), transfer 13 (sRGB curve),
     * matrix 0 (RGB), full range. A viewer that ignores cICP shows the file as sRGB, i.e. less saturated.
     */
    static bool insertCicp(std::vector<unsigned char>& png, color::Gamut gamut) {
        constexpr size_t kAfterIhdr = 8 + 25;   // signature + (length 4, type 4, data 13, crc 4)
        if (png.size() < kAfterIhdr) return false;
        unsigned char chunk[4 + 4 + 4 + 4];
        chunk[0] = 0; chunk[1] = 0; chunk[2] = 0; chunk[3] = 4;   // data length
        chunk[4] = 'c'; chunk[5] = 'I'; chunk[6] = 'C'; chunk[7] = 'P';
        chunk[8] = static_cast<unsigned char>(color::gamutPrimariesCode(gamut));
        chunk[9] = 13;   // transfer: sRGB curve
        chunk[10] = 0;   // matrix: RGB
        chunk[11] = 1;   // full range
        const uint32_t crc = pngCrc32(chunk + 4, 8);
        chunk[12] = static_cast<unsigned char>(crc >> 24); chunk[13] = static_cast<unsigned char>(crc >> 16);
        chunk[14] = static_cast<unsigned char>(crc >> 8);  chunk[15] = static_cast<unsigned char>(crc);
        png.insert(png.begin() + kAfterIhdr, chunk, chunk + sizeof(chunk));
        return true;
    }

    static bool writeBytes(const std::string& filename, const std::vector<unsigned char>& bytes) {
        std::ofstream f(filename, std::ios::binary);
        if (!f) return false;
        f.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
        return static_cast<bool>(f);
    }

    static void pngSink(void* ctx, void* data, int size) {
        auto* v = static_cast<std::vector<unsigned char>*>(ctx);
        const auto* p = static_cast<const unsigned char*>(data);
        v->insert(v->end(), p, p + size);
    }

    /// Writes an 8-bit RGB PNG tagged with a cICP chunk.
    static bool writePngWithCicp(const std::string& filename, int width, int height,
        const std::vector<unsigned char>& rgb8, color::Gamut gamut) {
        std::vector<unsigned char> png;
        if (!stbi_write_png_to_func(pngSink, &png, width, height, 3, rgb8.data(), width * 3)) return false;
        return insertCicp(png, gamut) && writeBytes(filename, png);
    }

    /**
     * @brief Writes a 16-bit RGB PNG (big-endian samples, 6 bytes per pixel), tagged with cICP unless the gamut is sRGB.
     *
     * stb_image_write has no 16-bit mode. A row of W pixels at 16 bits is 6 * W bytes: exactly a row of 2 * W pixels at 8 bits. So the data
     * is encoded as an 8-bit RGB image of twice the width with the PNG filter forced to "none" (the filters work on the bytes of a
     * whole pixel, which would be wrong here), and the IHDR is then patched to width W, 16 bits (and its CRC redone).
     */
    static bool writePng16(const std::string& filename, int width, int height,
        const std::vector<unsigned char>& data16, color::Gamut gamut) {
        std::vector<unsigned char> png;
        const int savedFilter = stbi_write_force_png_filter;
        stbi_write_force_png_filter = 0;
        const int ok = stbi_write_png_to_func(pngSink, &png, 2 * width, height, 3, data16.data(), 6 * width);
        stbi_write_force_png_filter = savedFilter;
        if (!ok || png.size() < 33) return false;
        for (int k = 0; k < 4; ++k) png[16 + k] = static_cast<unsigned char>(static_cast<uint32_t>(width) >> (8 * (3 - k)));   // width
        png[24] = 16;                                                                                                         // bit depth
        const uint32_t crc = pngCrc32(png.data() + 12, 17);                                                                   // type + data
        for (int k = 0; k < 4; ++k) png[29 + k] = static_cast<unsigned char>(crc >> (8 * (3 - k)));
        if (gamut != color::Gamut::SRGB && !insertCicp(png, gamut)) return false;
        return writeBytes(filename, png);
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

        if (ext == "pfm") {
            // --------------------------------------------------
            // PFM output: raw 32-bit float RGB, linear, no tone mapping. Unlike .hdr (RGBE: one shared exponent and an 8-bit mantissa per
            // pixel) it keeps the small channels of a strongly colored pixel (the blue of a candle-lit scene) and negative values (colors
            // outside the sRGB gamut), which is what a measurement needs. Rows are stored bottom to top, little endian ("-1" scale).
            // --------------------------------------------------
            std::vector<float> pfmData(static_cast<size_t>(m_width) * static_cast<size_t>(m_height) * 3);
            for (int y = 0; y < m_height; ++y) {
                for (int x = 0; x < m_width; ++x) {
                    Vector3 p = color::fromLinearSRGB(m_gamut, m_pixels[static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)] * s);
                    const size_t o = (static_cast<size_t>(m_height - 1 - y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)) * 3;
                    for (int c = 0; c < 3; ++c)
                        pfmData[o + c] = math::hasNonFinite(p[c]) ? 0.0f : static_cast<float>(p[c]);
                }
            }
            std::ofstream f(filename, std::ios::binary);
            f << "PF\n" << m_width << " " << m_height << "\n-1.0\n";
            f.write(reinterpret_cast<const char*>(pfmData.data()), static_cast<std::streamsize>(pfmData.size() * sizeof(float)));
            std::cout << "[Film] Saved PFM image: " << filename << std::endl;
            return;
        }

        if (ext == "hdr") {
            // --------------------------------------------------
            // HDR Output (always pack to float buffer safely)
            // --------------------------------------------------
            std::vector<float> hdrData(static_cast<size_t>(m_width) * static_cast<size_t>(m_height) * 3);

            const int numPixels = static_cast<int>(m_pixels.size());

            #pragma omp parallel for schedule(static)
            for (int i = 0; i < numPixels; ++i) {
                Vector3 p = color::fromLinearSRGB(m_gamut, m_pixels[i] * s);

                // NaN/Inf -> 0
                if (math::hasNonFinite(p[0])) p[0] = Real(0);
                if (math::hasNonFinite(p[1])) p[1] = Real(0);
                if (math::hasNonFinite(p[2])) p[2] = Real(0);

                hdrData[i * 3 + 0] = static_cast<float>(p[0]);
                hdrData[i * 3 + 1] = static_cast<float>(p[1]);
                hdrData[i * 3 + 2] = static_cast<float>(p[2]);
            }

            stbi_write_hdr(filename.c_str(), m_width, m_height, 3, hdrData.data());
            std::cout << "[Film] Saved HDR image: " << filename << std::endl;
            return;
        }

        // --------------------------------------------------
        // LDR Output (PNG / BMP / JPG)
        // --------------------------------------------------
        std::vector<unsigned char> outputData(static_cast<size_t>(m_width) * static_cast<size_t>(m_height) * 3);

        // TODO: pass this in from the UI / settings (fixed for now)
        const Real exposureStops = Real(0);
        const int numPixels = static_cast<int>(m_pixels.size());

        const bool want16 = m_png16 && ext == "png";
        std::vector<unsigned char> data16;
        if (want16) data16.resize(static_cast<size_t>(numPixels) * 6);

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < numPixels; ++i) {
            // Change of primaries first (exact, unclamped), then the display pipeline, which clips what is still negative:
            // exposure -> Reinhard -> sRGB curve -> clamp. The curve is the sRGB one for every gamut (as for Display P3).
            Vector3 pixel = color::fromLinearSRGB(m_gamut, m_pixels[i] * s);
            pixel = color::toDisplaySRGB_Reinhard(pixel, exposureStops);

            outputData[i * 3 + 0] = toByte(pixel[0]);
            outputData[i * 3 + 1] = toByte(pixel[1]);
            outputData[i * 3 + 2] = toByte(pixel[2]);
            if (want16) {
                for (int c = 0; c < 3; ++c) {
                    const unsigned v = static_cast<unsigned>(Real(65535.0) * rayt::math::saturate(pixel[c]) + Real(0.5));
                    data16[static_cast<size_t>(i) * 6 + static_cast<size_t>(c) * 2] = static_cast<unsigned char>(v >> 8);
                    data16[static_cast<size_t>(i) * 6 + static_cast<size_t>(c) * 2 + 1] = static_cast<unsigned char>(v & 255u);
                }
            }
        }

        if (ext == "png") {
            if (want16) {
                if (!writePng16(filename, m_width, m_height, data16, m_gamut)) {
                    std::cerr << "[Film] Error: could not write " << filename << std::endl;
                    return;
                }
            }
            else if (m_gamut == color::Gamut::SRGB)
                stbi_write_png(filename.c_str(), m_width, m_height, 3, outputData.data(), m_width * 3);
            else if (!writePngWithCicp(filename, m_width, m_height, outputData, m_gamut)) {
                std::cerr << "[Film] Error: could not write " << filename << std::endl;
                return;
            }
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

        std::cout << "[Film] Saved LDR image: " << filename << (want16 ? " (16 bit)" : "");
        if (m_gamut != color::Gamut::SRGB && ext != "png")
            std::cout << " (" << color::gamutName(m_gamut) << " values, but this format carries no color tag; use .png)";
        else if (m_gamut != color::Gamut::SRGB)
            std::cout << " (" << color::gamutName(m_gamut) << ", cICP-tagged)";
        std::cout << std::endl;
    }

    void Film::saveAlbedo(const std::string& filename) const {
        saveBuffer(filename, m_albedo, /*encodeSrgb=*/true, /*remapSigned=*/false);
    }

    void Film::saveNormal(const std::string& filename) const {
        saveBuffer(filename, m_normal, /*encodeSrgb=*/false, /*remapSigned=*/true);
    }

    void Film::saveBuffer(const std::string& filename, const std::vector<Vector3>& data,
        bool encodeSrgb, bool remapSigned) const {
        if (data.empty()) {
            std::cerr << "[Film] Error: AOVs are not enabled, nothing to save: " << filename << std::endl;
            return;
        }

        std::string ext = filename.substr(filename.find_last_of('.') + 1);
        std::transform(ext.begin(), ext.end(), ext.begin(),
            [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

        const size_t numPixels = data.size();

        if (ext == "hdr") {
            std::vector<float> hdrData(numPixels * 3);
            for (size_t i = 0; i < numPixels; ++i) {
                for (int c = 0; c < 3; ++c) {
                    const Real v = data[i][c];
                    hdrData[i * 3 + c] = math::hasNonFinite(v) ? 0.0f : static_cast<float>(v);
                }
            }
            stbi_write_hdr(filename.c_str(), m_width, m_height, 3, hdrData.data());
            std::cout << "[Film] Saved HDR image: " << filename << std::endl;
            return;
        }

        std::vector<unsigned char> out(numPixels * 3);
        for (size_t i = 0; i < numPixels; ++i) {
            for (int c = 0; c < 3; ++c) {
                Real v = data[i][c];
                if (math::hasNonFinite(v)) v = Real(0);
                if (remapSigned) v = v * Real(0.5) + Real(0.5);
                else if (encodeSrgb) v = color::linearToSRGB(rayt::math::saturate(v));
                out[i * 3 + c] = toByte(v);
            }
        }

        if (ext == "png") {
            stbi_write_png(filename.c_str(), m_width, m_height, 3, out.data(), m_width * 3);
        }
        else if (ext == "bmp") {
            stbi_write_bmp(filename.c_str(), m_width, m_height, 3, out.data());
        }
        else if (ext == "jpg") {
            stbi_write_jpg(filename.c_str(), m_width, m_height, 3, out.data(), 90);
        }
        else {
            std::cerr << "[Film] Error: Unsupported file extension: " << ext << std::endl;
            return;
        }

        std::cout << "[Film] Saved LDR image: " << filename << std::endl;
    }

} // namespace rayt