#include "pch.h"

#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <utility>
#include <cmath>

#include "Core/Types.hpp"
#include "IO/ImageLoader.hpp"
#include "Core/Image.hpp"

// Define STB_IMAGE_IMPLEMENTATION in only one source file (usually pch.cpp or here if not in pch)
// #define STB_IMAGE_IMPLEMENTATION 
#include "stb_image.h"

namespace rayt::io {

    /**
     * @brief Converts a string to lower case.
     * @param s The input string.
     * @return The lower case string.
     */
    static std::string toLower(std::string s) {
        for (char& c : s) {
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        }
        return s;
    }

    /**
     * @brief Extracts the file extension from a path.
     * @param path The file path.
     * @return The extension in lower case (without the dot), or empty if none found.
     * @todo Replace with std::filesystem::path::extension() for production use
     */
    static std::string getExt(const std::string& path) {
        auto p = path.find_last_of('.');
        if (p == std::string::npos) return "";
        return toLower(path.substr(p + 1));
    }

    /**
     * @brief Converts an sRGB color component to linear space.
     * * Applies the inverse gamma correction to convert non-linear sRGB values
     * to linear intensity values suitable for physically based rendering.
     *
     * @param c The sRGB color component (normalized 0.0 to 1.0).
     * @return The linear color component.
     */
    static float srgbToLinear(float c) {
        if (c <= 0.04045f)
            return c / 12.92f;
        return std::pow((c + 0.055f) / 1.055f, 2.4f);
    }

    /**
     * @brief Loads a High Dynamic Range (HDR) image.
     * * Reads the image as floating point data. No color space conversion
     * is applied as HDR images are assumed to be in linear space.
     *
     * @param filename The path to the HDR file.
     * @return A constructed Image object containing the pixel data.
     * @throws std::runtime_error If the file cannot be loaded.
     */
    Image loadHDR(const std::string& filename) {
        int w = 0, h = 0, n = 0;

        // Load as float RGB
        float* data = stbi_loadf(filename.c_str(), &w, &h, &n, 3);
        if (!data) {
            throw std::runtime_error(
                "Failed to load HDR: " + filename + " (" + stbi_failure_reason() + ")"
            );
        }

        std::vector<Vector3> pixels;
        pixels.resize(static_cast<size_t>(w) * h);

        for (int i = 0; i < w * h; ++i) {
            pixels[i] = Vector3(
                data[3 * i + 0],
                data[3 * i + 1],
                data[3 * i + 2]
            );
        }

        stbi_image_free(data);

        return Image(w, h, std::move(pixels));
    }

    /**
     * @brief Loads a Low Dynamic Range (LDR) image (png, jpg, etc.).
     * * Reads the image as 8-bit data, normalizes it to [0, 1], and converts
     * it from sRGB to linear space.
     *
     * @param filename The path to the image file.
     * @return A constructed Image object with linear color data.
     * @throws std::runtime_error If the file cannot be loaded.
     */
    Image loadLDR(const std::string& filename) {
        int w = 0, h = 0, n = 0;

        // Force load as RGB (3 channels)
        unsigned char* data = stbi_load(filename.c_str(), &w, &h, &n, 3);
        if (!data) {
            throw std::runtime_error(
                "Failed to load LDR: " + filename + " (" + stbi_failure_reason() + ")"
            );
        }

        std::vector<Vector3> pixels;
        pixels.resize(static_cast<size_t>(w) * h);

        for (int i = 0; i < w * h; ++i) {
            float r = data[3 * i + 0] / 255.0f;
            float g = data[3 * i + 1] / 255.0f;
            float b = data[3 * i + 2] / 255.0f;

            // Convert sRGB to Linear
            pixels[i] = Vector3(
                srgbToLinear(r),
                srgbToLinear(g),
                srgbToLinear(b)
            );
        }

        stbi_image_free(data);
        return Image(w, h, std::move(pixels));
    }

    /**
     * @brief Loads an image from a file, automatically detecting the format.
     * * Dispatches to loadHDR or loadLDR based on the file extension.
     * Supported HDR formats: .hdr
     * Supported LDR formats: .png, .jpg, .jpeg, .bmp, .tga
     *
     * @param filename The path to the image file.
     * @return The loaded Image.
     * @throws std::runtime_error If the format is unsupported or loading fails.
     */
    Image loadImage(const std::string& filename) {
        const std::string ext = getExt(filename);

        if (ext == "hdr") {
            return loadHDR(filename);
        }

        if (ext == "png" || ext == "jpg" || ext == "jpeg" ||
            ext == "bmp" || ext == "tga") {
            return loadLDR(filename);
        }

        throw std::runtime_error("Unsupported image format: " + ext);
    }

} // namespace rayt::io