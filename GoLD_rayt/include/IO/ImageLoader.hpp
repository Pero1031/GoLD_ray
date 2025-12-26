#pragma once

#include <string>
#include <stdexcept>

#include "Core/Image.hpp"

namespace rayt::io {

    /**
     * @brief Loads an image file, automatically determining the format/loader by extension.
     *
     * Dispatches to loadHDR() or loadLDR() based on the file extension.
     * Throws an exception if the format is unsupported or the file cannot be read.
     *
     * @param filename The path to the image file.
     * @return The loaded Image object.
     * @throws std::runtime_error If loading fails.
     */
    Image loadImage(const std::string& filename);

    /**
     * @brief Explicitly loads a High Dynamic Range (HDR) image.
     *
     * Reads the file as floating-point data. Assumes the data is already in linear space.
     *
     * @param filename The path to the HDR file (typically .hdr).
     * @return The loaded Image object containing linear float data.
     */
    Image loadHDR(const std::string& filename);

    /**
     * @brief Explicitly loads a Low Dynamic Range (LDR) image (PNG, JPG, etc.).
     *
     * Reads the file as 8-bit data and normalizes it.
     * IMPORTANT: This function automatically converts colors from sRGB to Linear space
     * for physically based rendering.
     *
     * @param filename The path to the LDR file.
     * @return The loaded Image object containing linear float data.
     */
    Image loadLDR(const std::string& filename);

} //namespace rayt::io