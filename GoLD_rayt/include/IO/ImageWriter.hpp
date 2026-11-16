/**
 * @file include/IO/ImageWriter.hpp
 * @brief Encoders of image files: PNG (8 / 16 bit, with a cICP tag for wide gamuts), BMP, JPG, PFM and Radiance HDR.
 *
 * They take values that are ready to store (bytes, 16-bit integers, floats), rows from the top of the picture down, and only encode;
 * what the values are (tone mapping, gamut, coverage) is the caller's business (Render/Film). Every function returns false when the
 * file could not be written.
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "Color/Gamut.hpp"

namespace rayt::io {

    /// .png, .bmp or .jpg (quality 90) by the extension, 8 bits per channel, 1 (gray) or 3 (RGB) channels. A PNG of 3 channels gets a cICP
    /// chunk when `tag` is not sRGB (the values are then in that gamut, sRGB curve); the other formats carry no tag.
    bool WriteImage8(const std::string& path, int width, int height, int channels, const std::vector<unsigned char>& data, color::Gamut tag = color::Gamut::SRGB);

    /// A 16-bit PNG (0..65535 per sample), 1 or 3 channels, with a cICP chunk as WriteImage8.
    bool WritePng16(const std::string& path, int width, int height, int channels, const std::vector<std::uint16_t>& data, color::Gamut tag = color::Gamut::SRGB);

    /// A PFM: one channel ("Pf") or three ("PF") floats, little endian, written bottom row first as the format wants. Non-finite values become 0.
    bool WritePfm(const std::string& path, int width, int height, int channels, const std::vector<float>& data);

    /// A Radiance .hdr (RGBE) of three float channels.
    bool WriteHdr(const std::string& path, int width, int height, const std::vector<float>& rgb);

} // namespace rayt::io
