#pragma once
/**
 * @file DngWriter.hpp
 * @brief A raw file in Adobe's Digital Negative format (DNG 1.4): the sensor's Bayer mosaic as recorded, for any raw developer.
 *
 * A DNG is a TIFF file: an 8-byte header ("II", 42, the offset of the first directory) and an image file directory (IFD) of tags, here
 * with the raw image itself in it (uncompressed, 16 bits, one strip). The tags a raw developer needs besides the size: the colour filter
 * pattern (CFARepeatPatternDim 2 x 2, CFAPattern R G G B), the black and white levels of the DN, ColorMatrix1 (CIE XYZ to the camera's
 * raw values under CalibrationIlluminant1, 21 = D65) and AsShotNeutral (the raw values of the white the camera was balanced to). From
 * them a developer subtracts black, divides by white - black, balances, demosaics and turns the camera's colours into XYZ.
 */

#include <cstdint>
#include <string>
#include <vector>

namespace rayt::io {

    struct DngInfo {
        int black = 0;                     ///< BlackLevel (DN of no light)
        int white = 65535;                 ///< WhiteLevel (the top DN)
        double colorMatrix[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };   ///< ColorMatrix1: XYZ -> camera raw, row by row
        double asShotNeutral[3] = { 1, 1, 1 };                     ///< AsShotNeutral: the white's raw values (largest 1)
        std::string model = "GoLD_rayt";   ///< UniqueCameraModel and Model
        std::string description;           ///< ImageDescription ("" = none)
    };

    /// Writes an RGGB mosaic (width x height DN, row by row from the top) as a DNG. False if the file cannot be written.
    bool WriteDng(const std::string& path, int width, int height, const std::vector<std::uint16_t>& mosaic, const DngInfo& info);

} // namespace rayt::io
