// Encoders of image files (see IO/ImageWriter.hpp). stb_image_write is compiled in IO/ImageIO.cpp.
#include "pch.h"
#include "IO/ImageWriter.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>

#include "stb_image_write.h"

namespace rayt::io {

    namespace {

        std::string LowerExtension(const std::string& path) {
            const size_t dot = path.find_last_of('.');
            std::string ext = (dot == std::string::npos) ? std::string() : path.substr(dot + 1);
            std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return ext;
        }

        /// CRC-32 of the PNG spec (polynomial 0xEDB88320), over a chunk's type and data.
        std::uint32_t PngCrc32(const unsigned char* data, size_t n) {
            std::uint32_t crc = 0xFFFFFFFFu;
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
        bool InsertCicp(std::vector<unsigned char>& png, color::Gamut gamut) {
            constexpr size_t kAfterIhdr = 8 + 25;   // signature + (length 4, type 4, data 13, crc 4)
            if (png.size() < kAfterIhdr) return false;
            unsigned char chunk[4 + 4 + 4 + 4];
            chunk[0] = 0; chunk[1] = 0; chunk[2] = 0; chunk[3] = 4;   // data length
            chunk[4] = 'c'; chunk[5] = 'I'; chunk[6] = 'C'; chunk[7] = 'P';
            chunk[8] = static_cast<unsigned char>(color::gamutPrimariesCode(gamut));
            chunk[9] = 13;   // transfer: sRGB curve
            chunk[10] = 0;   // matrix: RGB
            chunk[11] = 1;   // full range
            const std::uint32_t crc = PngCrc32(chunk + 4, 8);
            chunk[12] = static_cast<unsigned char>(crc >> 24); chunk[13] = static_cast<unsigned char>(crc >> 16);
            chunk[14] = static_cast<unsigned char>(crc >> 8);  chunk[15] = static_cast<unsigned char>(crc);
            png.insert(png.begin() + kAfterIhdr, chunk, chunk + sizeof(chunk));
            return true;
        }

        bool WriteBytes(const std::string& path, const std::vector<unsigned char>& bytes) {
            std::ofstream f(path, std::ios::binary);
            if (!f) return false;
            f.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
            return static_cast<bool>(f);
        }

        void PngSink(void* ctx, void* data, int size) {
            auto* v = static_cast<std::vector<unsigned char>*>(ctx);
            const auto* p = static_cast<const unsigned char*>(data);
            v->insert(v->end(), p, p + size);
        }

    } // namespace

    bool WriteImage8(const std::string& path, int width, int height, int channels, const std::vector<unsigned char>& data, color::Gamut tag) {
        const std::string ext = LowerExtension(path);
        if (ext == "png") {
            if (tag == color::Gamut::SRGB || channels != 3)
                return stbi_write_png(path.c_str(), width, height, channels, data.data(), width * channels) != 0;
            std::vector<unsigned char> png;
            if (!stbi_write_png_to_func(PngSink, &png, width, height, channels, data.data(), width * channels)) return false;
            return InsertCicp(png, tag) && WriteBytes(path, png);
        }
        if (ext == "bmp") return stbi_write_bmp(path.c_str(), width, height, channels, data.data()) != 0;
        if (ext == "jpg") return stbi_write_jpg(path.c_str(), width, height, channels, data.data(), 90) != 0;
        return false;
    }

    /**
     * stb_image_write has no 16-bit mode. A row of W pixels at 16 bits is 2 * channels * W bytes: exactly a row of 2 * W pixels at 8 bits.
     * So the big-endian samples are encoded as an 8-bit image of twice the width with the PNG filter forced to "none" (the filters work on
     * the bytes of a whole pixel, which would be wrong here), and the IHDR is then patched to width W, 16 bits (and its CRC redone).
     */
    bool WritePng16(const std::string& path, int width, int height, int channels, const std::vector<std::uint16_t>& data, color::Gamut tag) {
        std::vector<unsigned char> bytes(data.size() * 2);
        for (size_t i = 0; i < data.size(); ++i) {
            bytes[i * 2] = static_cast<unsigned char>(data[i] >> 8);
            bytes[i * 2 + 1] = static_cast<unsigned char>(data[i] & 0xFFu);
        }
        std::vector<unsigned char> png;
        const int savedFilter = stbi_write_force_png_filter;
        stbi_write_force_png_filter = 0;
        const int ok = stbi_write_png_to_func(PngSink, &png, 2 * width, height, channels, bytes.data(), 2 * channels * width);
        stbi_write_force_png_filter = savedFilter;
        if (!ok || png.size() < 33) return false;
        for (int k = 0; k < 4; ++k) png[16 + k] = static_cast<unsigned char>(static_cast<std::uint32_t>(width) >> (8 * (3 - k)));   // width
        png[24] = 16;                                                                                                              // bit depth
        const std::uint32_t crc = PngCrc32(png.data() + 12, 17);                                                                   // type + data
        for (int k = 0; k < 4; ++k) png[29 + k] = static_cast<unsigned char>(crc >> (8 * (3 - k)));
        if (tag != color::Gamut::SRGB && channels == 3 && !InsertCicp(png, tag)) return false;
        return WriteBytes(path, png);
    }

    bool WritePfm(const std::string& path, int width, int height, int channels, const std::vector<float>& data) {
        const size_t row = static_cast<size_t>(width) * static_cast<size_t>(channels);
        std::vector<float> rows(data.size());
        for (int y = 0; y < height; ++y)   // PFM rows go from the bottom up
            for (size_t k = 0; k < row; ++k) {
                const float v = data[static_cast<size_t>(y) * row + k];
                rows[static_cast<size_t>(height - 1 - y) * row + k] = std::isfinite(v) ? v : 0.0f;
            }
        std::ofstream f(path, std::ios::binary);
        if (!f) return false;
        f << (channels == 1 ? "Pf" : "PF") << "\n" << width << " " << height << "\n-1.0\n";
        f.write(reinterpret_cast<const char*>(rows.data()), static_cast<std::streamsize>(rows.size() * sizeof(float)));
        return static_cast<bool>(f);
    }

    bool WriteHdr(const std::string& path, int width, int height, const std::vector<float>& rgb) {
        return stbi_write_hdr(path.c_str(), width, height, 3, rgb.data()) != 0;
    }

} // namespace rayt::io
