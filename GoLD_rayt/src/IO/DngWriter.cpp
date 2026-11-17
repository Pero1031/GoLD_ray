// A raw file in the DNG format; see IO/DngWriter.hpp.
#include "pch.h"
#include "IO/DngWriter.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>

namespace rayt::io {

    namespace {

        enum : std::uint16_t { kByte = 1, kAscii = 2, kShort = 3, kLong = 4, kRational = 5, kSRational = 10 };

        struct Entry {
            std::uint16_t tag = 0, type = 0;
            std::uint32_t count = 0;
            std::vector<std::uint8_t> data;   ///< the value's bytes, little-endian
        };

        void Put16(std::vector<std::uint8_t>& b, std::uint32_t v) { b.push_back(std::uint8_t(v & 0xff)); b.push_back(std::uint8_t((v >> 8) & 0xff)); }
        void Put32(std::vector<std::uint8_t>& b, std::uint32_t v) { for (int k = 0; k < 4; ++k) b.push_back(std::uint8_t((v >> (8 * k)) & 0xff)); }

    } // namespace

    bool WriteDng(const std::string& path, int width, int height, const std::vector<std::uint16_t>& mosaic, const DngInfo& info) {
        if (width <= 0 || height <= 0 || mosaic.size() != static_cast<size_t>(width) * static_cast<size_t>(height)) return false;
        std::vector<Entry> entries;
        auto add = [&](std::uint16_t tag, std::uint16_t type, std::uint32_t count, std::vector<std::uint8_t> data) {
            entries.push_back(Entry{ tag, type, count, std::move(data) });
        };
        auto shorts = [&](std::uint16_t tag, std::initializer_list<std::uint32_t> v) { std::vector<std::uint8_t> d; for (auto x : v) Put16(d, x); add(tag, kShort, std::uint32_t(v.size()), d); };
        auto longs = [&](std::uint16_t tag, std::initializer_list<std::uint32_t> v) { std::vector<std::uint8_t> d; for (auto x : v) Put32(d, x); add(tag, kLong, std::uint32_t(v.size()), d); };
        auto bytes = [&](std::uint16_t tag, std::initializer_list<std::uint8_t> v) { add(tag, kByte, std::uint32_t(v.size()), std::vector<std::uint8_t>(v)); };
        auto ascii = [&](std::uint16_t tag, const std::string& s) { std::vector<std::uint8_t> d(s.begin(), s.end()); d.push_back(0); add(tag, kAscii, std::uint32_t(d.size()), d); };

        const std::uint32_t w = std::uint32_t(width), h = std::uint32_t(height), imageBytes = w * h * 2;
        longs(254, { 0 });                    // NewSubFileType: the main image
        longs(256, { w });                    // ImageWidth
        longs(257, { h });                    // ImageLength
        shorts(258, { 16 });                  // BitsPerSample
        shorts(259, { 1 });                   // Compression: none
        shorts(262, { 32803 });               // PhotometricInterpretation: CFA
        if (!info.description.empty()) ascii(270, info.description);   // ImageDescription
        ascii(271, "GoLD_rayt");              // Make
        ascii(272, info.model);               // Model
        longs(273, { 0 });                    // StripOffsets (set below)
        shorts(274, { 1 });                   // Orientation: row 0 at the top, column 0 on the left
        shorts(277, { 1 });                   // SamplesPerPixel
        longs(278, { h });                    // RowsPerStrip: one strip
        longs(279, { imageBytes });           // StripByteCounts
        shorts(284, { 1 });                   // PlanarConfiguration
        ascii(305, "GoLD_rayt");              // Software
        shorts(33421, { 2, 2 });              // CFARepeatPatternDim
        bytes(33422, { 0, 1, 1, 2 });         // CFAPattern: R G / G B
        bytes(50706, { 1, 4, 0, 0 });         // DNGVersion
        bytes(50707, { 1, 1, 0, 0 });         // DNGBackwardVersion
        ascii(50708, info.model);             // UniqueCameraModel
        bytes(50710, { 0, 1, 2 });            // CFAPlaneColor: red, green, blue
        shorts(50711, { 1 });                 // CFALayout: rectangular
        longs(50714, { std::uint32_t(info.black) });   // BlackLevel
        longs(50717, { std::uint32_t(info.white) });   // WhiteLevel
        auto matrix = [&](std::uint16_t tag, const double (&m)[9]) {   // nine signed rationals over 10000
            std::vector<std::uint8_t> d;
            for (const double v : m) { Put32(d, std::uint32_t(std::int32_t(std::lround(v * 10000.0)))); Put32(d, 10000); }
            add(tag, kSRational, 9, d);
        };
        matrix(50721, info.colorMatrixA);     // ColorMatrix1 (under illuminant A)
        matrix(50722, info.colorMatrix);      // ColorMatrix2 (under D65)
        {
            std::vector<std::uint8_t> d;      // AsShotNeutral: three rationals over 1000000
            for (const double v : info.asShotNeutral) { Put32(d, std::uint32_t(std::lround(std::max(v, 0.0) * 1000000.0))); Put32(d, 1000000); }
            add(50728, kRational, 3, d);
        }
        shorts(50778, { 17 });                // CalibrationIlluminant1: Standard light A
        shorts(50779, { 21 });                // CalibrationIlluminant2: D65
        matrix(50964, info.forwardMatrixA);   // ForwardMatrix1
        matrix(50965, info.forwardMatrix);    // ForwardMatrix2
        std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) { return a.tag < b.tag; });

        // Layout: header (8), the directory, the values longer than 4 bytes (each at an even offset), the image (even offset).
        const std::uint32_t ifdOffset = 8, ifdSize = 2 + 12 * std::uint32_t(entries.size()) + 4;
        std::uint32_t offset = ifdOffset + ifdSize;
        std::vector<std::uint32_t> valueOffset(entries.size(), 0);
        for (size_t i = 0; i < entries.size(); ++i)
            if (entries[i].data.size() > 4) { valueOffset[i] = offset; offset += std::uint32_t(entries[i].data.size()); offset += offset & 1; }
        const std::uint32_t imageOffset = offset;
        for (Entry& e : entries)
            if (e.tag == 273) { e.data.clear(); Put32(e.data, imageOffset); }

        std::vector<std::uint8_t> file;
        file.reserve(size_t(imageOffset) + imageBytes);
        file.push_back('I'); file.push_back('I'); Put16(file, 42); Put32(file, ifdOffset);
        Put16(file, std::uint32_t(entries.size()));
        for (size_t i = 0; i < entries.size(); ++i) {
            const Entry& e = entries[i];
            Put16(file, e.tag); Put16(file, e.type); Put32(file, e.count);
            if (e.data.size() > 4) Put32(file, valueOffset[i]);
            else { std::vector<std::uint8_t> v = e.data; v.resize(4, 0); file.insert(file.end(), v.begin(), v.end()); }   // left-justified
        }
        Put32(file, 0);   // no next directory
        for (size_t i = 0; i < entries.size(); ++i)
            if (entries[i].data.size() > 4) {
                file.insert(file.end(), entries[i].data.begin(), entries[i].data.end());
                if (file.size() & 1) file.push_back(0);
            }
        for (const std::uint16_t v : mosaic) Put16(file, v);

        std::ofstream out(path, std::ios::binary);
        if (!out) return false;
        out.write(reinterpret_cast<const char*>(file.data()), std::streamsize(file.size()));
        return bool(out);
    }

} // namespace rayt::io
