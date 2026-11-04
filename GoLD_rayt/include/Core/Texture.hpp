/**
 * @file Core/Texture.hpp
 * @brief 8-bit RGBA image texture with bilinear, repeating lookup (header-only).
 *
 * - Stored as RGBA8 (4 bytes per texel); colors are converted from sRGB to linear at lookup with a
 *   256-entry table, so a 4096 x 4096 texture costs 64 MB, not 256 MB.
 * - Texture coordinates follow the DirectX / MMD convention: (0, 0) is the TOP-left texel, v grows
 *   downward, and coordinates repeat (wrap) outside [0, 1].
 * - The file is read into memory first and decoded with stb_image, so paths with non-ASCII
 *   characters (CJK file names inside MMD models) work: stbi_load() would use the ANSI code page.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <vector>

#include "Core/Types.hpp"
#include "stb_image.h"   // declarations only; the implementation lives in src/IO/ImageIO.cpp

namespace rayt {

    class Texture2D {
    public:
        /// Loads PNG / JPG / BMP / TGA ... Returns nullptr if the file cannot be read or decoded.
        static std::shared_ptr<Texture2D> load(const std::filesystem::path& file) {
            std::ifstream f(file, std::ios::binary);
            if (!f) return nullptr;
            const std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
            if (bytes.empty()) return nullptr;

            int w = 0, h = 0, channels = 0;
            unsigned char* data = stbi_load_from_memory(bytes.data(), static_cast<int>(bytes.size()), &w, &h, &channels, 4);
            if (!data) return nullptr;

            auto tex = std::make_shared<Texture2D>();
            tex->m_w = w;
            tex->m_h = h;
            tex->m_rgba.assign(data, data + static_cast<size_t>(w) * static_cast<size_t>(h) * 4);
            stbi_image_free(data);

            // Does the texture contain any non-opaque texel? (decides whether an alpha test is needed)
            for (size_t i = 3; i < tex->m_rgba.size(); i += 4) {
                if (tex->m_rgba[i] < 255) { tex->m_hasAlpha = true; break; }
            }
            return tex;
        }

        /// Builds a texture from 8-bit sRGB RGBA texels (row 0 = top), e.g. a procedural pattern. nullptr if the size does not match.
        static std::shared_ptr<Texture2D> fromRgba8(int w, int h, std::vector<unsigned char> rgba) {
            if (w <= 0 || h <= 0 || rgba.size() != static_cast<size_t>(w) * static_cast<size_t>(h) * 4) return nullptr;
            auto tex = std::make_shared<Texture2D>();
            tex->m_w = w;
            tex->m_h = h;
            tex->m_rgba = std::move(rgba);
            for (size_t i = 3; i < tex->m_rgba.size(); i += 4) {
                if (tex->m_rgba[i] < 255) { tex->m_hasAlpha = true; break; }
            }
            return tex;
        }

        int width() const { return m_w; }
        int height() const { return m_h; }
        bool hasAlpha() const { return m_hasAlpha; }

        /// Linear RGB, bilinearly filtered.
        Vector3 rgb(Real u, Real v) const {
            Vector3 c(0.0);
            Real unused;
            sample(u, v, &c, unused);
            return c;
        }

        /// Alpha in [0, 1], bilinearly filtered.
        Real alpha(Real u, Real v) const {
            Vector3 unused;
            Real a = 1;
            sample(u, v, &unused, a, /*wantColor=*/false);
            return a;
        }

    private:
        static const std::array<float, 256>& srgbTable() {
            static const std::array<float, 256> table = [] {
                std::array<float, 256> t{};
                for (int i = 0; i < 256; ++i) {
                    const double c = i / 255.0;
                    t[static_cast<size_t>(i)] = static_cast<float>(c <= 0.04045 ? c / 12.92 : std::pow((c + 0.055) / 1.055, 2.4));
                }
                return t;
            }();
            return table;
        }

        void sample(Real u, Real v, Vector3* rgbOut, Real& alphaOut, bool wantColor = true) const {
            if (m_w <= 0 || m_h <= 0) { alphaOut = 1; return; }
            const Real x = (u - std::floor(u)) * Real(m_w) - Real(0.5);   // repeat, texel centers at +0.5
            const Real y = (v - std::floor(v)) * Real(m_h) - Real(0.5);
            const int x0 = static_cast<int>(std::floor(x)), y0 = static_cast<int>(std::floor(y));
            const Real fx = x - Real(x0), fy = y - Real(y0);

            const auto& lut = srgbTable();
            Real r = 0, g = 0, b = 0, a = 0;
            for (int dy = 0; dy < 2; ++dy) {
                for (int dx = 0; dx < 2; ++dx) {
                    const Real w = (dx ? fx : Real(1) - fx) * (dy ? fy : Real(1) - fy);
                    const int px = ((x0 + dx) % m_w + m_w) % m_w;
                    const int py = ((y0 + dy) % m_h + m_h) % m_h;
                    const unsigned char* t = &m_rgba[(static_cast<size_t>(py) * static_cast<size_t>(m_w) + static_cast<size_t>(px)) * 4];
                    if (wantColor) { r += w * lut[t[0]]; g += w * lut[t[1]]; b += w * lut[t[2]]; }
                    a += w * (t[3] / Real(255));
                }
            }
            if (wantColor && rgbOut) *rgbOut = Vector3(r, g, b);
            alphaOut = a;
        }

        int m_w = 0, m_h = 0;
        bool m_hasAlpha = false;
        std::vector<unsigned char> m_rgba;
    };

} // namespace rayt
