/*
* @file src/Render/Film.cpp
* @brief Film class implementation (pixel storage, accumulation, and the values of the written pictures).
*
* The Film stores per-pixel radiance in linear RGB and turns it into the values of the files (gamut, exposure, tone map, 8 / 16 bit);
* the files themselves are encoded by IO/ImageWriter.
*/

#include "pch.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Core/Math.hpp"
#include "Render/Film.hpp"
#include "Color/ColorSpace.hpp"
#include "Color/ColorTransform.hpp"
#include "IO/ImageWriter.hpp"

namespace rayt {

    namespace {

        /// A value in [0, 1] (saturated) to 8 bits, rounded to nearest.
        unsigned char ToByte(Real x) {
            return static_cast<unsigned char>(Real(255.0) * rayt::math::saturate(x) + Real(0.5));
        }

        /// A value in [0, 1] (saturated) to 16 bits, rounded to nearest.
        std::uint16_t To16(Real x) {
            return static_cast<std::uint16_t>(Real(65535.0) * rayt::math::saturate(x) + Real(0.5));
        }

        std::string LowerExtension(const std::string& filename) {
            const size_t dot = filename.find_last_of('.');
            std::string ext = (dot == std::string::npos) ? std::string() : filename.substr(dot + 1);
            std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return ext;
        }

        float Finite(Real v) { return math::hasNonFinite(v) ? 0.0f : static_cast<float>(v); }

    } // namespace

    /**
     * @brief Constructs a Film with the given resolution; all pixels are black.
     */
    Film::Film(int width, int height)
        : m_width(width), m_height(height) {
        m_pixels.resize(static_cast<size_t>(width) * static_cast<size_t>(height), Vector3(0.0));
    }

    /// Overwrites the radiance of a pixel ((0,0) is the top-left corner); outside the film: ignored.
    void Film::setPixel(int x, int y, const Vector3& radiance) {
        if (x < 0 || x >= m_width || y < 0 || y >= m_height) return;
        m_pixels[static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)] = radiance;
    }

    /// Adds radiance to a pixel (progressive rendering); outside the film: ignored.
    void Film::addPixel(int x, int y, const Vector3& radiance) {
        if (x < 0 || x >= m_width || y < 0 || y >= m_height) return;
        m_pixels[static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x)] += radiance;
    }

    /// Resets the film to black (restarting a progressive render).
    void Film::clear() {
        std::fill(m_pixels.begin(), m_pixels.end(), Vector3(0.0));
    }

    void Film::enableAOVs() {
        const size_t n = static_cast<size_t>(m_width) * static_cast<size_t>(m_height);
        m_albedo.assign(n, Vector3(0.0));
        m_normal.assign(n, Vector3(0.0));
    }

    void Film::enableStokes() {
        const size_t n = static_cast<size_t>(m_width) * static_cast<size_t>(m_height);
        for (auto& b : m_stokes) b.assign(n, Vector3(0.0));
    }

    void Film::setStokes(int x, int y, const Vector3 s[4]) {
        if (!hasStokes() || x < 0 || x >= m_width || y < 0 || y >= m_height) return;
        const size_t i = static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x);
        for (int k = 0; k < 4; ++k) m_stokes[k][i] = s[k];
    }

    void Film::saveStokes(const std::string& stem) const {
        if (!hasStokes()) return;
        for (int k = 0; k < 4; ++k) {
            const std::string filename = stem + "_S" + std::to_string(k) + ".pfm";
            std::vector<float> rgb(m_stokes[k].size() * 3);
            for (size_t i = 0; i < m_stokes[k].size(); ++i)
                for (int c = 0; c < 3; ++c) rgb[i * 3 + c] = Finite(m_stokes[k][i][c]);
            if (!io::WritePfm(filename, m_width, m_height, 3, rgb)) { std::cerr << "[Film] Error: could not write " << filename << std::endl; continue; }
            std::cout << "[Film] Saved Stokes image: " << filename << std::endl;
        }
    }

    void Film::polarizerImage(Real angleDeg, Film& out) const {
        if (!hasStokes() || out.m_width != m_width || out.m_height != m_height) return;
        const Real a = angleDeg * constants::PI / Real(180);
        const Real c2 = std::cos(Real(2) * a), s2 = std::sin(Real(2) * a);
        for (int y = 0; y < m_height; ++y)
            for (int x = 0; x < m_width; ++x) {
                const size_t i = static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x);
                out.m_pixels[i] = Real(0.5) * (m_stokes[0][i] + m_stokes[1][i] * c2 + m_stokes[2][i] * s2);
            }
    }

    void Film::saveLinearPolarization(const std::string& stem) const {
        if (!hasStokes()) return;
        const size_t n = static_cast<size_t>(m_width) * static_cast<size_t>(m_height);
        std::vector<Vector3> dolp(n), aolp(n);
        auto luma = [](const Vector3& c) { return Real(0.2126) * c.x + Real(0.7152) * c.y + Real(0.0722) * c.z; };
        for (size_t i = 0; i < n; ++i) {
            const Real s0 = luma(m_stokes[0][i]), s1 = luma(m_stokes[1][i]), s2 = luma(m_stokes[2][i]);
            const Real d = s0 > Real(1e-6) ? std::min(Real(1), std::sqrt(s1 * s1 + s2 * s2) / s0) : Real(0);
            dolp[i] = Vector3(d);
            // Hue = angle of polarization (0 .. 180 degrees round the color wheel), value = degree.
            const Real angle = Real(0.5) * std::atan2(s2, s1);          // -90 .. 90 degrees
            Real h = angle / constants::PI + Real(0.5);                  // 0 .. 1
            if (h >= Real(1)) h = Real(0);
            const Real hh = h * Real(6);
            const int sector = std::min(5, static_cast<int>(hh));
            const Real f = hh - Real(sector);
            const Real v = d, p = 0, q = v * (Real(1) - f), t = v * f;
            Vector3 rgb;
            switch (sector) {
            case 0: rgb = Vector3(v, t, p); break;
            case 1: rgb = Vector3(q, v, p); break;
            case 2: rgb = Vector3(p, v, t); break;
            case 3: rgb = Vector3(p, q, v); break;
            case 4: rgb = Vector3(t, p, v); break;
            default: rgb = Vector3(v, p, q); break;
            }
            aolp[i] = rgb;
        }
        saveBuffer(stem + "_dolp.png", dolp, /*encodeSrgb=*/false, /*remapSigned=*/false);
        saveBuffer(stem + "_aolp.png", aolp, /*encodeSrgb=*/false, /*remapSigned=*/false);
    }

    void Film::setAOV(int x, int y, const Vector3& albedo, const Vector3& normal) {
        if (!hasAOVs() || x < 0 || x >= m_width || y < 0 || y >= m_height) return;
        const size_t i = static_cast<size_t>(y) * static_cast<size_t>(m_width) + static_cast<size_t>(x);
        m_albedo[i] = albedo;
        m_normal[i] = normal;
    }

    void Film::save(const std::string& filename) const {
        save(filename, 1.0f);
    }

    /**
     * @brief Saves the film, every radiance multiplied by `scale` (e.g. 1 / spp of a progressive render). The extension picks the format:
     *
     * - .pfm / .hdr: linear floats in the output gamut, no exposure or tone map. A .pfm keeps what .hdr (RGBE: one shared exponent and an
     *   8-bit mantissa per pixel) loses: the small channels of a strongly colored pixel and negative values (colors outside sRGB).
     * - .png / .bmp / .jpg: the display values (setDisplay): primaries, exposure, Reinhard or a clip at white, sRGB curve, clamp; 8 bits,
     *   or 16 for a .png with setPng16. A .png in another gamut than sRGB is cICP-tagged; .bmp / .jpg carry no tag.
     */
    void Film::save(const std::string& filename, float scale) const {
        const std::string ext = LowerExtension(filename);
        const Real s = static_cast<Real>(scale);
        const int numPixels = static_cast<int>(m_pixels.size());

        if (ext == "pfm" || ext == "hdr") {
            std::vector<float> rgb(static_cast<size_t>(numPixels) * 3);
#pragma omp parallel for schedule(static)
            for (int i = 0; i < numPixels; ++i) {
                const Vector3 p = color::fromLinearSRGB(m_gamut, m_pixels[i] * s);
                for (int c = 0; c < 3; ++c) rgb[static_cast<size_t>(i) * 3 + c] = Finite(p[c]);
            }
            const bool ok = (ext == "pfm") ? io::WritePfm(filename, m_width, m_height, 3, rgb) : io::WriteHdr(filename, m_width, m_height, rgb);
            if (!ok) { std::cerr << "[Film] Error: could not write " << filename << std::endl; return; }
            std::cout << "[Film] Saved " << (ext == "pfm" ? "PFM" : "HDR") << " image: " << filename << std::endl;
            return;
        }
        if (ext != "png" && ext != "bmp" && ext != "jpg") {
            std::cerr << "[Film] Error: Unsupported file extension: " << ext << std::endl;
            return;
        }

        // --exposure / --tonemap of the CLI (the GUI applies its own exposure to its display copy).
        const bool clip = (m_toneMap == ToneMap::Clip);
        const bool want16 = m_png16 && ext == "png";
        std::vector<unsigned char> data8(want16 ? 0 : static_cast<size_t>(numPixels) * 3);
        std::vector<std::uint16_t> data16(want16 ? static_cast<size_t>(numPixels) * 3 : 0);
#pragma omp parallel for schedule(static)
        for (int i = 0; i < numPixels; ++i) {
            // Change of primaries first (exact, unclamped), then the display pipeline, which clips what is still negative:
            // exposure -> Reinhard (or a clip at white) -> sRGB curve -> clamp. The curve is the sRGB one for every gamut (as for Display P3).
            Vector3 pixel = color::fromLinearSRGB(m_gamut, m_pixels[i] * s);
            pixel = clip ? color::toDisplaySRGB_Clip(pixel, m_exposureStops) : color::toDisplaySRGB_Reinhard(pixel, m_exposureStops);
            for (int c = 0; c < 3; ++c) {
                if (want16) data16[static_cast<size_t>(i) * 3 + c] = To16(pixel[c]);
                else data8[static_cast<size_t>(i) * 3 + c] = ToByte(pixel[c]);
            }
        }
        const bool ok = want16 ? io::WritePng16(filename, m_width, m_height, 3, data16, m_gamut) : io::WriteImage8(filename, m_width, m_height, 3, data8, m_gamut);
        if (!ok) { std::cerr << "[Film] Error: could not write " << filename << std::endl; return; }

        std::cout << "[Film] Saved LDR image: " << filename << (want16 ? " (16 bit)" : "");
        if (m_gamut != color::Gamut::SRGB && ext != "png")
            std::cout << " (" << color::gamutName(m_gamut) << " values, but this format carries no color tag; use .png)";
        else if (m_gamut != color::Gamut::SRGB)
            std::cout << " (" << color::gamutName(m_gamut) << ", cICP-tagged)";
        std::cout << std::endl;

        if (clip) {
            const WhiteClip wc = whiteClip(scale);
            std::ostringstream line;   // its own stream: the format flags must not stay on std::cout
            line << "[Film] clip at white: " << std::fixed << std::setprecision(2) << Real(100) * wc.share << " % of the pixels have a channel above white";
            if (wc.brightest > Real(0))
                line << " (the brightest is " << std::showpos << std::log2(wc.brightest) << " stops from white)";
            std::cout << line.str() << std::endl;
        }
    }

    Film::WhiteClip Film::whiteClip(float scale) const {
        const Real gain = static_cast<Real>(scale) * std::pow(Real(2), m_exposureStops);
        size_t over = 0;
        WhiteClip wc;
        for (const Vector3& px : m_pixels) {
            const Vector3 p = color::fromLinearSRGB(m_gamut, px * gain);
            const Real m = std::max(p[0], std::max(p[1], p[2]));
            if (math::hasNonFinite(m)) continue;
            if (m > Real(1)) ++over;
            wc.brightest = std::max(wc.brightest, m);
        }
        wc.share = m_pixels.empty() ? Real(0) : Real(over) / Real(m_pixels.size());
        return wc;
    }

    void Film::saveMask(const std::string& filename) const {
        if (m_mask.empty()) { std::cerr << "[Film] Error: the mask is not enabled, nothing to save: " << filename << std::endl; return; }
        const std::string ext = LowerExtension(filename);
        bool ok = false;
        if (ext == "pfm") ok = io::WritePfm(filename, m_width, m_height, 1, m_mask);
        else if (ext == "png" && m_png16) {
            std::vector<std::uint16_t> v(m_mask.size());
            for (size_t i = 0; i < v.size(); ++i) v[i] = To16(Real(m_mask[i]));
            ok = io::WritePng16(filename, m_width, m_height, 1, v);
        }
        else if (ext == "png") {
            std::vector<unsigned char> v(m_mask.size());
            for (size_t i = 0; i < v.size(); ++i) v[i] = ToByte(Real(m_mask[i]));
            ok = io::WriteImage8(filename, m_width, m_height, 1, v);
        }
        else { std::cerr << "[Film] Error: a mask is written as .png or .pfm: " << filename << std::endl; return; }
        if (!ok) { std::cerr << "[Film] Error: could not write " << filename << std::endl; return; }
        std::cout << "[Film] Saved mask: " << filename << (ext == "png" && m_png16 ? " (16 bit)" : "") << std::endl;
    }

    void Film::saveDepth(const std::string& filename) const {
        if (m_depth.empty()) { std::cerr << "[Film] Error: the depth is not enabled, nothing to save: " << filename << std::endl; return; }
        if (!io::WritePfm(filename, m_width, m_height, 1, m_depth)) { std::cerr << "[Film] Error: could not write " << filename << std::endl; return; }
        std::cout << "[Film] Saved depth: " << filename << std::endl;
    }

    void Film::saveAlbedo(const std::string& filename) const {
        saveBuffer(filename, m_albedo, /*encodeSrgb=*/true, /*remapSigned=*/false);
    }

    void Film::saveNormal(const std::string& filename) const {
        saveBuffer(filename, m_normal, /*encodeSrgb=*/false, /*remapSigned=*/true);
    }

    void Film::saveBuffer(const std::string& filename, const std::vector<Vector3>& data, bool encodeSrgb, bool remapSigned) const {
        if (data.empty()) {
            std::cerr << "[Film] Error: AOVs are not enabled, nothing to save: " << filename << std::endl;
            return;
        }
        const std::string ext = LowerExtension(filename);
        if (ext == "hdr") {
            std::vector<float> rgb(data.size() * 3);
            for (size_t i = 0; i < data.size(); ++i)
                for (int c = 0; c < 3; ++c) rgb[i * 3 + c] = Finite(data[i][c]);
            if (!io::WriteHdr(filename, m_width, m_height, rgb)) { std::cerr << "[Film] Error: could not write " << filename << std::endl; return; }
            std::cout << "[Film] Saved HDR image: " << filename << std::endl;
            return;
        }
        if (ext != "png" && ext != "bmp" && ext != "jpg") {
            std::cerr << "[Film] Error: Unsupported file extension: " << ext << std::endl;
            return;
        }
        std::vector<unsigned char> out(data.size() * 3);
        for (size_t i = 0; i < data.size(); ++i)
            for (int c = 0; c < 3; ++c) {
                Real v = data[i][c];
                if (math::hasNonFinite(v)) v = Real(0);
                if (remapSigned) v = v * Real(0.5) + Real(0.5);
                else if (encodeSrgb) v = color::linearToSRGB(rayt::math::saturate(v));
                out[i * 3 + c] = ToByte(v);
            }
        if (!io::WriteImage8(filename, m_width, m_height, 3, out)) { std::cerr << "[Film] Error: could not write " << filename << std::endl; return; }
        std::cout << "[Film] Saved LDR image: " << filename << std::endl;
    }

} // namespace rayt
