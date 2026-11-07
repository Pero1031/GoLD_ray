/**
 * @file include/Render/Denoiser.hpp
 * @brief Intel Open Image Denoise (OIDN) wrapper: denoises the beauty buffer of a Film.
 *
 * The beauty image is passed as linear HDR (before tone mapping), together with the albedo and
 * normal AOVs when the film has them (Film::enableAOVs): they tell the network which edges and
 * textures are real. The result is a separate Film, so the raw image stays available (a denoised
 * image is a learned estimate, not an unbiased one: keep the raw image for numerical checks).
 *
 * Negative values: OIDN clips its input to >= 0. The film holds linear sRGB, and a color outside the sRGB triangle (a candle-lit
 * gold, a saturated LED) has a negative channel there, which the denoiser would turn into 0. DenoiseFilm therefore takes the
 * gamut of the output: the image is converted to that gamut's linear RGB (positive for colors inside it), denoised, and converted
 * back, so the result is again linear sRGB with its negative values (the denoiser's own smoothing of them is exact: the
 * conversion is linear). Only colors outside the working gamut itself are still clipped.
 *
 * Needs OpenImageDenoise.dll (and its device / TBB DLLs) next to the executable.
 */

#pragma once

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

#include <OpenImageDenoise/oidn.hpp>

#include "Core/Types.hpp"
#include "Render/Film.hpp"

namespace rayt {

    namespace detail {
        /// Linear sRGB pixels -> floats, in the linear RGB of `gamut` (sRGB: unchanged).
        inline std::vector<float> toFloat3(const Vector3* data, size_t pixels, color::Gamut gamut = color::Gamut::SRGB) {
            std::vector<float> out(pixels * 3);
            for (size_t i = 0; i < pixels; ++i) {
                Vector3 p = data[i];
                for (int c = 0; c < 3; ++c) if (!std::isfinite(p[c])) p[c] = Real(0);
                p = color::fromLinearSRGB(gamut, p);
                for (int c = 0; c < 3; ++c) out[i * 3 + static_cast<size_t>(c)] = static_cast<float>(p[c]);
            }
            return out;
        }
    } // namespace detail

    /**
     * @brief Denoises `in` into `out` (same size). Returns false and prints the reason on failure.
     *
     * @param cleanAux true if the albedo / normal AOVs are noise-free (true for pixel-averaged AOVs
     *                 at a reasonable spp); false lets OIDN prefilter them, which is safer at low spp.
     * @param workingGamut the gamut the image is denoised in (the output gamut): see the note on negative values above.
     */
    inline bool DenoiseFilm(const Film& in, Film& out, bool cleanAux = false, color::Gamut workingGamut = color::Gamut::SRGB) {
        const int w = in.width(), h = in.height();
        if (out.width() != w || out.height() != h) {
            std::cerr << "[Denoise] output film has a different size\n";
            return false;
        }
        const size_t pixels = static_cast<size_t>(w) * static_cast<size_t>(h);

        std::vector<float> color = detail::toFloat3(in.getData(), pixels, workingGamut);
        std::vector<float> albedo, normal;
        if (in.hasAOVs()) {
            albedo = detail::toFloat3(in.getAlbedoData(), pixels, workingGamut);   // inside sRGB, so inside [0,1] in a wider gamut too
            normal = detail::toFloat3(in.getNormalData(), pixels);
        }
        std::vector<float> result(pixels * 3);

        const auto t0 = std::chrono::steady_clock::now();
        oidn::DeviceRef device = oidn::newDevice(oidn::DeviceType::CPU);
        device.commit();

        oidn::FilterRef filter = device.newFilter("RT");
        filter.setImage("color", color.data(), oidn::Format::Float3, static_cast<size_t>(w), static_cast<size_t>(h));
        if (!albedo.empty()) {
            filter.setImage("albedo", albedo.data(), oidn::Format::Float3, static_cast<size_t>(w), static_cast<size_t>(h));
            filter.setImage("normal", normal.data(), oidn::Format::Float3, static_cast<size_t>(w), static_cast<size_t>(h));
            filter.set("cleanAux", cleanAux);
        }
        filter.setImage("output", result.data(), oidn::Format::Float3, static_cast<size_t>(w), static_cast<size_t>(h));
        filter.set("hdr", true);   // linear radiance, not limited to [0,1]
        filter.commit();
        filter.execute();

        const char* message = nullptr;
        if (device.getError(message) != oidn::Error::None) {
            std::cerr << "[Denoise] OIDN error: " << (message ? message : "unknown") << "\n";
            return false;
        }

        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                const size_t i = (static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x)) * 3;
                out.setPixel(x, y, color::toLinearSRGB(workingGamut, Vector3(result[i], result[i + 1], result[i + 2])));
            }
        }
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        std::cout << "[Denoise] OIDN: " << w << "x" << h << (albedo.empty() ? ", color only" : ", with albedo + normal")
            << (albedo.empty() ? "" : (cleanAux ? " (clean AOVs)" : " (AOVs prefiltered)"))
            << (workingGamut == color::Gamut::SRGB ? "" : std::string(", in ") + color::gamutName(workingGamut))
            << ", " << static_cast<int>(ms) << " ms\n";
        return true;
    }

} // namespace rayt
