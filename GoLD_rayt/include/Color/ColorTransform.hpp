/*
 * @file Color/ColorTransform.hpp
 * @brief Display transforms (exposure / tone mapping / OETF encoding).
 *
 * This header provides the *display pipeline* that maps linear RGB radiance
 * (typically HDR values from a path tracer) into a display-ready encoded RGB
 * (typically sRGB in [0,1]) suitable for 8-bit output.
 *
 * Design notes:
 * - This file is about *display mapping* (what you show on screen / save as PNG/JPG).
 * - Colorimetric transforms (spectral->XYZ->RGB, matrices, CMFs, etc.) live in ColorSpace.hpp.
 * - We intentionally keep these concerns separate:
 *     * ColorSpace.hpp: color science / spaces / CMFs / matrices / transfer functions
 *     * ColorTransform.hpp: exposure + tone mapping + (optional) encoding for display
 *
 * Assumptions:
 * - Spectrum is currently an RGB triplet (r,g,b) in linear space.
 * - If Spectrum becomes spectral in the future, this pipeline must be applied
 *   *after* spectral->RGB conversion.
 */

#pragma once

#include <cmath>

#include "Core/Constants.hpp"
#include "Core/Core.hpp"
#include "Core/Math.hpp"
#include "Color/ColorSpace.hpp"

namespace rayt::color {

    // -------------------------------------------------------------------------
    // Exposure
    // -------------------------------------------------------------------------

    /**
     * @brief Apply photographic exposure to a scalar value.
     *
     * Exposure is expressed in "stops" (EV steps). Each +1 stop doubles brightness:
     *   x' = x * 2^exposureStops
     *
     * Typical use:
     * - exposureStops = 0   : no change
     * - exposureStops = +1  : 2x brighter
     * - exposureStops = -1  : 2x darker
     *
     * @param x             Linear input value (e.g., radiance component)
     * @param exposureStops Exposure in stops (EV)
     * @return Exposure-adjusted value
     */
    inline Real applyExposure(Real x, Real exposureStops) {
        return x * std::pow(Real(2), exposureStops);
    }

    /**
     * @brief Apply photographic exposure to an RGB Spectrum.
     *
     * Component-wise exposure scaling using the same stop-based model:
     *   rgb' = rgb * 2^exposureStops
     *
     * @param x             Linear RGB spectrum
     * @param exposureStops Exposure in stops (EV)
     * @return Exposure-adjusted RGB spectrum
     */
    inline Spectrum applyExposure(const Spectrum& x, Real exposureStops) {
        return x * std::pow(Real(2), exposureStops);
    }

    // -------------------------------------------------------------------------
    // Tone Mapping
    // -------------------------------------------------------------------------

    /**
     * @brief Reinhard global tone mapping operator (scalar).
     *
     * Compresses HDR values into (0,1) smoothly:
     *   f(x) = x / (1 + x)
     *
     * Properties:
     * - f(0)=0, f(∞)->1
     * - Monotonic and stable
     * - Simple, but can desaturate highlights when used per-channel
     *
     * @param x Linear HDR value (>=0 recommended)
     * @return Tone-mapped value in [0,1)
     */
    inline Real toneMapReinhard(Real x) {
        return x / (Real(1) + x);
    }

    /**
     * @brief Reinhard tone mapping applied component-wise to RGB Spectrum.
     *
     * Note:
     * - This is a per-channel operator. It is simple and commonly used,
     *   but for more color-preserving behavior you may later switch to
     *   luminance-based tone mapping (operate on Y and rescale RGB).
     *
     * @param x Linear HDR RGB spectrum (>=0 recommended)
     * @return Tone-mapped RGB spectrum in [0,1)
     */
    inline Spectrum toneMapReinhard(const Spectrum& x) {
        return x / (x + Spectrum(Real(1)));
    }

    // -------------------------------------------------------------------------
    // Display pipeline
    // -------------------------------------------------------------------------

    /**
     * @brief Convert linear RGB (Spectrum) to display-ready *encoded* sRGB with exposure + Reinhard.
     *
     * Pipeline (per component):
     *  1) sanitize:
     *     - NaN/Inf -> 0
     *     - negative -> 0 (display pipeline assumes non-negative radiance)
     *  2) exposure:
     *     - multiply by 2^exposureStops
     *  3) tone map:
     *     - Reinhard compression into (0,1)
     *  4) encode:
     *     - apply sRGB OETF (linear -> sRGB encoded)
     *       (implementation lives in ColorSpace.hpp)
     *  5) clamp:
     *     - clamp to [0,1] for safe 8-bit conversion
     *
     * Important:
     * - Output is *encoded* sRGB, not linear RGB.
     *   Do NOT apply another gamma/OETF after calling this.
     *
     * Usage examples:
     * - UI preview path: convert to bytes and upload to an RGBA8 texture
     * - PNG/JPG saving: convert to bytes and write with stb_image_write
     * - HDR saving: do NOT call this; write linear values directly.
     *
     * @param linearRGB     Linear RGB input (Spectrum), typically HDR radiance or averaged radiance
     * @param exposureStops Exposure in stops (EV). Default 0.
     * @return Encoded sRGB in [0,1] per channel (ready for 8-bit quantization)
     */
    inline Spectrum toDisplaySRGB_Reinhard(const Spectrum& linearRGB, Real exposureStops = Real(0)) {
        auto sanitize = [](Real v) -> Real {
            if (rayt::math::hasNonFinite(v) || v < Real(0)) return Real(0);
            return v;
            };

        Spectrum x = linearRGB;
        x.r = sanitize(x.r);
        x.g = sanitize(x.g);
        x.b = sanitize(x.b);

        x = applyExposure(x, exposureStops);
        x = toneMapReinhard(x);

        // Encode to sRGB (canonical implementation in ColorSpace.hpp)
        x.r = linearToSRGB(x.r);
        x.g = linearToSRGB(x.g);
        x.b = linearToSRGB(x.b);

        // Clamp for 8-bit output safety
        x.r = rayt::math::saturate(x.r);
        x.g = rayt::math::saturate(x.g);
        x.b = rayt::math::saturate(x.b);

        return x;
    }

} //namespace rayt::color