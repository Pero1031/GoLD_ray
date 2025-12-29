#pragma once

#include "Core/Constants.hpp"
#include "Core/Core.hpp"

namespace rayt::color {

    /**
     * @brief Converts a linear RGB component to gamma-encoded space.
     *
     * This function applies a simple power-law gamma correction,
     * commonly approximated as gamma = 2.2.
     *
     * @note This is an approximation of sRGB and should only be used
     *       when exact sRGB transfer is not required.
     * @note Typically applied as the final step before image output.
     *
     * @param linearComponent Linear-space color component (>= 0).
     * @return Gamma-corrected color component.
     */
    inline Real linearToGamma(Real linearComponent) {
        if (linearComponent > 0) {
            return std::pow(linearComponent, 1.0 / 2.2);
        }
        return 0.0;
    }

    /**
     * @brief Applies photographic exposure adjustment.
     *
     * Exposure is applied in powers of two, following the convention:
     * x' = x * 2^exposure
     *
     * This matches the exposure model commonly used in HDR rendering
     * and physically based camera simulations.
     *
     * @param x Input linear radiance or color value.
     * @param exposure Exposure value in stops.
     * @return Exposure-adjusted value.
     */
    inline Real applyExposure(Real x, Real exposure) {
        return x * std::pow(2.0, exposure);
    }

    /**
     * @brief Reinhard global tone mapping operator.
     *
     * Compresses high dynamic range values into displayable range
     * using the classic Reinhard operator:
     * f(x) = x / (1 + x)
     *
     * @note This operator is simple, stable, and preserves relative
     *       contrast, but may desaturate highlights.
     *
     * @param x Input linear HDR value.
     * @return Tone-mapped value in [0, 1).
     */
    inline Real toneMapReinhard(Real x) {
        return x / (1.0 + x);
    }

    inline Spectrum toneMapReinhard(const Spectrum& x) {
        return x / (x + Spectrum(1.0));
    }

    /**
     * @brief Converts a linear RGB component to sRGB color space.
     *
     * Implements the standard sRGB transfer function (IEC 61966-2-1):
     * - Linear segment for low intensities
     * - Power-law segment for higher intensities
     *
     * @note This function assumes the input is in linear color space
     *       and typically follows tone mapping and exposure adjustment.
     *
     * @param x Linear-space color component.
     * @return sRGB-encoded color component.
     */
    inline Real linearToSRGB(Real x) {
        if (x <= 0.0031308)
            return 12.92 * x;
        return 1.055 * std::pow(x, 1.0 / 2.4) - 0.055;
    }

    inline Spectrum toDisplayGamma22(const Spectrum& linear) {
        Spectrum c = toneMapReinhard(linear);
        c.r = linearToGamma(c.r);
        c.g = linearToGamma(c.g);
        c.b = linearToGamma(c.b);
        return c;
    }

} //namespace rayt::color