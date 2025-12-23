#pragma once

/**
 * @file Spectrum.hpp
 * @brief Utilities for spectral data (color) manipulation and validation.
 * * Ensures numerical stability during light transport calculations by
 * filtering out invalid physical states (NaNs, Infinities, and negative energy).
 */

#include "Core/Types.hpp"
#include "Core/Assert.hpp"

#include <cmath>
#include <glm/glm.hpp>

namespace rayt {

    /**
     * @brief Checks if the spectrum contributes no energy to the scene.
     * Values less than or equal to zero are considered "black" (non-contributing).
     * @param s The spectral radiance/importance to check.
     */
    inline bool isBlack(const Spectrum& s) {
        // Simple threshold check; could be extended with an epsilon if needed for noise floor.
        return s.x <= 0.0 && s.y <= 0.0 && s.z <= 0.0;
    }

    /**
     * @brief Checks for invalid numerical states (NaN or Infinity).
     * @return True if any component is non-finite.
     */
    inline bool HasInvalidValues(const Spectrum& s) {
        // std::isfinite returns false for NaN, Inf, and -Inf.
        return !std::isfinite(s.x) || !std::isfinite(s.y) || !std::isfinite(s.z);
    }

    /**
     * @brief Ensures a spectrum is physically valid and numerically safe.
     * Clamps negative values to zero and recovers from NaNs/Infs by returning black.
     * * @note This is critical for robust Monte Carlo integration, as a single
     * invalid sample can otherwise corrupt an entire pixel or image.
     */
    inline Spectrum Sanitize(const Spectrum& s) {
#ifndef NDEBUG  // Strict validation in Debug builds
        Assert(!HasInvalidValues(s));
#endif
        // C++20 [[unlikely]]: Hints to the compiler that numerical failure is rare,
        // optimizing the fast path for valid radiance values.
        if (HasInvalidValues(s)) [[unlikely]] return Spectrum(0.0);  

        // Clamp negative values to zero to maintain physical energy conservation.
        return glm::max(s, Spectrum(0.0));
    }

} // namespace rayt