#pragma once

#include "Core/Types.hpp"
#include <cmath>
#include <glm/glm.hpp>

namespace rayt {

    // Checks if the spectrum contributes no energy.
    // Negative values are also considered "black" (non-contributing).
    inline bool isBlack(const Spectrum& s) {
        // (ごく小さい値を許容するためのイプシロン判定を入れても良い)
        return s.x <= 0.0 && s.y <= 0.0 && s.z <= 0.0;
    }

    // Checks for NaN (Not a Number) or Infinity.
    // std::isfinite returns true only if the value is normal, subnormal or zero.
    inline bool HasInvalidValues(const Spectrum& s) {
        return !std::isfinite(s.x) || !std::isfinite(s.y) || !std::isfinite(s.z);
    }

    // Clamps negative values to 0 and replaces NaNs/Infs with black.
    // Crucial for long-running Monte Carlo renders.
    inline Spectrum Sanitize(const Spectrum& s) {
#ifndef NDEBUG  // debug mode 
        Assert(!HasInvalidValues(s));
#endif
        if (HasInvalidValues(s)) [[unlikely]] return Spectrum(0.0);  // C++20: 異常値は「稀である」とコンパイラに伝える
        return glm::max(s, Spectrum(0.0));
    }

} // namespace rayt