#pragma once

#include "Core/Constants.hpp"
#include "Core/Core.hpp"

namespace renderer {

    using Real = rayt::Real;

    // Gamma Correction (Linear space -> sRGB space)
    // Standard gamma value is approximately 2.2.
    // Should be applied to the final pixel color before saving to an image.
    inline Real linearToGamma(Real linearComponent) {
        if (linearComponent > 0) {
            return std::pow(linearComponent, 1.0 / 2.2);
        }
        return 0.0;
    }
}