#pragma once

// for type def

// Include GLM here so all files have access to vector math.
#include <glm/glm.hpp>

namespace rayt {

    // Defines a type alias for floating-point numbers.
    // Using 'Real' allows easy switching between 'float' and 'double' precision.
    using Real = double;

    // ---------------------------------------------------------------------
   // Type Definitions (Precision Control)
   // ---------------------------------------------------------------------
   // Vector/Point Aliases
   // In strict PBRT, Points and Vectors are different classes. 
   // With GLM, we alias them for clarity, even if they are the same underlying type.
    using Vector3 = glm::vec<3, Real, glm::defaultp>;
    using Point3 = glm::vec<3, Real, glm::defaultp>;
    using Normal3 = glm::vec<3, Real, glm::defaultp>; // Normalized vector
    using Point2 = glm::vec<2, Real, glm::defaultp>;
    using UV = glm::vec<2, Real, glm::defaultp>;

    // Matrix Aliases
    using Matrix4x4 = glm::mat<4, 4, Real, glm::defaultp>;

    // ---------------------------------------------------------------------
    // Spectral Representation
    // ---------------------------------------------------------------------

    // Currently, we use RGB (vec3) as a placeholder.
    // TODO: Replace this with a 'SampledSpectrum' class later for full spectral rendering.
    // e.g., class SampledSpectrum { std::array<Real, 60> values; ... };
    using Spectrum = glm::vec<3, Real, glm::defaultp>;

} //namespace rayt