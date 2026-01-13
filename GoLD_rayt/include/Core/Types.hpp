/**
 * @file Core/Types.hpp
 * @brief Semantic type aliases for the ray tracing engine.
 *
 * - Real: Defined in Precision.hpp (defaults to double-precision).
 * - Vector types: Aliased to GLM equivalents for performance and reliability.
 * - Dependency: Including this header introduces a dependency on the GLM library.
 */

#pragma once

#include <glm/glm.hpp>

#include "Core/Precision.hpp"

namespace rayt {

    // ---------------------------------------------------------------------
    // Floating-Point Geometry Types
    // ---------------------------------------------------------------------
   
    // While GLM uses the same underlying type for these, we alias them 
    // to distinguish their semantic roles in rendering (e.g., transformations).
    using Vector3 = glm::vec<3, Real, glm::defaultp>;
    using Point3 = glm::vec<3, Real, glm::defaultp>;
    using Normal3 = glm::vec<3, Real, glm::defaultp>; // Note: Must be transformed by inverse-transpose.
    using Point2 = glm::vec<2, Real, glm::defaultp>;
    using UV = glm::vec<2, Real, glm::defaultp>;

    // ---------------------------------------------------------------------
    // Bandwidth-Efficient Types (Single Precision)
    // ---------------------------------------------------------------------
    
    // Used for high-volume data where double precision is not required, 
    // such as texture data, environment maps, or random samples.
    using Vector3f = glm::vec<3, Float, glm::defaultp>;
    using Point2f = glm::vec<2, Float, glm::defaultp>;

    // ---------------------------------------------------------------------
    // Matrix Types
    // ---------------------------------------------------------------------

    using Matrix4x4 = glm::mat<4, 4, Real, glm::defaultp>;

    // ---------------------------------------------------------------------
    // Color & Spectral Representation
    // ---------------------------------------------------------------------
    
    // TODO: Transition to a dedicated Spectrum class in Spectrum.hpp 
    // once full spectral rendering is implemented.
    // Currently used as a placeholder for RGB-based light transport.
    using Spectrum = glm::vec<3, Real, glm::defaultp>;

} //namespace rayt