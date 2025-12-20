#pragma once

// -------------------------------------------------------------------------
// Standard Library Includes
// -------------------------------------------------------------------------
#include <iostream>
#include <cassert>

// The following header is included by GLM, but included here explicitly for clarity.
#include <cmath>
#include <cstdlib>

// -------------------------------------------------------------------------
// Third-Party Library Includes (Linear Algebra)
// -------------------------------------------------------------------------
// Include GLM here so all files have access to vector math.
#include <glm/glm.hpp>

// -------------------------------------------------------------------------
// 1. Global Constants & Configuration
// -------------------------------------------------------------------------
#include "Constants.hpp" // custom constants (PI, etc.)

namespace rayt {

    // ---------------------------------------------------------------------
    // 2. Type Definitions (Precision Control)
    // ---------------------------------------------------------------------

    // Switch between 'double' and 'float' easily.
    // For spectral rendering research (complex refractive indices), 'double' is highly recommended.
    using Real = double;  // Define Real as an alias for double for easier type switching.

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
    // 3. Spectral Representation
    // ---------------------------------------------------------------------

    // Currently, we use RGB (vec3) as a placeholder.
    // TODO: Replace this with a 'SampledSpectrum' class later for full spectral rendering.
    // e.g., class SampledSpectrum { std::array<Real, 60> values; ... };
    using Spectrum = glm::vec<3, Real, glm::defaultp>;


    // ---------------------------------------------------------------------
    // 4. Forward Declarations
    // ---------------------------------------------------------------------
    // Declaring classes here allows us to use pointers/references in header files
    // without including the full definition, reducing compilation time and circular dependencies.

    class Scene;
    class Camera;
    class Sampler;
    class Integrator;

    struct Ray;
    struct Interaction;
    struct SurfaceInteraction;

    class Shape;
    class Primitive;
    class GeometricPrimitive;

    class Material;
    class BSDF;
    class BxDF;
    class Texture;

    // ---------------------------------------------------------------------
    // 5. Utility Macros / Functions
    // ---------------------------------------------------------------------
    // Custom assertion macro for debugging.
    //
    // - Useful for checking physical validity (e.g., energy conservation).
    // - Enabled only in Debug builds (when NDEBUG is NOT defined).
    // - Prints the failed expression, source file, and line number.
    // - Aborts execution immediately to catch invalid states early.
    // - Completely removed in Release builds (zero runtime cost).
    //
    // NOTE:
    //   Do NOT place expressions with side effects inside Assert(),
    //   because they will not be evaluated in Release builds.
    // -----------------------------------------------------------------------------

#ifndef NDEBUG
#define Assert(expr) \
        do { \
            if (!(expr)) { \
                std::cerr << "Assertion failed: " << #expr << " in " \
                        << __FILE__ << " line " << __LINE__ << std::endl; \
                std::abort(); \
            } \
        } while (0)
#else
#define Assert(expr) do { } while (0)   // Removed in Release build
#endif

    // =====================================================================
    // Numerical Validity Checks
    // =====================================================================

    // Check for Invalid Values (NaN or Infinity)
    //
    // Returns true if the spectrum contains any non-finite values.
    // This is more robust than just checking for NaN, as it also catches 
    // infinite values which can occur due to division by zero or extremely 
    // small probability densities (PDFs).
    //
    // - NaN (Not a Number): e.g., 0.0 / 0.0
    // - Infinity: e.g., 1.0 / 0.0
    //
    // Note: Keeping infinite values can cause artifacts like "fireflies" 
    // or destroy global averages (e.g., in auto-exposure).
    inline bool HasInvalidValues(const Spectrum & c) {
        return !std::isfinite(c.x) ||
               !std::isfinite(c.y) ||
               !std::isfinite(c.z);
    }

} // namespace rayt