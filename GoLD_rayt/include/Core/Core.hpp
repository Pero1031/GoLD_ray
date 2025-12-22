#pragma once

// Standard Library Includes
#include <iostream>
#include <cassert>
#include <algorithm> // for std::max

// The following header is included by GLM, but included here explicitly for clarity.
#include <cmath>
#include <cstdlib>

// Third-Party Library Includes (Linear Algebra)
// Include GLM here so all files have access to vector math.
#include <glm/glm.hpp>

// Global Constants & Configuration
//#include "Core/Constants.hpp" // custom constants (PI, etc.)
#include "Core/Types.hpp"  // Defines 'Real'

namespace rayt {

    // ---------------------------------------------------------------------
    // Type Definitions (Precision Control)
    // ---------------------------------------------------------------------

    // Switch between 'double' and 'float' easily.
    // For spectral rendering research (complex refractive indices), 'double' is highly recommended.
    // using Real = double;  // Define Real as an alias for double for easier type switching.

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
        if (HasInvalidValues(s)) {
            //std::cerr << "Invalid Spectrum detected\n";
            return Spectrum(0.0);
        }
        return glm::max(s, Spectrum(0.0));
    }

    // ---------------------------------------------------------------------
    // Forward Declarations
    // ---------------------------------------------------------------------
    // Declaring classes here allows us to use pointers/references in header files
    // without including the full definition, reducing compilation time and circular dependencies.

    // Scene & Geometry
    class Scene;
    class Shape;
    class Primitive;
    class GeometricPrimitive;
    class Aggregate; // BVH, etc.

    // Main Pipeline
    class Camera;
    class Sampler;
    class Integrator;
    class Film;      
    class Filter;    

    // Rays & Interactions
    struct Ray;
    struct RayDifferential; // Added: For texture filtering
    struct Interaction;
    struct SurfaceInteraction;

    // Lighting
    class Light;            
    class AreaLight;        
    class VisibilityTester; 

    // Shading & Materials
    class Material;
    class BSDF;
    class BxDF;
    class Texture;          // Note: If Texture is a template class later, remove this.
    class Fresnel;          
    class MicrofacetDistribution; 

    // Volumetrics
    class Medium;           
    struct MediumInterface; 
    class PhaseFunction;    

    // ---------------------------------------------------------------------
    // Utility Macros / Functions
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

} // namespace rayt