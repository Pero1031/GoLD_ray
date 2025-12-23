#pragma once

/**
 * @file Core.hpp
 * @brief Umbrella header for the rayt core library.
 * * This header provides a convenient single entry point for the engine's core types,
 * constants, and essential utilities.
 * * @attention
 * To maintain fast compilation and avoid circular dependencies, prefer including
 * individual headers (e.g., Types.hpp, Assert.hpp) in your internal library
 * implementations. This umbrella header is primarily intended for use in
 * top-level files like main.cpp.
 */

#include "Core/Types.hpp"           // Fundamental type aliases (Real, Point3, etc.)
#include "Core/Constants.hpp"       // Mathematical and physical constants
#include "Core/Forward.hpp"         // Forward declarations to break cycles
#include "Core/Assert.hpp"          // Debug-only validation macros
#include "Core/Math.hpp"            // Common math utilities
#include "Core/Ray.hpp"             // Ray and RayDifferential structures
#include "Core/SpectrumUtils.hpp"   // Color/Spectral safety routines

namespace rayt {

	// Global core-level configuration or engine-wide utilities can be placed here.

} // namespace rayt