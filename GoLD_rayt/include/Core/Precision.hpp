/**
 * @file Core/Precision.hpp
 * @brief Project-wide floating-point precision configuration.
 *
 * Defines the canonical scalar types used across the renderer:
 * - Real  : primary high-precision scalar (default: double)
 * - Float : secondary scalar for bandwidth-sensitive data (default: float)
 *
 * @note This header is intentionally dependency-free to keep compilation fast
 *       and to avoid circular includes. It should be safe to include anywhere.
 * @warning Changing these aliases affects the entire ABI/layout of many types.
 *          Prefer a full rebuild after modification.
 *
 * It also pulls in Core/RenderMode.hpp, the compile-time RGB / spectral switch, so that every
 * file that sees Real also sees the same rendering mode.
 */

#pragma once

#include "Core/RenderMode.hpp"   // RAYT_ENABLE_SPECTRAL -> RAYT_SPECTRAL_MODE

namespace rayt {

    // Primary scalar type used for geometry, transport, and most computations.
    using Real = double;

    // Secondary scalar type used for images, texture storage, and other
    // bandwidth-sensitive data where double precision is unnecessary.
    using Float = float;

} // namespace rayt