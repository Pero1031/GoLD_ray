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
 */

#pragma once

namespace rayt {

    // Primary scalar type used for geometry, transport, and most computations.
    using Real = double;

    // Secondary scalar type used for images, texture storage, and other
    // bandwidth-sensitive data where double precision is unnecessary.
    using Float = float;

} // namespace rayt