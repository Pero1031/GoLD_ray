#pragma once

#include <cstdlib>
#include <iostream>

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
// ---------------------------------------------------------------------

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