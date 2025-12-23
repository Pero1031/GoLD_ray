#pragma once

/**
 * @file Assert.hpp
 * @brief Diagnostic macros for state validation and debugging.
 * * Provides a lightweight assertion mechanism to ensure internal consistency
 * and physical correctness (e.g., energy conservation, non-negative PDFs)
 * without impacting performance in production builds.
 */

#include <cstdlib>
#include <iostream>

 // ---------------------------------------------------------------------
 // Utility Macros / Functions
 // ---------------------------------------------------------------------
 // Custom assertion macro for debugging.
 //
 // - Purpose: Validate logical assumptions and physical invariants.
 // - Build Target: Active in Debug builds; stripped in Release (zero overhead).
 // - Behavior: Prints a diagnostic message (expression, file, line) to 
 //   stderr and calls std::abort() upon failure.
 //
 // WARNING:
 //   Avoid expressions with side effects (e.g., variable increments) 
 //   inside Assert(), as they will be omitted in Release builds.
 // ---------------------------------------------------------------------

#ifndef NDEBUG
/**
 * @brief Asserts that the given expression is true.
 * @param expr The boolean expression to evaluate.
 */
#define Assert(expr) \
        do { \
            if (!(expr)) { \
                std::cerr << "Assertion failed: " << #expr << " in " \
                        << __FILE__ << " line " << __LINE__ << std::endl; \
                std::abort(); \
            } \
        } while (0)
#else
/**
 * @brief NOP (No-Operation) for production builds.
 */
#define Assert(expr) do { } while (0)   // Removed in Release build
#endif