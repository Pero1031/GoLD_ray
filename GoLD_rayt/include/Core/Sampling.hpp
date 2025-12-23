#pragma once

/**
 * @file Sampling.hpp
 * @brief Random number generation and geometric sampling routines.
 * * Provides utilities for generating samples on various manifolds (spheres, disks, hemispheres).
 * Includes importance sampling methods for reducing variance in Monte Carlo integration.
 */

#include <random>
#include <numbers>
#include <algorithm>
#include <cmath>

// --- Math / Geometry ---
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

#include "Core/Constants.hpp"   // Include the constants header if available
#include "Core/Types.hpp"  // to use Real = double
#include "Core/Math.hpp"

namespace rayt::sampling {

    // -------------------------------------------------------------------------
    // Random Number Generation (RNG)
    // -------------------------------------------------------------------------

    /**
     * @brief Generates a random real number in the range [0, 1).
     * * Uses the Mersenne Twister engine with 'thread_local' storage to ensure
     * high-performance, thread-safe parallel rendering without mutex contention.
     */
    inline Real Random() {
        // static thread_local std::mt19937 generator(std::random_device{}());

        // Fixed seed for deterministic debugging
        static thread_local std::mt19937 generator(12345);
        static thread_local std::uniform_real_distribution<Real> distribution(0.0, 1.0);
        return distribution(generator);
    }

    /**
     * @brief Returns a random double in the range [min, max).
     */
    inline Real Random(Real min, Real max) {
        return min + (max - min) * Random();
    }

    /**
     * @brief Generates a 2D random sample in the unit square [0, 1)^2.
     */
    inline Point2 Random2D() {
        return Point2(Random(), Random());
    }

    // -------------------------------------------------------------------------
    // Geometric Sampling
    // -------------------------------------------------------------------------

    /**
     * @brief Returns a random point inside a unit sphere (radius 1).
     * * Uses the rejection sampling method. Suitable for simple diffuse scattering
     * and rough reflection approximations.
     */
    inline Vector3 randomInUnitSphere() {
        while (true) {
            // Generate a random vector in [-1, 1]^3 cube
            Vector3 p(Random(Real(-1), Real(1)),
                Random(Real(-1), Real(1)),
                Random(Real(-1), Real(1)));

            // If inside sphere, return it
            if (glm::dot(p, p) < Real(1)) {
                return p;
            }
        }
    }

    /**
     * @brief Maps a 2D sample from a unit square to a unit disk.
     * * Uses the square-root transform to ensure a uniform distribution across
     * the disk's area (prevents samples from clustering at the center).
     * @param u Random numbers in [0, 1)^2.
     */
    inline Vector3 UniformSampleDisk(const Point2& u) {
        Real r = std::sqrt(u.x);
        Real theta = constants::TWO_PI * u.y;

        return Vector3(r * std::cos(theta), r * std::sin(theta), 0.0);
    }

    /**
     * @brief Uniformly samples a direction on the unit sphere.
     * * Based on Archimedes' Hat-Box Theorem.
     * @param u Random numbers in [0, 1)^2.
     */
    inline Vector3 UniformSampleSphere(const Point2& u) {
        // アルキメデスのハットボックス定理に基づく変換
        Real z = 1.0 - 2.0 * u.x; // z in [-1, 1]
        Real r = math::safe_sqrt(1.0 - z * z);
        Real phi = constants::TWO_PI * u.y;

        return Vector3(r * std::cos(phi), r * std::sin(phi), z);
    }

    inline Vector3 CosineSampleHemisphere() {
        auto r1 = Random();
        auto r2 = Random();
        auto z = std::sqrt(1.0 - r2);
        auto phi = constants::TWO_PI * r1;
        auto x = std::cos(phi) * std::sqrt(r2);
        auto y = std::sin(phi) * std::sqrt(r2);
        return Vector3(x, y, z);
    }

    /**
     * @brief Cosine-weighted sampling on the unit hemisphere.
     * * This method is essential for Lambertian (diffuse) reflection, as it
     * automatically accounts for the cosine term in the Rendering Equation.
     * * @param u Random numbers in [0, 1)^2.
     * @return A direction vector in the local hemisphere (up = +Z).
     */
    inline Vector3 CosineSampleHemisphere(const Point2& u) {
        Real r1 = u.x;
        Real r2 = u.y;
        Real z = std::sqrt(std::max(Real(0.0), Real(1.0) - r2)); // cos(theta)
        Real phi = constants::TWO_PI * r1;
        Real sinTheta = std::sqrt(r2);

        Real x = std::cos(phi) * sinTheta; // sin(theta)*cos(phi)
        Real y = std::sin(phi) * sinTheta; // sin(theta)*sin(phi)
        return Vector3(x, y, z);
    }
}