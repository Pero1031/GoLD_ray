/**
 * @file Core/Sampling.hpp
 * @brief Random number generation and geometric sampling routines.
 * * Provides utilities for generating samples on various manifolds (spheres, disks, hemispheres).
 * Includes importance sampling methods for reducing variance in Monte Carlo integration.
 */

#pragma once

#include <random>
#include <numbers>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>

#include <glm/gtc/constants.hpp>

#include "Core/Constants.hpp"   
#include "Core/Types.hpp"  
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
    namespace detail {
        /// This thread's generator. Until SeedStream() is called on the thread it starts from the operating system's entropy.
        inline std::mt19937& Generator() {
            static thread_local std::mt19937 generator(std::random_device{}());
            return generator;
        }
        inline std::atomic<std::uint64_t>& RunSeed() { static std::atomic<std::uint64_t> seed{ 0 }; return seed; }
        inline std::atomic<bool>& HasRunSeed() { static std::atomic<bool> has{ false }; return has; }
        /// The finalizer of splitmix64: turns nearby numbers into unrelated ones.
        inline std::uint64_t Mix(std::uint64_t z) {
            z += 0x9E3779B97F4A7C15ULL;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
            return z ^ (z >> 31);
        }
    }

    inline Real Random() {
        static thread_local std::uniform_real_distribution<Real> distribution(0.0, 1.0);
        return distribution(detail::Generator());
    }

    /**
     * @brief The seed of the whole run (--seed). Without it the random numbers differ from run to run.
     *
     * The seed alone does not fix anything: the threads share the work in a different order every time. A piece of work
     * that must repeat calls SeedStream() with its own number when it starts (the path tracer: once per image row).
     */
    inline void SetRunSeed(std::uint64_t seed) {
        detail::RunSeed().store(seed);
        detail::HasRunSeed().store(true);
    }

    /**
     * @brief Restarts this thread's random numbers for one piece of work: the stream is a function of the run seed and of
     *        (domain, a, b) only, not of the thread or of what ran before. Does nothing when no run seed is set.
     * @param domain which loop asks (1 = image rows of the path tracer, ...), so that two loops never share a stream
     */
    inline void SeedStream(std::uint64_t domain, std::uint64_t a = 0, std::uint64_t b = 0) {
        if (!detail::HasRunSeed().load(std::memory_order_relaxed)) return;
        std::uint64_t h = detail::Mix(detail::RunSeed().load(std::memory_order_relaxed) ^ detail::Mix(domain));
        h = detail::Mix(h ^ detail::Mix(a));
        h = detail::Mix(h ^ detail::Mix(b));
        std::seed_seq seq{ static_cast<std::uint32_t>(h), static_cast<std::uint32_t>(h >> 32) };
        detail::Generator().seed(seq);
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
        // Transformation based on Archimedes' hat-box theorem
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