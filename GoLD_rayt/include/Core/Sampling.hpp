#pragma once

#include <random>
#include <numbers>
#include <algorithm>
#include <cmath>

// --- Math / Geometry ---
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

#include "Core/Constants.hpp"   // Include the constants header if available
#include "Core/Core.hpp"
#include "Core/Types.hpp"  // to use Real = double

namespace rayt::sampling {

	using Real = rayt::Real;

    // -------------------------------------------------------------------------
    // Random Number Generation (RNG) and Sampling Helper
    // -------------------------------------------------------------------------
    // Uses modern C++ <random> which provides better distribution properties.
    // 'thread_local' ensures thread safety during parallel rendering.
    // No mutex/lock is required even with OpenMP parallelization, ensuring high performance without thread contention.
    inline Real Random() {
        // Initialize generator once per thread
        // static thread_local std::mt19937 generator(std::random_device{}());

        // fix for debug 
        static thread_local std::mt19937 generator(12345);
        static thread_local std::uniform_real_distribution<Real> distribution(0.0, 1.0);
        return distribution(generator);
    }

    // Utility for Scattering
    // Returns a random float in range [min, max)
    inline Real Random(Real min, Real max) {
        return min + (max - min) * Random();
    }

    // 2D random
    inline Point2 Random2D() {
        return Point2(Random(), Random());
    }

    // Returns a random point inside a unit sphere (radius 1).
    // Used for diffuse scattering and rough reflection.
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

    // Square to Disk mapping (Polar coordinates)
    // Input: u (random numbers in [0,1])   
    // Output: Point inside unit disk (z=0)
    inline Vector3 UniformSampleDisk(const Point2& u) {
        // r = sqrt(u1) にすることで、中心に集まらず一様に分布します
        Real r = std::sqrt(u.x);
        Real theta = constants::TWO_PI * u.y;

        return Vector3(r * std::cos(theta), r * std::sin(theta), 0.0);
    }

    // Uniformly sample a direction on the unit sphere
    // Input: u (random numbers in [0,1])
    // Output: Unit vector
    inline Vector3 UniformSampleSphere(const Point2& u) {
        // アルキメデスのハットボックス定理に基づく変換
        Real z = 1.0 - 2.0 * u.x; // z in [-1, 1]
        Real r = math::safe_sqrt(1.0 - z * z);
        Real phi = constants::TWO_PI * u.y;

        return Vector3(r * std::cos(phi), r * std::sin(phi), z);
    }

    // 古い古典的な関数なのでいずれ消す-------------------------------------
    /**
     * @brief Generates a random point inside a unit disk (radius 1) on the XY plane.
     * Z-coordinate is always 0.
     * Used for camera aperture sampling (Depth of Field).
     
    inline Vector3 randomInUnitDisk() {
        while (true) {
            // Generate x, y in [-1, 1]
            Vector3 p(Random(Real(-1.0), Real(1.0)), Random(Real(-1.0), Real(1.0)), Real(0.0));
            // Check if the point is within the circle (x^2 + y^2 < 1)
            if (glm::dot(p, p) >= Real(1.0)) continue; // Retry if outside the disk
            return p;
        }
    }

    // Returns a random unit vector (on the sphere surface).
    // Used for ideal diffuse reflection (Lambertian).
    inline Vector3 randomUnitVector() {
        while (true) {
            Vector3 v = randomInUnitSphere();
            Real len2 = glm::dot(v, v);
            if (len2 > Real(1e-12)) return v / std::sqrt(len2);
        }
    }*/

    inline Vector3 CosineSampleHemisphere() {
        auto r1 = Random();
        auto r2 = Random();
        auto z = std::sqrt(1.0 - r2);
        auto phi = constants::TWO_PI * r1;
        auto x = std::cos(phi) * std::sqrt(r2);
        auto y = std::sin(phi) * std::sqrt(r2);
        return Vector3(x, y, z);
    }

    // 外部から乱数(u)を受け取ってサンプリングするバージョン（PBRに必須）
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