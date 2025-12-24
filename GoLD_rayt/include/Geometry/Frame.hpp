#pragma once

/**
 * @file Frame.hpp
 * @brief Utilities for Orthonormal Basis (ONB) and coordinate space transformations.
 * * Provides routines to construct and manage local coordinate systems (TBN frames).
 * This is essential for transforming microfacet samples or BRDF directions
 * from local tangent space to world space and vice-versa.
 */

#include <numbers>
#include <algorithm>
#include <cmath>
#include <utility>

// --- Math / Geometry ---
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

#include "Core/Constants.hpp"   // Include the constants header if available
#include "Core/Core.hpp"
#include "Core/Math.hpp"
#include "Core/Types.hpp"  // to use Real = double

// localのx= tangent(u), y= bitangent(v), z= normal(w)(right hand)

namespace rayt::frame {

    // -------------------------------------------------------------------------
    // Coordinate System (Orthonormal Basis) Construction
    // -------------------------------------------------------------------------

    /**
     * @brief Builds an orthonormal basis (T, B, N) from a given unit normal N.
     * * Implements the branchless construction method by Duff et al. (Building an
     * Orthonormal Basis, Revisited), which is a robust improvement over Frisvad's method.
     * * @param N The unit normal vector (serves as the z-axis of the local frame). |N| == 1.
     * @param T (Output) The computed tangent vector (x-axis).
     * @param B (Output) The computed bitangent vector (y-axis).
     */
    inline void makeOrthonormalBasis(const Vector3& N, Vector3& T, Vector3& B) {
        // Robust branchless method to avoid singularities
        Real sign = std::copysign(Real(1.0), N.z);
        Real a = -Real(1.0) / (sign + N.z);
        Real b = N.x * N.y * a;

        T = Vector3(Real(1.0) + sign * N.x * N.x * a, sign * b, -sign * N.x);
        B = Vector3(b, sign + N.y * N.y * a, -N.y);
    }

    /**
     * @brief Transforms a local direction vector to world space relative to normal n.
     * * @param n The reference world-space normal.
     * @param v The local-space vector (where z is the normal direction).
     * @return Vector3 The world-space equivalent of v.
     */
    inline Vector3 localToWorld(const Vector3& n, const Vector3& v) {
        Vector3 nn = glm::normalize(n);
        Vector3 T = Vector3(0.0);
        Vector3 B = Vector3(0.0);

        makeOrthonormalBasis(nn, T, B);
        return v.x * T + v.y * B + v.z * nn;
    }

    /**
     * @brief Orthonormal Basis (ONB) structure.
     * * Represents a local coordinate system defined by three orthogonal unit vectors.
     * Convention: s = Tangent, t = Bitangent, n = Normal (Right-handed system).
     */
    struct Frame {

        Vector3 s = Vector3(0.0);   // tangent local X
        Vector3 t = Vector3(0.0);   // bitangent local Y
        Vector3 n = Vector3(0.0);   // normal local Z

        Frame() = default;

        /**
         * @brief Constructs an ONB from a single normal vector.
         */
        Frame(const Vector3& nn) {
            buildFromW(nn);
        }

        /**
         * @brief Builds the basis using the provided normal as the 'w' (Z) axis.
         */
        void buildFromW(const Vector3& nn) {
            this->n = glm::normalize(nn);
            makeOrthonormalBasis(this->n, s, t);
        }

        /**
         * @brief Builds the basis aligned with a specific normal and tangent.
         * * Useful for Anisotropic materials or Normal Mapping where the
         * orientation of the UV coordinates must be preserved.
         * * Uses Gram-Schmidt orthogonalization to ensure a perfect basis.
         */
        void buildFromNormalAndTangent(const Vector3& nn, const Vector3& tangent) {
            this->n = glm::normalize(nn);

            // Projection of tangent onto the plane perpendicular to w
            Vector3 tt = tangent - this->n * glm::dot(this->n, tangent);  
            if (glm::dot(tt, tt) < 1e-12) {      // Fallback for degenerate cases
                makeOrthonormalBasis(this->n, s, this->t); 
                return;
            }

            s = glm::normalize(tt);
            // Right-handed system: t = n x s
            this->t = glm::cross(this->n, s);

            // Conditional renormalization (cheap in the common case)
            Real len2 = glm::dot(this->t, this->t);
            if (len2 > Real(0.0)) {
                // Only fix when the deviation is noticeable
                if (std::abs(len2 - Real(1.0)) > Real(1e-3)) {
                    this->t *= Real(1.0) / std::sqrt(len2);
                }
            }
        }

        /**
         * @brief Transforms a vector from local space to world space.
         */
        Vector3 localToWorld(const Vector3& a) const {
            return a.x * s + a.y * t + a.z * n;
        }

        /**
         * @brief Transforms a vector from world space to local space.
         */
        Vector3 worldToLocal(const Vector3& a) const {
            return Vector3(glm::dot(a, s), glm::dot(a, t), glm::dot(a, n));
        }

        // ---------------------------------------------------------------------
        // Trigonometric Helpers (Assumes 'v' is in LOCAL space!)
        // ---------------------------------------------------------------------

        // Cosine of the angle between v and the normal (Z-axis)
        static Real cosTheta(const Vector3& v) { return v.z; }

        static Real cosTheta2(const Vector3& v) { return v.z * v.z; }

        // Absolute Cosine (useful for two-sided materials)
        static Real absCosTheta(const Vector3& v) { return std::abs(v.z); }

        static Real sinTheta2(const Vector3& v) {
            return std::max(Real(0.0), Real(1.0) - cosTheta2(v));
        }

        static Real sinTheta(const Vector3& v) {
            return std::sqrt(sinTheta2(v));
        }

        static Real tanTheta(const Vector3& v) {
            Real temp = Real(1.0) - v.z * v.z;
            if (temp <= 0.0) return 0.0;
            return std::sqrt(temp) / v.z;
        }

        static Real tanTheta2(const Vector3& v) {
            Real temp = Real(1.0) - v.z * v.z;
            if (temp <= 0.0) return 0.0;
            return temp / (v.z * v.z);
        }

        // Azimuthal angle helpers (Phi)
        static Real sinPhi(const Vector3& v) {
            Real sinThetaSq = sinTheta2(v);
            if (sinThetaSq <= 1e-8) return 0.0; // Prevent div by zero
            return std::clamp(v.y * std::sqrt(1.0 / sinThetaSq), -1.0, 1.0);
        }

        static Real cosPhi(const Vector3& v) {
            Real sinThetaSq = sinTheta2(v);
            if (sinThetaSq <= 1e-8) return 1.0;
            return std::clamp(v.x * std::sqrt(1.0 / sinThetaSq), -1.0, 1.0);
        }

        // Returns {cosPhi, sinPhi} optimized
        static std::pair<Real, Real> sincosPhi(const Vector3& v) {
            Real sinThetaSq = sinTheta2(v);
            if (sinThetaSq <= 1e-8) return { 1.0, 0.0 };

            Real invSinTheta = std::sqrt(1.0 / sinThetaSq);
            return {
                std::clamp(v.x * invSinTheta, -1.0, 1.0), // cos
                std::clamp(v.y * invSinTheta, -1.0, 1.0)  // sin
            };
        }

        // Equality checks
        bool operator==(const Frame& f) const { return s == f.s && t == f.t && n == f.n; }
        bool operator!=(const Frame& f) const { return !(*this == f); }

    };

} // namespace rayt::frame