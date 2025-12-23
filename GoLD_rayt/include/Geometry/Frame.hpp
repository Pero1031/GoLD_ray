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
     * Convention: u = Tangent, v = Bitangent, w = Normal (Right-handed system).
     */
    struct Frame {

        Vector3 s = Vector3(0.0);   // tangent local X
        Vector3 v = Vector3(0.0);   // bitangent local Y
        Vector3 w = Vector3(0.0);   // normal local Z

        Frame() {}

        /**
         * @brief Constructs an ONB from a single normal vector.
         */
        Frame(const Vector3& n) {
            buildFromW(n);
        }

        /**
         * @brief Builds the basis using the provided normal as the 'w' (Z) axis.
         */
        void buildFromW(const Vector3& n) {
            w = glm::normalize(n);
            makeOrthonormalBasis(w, s, v);
        }

        /**
         * @brief Builds the basis aligned with a specific normal and tangent.
         * * Useful for Anisotropic materials or Normal Mapping where the
         * orientation of the UV coordinates must be preserved.
         * * Uses Gram-Schmidt orthogonalization to ensure a perfect basis.
         */
        void buildFromNormalAndTangent(const Vector3& n, const Vector3& tangent) {
            w = glm::normalize(n);

            // Projection of tangent onto the plane perpendicular to w
            Vector3 t = tangent - w * glm::dot(w, tangent);  
            if (glm::dot(t, t) < 1e-12) {      // Fallback for degenerate cases
                makeOrthonormalBasis(w, s, v); 
                return;
            }

            s = glm::normalize(t);
            v = glm::cross(w, s);  // Right-handed system
        }

        /**
         * @brief Transforms a vector from local space to world space.
         */
        Vector3 localToWorld(const Vector3& a) const {
            return a.x * s + a.y * v + a.z * w;
        }

        /**
         * @brief Transforms a vector from world space to local space.
         */
        Vector3 worldToLocal(const Vector3& a) const {
            return Vector3(glm::dot(a, s), glm::dot(a, v), glm::dot(a, w));
        }
    };

} // namespace rayt::frame