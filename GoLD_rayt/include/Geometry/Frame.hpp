#pragma once

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

namespace rayt::frame {

    // -------------------------------------------------------------------------
    // Coordinate System (Orthonormal Basis)
    // -------------------------------------------------------------------------

    /**
     * @brief Builds an orthonormal basis (T, B, N) from a given normal N.
     * * A robust, branchless orthonormal basis construction (Frisvad-style).
     * Essential for transforming sampled vectors from local space to world space.
     * * @param N  The unit normal vector (z-axis of the local frame).
     * @param T  (Output) The computed tangent vector.
     * @param B  (Output) The computed bitangent vector.
     */
    inline void makeOrthonormalBasis(const Vector3& N, Vector3& T, Vector3& B) {
        Real sign = std::copysign(Real(1.0), N.z);
        Real a = -Real(1.0) / (sign + N.z);
        Real b = N.x * N.y * a;

        T = Vector3(Real(1.0) + sign * N.x * N.x * a, sign * b, -sign * N.x);
        B = Vector3(b, sign + N.y * N.y * a, -N.y);
    }

    // 法線 n を基準に、ローカルベクトル v をワールド座標へ変換
    inline Vector3 localToWorld(const Vector3& n, const Vector3& v) {
        // 簡易的なONB生成 (Frisvad's methodなど)
        Vector3 a = (std::abs(n.x) > 0.9) ? Vector3(0, 1, 0) : Vector3(1, 0, 0);
        Vector3 b = glm::normalize(glm::cross(n, a)); // 接線
        Vector3 c = glm::cross(n, b);                 // 従法線
        return v.x * b + v.y * c + v.z * n;
    }

    struct Onb {

        Vector3 u = Vector3(0.0);
        Vector3 v = Vector3(0.0);
        Vector3 w = Vector3(0.0);

        Onb() {}

        // 法線 n から基底を作るコンストラクタ
        Onb(const Vector3& n) {
            buildFromW(n);
        }

        void buildFromW(const Vector3& n) {
            w = glm::normalize(n);
            makeOrthonormalBasis(w, u, v);
        }

        // 異方性反射や法線マップで、「UVの向き」に合わせたいときに使います。
        void buildFromUW(const Vector3& n, const Vector3& tangent) {
            w = glm::normalize(n);
            // 接線(tangent)をnに対して直交化する (Gram-Schmidt)
            u = glm::normalize(tangent - w * glm::dot(w, tangent));
            v = glm::cross(w, u);
        }

        // ローカル座標 (v) を ワールド座標に変換
        Vector3 localToWorld(const Vector3& a) const {
            return a.x * u + a.y * v + a.z * w;
        }

        // ワールド座標 (a) を ローカル座標に変換
        Vector3 worldToLocal(const Vector3& a) const {
            return Vector3(glm::dot(a, u), glm::dot(a, v), glm::dot(a, w));
        }
    };

}