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

// localのx= tangent(u), y= bitangent(v), z= normal(w)(右手系)

namespace rayt::frame {

    // -------------------------------------------------------------------------
    // Coordinate System (Orthonormal Basis)
    // -------------------------------------------------------------------------

    /**
     * @brief Builds an orthonormal basis (T, B, N) from a given normal N.
     * * A robust, branchless orthonormal basis construction (Frisvad-style).
     * Essential for transforming sampled vectors from local space to world space.
     * * @param N  The unit normal vector (z-axis of the local frame).|N| == 1
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

    // Transform the local vector v to world coordinates based on the normal n.
    inline Vector3 localToWorld(const Vector3& n, const Vector3& v) {
        // Frisvad's method
        Vector3 nn = glm::normalize(n);
        Vector3 T = Vector3(0.0);
        Vector3 B = Vector3(0.0);

        makeOrthonormalBasis(nn, T, B);
        return v.x * T + v.y * B + v.z * nn;
    }

    struct Onb {

        Vector3 u = Vector3(0.0);   // tangent
        Vector3 v = Vector3(0.0);   // bitangent
        Vector3 w = Vector3(0.0);   // normal

        Onb() {}

        // 法線 n から基底を作るコンストラクタ
        Onb(const Vector3& n) {
            buildFromW(n);
        }

        void buildFromW(const Vector3& n) {
            w = glm::normalize(n);
            makeOrthonormalBasis(w, u, v);
        }

        // 異方性反射や法線マップで、「UVの向き」に合わせたいとき
        void buildFromNormalAndTangent(const Vector3& n, const Vector3& tangent) {
            w = glm::normalize(n);

            Vector3 t = tangent - w * glm::dot(w, tangent);  
            if (glm::dot(t, t) < 1e-12) {      // 退化対策
                makeOrthonormalBasis(w, u, v); // フォールバック
                return;
            }

            // 接線(tangent)をnに対して直交化する (Gram-Schmidt)
            u = glm::normalize(t);
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