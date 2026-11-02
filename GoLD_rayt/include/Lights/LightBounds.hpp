/**
 * @file Lights/LightBounds.hpp
 * @brief Spatial and directional emission bounds of a light, and the importance of the light
 *        for a shading point (PBRT-v4 "LightBounds::Importance").
 *
 * A light is summarized by:
 *  - bounds    : a box around everything that emits,
 *  - w, cosThetaO : a cone of surface normals (axis w, half angle thetaO),
 *  - cosThetaE : emission spreads at most thetaE from each normal (pi/2 for a diffuse emitter),
 *  - phi       : emitted power (relative units are enough),
 *  - twoSided  : emits on both sides of the surface.
 *
 * importance(p, n) bounds how much the light can contribute at point p with normal n:
 *     phi * cos(theta') / d^2  *  cos(theta'_i)
 * where d is the distance to the box center (clamped), theta' the smallest possible angle
 * between the emission cone and the direction to p, and theta'_i the smallest possible angle
 * between the receiver normal and the light. It is 0 if the light cannot reach p at all.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include "Core/Types.hpp"
#include "Core/Math.hpp"
#include "Core/AABB.hpp"

namespace rayt {

    struct LightBounds {
        AABB bounds;
        Vector3 w = Vector3(0, 0, 1);   ///< axis of the cone of emitting normals
        Real phi = 0;                   ///< emitted power
        Real cosThetaO = -1;            ///< cos of the normal cone's half angle (-1: all directions)
        Real cosThetaE = 0;             ///< cos of the emission spread around a normal (0: hemisphere)
        bool twoSided = false;

        /**
         * @param p reference point
         * @param n shading normal at p, or (0,0,0) for no orientation (e.g. a medium point)
         */
        Real importance(const Point3& p, const Vector3& n) const {
            const Point3 pc = bounds.center();
            Real d2 = glm::dot(p - pc, p - pc);
            d2 = std::max(d2, glm::length(bounds.max - bounds.min) / Real(2));   // PBRT's clamp

            const Vector3 wi = (d2 > Real(0)) ? glm::normalize(p - pc) : Vector3(0, 0, 1);
            Real cosThetaW = glm::dot(w, wi);
            if (twoSided) cosThetaW = std::abs(cosThetaW);
            const Real sinThetaW = math::safe_sqrt(Real(1) - cosThetaW * cosThetaW);

            // Angle subtended by the box as seen from p.
            const Real cosThetaB = boundSubtendedCos(p);
            const Real sinThetaB = math::safe_sqrt(Real(1) - cosThetaB * cosThetaB);

            // cos(theta') = cos(max(0, thetaW - thetaO - thetaB)); bail out beyond thetaE.
            const Real sinThetaO = math::safe_sqrt(Real(1) - cosThetaO * cosThetaO);
            const Real cosThetaX = cosSubClamped(sinThetaW, cosThetaW, sinThetaO, cosThetaO);
            const Real sinThetaX = sinSubClamped(sinThetaW, cosThetaW, sinThetaO, cosThetaO);
            const Real cosThetaP = cosSubClamped(sinThetaX, cosThetaX, sinThetaB, cosThetaB);
            if (cosThetaP <= cosThetaE) return Real(0);

            Real imp = phi * cosThetaP / d2;

            if (n != Vector3(0.0)) {
                const Real cosThetaI = std::abs(glm::dot(wi, n));
                const Real sinThetaI = math::safe_sqrt(Real(1) - cosThetaI * cosThetaI);
                imp *= cosSubClamped(sinThetaI, cosThetaI, sinThetaB, cosThetaB);
            }
            return std::max(imp, Real(0));
        }

    private:
        /// cos of the half angle of the cone from p that encloses the box's bounding sphere (-1 if p is inside).
        Real boundSubtendedCos(const Point3& p) const {
            const Point3 c = bounds.center();
            const Real r2 = glm::dot(bounds.max - c, bounds.max - c);
            const Real d2 = glm::dot(p - c, p - c);
            if (d2 < r2) return Real(-1);
            return math::safe_sqrt(Real(1) - r2 / d2);
        }

        // cos(max(0, a - b)) and sin(max(0, a - b)) from sines/cosines (a - b clamped at 0).
        static Real cosSubClamped(Real sinA, Real cosA, Real sinB, Real cosB) {
            if (cosA > cosB) return Real(1);
            return cosA * cosB + sinA * sinB;
        }
        static Real sinSubClamped(Real sinA, Real cosA, Real sinB, Real cosB) {
            if (cosA > cosB) return Real(0);
            return sinA * cosB - cosA * sinB;
        }
    };

} // namespace rayt
