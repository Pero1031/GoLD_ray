/**
 * @file Microfacet/GGX.cpp
 * @brief GGX (Trowbridge-Reitz) Microfacet Distribution implementation.
 * * This file contains the full implementation previously located in GGX.hpp.
 * * Behavior is intentionally kept identical to the original header-only version.
 */

#include "pch.h"

#include "Microfacet/GGX.hpp"

#include <cmath>        // std::cos, std::sin, std::abs
#include <limits>       // std::numeric_limits
#include <algorithm>    // std::max

namespace rayt {

    // -------------------------------------------------------------------------
    // Constructor
    // -------------------------------------------------------------------------
    GGXDistribution::GGXDistribution(Real ax, Real ay)
        : MicrofacetDistribution(ax, ay) {
        // NOTE:
        // alpha_x / alpha_y are clamped in MicrofacetDistribution's constructor
        // (via constants::ALPHA_EPSILON) to avoid degeneracy.
    }

    // -------------------------------------------------------------------------
    // D(wh): GGX Normal Distribution Function (NDF)
    // -------------------------------------------------------------------------
    Real GGXDistribution::D(const Vector3& wh) const {

        // Same guard as the original:
        // only upper hemisphere micro-normals contribute.
        if (wh.z <= Real(0.0)) return Real(0.0);

        // Precompute squared anisotropic roughness.
        Real ax2 = math::sqr(alpha_x);
        Real ay2 = math::sqr(alpha_y);

        // Optimized algebraic form:
        // e = (wx^2/ax^2) + (wy^2/ay^2) + wz^2
        // D = 1 / (PI * ax * ay * e^2)
        Real e =
            math::sqr(wh.x) / ax2 +
            math::sqr(wh.y) / ay2 +
            math::sqr(wh.z);

        return Real(1.0) / (constants::PI * alpha_x * alpha_y * math::sqr(e));
    }

    // -------------------------------------------------------------------------
    // lambda(w): Smith Lambda for GGX (anisotropic)
    // -------------------------------------------------------------------------
    Real GGXDistribution::lambda(const Vector3& w) const {

        // absCosTheta = |cos(theta)| in local space, where cos(theta)=w.z
        Real absCosTheta = std::abs(w.z);

        // Same handling as the original:
        // grazing/invalid -> +inf
        if (absCosTheta <= Real(0.0))
            return std::numeric_limits<Real>::infinity();

        // sin^2(theta) = 1 - cos^2(theta)
        Real sin2Theta = std::max(Real(0.0), Real(1.0) - math::sqr(absCosTheta));

        // If sin^2 is ~0, direction is near normal -> lambda=0
        if (sin2Theta <= Real(0.0))
            return Real(0.0);

        // tan^2(theta) = sin^2 / cos^2
        Real tan2Theta = sin2Theta * math::safe_recip(math::sqr(absCosTheta));

        // For anisotropic GGX, effective alpha depends on azimuth of w.
        Real wxy2 = math::sqr(w.x) + math::sqr(w.y);
        Real ax2 = math::sqr(alpha_x);
        Real ay2 = math::sqr(alpha_y);

        // Compute effective alpha^2 for the azimuthal direction of w:
        // alpha^2 = (wx^2*ax^2 + wy^2*ay^2) / (wx^2 + wy^2)
        // Fallback to ax^2 when (wx,wy) ~ 0.
        Real alpha2 = (wxy2 > Real(0.0))
            ? (math::sqr(w.x) * ax2 + math::sqr(w.y) * ay2) / wxy2
            : ax2;

        // Smith GGX lambda:
        // lambda = 0.5 * ( sqrt(1 + alpha^2 * tan^2(theta)) - 1 )
        return Real(0.5) * (math::safe_sqrt(Real(1.0) + alpha2 * tan2Theta) - Real(1.0));
    }

    // -------------------------------------------------------------------------
    // sample_wh(wo, u): VNDF sampling (Heitz 2018) for anisotropic GGX
    // -------------------------------------------------------------------------
    Vector3 GGXDistribution::sample_wh(const Vector3& wo, const Point2& u) const {

        // Same guard as the original:
        // if wo is not in upper hemisphere, return (0,0,1).
        if (wo.z <= Real(0.0)) return Vector3(0, 0, 1);

        // 1) Stretch view vector by anisotropic roughness and normalize.
        //    This maps anisotropic GGX to an isotropic configuration.
        Vector3 Vh = glm::normalize(Vector3(alpha_x * wo.x, alpha_y * wo.y, wo.z));

        // 2) Build an orthonormal basis (T1, T2, Vh) around Vh.
        //    Handle near-pole case without instability.
        Real lenSq = math::sqr(Vh.x) + math::sqr(Vh.y);
        Vector3 T1 = lenSq > 0
            ? Vector3(-Vh.y, Vh.x, 0) * math::safe_recip(math::safe_sqrt(lenSq))
            : Vector3(1, 0, 0);
        Vector3 T2 = glm::cross(Vh, T1);

        // 3) Sample the unit disk:
        //    r = sqrt(u.x), phi = 2πu.y
        Real r = math::safe_sqrt(u.x);
        Real phi = Real(2.0) * constants::PI * u.y;
        Real t1 = r * std::cos(phi);
        Real t2 = r * std::sin(phi);

        // 4) Warp t2 to fit the projection of the visible hemisphere.
        Real s = Real(0.5) * (Real(1.0) + Vh.z);
        t2 = (Real(1.0) - s) * math::safe_sqrt(Real(1.0) - math::sqr(t1)) + s * t2;

        // 5) Reconstruct Nh in stretched space.
        //    z component ensures Nh is on the hemisphere.
        Vector3 Nh = t1 * T1 + t2 * T2 +
            math::safe_sqrt(std::max(Real(0.0), Real(1.0) - math::sqr(t1) - math::sqr(t2))) * Vh;

        // 6) Unstretch and renormalize back to original space.
        //    NOTE: Keep Nh.z >= 0 (same as original).
        Vector3 ne = Vector3(
            alpha_x * Nh.x,
            alpha_y * Nh.y,
            std::max(Real(0.0), Nh.z)
        );

        return glm::normalize(ne);
    }

    // -------------------------------------------------------------------------
    // pdf(wo, wh): PDF of VNDF-sampled wh (wrt solid angle of wh)
    // -------------------------------------------------------------------------
    Real GGXDistribution::pdf(const Vector3& wo, const Vector3& wh) const {

        // Same guard as the original:
        // only upper hemisphere directions contribute.
        if (wo.z <= Real(0.0) || wh.z <= Real(0.0))
            return Real(0.0);

        // VNDF PDF:
        // pdf = G1(wo) * |wo·wh| * D(wh) / wo.z
        // (wrt solid angle of wh)
        return G1(wo) * std::abs(glm::dot(wo, wh)) * D(wh) * math::safe_recip(wo.z);
    }

} // namespace rayt