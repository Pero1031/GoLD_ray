/**
 * @file Microfacet/GGX.cpp
 * @brief GGX (Trowbridge-Reitz) Microfacet Distribution implementation.
 *
 * The formulas live in Kernel/Microfacet.hpp, which the GPU renderer compiles too (docs/GPU.md): these methods call them in
 * double. They are the code that was here before, in the same order of operations, so the results did not change.
 */

#include "pch.h"

#include "Microfacet/GGX.hpp"
#include "Kernel/Microfacet.hpp"

namespace rayt {

    // alpha_x / alpha_y are clamped in MicrofacetDistribution's constructor (constants::ALPHA_EPSILON) to avoid degeneracy.
    GGXDistribution::GGXDistribution(Real ax, Real ay)
        : MicrofacetDistribution(ax, ay) {
    }

    // D(wh) = 1 / (pi ax ay e^2), e = wx^2 / ax^2 + wy^2 / ay^2 + wz^2; 0 below the surface.
    Real GGXDistribution::D(const Vector3& wh) const {
        return kernel::ggxD(wh, alpha_x, alpha_y);
    }

    // Smith lambda = 0.5 (sqrt(cos^2 + alpha^2 sin^2) / |cos| - 1), alpha along the azimuth of w; +inf at grazing.
    // (Dividing by |cos| only: tan^2 through math::safe_recip capped lambda near ~500 for |cos| < 1e-3, which the multiple-bounce
    // walks (Microfacet/MultiBounce.hpp) do reach; the projected area lambda |cos| then came out far too small.)
    Real GGXDistribution::lambda(const Vector3& w) const {
        return kernel::ggxLambda(w, alpha_x, alpha_y);
    }

    // VNDF sampling (Heitz 2018): stretch wo, sample the projected disk of the visible hemisphere, unstretch.
    Vector3 GGXDistribution::sample_wh(const Vector3& wo, const Point2& u) const {
        return kernel::ggxSampleVisibleNormal(wo, alpha_x, alpha_y, u.x, u.y);
    }

    // pdf of the VNDF-sampled wh (solid angle of wh): G1(wo) |wo . wh| D(wh) / wo.z.
    Real GGXDistribution::pdf(const Vector3& wo, const Vector3& wh) const {
        return kernel::ggxPdfVisibleNormal(wo, wh, alpha_x, alpha_y);
    }

} // namespace rayt
