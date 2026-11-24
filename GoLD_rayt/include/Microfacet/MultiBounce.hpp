/**
 * @file Microfacet/MultiBounce.hpp
 * @brief Light that bounces more than once between the microfacets of a rough conductor (docs/MULTIBOUNCE.md).
 *
 * Two estimators of the multiple-bounce reflectance rho(wi, wo) = f(wi, wo) cos(theta_o) of a Smith GGX microsurface:
 *
 * - heitzRho: the random walk of Heitz, Hanika, d'Eon and Dachsbacher 2016 ("Multiple-scattering microfacet BSDFs with the
 *   Smith model"). The walk keeps the height of each hit on the microsurface; it is unbiased and serves as the reference.
 * - zhangRho: the position-free path formulation of Zhang, Funatomi, Fujimura, Mukaigawa and Kubo 2026 ("Successive height
 *   preintegration for a height-interdependent path formulation of multiple bounces in Smith microfacet BRDFs"). The heights are
 *   integrated out in closed form (their supplement, Eq. 5): no heights are traced, the noise is lower, the bias is small.
 *
 * Both are Monte Carlo estimators: every call draws a new path on the microsurface and returns a different value whose mean is
 * rho. The single bounce is the ordinary GGX term with the height-correlated Smith G2 in both models; minOrder = 2 leaves it
 * out so that a material can add it analytically (less noise) and keep only the multiple bounces stochastic.
 *
 * All directions are in the local frame of the surface (normal = +z) and point away from the surface, except the propagation
 * directions d of a path, which point along the light: d0 = -wi goes down into the surface.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Core/Math.hpp"
#include "Microfacet/GGX.hpp"

namespace rayt::multibounce {

    /// How a rough conductor accounts for light that bounces more than once between its microfacets.
    enum class Model {
        Off,          ///< single bounce only (energy is lost on rough surfaces)
        KullaConty,   ///< a deterministic compensation lobe (Kulla and Conty 2017): energy-conserving in a white furnace, colour approximate
        Heitz,        ///< the random walk of Heitz et al. 2016 (unbiased, noisy)
        Zhang         ///< the position-free formulation of Zhang et al. 2026 (low noise, tiny bias)
    };

    /// The longest path on the microsurface counted by zhangRho (the papers' comparisons use 10).
    constexpr int kMaxBounces = 10;

    /// Projected area of the microsurface seen from w, per unit area of the macrosurface: (1 + Lambda) cos for w above the
    /// surface, Lambda |cos| for w below it (only the facets that face down-going light can be seen from below).
    Real projectedArea(const GGXDistribution& dist, const Vector3& w);

    /// Samples a microfacet normal visible from w, from either hemisphere: pdf(m) = max(0, w.m) D(m) / projectedArea(w).
    /// Spherical caps (Dupuy and Benyoub 2023) in the space stretched to alpha = 1.
    Vector3 sampleVisibleNormal(Real alphaX, Real alphaY, const Vector3& w, Real u1, Real u2);

    /// The position-free segment term of Zhang et al. 2026 (supplement, Eq. 5) for a path of n bounces with the propagation
    /// directions e[0..n]: e[0] goes down (the incoming light), e[n] goes up (the light that leaves); lam[j] = Lambda(e[j]) and
    /// up[j] = e[j].z > 0. The term depends only on Lambda of each direction and on where the path first turns upward.
    Real zhangSegment(const Real* lam, const bool* up, int n);

    /// The vertex term of a bounce from the propagation direction dIn to dOut (Zhang et al., Eq. 16):
    /// F((-dIn).h) D(h) / (4 |dIn.z|) with h = normalize(-dIn + dOut).
    template <class FresnelFn>
    Spectrum vertexTerm(const GGXDistribution& dist, const Vector3& dIn, const Vector3& dOut, FresnelFn&& fresnel) {
        const Vector3 w = -dIn;
        Vector3 h = w + dOut;
        const Real len = glm::length(h);
        if (!(len > Real(0))) return Spectrum(0.0);
        h /= len;
        const Real c = glm::dot(w, h);
        if (!(c > Real(0)) || !(h.z > Real(0)) || !(std::abs(w.z) > Real(0))) return Spectrum(0.0);
        return fresnel(c) * (dist.D(h) / (Real(4) * std::abs(w.z)));
    }

    /**
     * @brief Zhang et al. 2026: one sample of rho(wi, wo) with the position-free path formulation (without their path type-aware
     *        evaluation, PTABE, which adds a further empirical bias).
     *
     * The path's inner directions are drawn from the visible normals of the previous direction; the last bounce is connected to
     * wo. Russian roulette from the bounce `rrDepth` on, with the segment term as the survival probability (their Alg. 1). It is
     * off by default (rrDepth = kMaxBounces): the segment term falls quickly along a path, so the few paths that survive carry
     * large weights and the estimate becomes heavy-tailed (furnace tests in docs/MULTIBOUNCE.md); ten bounces cost little.
     * @param fresnel Spectrum fresnel(Real cosThetaD): the Fresnel reflectance at a microfacet, per wavelength / channel.
     * @param random  Real random(): uniform numbers in [0, 1).
     * @param minOrder The shortest path counted (1: all; 2: multiple bounces only).
     */
    template <class FresnelFn, class RandomFn>
    Spectrum zhangRho(const GGXDistribution& dist, Real alphaX, Real alphaY, const Vector3& wi, const Vector3& wo,
        FresnelFn&& fresnel, RandomFn&& random, int minOrder = 1, int maxBounces = kMaxBounces, int rrDepth = kMaxBounces) {

        Spectrum f(0.0);
        if (!(wi.z > Real(0)) || !(wo.z > Real(0))) return f;
        maxBounces = std::clamp(maxBounces, 1, kMaxBounces);
        Vector3 d[kMaxBounces + 1];
        Real lam[kMaxBounces + 2];
        bool up[kMaxBounces + 2];
        const Real lamO = dist.lambda(wo);
        d[0] = -wi;
        lam[0] = dist.lambda(wi);
        up[0] = false;
        // segment term of the path d[0..k] followed by wo (k + 1 bounces)
        auto segment = [&](int k) {
            lam[k + 1] = lamO;
            up[k + 1] = true;
            return zhangSegment(lam, up, k + 1);
        };
        if (minOrder <= 1) f += vertexTerm(dist, d[0], wo, fresnel) * segment(0);

        Spectrum weight(1.0);
        for (int k = 1; k < maxBounces; ++k) {
            const Vector3 w = -d[k - 1];
            const Real u1 = random(), u2 = random();
            const Vector3 m = sampleVisibleNormal(alphaX, alphaY, w, u1, u2);
            const Real c = glm::dot(w, m);
            if (!(c > Real(0))) break;
            const Vector3 dk = d[k - 1] + Real(2) * c * m;   // mirror reflection of d[k-1] about m
            if (!(std::abs(dk.z) > Real(1e-9))) break;
            // vertex term / pdf of the sampled direction = F (projected area / |cos|)
            weight = weight * fresnel(c) * (projectedArea(dist, w) / std::abs(w.z));
            d[k] = dk;
            lam[k] = dist.lambda(dk);
            up[k] = dk.z > Real(0);
            const Real s = segment(k);
            if (k >= rrDepth) {
                const Real q = std::min(Real(1), s);
                if (!(q > Real(0)) || random() >= q) break;
                weight = weight / q;
            }
            if (k + 1 >= minOrder) f += weight * vertexTerm(dist, dk, wo, fresnel) * s;
        }
        return f;
    }

    /**
     * @brief Heitz et al. 2016: one sample of rho(wi, wo) by a random walk on the microsurface that keeps the height of each hit.
     *
     * Heights are handled through C1(h), the fraction of the microsurface below h (any height distribution gives the same
     * result with the Smith model). At each hit the light that leaves towards wo is added (next-event estimation, with the
     * probability C1(h)^Lambda(wo) that nothing blocks it); the walk then reflects off a visible normal and goes on until it leaves.
     * The estimate has a heavy tail: light that goes nearly straight up is rarely blocked, but when it is, the facets that face it
     * have a tiny projected area and the term above is large (the variance diverges logarithmically; the mean is fine).
     * @param maxOrder The longest path followed (paths end on their own; this bounds the cost for F = 1).
     */
    template <class FresnelFn, class RandomFn>
    Spectrum heitzRho(const GGXDistribution& dist, Real alphaX, Real alphaY, const Vector3& wi, const Vector3& wo,
        FresnelFn&& fresnel, RandomFn&& random, int minOrder = 1, int maxOrder = 64) {

        Spectrum sum(0.0);
        if (!(wi.z > Real(0)) || !(wo.z > Real(0))) return sum;
        const Real lamO = dist.lambda(wo);
        Spectrum energy(1.0);
        Vector3 wr = -wi;     // propagation direction
        Real c1 = 1;          // C1 of the current height (1 = above the surface)
        for (int order = 1; order <= maxOrder; ++order) {
            // the next hit along wr, or the light leaves
            const Real u = random();
            const Real lw = dist.lambda(wr);
            if (wr.z > Real(0)) {
                if (!(lw > Real(0))) break;                        // straight up: nothing above
                if (u > Real(1) - std::pow(c1, lw)) break;         // not blocked: leaves the surface
                c1 = c1 / std::pow(Real(1) - u, Real(1) / lw);
            }
            else {
                c1 = c1 * std::pow(Real(1) - u, Real(1) / (Real(1) + lw));
            }
            if (!(c1 > Real(0))) break;
            // the light that leaves this hit towards wo
            if (order >= minOrder) {
                const Vector3 w = -wr;
                Vector3 h = w + wo;
                const Real len = glm::length(h);
                if (len > Real(0)) {
                    h /= len;
                    const Real c = glm::dot(w, h);
                    if (c > Real(0) && h.z > Real(0))
                        sum += energy * fresnel(c) * (dist.D(h) / (Real(4) * projectedArea(dist, w)) * std::pow(c1, lamO));
                }
            }
            // reflect off a normal visible from the arriving light
            const Vector3 w = -wr;
            const Real u1 = random(), u2 = random();
            const Vector3 m = sampleVisibleNormal(alphaX, alphaY, w, u1, u2);
            const Real c = glm::dot(w, m);
            if (!(c > Real(0))) break;
            energy = energy * fresnel(c);
            wr = wr + Real(2) * c * m;
        }
        return sum;
    }

} // namespace rayt::multibounce
