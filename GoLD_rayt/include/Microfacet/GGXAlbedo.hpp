/**
 * @file Microfacet/GGXAlbedo.hpp
 * @brief Directional albedo E(mu, alpha) of the single-scattering GGX BRDF (Fresnel = 1).
 *
 * E(mu_o, alpha) = Integral f_ss(wo, wi) cos(theta_i) dwi, the fraction of the light that leaves a
 * rough mirror after ONE microfacet bounce. 1 - E is the energy the single-scattering model
 * loses (the missing multiple scattering), which Kulla-Conty compensation adds back.
 *
 * Computed once at first use, with the SAME distribution and masking functions the BRDF uses
 * (GGXDistribution, height-correlated Smith). For F = 1 the VNDF-sampled estimator is
 *     E(wo) = E_{wh ~ VNDF}[ G2(wo, wi) / G1(wo) ]
 * which has very low variance; a 32 x 32 stratified set is used per table cell.
 *
 * Table axes: mu = cos(theta_o) uniform in [0, 1], and a = sqrt(alpha) uniform in [0, 1]
 * (denser at small alpha, where E changes fastest).
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Microfacet/GGX.hpp"

namespace rayt {

    class GGXAlbedoTable {
    public:
        static constexpr int kMu = 32;
        static constexpr int kAlpha = 32;

        static const GGXAlbedoTable& get() {
            static const GGXAlbedoTable table;   // thread-safe initialization
            return table;
        }

        /// E(mu, alpha), bilinear in (mu, sqrt(alpha)).
        Real E(Real mu, Real alpha) const {
            const Real x = std::clamp(mu, Real(0), Real(1)) * Real(kMu - 1);
            const Real y = std::clamp(std::sqrt(std::max(alpha, Real(0))), Real(0), Real(1)) * Real(kAlpha - 1);
            const int i = std::min(static_cast<int>(x), kMu - 2);
            const int j = std::min(static_cast<int>(y), kAlpha - 2);
            const Real fx = x - Real(i), fy = y - Real(j);
            const Real e00 = m_E[idx(j, i)], e10 = m_E[idx(j, i + 1)];
            const Real e01 = m_E[idx(j + 1, i)], e11 = m_E[idx(j + 1, i + 1)];
            return (e00 * (1 - fx) + e10 * fx) * (1 - fy) + (e01 * (1 - fx) + e11 * fx) * fy;
        }

        /// E at the kMu mu-nodes for one alpha (a row to interpolate in mu quickly) and the cosine-weighted
        /// average E_avg = 2 Integral E(mu) mu dmu.
        void row(Real alpha, std::array<Real, kMu>& out, Real& eAvg) const {
            Real sum = 0;
            for (int i = 0; i < kMu; ++i) {
                const Real mu = Real(i) / Real(kMu - 1);
                out[static_cast<size_t>(i)] = E(mu, alpha);
                const Real w = (i == 0 || i == kMu - 1) ? Real(0.5) : Real(1);   // trapezoid
                sum += w * out[static_cast<size_t>(i)] * mu;
            }
            eAvg = Real(2) * sum / Real(kMu - 1);
        }

    private:
        GGXAlbedoTable() {
            constexpr int kSide = 32;   // kSide x kSide stratified samples per cell
            for (int j = 0; j < kAlpha; ++j) {
                const Real a = Real(j) / Real(kAlpha - 1);
                const Real alpha = std::max(a * a, Real(1e-4));
                const GGXDistribution dist(alpha, alpha);
                for (int i = 0; i < kMu; ++i) {
                    const Real mu = std::max(Real(i) / Real(kMu - 1), Real(1e-3));
                    const Vector3 wo(std::sqrt(std::max(Real(0), 1 - mu * mu)), 0, mu);
                    const Real g1 = dist.G1(wo);
                    Real sum = 0;
                    for (int sy = 0; sy < kSide; ++sy) {
                        for (int sx = 0; sx < kSide; ++sx) {
                            const Point2 u((Real(sx) + Real(0.5)) / kSide, (Real(sy) + Real(0.5)) / kSide);
                            const Vector3 wh = dist.sample_wh(wo, u);
                            const Vector3 wi = wh * (Real(2) * glm::dot(wo, wh)) - wo;
                            if (wi.z > 0) sum += dist.G(wo, wi) / std::max(g1, Real(1e-12));
                        }
                    }
                    m_E[idx(j, i)] = std::min(sum / Real(kSide * kSide), Real(1));
                }
            }
        }

        static size_t idx(int alphaIdx, int muIdx) { return static_cast<size_t>(alphaIdx * kMu + muIdx); }

        std::array<Real, kMu * kAlpha> m_E{};
    };

} // namespace rayt
