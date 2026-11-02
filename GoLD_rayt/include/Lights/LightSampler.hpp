/**
 * @file Lights/LightSampler.hpp
 * @brief Light selection strategies for next event estimation (PBRT-v4 style).
 *
 * A sampler picks ONE light per shading point and reports the probability p of that choice.
 * The integrator divides by p, and the MIS weights of both strategies use the combined density
 * pmf(ctx, light) * light.pdfLi(), so every sampler must return the same pmf from sample() and
 * pmf() for the same context.
 *
 * Samplers (all include the environment light in the list):
 *  - UniformLightSampler : every light 1/N.
 *  - PowerLightSampler   : bounded lights proportional to emitted power (alias table, O(1)),
 *                          independent of the shading point (PBRT "power").
 *  - SpatialLightSampler : bounded lights proportional to LightBounds::importance(p, n), i.e.
 *                          power x distance falloff x orientation cones of the light and of the
 *                          receiver (PBRT's BVH importance, evaluated for every light = PBRT's
 *                          "exhaustive" sampler; O(N) per query, fine up to a few dozen lights).
 *
 * Infinite / unbounded lights (environment maps) have no meaningful power, so, as in PBRT's BVH
 * sampler, they share a fixed probability: each unbounded light gets 1/(nUnbounded + 1) and
 * the bounded lights together get the remaining 1/(nUnbounded + 1).
 */

#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

#include "Core/Types.hpp"
#include "Lights/Light.hpp"

namespace rayt {

    /// The chosen light and the probability with which it was chosen.
    struct SampledLight {
        const Light* light = nullptr;
        Real p = 0;
    };

    class LightSampler {
    public:
        virtual ~LightSampler() = default;

        /// Picks a light for the shading point `ctx` using the uniform random number u in [0,1).
        virtual std::optional<SampledLight> sample(const LightSampleContext& ctx, Real u) const = 0;

        /// Probability that sample(ctx, .) returns `light` (0 for lights it does not know).
        virtual Real pmf(const LightSampleContext& ctx, const Light* light) const = 0;

        virtual size_t numLights() const = 0;
        bool hasLights() const { return numLights() > 0; }

        /// True for a sampler that does not choose: the integrator samples every light i < numLights()
        /// returned by lightAt(i) with pmf(ctx, light) > 0 (probability 1 each) at every shading point.
        virtual bool samplesAll() const { return false; }
        virtual const Light* lightAt(size_t /*i*/) const { return nullptr; }
    };

    // ---------------------------------------------------------------------------------------
    // Uniform
    // ---------------------------------------------------------------------------------------
    class UniformLightSampler : public LightSampler {
    public:
        explicit UniformLightSampler(const std::vector<std::shared_ptr<Light>>& lights) {
            for (const auto& l : lights) m_lights.push_back(l.get());
        }

        std::optional<SampledLight> sample(const LightSampleContext&, Real u) const override {
            if (m_lights.empty()) return std::nullopt;
            const int n = static_cast<int>(m_lights.size());
            const int i = std::min(static_cast<int>(u * Real(n)), n - 1);
            return SampledLight{ m_lights[static_cast<size_t>(i)], Real(1) / Real(n) };
        }

        Real pmf(const LightSampleContext&, const Light* light) const override {
            if (std::find(m_lights.begin(), m_lights.end(), light) == m_lights.end()) return Real(0);
            return Real(1) / Real(m_lights.size());
        }

        size_t numLights() const override { return m_lights.size(); }

    private:
        std::vector<const Light*> m_lights;
    };

    // ---------------------------------------------------------------------------------------
    // Shared bookkeeping: bounded lights (with LightBounds) and unbounded ones
    // ---------------------------------------------------------------------------------------
    class BoundedLightSet {
    public:
        explicit BoundedLightSet(const std::vector<std::shared_ptr<Light>>& lights) {
            for (const auto& sp : lights) {
                const Light* l = sp.get();
                if (const auto b = l->bounds()) {
                    m_index[l] = static_cast<int>(m_bounded.size());
                    m_bounded.push_back(l);
                    m_bounds.push_back(*b);
                }
                else {
                    m_index[l] = -1 - static_cast<int>(m_unbounded.size());
                    m_unbounded.push_back(l);
                }
            }
        }

        size_t numBounded() const { return m_bounded.size(); }
        size_t numUnbounded() const { return m_unbounded.size(); }
        const Light* bounded(size_t i) const { return m_bounded[i]; }
        const Light* unbounded(size_t i) const { return m_unbounded[i]; }
        const LightBounds& bounds(size_t i) const { return m_bounds[i]; }

        /// >= 0: index among bounded lights; < 0: -1 - index among unbounded; INT_MIN: unknown.
        int find(const Light* l) const {
            const auto it = m_index.find(l);
            return it == m_index.end() ? std::numeric_limits<int>::min() : it->second;
        }

        /// Probability of the unbounded group; `boundedUsable` = whether bounded lights can be picked at all.
        Real pUnbounded(bool boundedUsable) const {
            const Real nu = Real(m_unbounded.size());
            if (m_unbounded.empty()) return Real(0);
            return boundedUsable ? nu / (nu + Real(1)) : Real(1);
        }

        /**
         * Probability of the unbounded group when the bounded lights have a combined importance
         * `boundedTotal` at the shading point: the unbounded lights' irradiance against the bounded
         * lights' importance (both are about "L x solid angle"). A dim sky then no longer takes half
         * of all samples. Each group keeps at least kMinGroupShare so a wrong estimate cannot starve it;
         * lights that cannot estimate their strength give the fixed share of pUnbounded().
         */
        Real pUnboundedByImportance(Real boundedTotal) const {
            if (m_unbounded.empty()) return Real(0);
            if (boundedTotal <= Real(0)) return Real(1);
            Real eu = 0;
            for (const Light* l : m_unbounded) eu += l->irradianceEstimate();
            if (eu <= Real(0)) return pUnbounded(true);
            return std::clamp(eu / (eu + boundedTotal), kMinGroupShare, Real(1) - kMinGroupShare);
        }

        static constexpr Real kMinGroupShare = Real(0.05);

    private:
        std::vector<const Light*> m_bounded, m_unbounded;
        std::vector<LightBounds> m_bounds;
        std::unordered_map<const Light*, int> m_index;
    };

    // ---------------------------------------------------------------------------------------
    // Alias table (Vose): O(1) sampling from a discrete distribution
    // ---------------------------------------------------------------------------------------
    class AliasTable {
    public:
        AliasTable() = default;

        explicit AliasTable(const std::vector<Real>& weights) {
            const size_t n = weights.size();
            m_pmf.assign(n, Real(0));
            m_q.assign(n, Real(1));
            m_alias.assign(n, 0);
            if (n == 0) return;

            Real sum = 0;
            for (Real w : weights) sum += std::max(w, Real(0));
            for (size_t i = 0; i < n; ++i)
                m_pmf[i] = sum > Real(0) ? std::max(weights[i], Real(0)) / sum : Real(1) / Real(n);

            std::vector<Real> scaled(n);
            std::vector<int> small, large;
            for (size_t i = 0; i < n; ++i) {
                scaled[i] = m_pmf[i] * Real(n);
                (scaled[i] < Real(1) ? small : large).push_back(static_cast<int>(i));
            }
            while (!small.empty() && !large.empty()) {
                const int s = small.back(); small.pop_back();
                const int l = large.back(); large.pop_back();
                m_q[static_cast<size_t>(s)] = scaled[static_cast<size_t>(s)];
                m_alias[static_cast<size_t>(s)] = l;
                scaled[static_cast<size_t>(l)] += scaled[static_cast<size_t>(s)] - Real(1);
                (scaled[static_cast<size_t>(l)] < Real(1) ? small : large).push_back(l);
            }
            // Leftovers have probability 1 (numerical round-off).
            for (int i : small) m_q[static_cast<size_t>(i)] = Real(1);
            for (int i : large) m_q[static_cast<size_t>(i)] = Real(1);
        }

        size_t size() const { return m_pmf.size(); }
        Real pmf(size_t i) const { return m_pmf[i]; }

        int sample(Real u) const {
            const int n = static_cast<int>(m_pmf.size());
            const Real x = u * Real(n);
            const int i = std::min(static_cast<int>(x), n - 1);
            const Real frac = x - Real(i);
            return frac < m_q[static_cast<size_t>(i)] ? i : m_alias[static_cast<size_t>(i)];
        }

    private:
        std::vector<Real> m_pmf, m_q;
        std::vector<int> m_alias;
    };

    // ---------------------------------------------------------------------------------------
    // Power
    // ---------------------------------------------------------------------------------------
    class PowerLightSampler : public LightSampler {
    public:
        explicit PowerLightSampler(const std::vector<std::shared_ptr<Light>>& lights)
            : m_set(lights)
        {
            std::vector<Real> phi;
            for (size_t i = 0; i < m_set.numBounded(); ++i) phi.push_back(m_set.bounds(i).phi);
            m_table = AliasTable(phi);   // all-zero power falls back to uniform inside AliasTable
            m_pUnbounded = m_set.pUnbounded(m_set.numBounded() > 0);
        }

        std::optional<SampledLight> sample(const LightSampleContext&, Real u) const override {
            if (m_set.numUnbounded() > 0 && u < m_pUnbounded) {
                const size_t n = m_set.numUnbounded();
                const size_t i = std::min(static_cast<size_t>(u / m_pUnbounded * Real(n)), n - 1);
                return SampledLight{ m_set.unbounded(i), m_pUnbounded / Real(n) };
            }
            if (m_set.numBounded() == 0) return std::nullopt;
            const Real u2 = (u - m_pUnbounded) / (Real(1) - m_pUnbounded);
            const int i = m_table.sample(std::clamp(u2, Real(0), Real(0.99999999)));
            return SampledLight{ m_set.bounded(static_cast<size_t>(i)),
                (Real(1) - m_pUnbounded) * m_table.pmf(static_cast<size_t>(i)) };
        }

        Real pmf(const LightSampleContext&, const Light* light) const override {
            const int idx = m_set.find(light);
            if (idx >= 0) return (Real(1) - m_pUnbounded) * m_table.pmf(static_cast<size_t>(idx));
            if (idx == std::numeric_limits<int>::min()) return Real(0);
            return m_pUnbounded / Real(m_set.numUnbounded());
        }

        size_t numLights() const override { return m_set.numBounded() + m_set.numUnbounded(); }

    private:
        BoundedLightSet m_set;
        AliasTable m_table;
        Real m_pUnbounded = 0;
    };

    // ---------------------------------------------------------------------------------------
    // Spatial (importance of every light at the shading point)
    // ---------------------------------------------------------------------------------------
    class SpatialLightSampler : public LightSampler {
    public:
        explicit SpatialLightSampler(const std::vector<std::shared_ptr<Light>>& lights)
            : m_set(lights) {}

        std::optional<SampledLight> sample(const LightSampleContext& ctx, Real u) const override {
            const size_t nb = m_set.numBounded(), nu = m_set.numUnbounded();
            std::array<Real, kStack> stackBuf;
            std::vector<Real> heapBuf;
            Real* imp = stackBuf.data();
            if (nb > kStack) { heapBuf.resize(nb); imp = heapBuf.data(); }
            const Real total = importances(ctx, imp);

            const Real pU = m_set.pUnboundedByImportance(total);
            if (nu > 0 && u < pU) {
                const size_t i = std::min(static_cast<size_t>(u / pU * Real(nu)), nu - 1);
                return SampledLight{ m_set.unbounded(i), pU / Real(nu) };
            }
            if (total <= Real(0)) return std::nullopt;   // no bounded light can reach this point

            const Real target = std::clamp((u - pU) / (Real(1) - pU), Real(0), Real(0.99999999)) * total;
            Real acc = 0;
            size_t chosen = nb;
            for (size_t i = 0; i < nb; ++i) {
                if (imp[i] <= Real(0)) continue;
                acc += imp[i];
                chosen = i;                       // remember the last usable light (round-off safe)
                if (target < acc) break;
            }
            if (chosen == nb) return std::nullopt;
            return SampledLight{ m_set.bounded(chosen), (Real(1) - pU) * imp[chosen] / total };
        }

        Real pmf(const LightSampleContext& ctx, const Light* light) const override {
            const int idx = m_set.find(light);
            if (idx == std::numeric_limits<int>::min()) return Real(0);

            const size_t nb = m_set.numBounded();
            std::array<Real, kStack> stackBuf;
            std::vector<Real> heapBuf;
            Real* imp = stackBuf.data();
            if (nb > kStack) { heapBuf.resize(nb); imp = heapBuf.data(); }
            const Real total = importances(ctx, imp);
            const Real pU = m_set.pUnboundedByImportance(total);

            if (idx < 0) return pU / Real(m_set.numUnbounded());
            if (total <= Real(0)) return Real(0);
            return (Real(1) - pU) * imp[idx] / total;
        }

        size_t numLights() const override { return m_set.numBounded() + m_set.numUnbounded(); }

    private:
        static constexpr size_t kStack = 64;

        Real importances(const LightSampleContext& ctx, Real* out) const {
            Real total = 0;
            for (size_t i = 0; i < m_set.numBounded(); ++i) {
                out[i] = m_set.bounds(i).importance(ctx.p, ctx.gn);
                total += out[i];
            }
            return total;
        }

        BoundedLightSet m_set;
    };

    // ---------------------------------------------------------------------------------------
    // Exhaustive (every light that can reach the point, probability 1 each)
    // ---------------------------------------------------------------------------------------
    /**
     * With few lights, one shadow ray per light is often cheaper than the noise of picking one at
     * random: the variance of "choose light i with probability p_i" is V_i / p_i, with p_i = 1 it is V_i.
     * A bounded light whose LightBounds::importance is 0 cannot reach the point and is skipped
     * (pmf 0); everything else, the environment light included, has pmf 1. O(N) shadow rays per
     * shading point, so only for a handful of lights.
     */
    class ExhaustiveLightSampler : public LightSampler {
    public:
        explicit ExhaustiveLightSampler(const std::vector<std::shared_ptr<Light>>& lights)
            : m_set(lights)
        {
            for (size_t i = 0; i < m_set.numBounded(); ++i) m_all.push_back(m_set.bounded(i));
            for (size_t i = 0; i < m_set.numUnbounded(); ++i) m_all.push_back(m_set.unbounded(i));
        }

        // Not used (samplesAll() is true); returns the first reachable light for completeness.
        std::optional<SampledLight> sample(const LightSampleContext& ctx, Real) const override {
            for (const Light* l : m_all) if (pmf(ctx, l) > Real(0)) return SampledLight{ l, Real(1) };
            return std::nullopt;
        }

        Real pmf(const LightSampleContext& ctx, const Light* light) const override {
            const int idx = m_set.find(light);
            if (idx == std::numeric_limits<int>::min()) return Real(0);
            if (idx < 0) return Real(1);
            return m_set.bounds(static_cast<size_t>(idx)).importance(ctx.p, ctx.gn) > Real(0) ? Real(1) : Real(0);
        }

        size_t numLights() const override { return m_all.size(); }
        bool samplesAll() const override { return true; }
        const Light* lightAt(size_t i) const override { return m_all[i]; }

    private:
        BoundedLightSet m_set;
        std::vector<const Light*> m_all;
    };

} // namespace rayt
