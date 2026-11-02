/**
 * @file Color/SpectralCurve.hpp
 * @brief A tabulated function of wavelength (nm) with linear interpolation.
 *
 * Generic container for measured / generated spectral data (light source SPDs, reflectances, ...).
 * Fill it with add() (any order), call finalize() once, then evaluate().
 *
 * Out-of-range policy is chosen per query: SPDs and reflectances should be 0 outside the
 * measured range (Extrapolate::Zero), slowly varying optical constants are better clamped to the
 * end values (Extrapolate::Clamp).
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "Core/Precision.hpp"

namespace rayt::color {

    enum class Extrapolate { Zero, Clamp };

    class SpectralCurve {
    public:
        struct Sample {
            Real lambda_nm = Real(0);
            Real value = Real(0);
        };

        void add(Real lambda_nm, Real value) { m_samples.push_back({ lambda_nm, value }); }

        /// Sorts by wavelength and averages samples that share a wavelength. Call after the last add().
        void finalize() {
            std::sort(m_samples.begin(), m_samples.end(),
                [](const Sample& a, const Sample& b) { return a.lambda_nm < b.lambda_nm; });
            std::vector<Sample> merged;
            size_t i = 0;
            while (i < m_samples.size()) {
                size_t j = i;
                Real sum = Real(0);
                while (j < m_samples.size() && std::abs(m_samples[j].lambda_nm - m_samples[i].lambda_nm) < Real(1e-9)) {
                    sum += m_samples[j].value;
                    ++j;
                }
                merged.push_back({ m_samples[i].lambda_nm, sum / Real(j - i) });
                i = j;
            }
            m_samples.swap(merged);
        }

        bool empty() const { return m_samples.empty(); }
        size_t size() const { return m_samples.size(); }
        Real firstNm() const { return m_samples.empty() ? Real(0) : m_samples.front().lambda_nm; }
        Real lastNm() const { return m_samples.empty() ? Real(0) : m_samples.back().lambda_nm; }
        const std::vector<Sample>& samples() const { return m_samples; }

        /// Linear interpolation. Requires finalize().
        Real evaluate(Real lambda_nm, Extrapolate outside = Extrapolate::Zero) const {
            if (m_samples.empty()) return Real(0);
            if (lambda_nm < m_samples.front().lambda_nm)
                return outside == Extrapolate::Clamp ? m_samples.front().value : Real(0);
            if (lambda_nm > m_samples.back().lambda_nm)
                return outside == Extrapolate::Clamp ? m_samples.back().value : Real(0);

            const auto it = std::lower_bound(m_samples.begin(), m_samples.end(), lambda_nm,
                [](const Sample& s, Real l) { return s.lambda_nm < l; });
            if (it == m_samples.begin()) return it->value;
            const Sample& b = *it;
            const Sample& a = *(it - 1);
            const Real t = (lambda_nm - a.lambda_nm) / (b.lambda_nm - a.lambda_nm);
            return a.value + (b.value - a.value) * t;
        }

    private:
        std::vector<Sample> m_samples;
    };

} // namespace rayt::color
