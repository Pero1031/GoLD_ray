// include/Color/ComplexIorTable.hpp
#pragma once

#include <vector>
#include <utility>
#include <algorithm>

#include "Core/Precision.hpp"

namespace rayt::color {

    /**
     * @brief Tabulated complex IOR (eta, k) as a function of wavelength.
     *
     * Unit:
     * - Wavelength is stored in nanometers (nm).
     *
     * This is a data container + interpolation accessor.
     * File parsing is handled by IO loaders.
     */
    class ComplexIorTable {
    public:
        struct Sample {
            Real lambda_nm = Real(0);
            Real eta = Real(1);
            Real k = Real(0);

            bool operator<(const Sample& o) const { return lambda_nm < o.lambda_nm; }
        };

        bool empty() const { return samples_.empty(); }
        const std::vector<Sample>& samples() const { return samples_; }
        std::vector<Sample>& samples() { return samples_; }

        /// Linear interpolation. Returns (eta, k).
        std::pair<Real, Real> evaluateNK(Real lambda_nm) const {
            if (samples_.empty()) return { Real(1), Real(0) };

            if (lambda_nm <= samples_.front().lambda_nm)
                return { samples_.front().eta, samples_.front().k };
            if (lambda_nm >= samples_.back().lambda_nm)
                return { samples_.back().eta, samples_.back().k };

            Sample key{ lambda_nm, Real(0), Real(0) };
            auto it = std::lower_bound(samples_.begin(), samples_.end(), key);

            const Sample& b = *it;
            const Sample& a = *(it - 1);

            Real t = (lambda_nm - a.lambda_nm) / (b.lambda_nm - a.lambda_nm);
            Real eta = a.eta + (b.eta - a.eta) * t;
            Real k = a.k + (b.k - a.k) * t;
            return { eta, k };
        }

    private:
        std::vector<Sample> samples_;
    };

} // namespace rayt::color