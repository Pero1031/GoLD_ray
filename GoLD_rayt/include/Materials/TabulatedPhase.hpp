/**
 * @file Materials/TabulatedPhase.hpp
 * @brief A phase function given as a table p(theta) (for example Mie scattering of cloud droplets).
 *
 * Mie scattering by cloud droplets (radius ~ 10 um, a few times the wavelength x 100) is far from Henyey-Greenstein: a diffraction
 * peak in the forward direction that is about 100 times higher than a Henyey-Greenstein lobe of the same mean cosine and only
 * a degree wide, a rainbow-like bump near 140 degrees and a small peak straight backward (the glory). The table keeps all of it.
 *
 * Eval and sampling are exactly consistent, which the unbiased walk needs. The nodes are equally spaced in theta, f_i = 2 pi
 * sin(theta_i) p_i is the density in theta at the node, and f is taken to be linear between nodes:
 *   - sampling: pick the bin by the cumulative distribution at the nodes, then invert the quadratic inside the bin;
 *   - eval: f(theta) / (2 pi sin(theta)) from the same piecewise linear f.
 * The table is renormalized so that f integrates to exactly 1.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"

namespace rayt {

    class TabulatedPhase {
    public:
        /// @param thetaStep  spacing of the nodes in radians
        /// @param p          density per steradian at theta = 0, thetaStep, 2 thetaStep, ... pi (at least 3 values; the scale does not matter)
        TabulatedPhase(double thetaStep, const std::vector<double>& p)
            : m_step(thetaStep), m_n(int(p.size()))
        {
            const double pi = 3.14159265358979323846;
            m_f.resize(m_n);
            for (int i = 0; i < m_n; ++i) m_f[i] = 2.0 * pi * std::sin(i * m_step) * std::max(p[i], 0.0);
            m_cdf.assign(m_n, 0.0);
            for (int i = 1; i < m_n; ++i) m_cdf[i] = m_cdf[i - 1] + 0.5 * (m_f[i - 1] + m_f[i]) * m_step;
            const double total = m_cdf[m_n - 1];
            for (int i = 0; i < m_n; ++i) { m_f[i] /= total; m_cdf[i] /= total; }
            m_cdf[m_n - 1] = 1.0;
            // The mean cosine, by integrating f(theta) cos(theta) over each linear piece with the trapezoid rule on a fine split.
            double g = 0.0;
            for (int i = 0; i + 1 < m_n; ++i) {
                const int sub = 8;
                for (int k = 0; k < sub; ++k) {
                    const double t0 = double(k) / sub, t1 = double(k + 1) / sub;
                    const double a0 = (i + t0) * m_step, a1 = (i + t1) * m_step;
                    const double f0 = m_f[i] + (m_f[i + 1] - m_f[i]) * t0, f1 = m_f[i] + (m_f[i + 1] - m_f[i]) * t1;
                    g += 0.5 * (f0 * std::cos(a0) + f1 * std::cos(a1)) * (a1 - a0);
                }
            }
            m_meanCosine = g;
        }

        /// Reads a CSV with the lines "theta_deg, value" (equal spacing from 0 to 180 degrees; lines starting with # are comments).
        /// Returns nullptr (and fills `error`) if the file cannot be used.
        static std::shared_ptr<const TabulatedPhase> fromCsv(const std::string& path, std::string* error = nullptr) {
            std::ifstream in(path);
            if (!in) { if (error) *error = "cannot open the file"; return nullptr; }
            std::vector<double> theta, value;
            std::string line;
            while (std::getline(in, line)) {
                if (line.empty() || line[0] == '#') continue;
                std::replace(line.begin(), line.end(), ',', ' ');
                std::istringstream ss(line);
                double t, v;
                if (ss >> t >> v) { theta.push_back(t); value.push_back(v); }
            }
            if (theta.size() < 3) { if (error) *error = "needs at least 3 lines of \"theta_deg, value\""; return nullptr; }
            const double step = (theta.back() - theta.front()) / double(theta.size() - 1);
            if (std::abs(theta.front()) > 1e-9 || std::abs(theta.back() - 180.0) > 1e-6 || !(step > 0)) {
                if (error) *error = "the angles must run from 0 to 180 degrees";
                return nullptr;
            }
            for (size_t i = 1; i < theta.size(); ++i)
                if (std::abs(theta[i] - theta.front() - double(i) * step) > 1e-6 * std::max(1.0, step)) {
                    if (error) *error = "the angles must be equally spaced";
                    return nullptr;
                }
            return std::make_shared<TabulatedPhase>(step * 3.14159265358979323846 / 180.0, value);
        }

        /// A table of an analytic phase function p(cosTheta) (a test: the table of Henyey-Greenstein must act like Henyey-Greenstein).
        static std::shared_ptr<const TabulatedPhase> fromFunction(const std::function<double(double)>& pOfCos, int nodes) {
            const double step = 3.14159265358979323846 / double(nodes - 1);
            std::vector<double> p(nodes);
            for (int i = 0; i < nodes; ++i) p[i] = pOfCos(std::cos(i * step));
            return std::make_shared<TabulatedPhase>(step, p);
        }

        Real meanCosine() const { return Real(m_meanCosine); }

        /// Density per steradian for the angle with cosine cosTheta between the old and the new propagation direction.
        Real eval(Real cosTheta) const {
            const double c = std::clamp(double(cosTheta), -1.0, 1.0);
            const double theta = std::acos(c);
            const double x = theta / m_step;
            const int i = std::min(int(x), m_n - 2);
            const double t = x - i;
            const double f = m_f[i] + (m_f[i + 1] - m_f[i]) * t;
            const double s = std::max(std::sqrt(std::max(1.0 - c * c, 0.0)), 1e-9);
            return Real(f / (2.0 * 3.14159265358979323846 * s));
        }

        /// cos(theta) drawn with the density eval() describes (u in [0, 1)); the azimuth is uniform and is left to the caller.
        Real sampleCosine(Real u) const {
            const double uu = std::clamp(double(u), 0.0, 1.0 - 1e-12);
            // The last node whose cumulative value is <= uu
            const int i = std::clamp(int(std::upper_bound(m_cdf.begin(), m_cdf.end(), uu) - m_cdf.begin()) - 1, 0, m_n - 2);
            const double r = uu - m_cdf[i];
            const double f0 = m_f[i], f1 = m_f[i + 1];
            const double k = (f1 - f0) / m_step;                 // slope of f in the bin: f(s) = f0 + k s, integral f0 s + k s^2 / 2 = r
            const double disc = std::max(f0 * f0 + 2.0 * k * r, 0.0);
            const double denom = f0 + std::sqrt(disc);
            double s = denom > 0.0 ? 2.0 * r / denom : 0.0;      // the stable form of the quadratic's root
            s = std::clamp(s, 0.0, m_step);
            return Real(std::cos(i * m_step + s));
        }

    private:
        double m_step;
        int m_n;
        std::vector<double> m_f;     // density in theta at the nodes (integrates to 1)
        std::vector<double> m_cdf;   // cumulative distribution at the nodes
        double m_meanCosine = 0.0;
    };

} // namespace rayt
