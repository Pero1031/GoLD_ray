/**
 * @file Materials/SubsurfaceTable.hpp
 * @brief A tabulated radial profile for subsurface scattering (e.g. computed for a layered skin).
 *
 * The analytic profile of Materials/Subsurface.hpp (Burley's normalized diffusion) has one shape
 * for every material. A ChannelTable holds any radial profile Sp(r) as a histogram: the density of
 * exit points per unit area in log-spaced annuli, constant inside each annulus. That makes the
 * three things the probe-ray sampler needs consistent by construction:
 *   - profile(r): the area density (integrates to 1 over the plane, like sss::profile),
 *   - sampleRadius(u): the exact inverse of cdf() for this piecewise constant density,
 *   - albedo(): the total diffuse reflectance (the integral of the unnormalized profile).
 *
 * A LayeredTable bundles one ChannelTable per wavelength of the simulation and the three RGB
 * tables derived from them (see Materials/LayeredSkin.hpp).
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Core/SampledWavelengths.hpp"

namespace rayt::sss {

    class ChannelTable {
    public:
        /// @param edges    n + 1 increasing radii in meters
        /// @param density  n values: exit points per m^2 for one unit of incident power (unnormalized)
        ChannelTable(std::vector<double> edges, std::vector<double> density)
            : m_edges(std::move(edges))
        {
            const size_t n = density.size();
            m_pdf.assign(n, 0.0);
            m_cum.assign(n + 1, 0.0);
            double mass = 0;
            for (size_t k = 0; k < n; ++k) {
                const double area = area_k(k);
                mass += std::max(density[k], 0.0) * area;
                m_cum[k + 1] = mass;
            }
            m_albedo = mass;
            if (mass > 0) {
                for (size_t k = 0; k < n; ++k) m_pdf[k] = std::max(density[k], 0.0) / mass;
                for (double& c : m_cum) c /= mass;
            }
        }

        /// Total diffuse reflectance (power that leaves per unit incident power).
        Real albedo() const { return static_cast<Real>(m_albedo); }

        /// Area density of exit points (integrates to 1 over the plane). Radii below the first edge use the first annulus.
        Real profile(Real r) const {
            if (m_pdf.empty() || r >= m_edges.back()) return Real(0);
            const size_t k = (r <= m_edges.front()) ? 0 : binOf(r);
            return static_cast<Real>(m_pdf[k]);
        }

        /// Fraction of the profile inside radius r.
        Real cdf(Real r) const {
            if (m_pdf.empty() || r <= m_edges.front()) return Real(0);
            if (r >= m_edges.back()) return Real(1);
            const size_t k = binOf(r);
            const double a0 = m_edges[k] * m_edges[k], a1 = m_edges[k + 1] * m_edges[k + 1];
            const double t = (double(r) * double(r) - a0) / (a1 - a0);
            return static_cast<Real>(m_cum[k] + t * (m_cum[k + 1] - m_cum[k]));
        }

        /// Inverse of cdf(): radius for the uniform number u in [0, 1).
        Real sampleRadius(Real u) const {
            if (m_pdf.empty()) return Real(0);
            const auto it = std::upper_bound(m_cum.begin(), m_cum.end(), static_cast<double>(u));
            size_t k = (it == m_cum.begin()) ? 0 : static_cast<size_t>(it - m_cum.begin()) - 1;
            k = std::min(k, m_pdf.size() - 1);
            const double span = m_cum[k + 1] - m_cum[k];
            const double t = span > 0 ? std::clamp((double(u) - m_cum[k]) / span, 0.0, 1.0) : 0.5;
            const double a0 = m_edges[k] * m_edges[k], a1 = m_edges[k + 1] * m_edges[k + 1];
            return static_cast<Real>(std::sqrt(a0 + t * (a1 - a0)));
        }

        /// Mean exit distance (m): the radial first moment of the profile.
        Real meanRadius() const {
            double sum = 0;
            for (size_t k = 0; k < m_pdf.size(); ++k) {
                const double r0 = m_edges[k], r1 = m_edges[k + 1];
                const double mass = m_cum[k + 1] - m_cum[k];
                sum += mass * (2.0 / 3.0) * (r1 * r1 * r1 - r0 * r0 * r0) / (r1 * r1 - r0 * r0);
            }
            return static_cast<Real>(sum);
        }

    private:
        double area_k(size_t k) const { return constants::PI * (m_edges[k + 1] * m_edges[k + 1] - m_edges[k] * m_edges[k]); }

        size_t binOf(Real r) const {
            const auto it = std::upper_bound(m_edges.begin(), m_edges.end(), static_cast<double>(r));
            const size_t k = (it == m_edges.begin()) ? 0 : static_cast<size_t>(it - m_edges.begin()) - 1;
            return std::min(k, m_pdf.size() - 1);
        }

        std::vector<double> m_edges;
        std::vector<double> m_pdf;   // normalized area density per annulus
        std::vector<double> m_cum;   // cumulative mass at the edges (0 .. 1)
        double m_albedo = 0;
    };

    /// Profiles of a layered medium: one per simulated wavelength, and the three RGB tables.
    struct LayeredTable {
        std::vector<Real> wavelengths;                              ///< nm, evenly spaced
        std::vector<std::shared_ptr<const ChannelTable>> perWavelength;
        std::array<std::shared_ptr<const ChannelTable>, 3> rgb;     ///< projected under D65 with the film sensor

        /// Table for channel i of the current path: RGB mode uses the RGB tables, spectral mode the
        /// table of the nearest simulated wavelength.
        const ChannelTable* forChannel(int i, const SampledWavelengths& lambda) const {
#if defined(RAYT_SPECTRAL_MODE)
            const Real step = wavelengths.size() > 1 ? wavelengths[1] - wavelengths[0] : Real(1);
            const int idx = static_cast<int>(std::lround((lambda[i] - wavelengths.front()) / step));
            return perWavelength[static_cast<size_t>(std::clamp(idx, 0, static_cast<int>(perWavelength.size()) - 1))].get();
#else
            (void)lambda;
            return rgb[static_cast<size_t>(std::min(i, 2))].get();
#endif
        }
    };

} // namespace rayt::sss
