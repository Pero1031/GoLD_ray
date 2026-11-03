/**
 * @file Materials/LayeredSkin.hpp
 * @brief Radial reflectance profile of a layered skin, by Monte Carlo, and the tables built from it.
 *
 * Skin is modeled as a planar stack (Jacques-style chromophore model):
 *   - epidermis (thickness ~0.1 mm): absorption = melanin (fraction f_mel) + a baseline;
 *   - dermis (semi-infinite):        absorption = blood (volume fraction B, oxygen saturation S) + a baseline;
 *   - reduced scattering mu_s'(l) = a (l / 500 nm)^-b in both layers.
 * A pencil beam enters at the origin; photons random-walk (isotropic scattering with mu_s', the
 * similarity relation for anisotropic scattering) and every photon that leaves the top surface adds
 * its weight to a log-spaced annulus around the entry point. The top boundary and the layer
 * boundaries are index matched: the Fresnel entry / exit factors are not part of the profile, the
 * renderer applies them (CoatedDiffuse and SubsurfaceExit).
 *
 * The result Rd(r, lambda) is stored per wavelength and projected to RGB (D65, film sensor), which is
 * the profile that the RGB build uses; the spectral build uses the nearest wavelength.
 *
 * The numbers of the chromophore model are typical literature values written from memory (see
 * SkinModel); the spectra of hemoglobin are Prahl's table at 20 nm steps (omlc.org).
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "Core/Constants.hpp"
#include "Color/ColorSpace.hpp"
#include "Materials/SubsurfaceTable.hpp"

namespace rayt::skin {

    // ------------------------------------------------------------------------------------------
    // Optical properties
    // ------------------------------------------------------------------------------------------

    /// One homogeneous layer at one wavelength. The last layer of a stack is semi-infinite.
    struct Layer {
        double thickness = 0;   ///< m
        double mua = 0;         ///< absorption coefficient, 1/m
        double musp = 0;        ///< reduced scattering coefficient, 1/m
    };

    struct SkinModel {
        double melanin = 0.04;               ///< melanosome volume fraction of the epidermis (light ~0.02, dark ~0.25)
        double blood = 0.02;                 ///< blood volume fraction of the dermis
        double oxygen = 0.75;                ///< oxygen saturation of the blood
        double epidermis = 0.1e-3;           ///< epidermis thickness (m)
        double musp500Epidermis = 5000;      ///< mu_s' at 500 nm (1/m)
        double musp500Dermis = 4500;
        double scatterExponent = 1.3;        ///< b in a (l / 500 nm)^-b
        bool singleLayer = false;            ///< dermis only (a homogeneous semi-infinite medium)
    };

    namespace detail {
        inline constexpr double kHbFirstNm = 400, kHbStepNm = 20;
        // Molar extinction coefficient (cm^-1 / M) of HbO2 and Hb, 400 .. 740 nm every 20 nm (Prahl, omlc.org).
        inline constexpr double kHbO2[18] = { 266232, 480360, 102580, 44480, 26629.2, 20932.8, 24202.4, 53236, 32613.2,
            50104, 3200, 942, 442, 319.6, 277.6, 290, 348, 446 };
        inline constexpr double kHb[18] = { 223296, 407560, 413280, 23388.8, 14550, 20862, 31589.6, 46592, 53788,
            37020, 14677.2, 6509.6, 4345.2, 3226.56, 2407.92, 1794.28, 1325.88, 1115.88 };

        inline double lerpTable(const double* t, double nm) {
            const double x = std::clamp((nm - kHbFirstNm) / kHbStepNm, 0.0, 17.0);
            const int i = std::min(static_cast<int>(x), 16);
            const double f = x - i;
            return t[i] * (1 - f) + t[i + 1] * f;
        }
    }

    /// Baseline absorption of bloodless tissue (1/m), after Jacques: 0.244 + 85.3 exp(-(l - 154) / 66.2) cm^-1.
    inline double muaBaseline(double nm) { return (0.244 + 85.3 * std::exp(-(nm - 154.0) / 66.2)) * 100.0; }

    /// Absorption of pure melanosomes (1/m): 6.6e11 l^-3.33 cm^-1 (l in nm).
    inline double muaMelanin(double nm) { return 6.6e11 * std::pow(nm, -3.33) * 100.0; }

    /// Absorption of whole blood (1/m) at 150 g/L hemoglobin: mu_a = 2.303 * eps * 150 / 64500 cm^-1.
    inline double muaBlood(double nm, double oxygen) {
        const double k = 2.303 * 150.0 / 64500.0 * 100.0;
        return k * (oxygen * detail::lerpTable(detail::kHbO2, nm) + (1 - oxygen) * detail::lerpTable(detail::kHb, nm));
    }

    inline std::vector<Layer> layersAt(const SkinModel& m, double nm) {
        const double scatter = std::pow(nm / 500.0, -m.scatterExponent);
        Layer dermis;
        dermis.thickness = 0;
        dermis.mua = m.blood * muaBlood(nm, m.oxygen) + (1 - m.blood) * muaBaseline(nm);
        dermis.musp = m.musp500Dermis * scatter;
        if (m.singleLayer) return { dermis };
        Layer epi;
        epi.thickness = m.epidermis;
        epi.mua = m.melanin * muaMelanin(nm) + (1 - m.melanin) * muaBaseline(nm);
        epi.musp = m.musp500Epidermis * scatter;
        return { epi, dermis };
    }

    // ------------------------------------------------------------------------------------------
    // Monte Carlo
    // ------------------------------------------------------------------------------------------

    struct Rng {
        uint64_t s;
        explicit Rng(uint64_t seed) : s(seed) {}
        double next() {   // splitmix64 -> [0, 1)
            uint64_t z = (s += 0x9E3779B97F4A7C15ull);
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
            z ^= z >> 31;
            return (z >> 11) * (1.0 / 9007199254740992.0);
        }
    };

    struct RadialResult {
        std::vector<double> edges;         ///< nBins + 1 radii (m)
        std::vector<double> density;       ///< exit weight per m^2 per unit incident power
        double reflected = 0;              ///< total weight that left the top surface
        double lostOutside = 0;            ///< left the top surface beyond the last edge
        double absorbed = 0;               ///< absorbed (including roulette losses)
    };

    inline RadialResult simulateRadialProfile(const std::vector<Layer>& layers, int nBins, double rMin, double rMax,
        long long nPhotons, uint64_t seed)
    {
        const int nLayers = static_cast<int>(layers.size());
        std::vector<double> top(nLayers), bottom(nLayers);
        double z0 = 0;
        for (int i = 0; i < nLayers; ++i) {
            top[i] = z0;
            z0 += layers[i].thickness;
            bottom[i] = (i == nLayers - 1) ? 1e300 : z0;
        }
        std::vector<double> mut(nLayers), albedo(nLayers);
        for (int i = 0; i < nLayers; ++i) {
            mut[i] = layers[i].mua + layers[i].musp;
            albedo[i] = layers[i].musp / mut[i];
        }

        RadialResult res;
        res.edges.resize(nBins + 1);
        const double ratio = std::pow(rMax / rMin, 1.0 / nBins);
        for (int k = 0; k <= nBins; ++k) res.edges[k] = rMin * std::pow(ratio, k);
        std::vector<double> hist(nBins, 0.0);
        double lost = 0, absorbed = 0;
        const double logRatio = std::log(ratio);

#if defined(_OPENMP)
#pragma omp parallel
#endif
        {
            std::vector<double> h(nBins, 0.0);
            double lostLocal = 0, absLocal = 0;
#if defined(_OPENMP)
            const int tid = omp_get_thread_num(), nThreads = omp_get_num_threads();
#else
            const int tid = 0, nThreads = 1;
#endif
            Rng rng(seed + 0x1234567ull * (tid + 1));
            const long long first = nPhotons * tid / nThreads, last = nPhotons * (tid + 1) / nThreads;
            for (long long p = first; p < last; ++p) {
                double x = 0, y = 0, z = 0, ux = 0, uy = 0, uz = 1, w = 1;
                int layer = 0;
                for (long long step = 0; step < 100000000ll; ++step) {
                    const double s = -std::log(1.0 - rng.next()) / mut[layer];
                    double dist = 1e300;
                    if (uz > 0 && bottom[layer] < 1e299) dist = (bottom[layer] - z) / uz;
                    else if (uz < 0) dist = (top[layer] - z) / uz;
                    if (s < dist) {
                        x += ux * s; y += uy * s; z += uz * s;
                        absLocal += w * (1 - albedo[layer]);
                        w *= albedo[layer];
                        const double c = 2 * rng.next() - 1, sn = std::sqrt(std::max(0.0, 1 - c * c));
                        const double phi = 2 * constants::PI * rng.next();
                        ux = sn * std::cos(phi); uy = sn * std::sin(phi); uz = c;
                        if (w < 1e-4) {
                            if (rng.next() < 0.1) w /= 0.1;
                            else { absLocal += w; break; }
                        }
                    }
                    else {
                        x += ux * dist; y += uy * dist; z += uz * dist;
                        if (uz < 0) {
                            if (layer == 0) {   // leaves through the top surface
                                const double r = std::sqrt(x * x + y * y);
                                if (r >= rMax) lostLocal += w;
                                else {
                                    int k = r <= rMin ? 0 : static_cast<int>(std::log(r / rMin) / logRatio);
                                    h[std::min(k, nBins - 1)] += w;
                                }
                                break;
                            }
                            --layer;
                        }
                        else ++layer;
                    }
                }
            }
#if defined(_OPENMP)
#pragma omp critical
#endif
            {
                for (int k = 0; k < nBins; ++k) hist[k] += h[k];
                lost += lostLocal;
                absorbed += absLocal;
            }
        }

        res.density.resize(nBins);
        for (int k = 0; k < nBins; ++k) {
            const double area = constants::PI * (res.edges[k + 1] * res.edges[k + 1] - res.edges[k] * res.edges[k]);
            res.density[k] = hist[k] / double(nPhotons) / area;
            res.reflected += hist[k] / double(nPhotons);
        }
        res.lostOutside = lost / double(nPhotons);
        res.absorbed = absorbed / double(nPhotons);
        return res;
    }

    // ------------------------------------------------------------------------------------------
    // Tables
    // ------------------------------------------------------------------------------------------

    /// Builds the per-wavelength tables and the RGB tables (D65 projection) from the simulated densities.
    inline std::shared_ptr<sss::LayeredTable> buildLayeredTable(const std::vector<double>& wavelengths,
        const std::vector<double>& edges, const std::vector<std::vector<double>>& density /* [wavelength][bin] */)
    {
        auto table = std::make_shared<sss::LayeredTable>();
        for (size_t i = 0; i < wavelengths.size(); ++i) {
            table->wavelengths.push_back(static_cast<Real>(wavelengths[i]));
            table->perWavelength.push_back(std::make_shared<sss::ChannelTable>(edges, density[i]));
        }
        const size_t nBins = edges.size() - 1;
        std::array<std::vector<double>, 3> rgbDensity;
        for (auto& v : rgbDensity) v.assign(nBins, 0.0);
        for (size_t k = 0; k < nBins; ++k) {
            const auto fn = [&](Real nm) {
                const double x = std::clamp((double(nm) - wavelengths.front()) / (wavelengths[1] - wavelengths[0]), 0.0, double(wavelengths.size() - 1));
                const size_t i = std::min(static_cast<size_t>(x), wavelengths.size() - 2);
                const double f = x - double(i);
                return density[i][k] * (1 - f) + density[i + 1][k] * f;
            };
            const Vector3 c = color::spectrumToRGB_D65(fn);
            for (int ch = 0; ch < 3; ++ch) rgbDensity[ch][k] = std::max(0.0, double(c[ch]));
        }
        for (int ch = 0; ch < 3; ++ch) table->rgb[ch] = std::make_shared<sss::ChannelTable>(edges, rgbDensity[ch]);
        return table;
    }

    inline bool saveTableCsv(const std::string& path, const std::vector<double>& wavelengths,
        const std::vector<double>& edges, const std::vector<std::vector<double>>& density, const std::string& note)
    {
        std::ofstream f(path);
        if (!f) return false;
        f << "# layered skin radial profile: exit density (1/m^2 per unit incident power), index matched boundaries\n";
        f << "# " << note << "\n";
        f << "r_lo_m,r_hi_m";
        for (double nm : wavelengths) f << ",D_" << nm;
        f << "\n";
        f.precision(9);
        for (size_t k = 0; k + 1 < edges.size(); ++k) {
            f << edges[k] << "," << edges[k + 1];
            for (size_t i = 0; i < wavelengths.size(); ++i) f << "," << density[i][k];
            f << "\n";
        }
        return true;
    }

    inline std::shared_ptr<sss::LayeredTable> loadTableCsv(const std::string& path, std::string* error = nullptr) {
        std::ifstream f(path);
        if (!f) { if (error) *error = "cannot open " + path; return nullptr; }
        std::vector<double> wavelengths, edges;
        std::vector<std::vector<double>> density;
        std::string line;
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::vector<std::string> cells;
            std::stringstream ss(line);
            for (std::string c; std::getline(ss, c, ',');) cells.push_back(c);
            if (cells.size() < 4) continue;
            if (cells[0] == "r_lo_m") {
                for (size_t i = 2; i < cells.size(); ++i) wavelengths.push_back(std::stod(cells[i].substr(2)));
                density.assign(wavelengths.size(), {});
                continue;
            }
            if (wavelengths.empty()) continue;
            if (edges.empty()) edges.push_back(std::stod(cells[0]));
            edges.push_back(std::stod(cells[1]));
            for (size_t i = 0; i < wavelengths.size(); ++i) density[i].push_back(std::stod(cells[2 + i]));
        }
        if (wavelengths.size() < 2 || edges.size() < 3) { if (error) *error = "no table found in " + path; return nullptr; }
        return buildLayeredTable(wavelengths, edges, density);
    }

    /// Scatter distance l (m) of Burley's profile that has the same albedo A and mean radius as the table channel.
    /// (Burley's mean radius is 2.5 d with d = l / shape(A).)
    inline Real matchedScatterDistance(const sss::ChannelTable& t) {
        const Real d = t.meanRadius() / Real(2.5);
        const Real a = std::clamp(t.albedo(), Real(0), Real(1));
        const Real e = std::abs(a - Real(0.8));
        const Real shape = Real(1.85) - a + Real(7) * e * e * e;
        return d * shape;
    }

} // namespace rayt::skin
