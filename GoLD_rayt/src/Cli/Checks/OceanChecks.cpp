// Checks of the ocean surface (--scene ocean): the wave spectrum, the sampled set of waves, the dispersion relation, the mesh.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Geometry/ProceduralMesh.hpp"

namespace rayt::cli {

int RunOceanTests(const CliOptions& opt) {
    bool allPassed = true;
    (void)opt;
    const double g = 9.81;

    // ---- 1. The set of waves has the variance of the spectrum it was drawn from ------------------------------------------------
    for (const auto& [U, gamma] : { std::pair<double, double>{ 6.0, 1.0 }, { 10.0, 1.0 }, { 10.0, 3.3 } }) {
        OceanParams P;
        P.windSpeed = Real(U); P.peakEnhance = Real(gamma);
        const auto waves = MakeOceanWaves(P);
        const OceanStats s = OceanSpectrumStats(P, waves);
        const double dh = std::abs(s.heightVarianceSet / s.heightVarianceSpectrum - 1.0);
        const double ds = std::abs(s.slopeVarianceSet / s.slopeVarianceSpectrum - 1.0);
        const bool ok = dh < 0.04 && ds < 0.04;
        std::cout << "[Ocean]   wind " << U << " m/s, gamma " << gamma << ": " << waves.size() << " waves, height variance set " << s.heightVarianceSet
            << " spectrum " << s.heightVarianceSpectrum << " (" << dh * 100 << " %), slope variance set " << s.slopeVarianceSet << " spectrum "
            << s.slopeVarianceSpectrum << " (" << ds * 100 << " %)" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    // ---- 2. Against the closed forms --------------------------------------------------------------------------------------------
    {
        // Pierson-Moskowitz: the variance of the heights of the whole spectrum is 0.0027 U^4 / g^2 (significant wave height 0.21 U^2 / g).
        const double U = 10.0;
        OceanParams P; P.windSpeed = Real(U);
        const OceanStats s = OceanSpectrumStats(P, MakeOceanWaves(P));
        const double full = 8.1e-3 / (4.0 * 0.74) * std::pow(U, 4) / (g * g);
        const double captured = s.heightVarianceSpectrum / full;
        const bool ok = captured > 0.93 && captured < 1.01;
        std::cout << "[Ocean]   Pierson-Moskowitz at 10 m/s: the height variance of the whole spectrum is " << full << " m^2 (Hs = " << 4 * std::sqrt(full)
            << " m, 0.21 U^2/g = " << 0.21 * U * U / g << " m); the waves kept (k > peak / 3, wavelength > " << P.minWavelength << " m) hold " << captured * 100 << " %"
            << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
        // Cox & Munk (1954): the mean-square slope of the sea surface is about 0.003 + 0.00512 U (U at 12.5 m). The Pierson-Moskowitz spectrum
        // saturates the slope (about alpha / 2 ln(k_max / k_peak), nearly independent of the wind), the real sea does not: the short waves
        // grow with the wind. So the two agree at 6 m/s and drift apart (to a factor 2.4 at 14 m/s); the far roughness of the water material
        // is taken from Cox-Munk for that reason. Here: the same order of magnitude.
        for (double u : { 6.0, 10.0, 14.0 }) {
            OceanParams Q; Q.windSpeed = Real(u);
            const OceanStats t = OceanSpectrumStats(Q, MakeOceanWaves(Q));
            const double cm = 0.003 + 0.00512 * u;
            const bool close = t.slopeVarianceSpectrum > cm / 2.6 && t.slopeVarianceSpectrum < cm * 2.6;
            std::cout << "[Ocean]   wind " << u << " m/s: mean-square slope " << t.slopeVarianceSpectrum << " (to " << Q.minWavelength << " m waves), Cox-Munk " << cm << (close ? "   OK (within a factor 2.6)" : "   FAIL") << "\n";
            allPassed = close && allPassed;
        }
    }

    // ---- 3. The heights of the sum, at random places, have that variance and zero mean -----------------------------------------
    {
        OceanParams P; P.windSpeed = Real(8);
        const auto waves = MakeOceanWaves(P);
        const OceanStats s = OceanSpectrumStats(P, waves);
        std::mt19937 rng(5);
        std::uniform_real_distribution<double> U(-3000.0, 3000.0);
        double sum = 0, sum2 = 0;
        const int n = 40000;
        for (int i = 0; i < n; ++i) {
            const double h = OceanHeightAt(waves, Real(0), Real(U(rng)), Real(U(rng)));
            sum += h; sum2 += h * h;
        }
        const double mean = sum / n, var = sum2 / n - mean * mean;
        const bool ok = std::abs(var / s.heightVarianceSet - 1.0) < 0.12 && std::abs(mean) < 0.1 * std::sqrt(var);
        std::cout << "[Ocean]   heights at " << n << " random places: mean " << mean << " m, variance " << var << " (the waves' " << s.heightVarianceSet << ")" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    // ---- 4. Deep-water dispersion: omega^2 = g k, and a wave moves with its phase speed sqrt(g / k) ----------------------------
    {
        OceanParams P;
        const auto waves = MakeOceanWaves(P);
        double worstOmega = 0, worstMove = 0;
        for (size_t i = 0; i < waves.size(); i += 17) {
            const OceanWave& w = waves[i];
            worstOmega = std::max(worstOmega, std::abs(double(w.omega) * w.omega / (g * w.k) - 1.0));
            const std::vector<OceanWave> one = { w };
            const double t = 3.7, c = std::sqrt(g / w.k);
            const double x = 12.3, z = -4.5;
            const double h0 = OceanHeightAt(one, Real(0), Real(x), Real(z));
            const double h1 = OceanHeightAt(one, Real(t), Real(x + c * t * w.kx / w.k), Real(z + c * t * w.kz / w.k));   // the crest has moved by c t along k
            worstMove = std::max(worstMove, std::abs(h0 - h1) / w.a);
        }
        const bool ok = worstOmega < 1e-5 && worstMove < 1e-3;
        std::cout << "[Ocean]   dispersion: worst omega^2 / (g k) - 1 = " << worstOmega << ", worst height change after moving with the phase speed (relative to the amplitude) = " << worstMove << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    // ---- 5. The directions follow cos^(2s)(theta / 2) around the wind -----------------------------------------------------------
    {
        OceanParams P; P.windDeg = Real(30); P.spread = Real(4);
        const auto waves = MakeOceanWaves(P);
        double within = 0, total = 0;
        const double wind = 30.0 * constants::PI / 180.0;
        for (const OceanWave& w : waves) {
            const double a2 = double(w.a) * w.a;
            total += a2;
            if (std::cos(std::atan2(double(w.kz), double(w.kx)) - wind) > std::cos(constants::PI / 4)) within += a2;   // within 45 deg of the wind
        }
        // The same fraction from the density: Integral over |theta| < 45 deg of cos^8(theta / 2) divided by the integral over all directions.
        double num = 0, den = 0;
        const int n = 20000;
        for (int i = 0; i < n; ++i) {
            const double th = -constants::PI + (i + 0.5) * 2 * constants::PI / n;
            const double d = std::pow(std::max(std::cos(0.5 * th), 0.0), 8.0);
            den += d;
            if (std::abs(th) < constants::PI / 4) num += d;
        }
        const bool ok = std::abs(within / total - num / den) < 0.03;
        std::cout << "[Ocean]   energy within 45 degrees of the wind: " << within / total << " (density: " << num / den << ")" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    // ---- 6. The mesh -------------------------------------------------------------------------------------------------------------
    {
        OceanParams P; P.rings = 80; P.sectors = 160; P.rInner = Real(1); P.rOuter = Real(600);
        const TriangleMesh m = MakeOceanMesh(P);
        bool finite = true, up = true;
        Real maxH = 0;
        for (size_t i = 0; i < m.p.size(); ++i) {
            if (math::hasNonFinite(m.p[i]) || math::hasNonFinite(m.n[i])) finite = false;
            if (m.n[i].y <= Real(0)) up = false;
            maxH = std::max(maxH, std::abs(m.p[i].y));
        }
        const bool ok = finite && up && m.numTriangles() > 10000 && maxH < Real(30);
        std::cout << "[Ocean]   mesh: " << m.p.size() << " vertices, " << m.numTriangles() << " triangles, largest height " << maxH << " m, all normals face up and are finite" << (ok ? "   OK" : "   FAIL") << "\n";
        allPassed = ok && allPassed;
    }

    std::cout << "[Ocean] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
