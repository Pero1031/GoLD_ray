/**
 * @file IO/FadingSeries.hpp
 * @brief A measured reflectance that changes with the time of exposure: R(lambda, t) from a series of measurements.
 *
 * Layout of a series folder (what tools/data/fading_refangle.py writes, see assets/Fading_measured):
 *   <dir>/hours.csv          "round,hours" (the exposure time of each measurement round, increasing)
 *   <dir>/r<round>/<color>.csv   "wavelength_nm, reflectance" of one patch in that round
 *
 * At a time between two rounds the curve is the linear blend of the two measured curves, wavelength by wavelength.
 * Outside the measured times the curve is held at the first or last round (nothing is extrapolated).
 * <dir>/km/<color>.csv, when present, holds a Kubelka-Munk fading law fitted to the rounds (FadingKm below): it answers any time.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Core/Types.hpp"
#include "Color/SpectralCurve.hpp"
#include "IO/SpectralCurveLoader.hpp"

namespace rayt::io {

    struct FadingRound {
        int round = 0;
        Real hours = Real(0);
    };

    /// The curve at one time and how it was made.
    struct FadingCurve {
        color::SpectralCurve curve;
        int roundLo = 0, roundHi = 0;   ///< the rounds that were blended (equal when the time sits on one round or outside the series)
        Real blend = Real(0);           ///< 0 = roundLo, 1 = roundHi
        Real firstHours = Real(0), lastHours = Real(0);
        bool clamped = false;           ///< the time lay outside the measured range and was held at its end
    };

    /// (1 - t) a + t b at every wavelength of either curve (the ends are held where one curve is shorter).
    inline color::SpectralCurve blendCurves(const color::SpectralCurve& a, const color::SpectralCurve& b, Real t) {
        color::SpectralCurve out;
        for (const auto* c : { &a, &b })
            for (const auto& s : c->samples())
                out.add(s.lambda_nm, (Real(1) - t) * a.evaluate(s.lambda_nm, color::Extrapolate::Clamp) + t * b.evaluate(s.lambda_nm, color::Extrapolate::Clamp));
        out.finalize();
        return out;
    }

    /// Reads "round,hours" lines (a header, comments and blank lines are skipped).
    inline bool loadFadingHours(const std::string& file, std::vector<FadingRound>& rounds, std::string* err = nullptr) {
        std::ifstream in(file);
        if (!in) { if (err) *err = "cannot open " + file; return false; }
        rounds.clear();
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            const size_t comma = line.find(',');
            if (comma == std::string::npos) continue;
            try {
                size_t i0 = 0, i1 = 0;
                const int r = std::stoi(line.substr(0, comma), &i0);
                const double h = std::stod(line.substr(comma + 1), &i1);
                if (i0 == 0 || i1 == 0) continue;
                rounds.push_back({ r, Real(h) });
            }
            catch (...) { continue; }   // the header
        }
        if (rounds.empty()) { if (err) *err = "no \"round,hours\" rows in " + file; return false; }
        for (size_t i = 1; i < rounds.size(); ++i)
            if (!(rounds[i].hours > rounds[i - 1].hours)) { if (err) *err = "the hours must increase from row to row in " + file; return false; }
        return true;
    }

    /// All the rounds of one patch: the measured curves and their exposure hours.
    struct FadingSeries {
        std::vector<FadingRound> rounds;
        std::vector<color::SpectralCurve> curves;   ///< one per round

        Real firstHours() const { return rounds.front().hours; }
        Real lastHours() const { return rounds.back().hours; }

        /// The two rounds around `hours` (held at the ends: nothing is extrapolated) and the weight of the second.
        void bracket(Real hours, size_t& lo, size_t& hi, Real& blend) const {
            const Real h = std::clamp(hours, firstHours(), lastHours());
            hi = 0;
            while (hi + 1 < rounds.size() && rounds[hi].hours < h) ++hi;   // the first round at or after h
            lo = (hi > 0 && rounds[hi].hours > h) ? hi - 1 : hi;
            blend = lo == hi ? Real(0) : (h - rounds[lo].hours) / (rounds[hi].hours - rounds[lo].hours);
        }
    };

    /// Reads hours.csv and the curve of patch `color` in every round of the series in `dir`.
    inline bool loadFadingSeries(const std::string& dir, const std::string& color, FadingSeries& out, std::string* err = nullptr) {
        out = FadingSeries();
        if (!loadFadingHours(dir + "/hours.csv", out.rounds, err)) return false;
        out.curves.resize(out.rounds.size());
        for (size_t i = 0; i < out.rounds.size(); ++i)
            if (!loadReflectanceCsv(dir + "/r" + std::to_string(out.rounds[i].round) + "/" + color + ".csv", out.curves[i], err)) return false;
        return true;
    }

    /// The reflectance of patch `color` of the series in `dir` after `hours` of exposure.
    inline bool loadFadingCurve(const std::string& dir, const std::string& color, Real hours, FadingCurve& out, std::string* err = nullptr) {
        FadingSeries series;
        if (!loadFadingSeries(dir, color, series, err)) return false;
        out = FadingCurve();
        out.firstHours = series.firstHours();
        out.lastHours = series.lastHours();
        out.clamped = hours < out.firstHours || hours > out.lastHours;
        size_t lo = 0, hi = 0;
        series.bracket(hours, lo, hi, out.blend);
        out.roundLo = series.rounds[lo].round;
        out.roundHi = series.rounds[hi].round;
        out.curve = lo == hi ? series.curves[lo] : blendCurves(series.curves[lo], series.curves[hi], out.blend);
        return true;
    }

    /// The fading law fitted to one patch of a series (tools/data/fading_km_fit.py writes <dir>/km/<color>.csv):
    ///     Q = K/S of R / scale,   Q(t) = Qinf + (Q0 - Qinf) exp(-k t),   R(t) = scale * (1 + Q - sqrt(Q^2 + 2 Q))
    /// per wavelength. Unlike the blend of the rounds it answers any time: after the last round (extrapolated) and before the first
    /// (negative hours: the colorant as it was before the exposure).
    struct FadingKm {
        std::vector<Real> lambda, r0, rInf, k;   ///< nm, reflectance at 0 h, reflectance when the colorant is gone, rate per hour
        Real scale = Real(1);                    ///< > 1 for reflectance factors above 1 (K/S needs R / scale < 1)

        static Real ks(Real r) { r = std::clamp(r, Real(1e-3), Real(0.999)); return (Real(1) - r) * (Real(1) - r) / (Real(2) * r); }
        /// 1 + q - sqrt(q^2 + 2 q), written so that a huge K/S gives a reflectance near 0 instead of inf - inf.
        static Real ksInv(Real q) { return Real(1) / (Real(1) + q + std::sqrt(q * q + Real(2) * q)); }

        /// R(lambda_i, hours) for sample i (the exponent is capped: far back in time the colorant only gets darker, never infinite).
        Real reflectance(size_t i, Real hours) const {
            const Real q0 = ks(r0[i] / scale), qInf = ks(rInf[i] / scale);
            const Real q = std::max(qInf + (q0 - qInf) * std::exp(std::min(-k[i] * hours, Real(50))), Real(0));
            return scale * ksInv(q);
        }

        color::SpectralCurve curveAt(Real hours) const {
            color::SpectralCurve c;
            for (size_t i = 0; i < lambda.size(); ++i) c.add(lambda[i], reflectance(i, hours));
            c.finalize();
            return c;
        }

        /// The law sampled at n evenly spaced times from h0 to h1 (a series whose blend of neighbouring rounds stays within a
        /// fraction of a percent of the law), for surfaces where the time varies from point to point.
        FadingSeries sampled(Real h0, Real h1, int n = 65) const {
            FadingSeries s;
            if (!(h1 > h0)) n = 1;
            for (int i = 0; i < n; ++i) {
                const Real h = n == 1 ? h0 : h0 + (h1 - h0) * Real(i) / Real(n - 1);
                s.rounds.push_back({ i + 1, h });
                s.curves.push_back(curveAt(h));
            }
            return s;
        }
    };

    /// Reads <dir>/km/<color>.csv ("wavelength_nm,R0,Rinf,k_per_hour,scale"; comments and the header are skipped).
    inline bool loadFadingKm(const std::string& dir, const std::string& color, FadingKm& out, std::string* err = nullptr) {
        const std::string file = dir + "/km/" + color + ".csv";
        std::ifstream in(file);
        if (!in) { if (err) *err = "cannot open " + file + " (make it with tools/data/fading_km_fit.py " + dir + ")"; return false; }
        out = FadingKm();
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);
            std::string cell;
            double v[5];
            int n = 0;
            try {
                while (n < 5 && std::getline(ss, cell, ',')) { size_t used = 0; v[n] = std::stod(cell, &used); if (used == 0) break; ++n; }
            }
            catch (...) { continue; }   // the header
            if (n < 5) continue;
            if (!out.lambda.empty() && !(Real(v[0]) > out.lambda.back())) { if (err) *err = "the wavelengths must increase in " + file; return false; }
            if (!(v[3] >= 0) || !(v[4] >= 1)) { if (err) *err = "k must not be negative and the scale must be at least 1 in " + file; return false; }
            out.lambda.push_back(Real(v[0])); out.r0.push_back(Real(v[1])); out.rInf.push_back(Real(v[2])); out.k.push_back(Real(v[3]));
            out.scale = Real(v[4]);
        }
        if (out.lambda.size() < 2) { if (err) *err = "fewer than two rows \"wavelength_nm,R0,Rinf,k_per_hour,scale\" in " + file; return false; }
        return true;
    }

    /// Hours of the fading experiment's pseudo-sunlight (SXL-500V2 at full power, 371.76 W/m^2) that do the damage of one day of
    /// direct sunlight outdoors: 1 day in the device = 2.911 days outdoors (Kataoka's thesis, 2025, section 6.1: the annual direct
    /// irradiation at Tsukuba in 2022, 2692.77 MJ/m^2, against the device, ultraviolet weighted 0.95 and visible light 0.05).
    /// Heat, humidity, weather and UV-absorbing coatings are left out, so this is an order of magnitude, not a prediction.
    constexpr Real kOutdoorDaysPerDeviceDay = Real(2.911);
    inline Real outdoorDaysToDeviceHours(Real days) { return days * Real(24) / kOutdoorDaysPerDeviceDay; }

} // namespace rayt::io
