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
 */

#pragma once

#include <algorithm>
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

} // namespace rayt::io
