// src/IO/ComplexIorLoader.cpp

#include "pch.h"

#include "IO/ComplexIorLoader.hpp"

#include <vector>
#include <algorithm>
#include <cctype>

#include "Core/Precision.hpp"

namespace rayt::io {

    namespace {

        enum class BlockMode { None, N, K };

        struct WlVal {
            Real lambda_nm;
            Real v;
            bool operator<(const WlVal& o) const { return lambda_nm < o.lambda_nm; }
        };

        static inline bool iequals(const std::string& a, const std::string& b) {
            if (a.size() != b.size()) return false;
            for (size_t i = 0; i < a.size(); ++i) {
                if (std::tolower((unsigned char)a[i]) != std::tolower((unsigned char)b[i])) return false;
            }
            return true;
        }

        static inline bool parseReal(const std::string& s, Real& out) {
            // Very small helper: use std::stod but keep it isolated.
            // (In production you might want from_chars; stod is fine for IO.)
            try {
                size_t idx = 0;
                double d = std::stod(s, &idx);
                if (idx == 0) return false;
                out = (Real)d;
                return true;
            }
            catch (...) {
                return false;
            }
        }

        static inline Real toNmHeuristic(Real wl) {
            // RefractiveIndex.info often uses micrometers (um): ~0.18 .. 2.0
            // nm would be ~180 .. 2000. We detect by magnitude.
            if (wl < Real(50)) {
                return wl * Real(1000); // um -> nm
            }
            return wl; // assume already nm
        }

        static inline void sortAndUnique(std::vector<Real>& xs) {
            std::sort(xs.begin(), xs.end());
            xs.erase(std::unique(xs.begin(), xs.end(),
                [](Real a, Real b) { return std::abs(a - b) < Real(1e-12); }), xs.end());
        }

        static Real evalLinear(const std::vector<WlVal>& data, Real lambda_nm, Real defaultValue) {
            if (data.empty()) return defaultValue;

            if (lambda_nm <= data.front().lambda_nm) return data.front().v;
            if (lambda_nm >= data.back().lambda_nm)  return data.back().v;

            WlVal key{ lambda_nm, Real(0) };
            auto it = std::lower_bound(data.begin(), data.end(), key);
            const WlVal& b = *it;
            const WlVal& a = *(it - 1);

            Real t = (lambda_nm - a.lambda_nm) / (b.lambda_nm - a.lambda_nm);
            return a.v + (b.v - a.v) * t;
        }

    } // namespace

    bool loadComplexIor_RefractiveIndexInfoCsv(
        const std::string& filename,
        rayt::color::ComplexIorTable& outTable,
        std::string* outError,
        const CsvOptions& csvOpt)
    {
        outTable.samples().clear();

        CsvOptions opt = csvOpt;
        opt.delimiter = ',';        // RefractiveIndex.info uses comma
        opt.allow_tabs = true;
        opt.allow_spaces = false;

        CsvTable table;
        std::string err;
        if (!readDelimitedFile(filename, table, &err, opt)) {
            if (outError) *outError = err;
            return false;
        }

        BlockMode mode = BlockMode::None;
        std::vector<WlVal> nData;
        std::vector<WlVal> kData;

        for (const auto& row : table) {
            if (row.size() < 2) continue;

            // Header switch: "wl,n" or "wl,k"
            if (iequals(row[0], "wl") && iequals(row[1], "n")) { mode = BlockMode::N; continue; }
            if (iequals(row[0], "wl") && iequals(row[1], "k")) { mode = BlockMode::K; continue; }

            Real wl = 0, val = 0;
            if (!parseReal(row[0], wl) || !parseReal(row[1], val)) {
                // Non-numeric line (or junk) -> ignore
                continue;
            }

            Real wl_nm = toNmHeuristic(wl);

            if (mode == BlockMode::N) {
                nData.push_back(WlVal{ wl_nm, val });
            }
            else if (mode == BlockMode::K) {
                kData.push_back(WlVal{ wl_nm, val });
            }
            else {
                // If no header was seen yet, ignore.
                continue;
            }
        }

        if (nData.empty()) {
            if (outError) *outError = "No 'wl,n' block data found in: " + filename;
            return false;
        }
        if (kData.empty()) {
            // kが無い場合も、金属以外で“etaだけ”みたいなケースがあり得るので
            // ここは失敗にせず k=0 として扱う方針もありですわ。
            // 今回は警告メッセージを返しつつ続行します。
            if (outError) *outError = "Warning: No 'wl,k' block data found. k will be set to 0.";
        }

        std::sort(nData.begin(), nData.end());
        std::sort(kData.begin(), kData.end());

        // Build union wavelength grid from both blocks
        std::vector<Real> wlGrid;
        wlGrid.reserve(nData.size() + kData.size());
        for (auto& p : nData) wlGrid.push_back(p.lambda_nm);
        for (auto& p : kData) wlGrid.push_back(p.lambda_nm);
        sortAndUnique(wlGrid);

        auto& dst = outTable.samples();
        dst.reserve(wlGrid.size());

        for (Real wl_nm : wlGrid) {
            Real eta = evalLinear(nData, wl_nm, Real(1));
            Real k = evalLinear(kData, wl_nm, Real(0));
            dst.push_back(rayt::color::ComplexIorTable::Sample{ wl_nm, eta, k });
        }

        // Final sanity: ensure sorted
        std::sort(dst.begin(), dst.end());
        return true;
    }

} // namespace rayt::io