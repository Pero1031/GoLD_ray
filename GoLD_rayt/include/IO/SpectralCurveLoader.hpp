/**
 * @file IO/SpectralCurveLoader.hpp
 * @brief Loaders for two-column "wavelength, value" tables (light source SPD, reflectance, ...).
 *
 * Accepted: comma / tab / space separated text. Rows whose first two fields are not numbers
 * (headers such as "wavelength,power", comments) are skipped.
 *
 * Wavelength units (WavelengthUnit):
 *  - Auto : below 50 -> micrometers (RefractiveIndex.info style), otherwise nm
 *  - Nm / Um : as stated
 *  - Ev   : photon energy in eV (lambda_nm = 1239.84 / E), common in metal optical-constant tables
 */

#pragma once

#include <string>

#include "IO/CsvReader.hpp"
#include "Color/SpectralCurve.hpp"

namespace rayt::io {

    enum class WavelengthUnit { Auto, Nm, Um, Ev };

    /// Converts a table's wavelength-axis value to nm.
    inline Real toNm(Real x, WavelengthUnit unit) {
        switch (unit) {
        case WavelengthUnit::Nm: return x;
        case WavelengthUnit::Um: return x * Real(1000);
        case WavelengthUnit::Ev: return x > Real(0) ? Real(1239.841984) / x : Real(0);
        default:                 return x < Real(50) ? x * Real(1000) : x;
        }
    }

    inline bool loadSpectralCurveCsv(const std::string& filename,
        rayt::color::SpectralCurve& out, std::string* outError = nullptr,
        WavelengthUnit unit = WavelengthUnit::Auto)
    {
        CsvOptions opt;
        opt.delimiter = ',';
        opt.allow_tabs = true;
        opt.allow_spaces = true;

        CsvTable table;
        std::string err;
        if (!readDelimitedFile(filename, table, &err, opt)) {
            if (outError) *outError = err;
            return false;
        }

        out = rayt::color::SpectralCurve();
        for (const CsvRow& row : table) {
            if (row.size() < 2) continue;
            try {
                size_t i0 = 0, i1 = 0;
                const double wl = std::stod(row[0], &i0);
                const double v = std::stod(row[1], &i1);
                if (i0 == 0 || i1 == 0) continue;
                const Real nm = toNm(Real(wl), unit);
                if (nm > Real(0)) out.add(nm, Real(v));
            }
            catch (...) { continue; }   // header / junk line
        }
        out.finalize();

        if (out.size() < 2) {
            if (outError) *outError = "Fewer than two numeric (wavelength, value) rows in: " + filename;
            return false;
        }
        return true;
    }

    /**
     * @brief Loads a reflectance table R(lambda). Values are fractions (0..1); a table whose largest
     *        value exceeds 1.5 is taken to be in percent and divided by 100.
     */
    inline bool loadReflectanceCsv(const std::string& filename,
        rayt::color::SpectralCurve& out, std::string* outError = nullptr,
        WavelengthUnit unit = WavelengthUnit::Auto)
    {
        if (!loadSpectralCurveCsv(filename, out, outError, unit)) return false;

        Real maxValue = Real(0);
        for (const auto& s : out.samples()) maxValue = std::max(maxValue, s.value);
        if (maxValue > Real(1.5)) {
            rayt::color::SpectralCurve scaled;
            for (const auto& s : out.samples()) scaled.add(s.lambda_nm, s.value / Real(100));
            scaled.finalize();
            out = scaled;
        }
        return true;
    }

} // namespace rayt::io
