/**
 * @file IO/SpectralCurveLoader.hpp
 * @brief Loads a two-column "wavelength, value" table (light source SPD, reflectance, ...).
 *
 * Accepted: comma / tab / space separated text. Rows whose first two fields are not numbers
 * (headers such as "wavelength,power", comments) are skipped. Wavelengths below 50 are taken as
 * micrometers (RefractiveIndex.info style) and converted to nm, otherwise nm is assumed.
 */

#pragma once

#include <string>

#include "IO/CsvReader.hpp"
#include "Color/SpectralCurve.hpp"

namespace rayt::io {

    inline bool loadSpectralCurveCsv(const std::string& filename,
        rayt::color::SpectralCurve& out, std::string* outError = nullptr)
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
                out.add(Real(wl < 50.0 ? wl * 1000.0 : wl), Real(v));
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

} // namespace rayt::io
