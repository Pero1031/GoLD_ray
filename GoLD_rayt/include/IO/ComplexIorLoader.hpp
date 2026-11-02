// include/IO/ComplexIorLoader.hpp

#pragma once

#include <string>

#include "IO/CsvReader.hpp"
#include "IO/SpectralCurveLoader.hpp"
#include "Color/ComplexIorTable.hpp"

namespace rayt::io {

    /**
     * @brief Loads complex IOR data (eta, k) from RefractiveIndex.info-style CSV.
     *
     * Supported format (two blocks in one file):
     *   wl,n
     *   0.1879,1.28
     *   ...
     *
     *   wl,k
     *   0.1879,1.188
     *   ...
     *
     * Notes:
     * - wl is often in micrometers (um) on RefractiveIndex.info; this loader converts it to nm.
     * - The resulting table stores wavelengths in nm.
     */
    bool loadComplexIor_RefractiveIndexInfoCsv(
        const std::string& filename,
        rayt::color::ComplexIorTable& outTable,
        std::string* outError = nullptr,
        const CsvOptions& csvOpt = {}
    );

    /**
     * @brief Loads a table with one row per wavelength: "wavelength, n, k" (comma / tab / space
     *        separated; header and comment lines are skipped).
     * @param unit Unit of the wavelength column (Auto: <50 -> um, else nm; Ev for photon energy).
     */
    bool loadComplexIor_ThreeColumn(
        const std::string& filename,
        rayt::color::ComplexIorTable& outTable,
        WavelengthUnit unit = WavelengthUnit::Auto,
        std::string* outError = nullptr
    );

    /// Loads n(l) and k(l) from two separate two-column files.
    inline bool loadComplexIor_SeparateFiles(
        const std::string& nFile,
        const std::string& kFile,
        rayt::color::ComplexIorTable& outTable,
        WavelengthUnit unit = WavelengthUnit::Auto,
        std::string* outError = nullptr)
    {
        rayt::color::SpectralCurve n, k;
        if (!loadSpectralCurveCsv(nFile, n, outError, unit)) return false;
        if (!loadSpectralCurveCsv(kFile, k, outError, unit)) return false;
        outTable = rayt::color::ComplexIorTable::fromCurves(n, k);
        return true;
    }

    /// Detects the format (RefractiveIndex.info "wl,n / wl,k" blocks, otherwise 3 columns) and loads it.
    bool loadComplexIor(
        const std::string& filename,
        rayt::color::ComplexIorTable& outTable,
        WavelengthUnit unit = WavelengthUnit::Auto,
        std::string* outError = nullptr
    );

} // namespace rayt::io