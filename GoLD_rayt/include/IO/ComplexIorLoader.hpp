// include/IO/ComplexIorLoader.hpp

#pragma once

#include <string>

#include "IO/CsvReader.hpp"
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

} // namespace rayt::io