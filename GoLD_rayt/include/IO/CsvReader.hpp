/**
 * @file IO/CsvReader.hpp
 * @brief Lightweight delimited-text reader for CSV/TSV-like files.
 *
 * Provides a small utility to read a text table (CSV/TSV/space-delimited) into
 * rows of string fields. This module intentionally does not impose any semantic
 * meaning (e.g., IOR, SPD). Domain-specific parsers should be implemented on top.
 *
 * See also:
 * - IO/ComplexIorLoader (domain-specific)
 * - IO/SpectrumTableLoader (domain-specific)
 */

#pragma once

#include <string>
#include <vector>

namespace rayt::io {

    struct CsvOptions {
        char delimiter = ',';          ///< Primary delimiter (e.g., ',' for CSV, '\t' for TSV)
        bool allow_tabs = true;        ///< Treat tabs as delimiters as well
        bool allow_spaces = false;     ///< Treat spaces as delimiters (useful for whitespace tables)
        bool skip_empty_lines = true;  ///< Skip empty lines
        bool trim_whitespace = true;   ///< Trim whitespace around fields
        bool allow_comments = true;    ///< Skip comment lines
        char comment_char = '#';       ///< Comment prefix
    };

    using CsvRow = std::vector<std::string>;
    using CsvTable = std::vector<CsvRow>;

    /**
     * @brief Reads a delimited text file into a table of string fields.
     * @param filename Path to file.
     * @param outTable Output rows.
     * @param outError Optional error message.
     * @param opt Reader options.
     * @return true on success.
     */
    bool readDelimitedFile(
        const std::string& filename,
        CsvTable& outTable,
        std::string* outError = nullptr,
        const CsvOptions& opt = {}
    );

} // namespace rayt::io