/**
 * @file IO/CsvReader.cpp
 * @brief Implementation of the lightweight CSV/TSV reader.
 */

#include "IO/CsvReader.hpp"

#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>

namespace rayt::io {

    /**
     * @brief Trims leading and trailing whitespace from a string in-place.
     * @param s The string to modify.
     */
    static inline void trimInPlace(std::string& s) {
        // Lambda to check for whitespace
        auto is_ws = [](unsigned char c) { return std::isspace(c) != 0; };

        // Find the first non-whitespace character
        auto it1 = std::find_if_not(s.begin(), s.end(), is_ws);

        // Find the last non-whitespace character (searching in reverse)
        auto it2 = std::find_if_not(s.rbegin(), s.rend(), is_ws).base();

        // If the string is all whitespace or empty
        if (it1 >= it2) { s.clear(); return; }

        // Update the string content to the trimmed range
        s.assign(it1, it2);
    }

    /**
     * @brief Checks if a line is empty or consists solely of whitespace.
     */
    static inline bool isLineEmptyOrWs(const std::string& s) {
        for (unsigned char c : s) {
            if (!std::isspace(c)) return false;
        }
        return true;
    }

    bool readDelimitedFile(
        const std::string& filename,
        CsvTable& outTable,
        std::string* outError,
        const CsvOptions& opt)
    {
        // Clear existing data in the output table
        outTable.clear();

        std::ifstream file(filename);
        if (!file.is_open()) {
            if (outError) *outError = "Could not open file: " + filename;
            return false;
        }

        std::string line;
        // Read the file line by line
        while (std::getline(file, line)) {

            // Skip empty lines if requested
            if (opt.skip_empty_lines && (line.empty() || isLineEmptyOrWs(line))) continue;

            // Handle comments: check if the line starts with the comment character
            if (opt.allow_comments) {
                std::string tmp = line;
                // Trim logic is needed because comments might be indented
                if (opt.trim_whitespace) trimInPlace(tmp);
                if (!tmp.empty() && tmp[0] == opt.comment_char) continue;  // Skip this line
            }

            CsvRow row;
            row.reserve(8);    // Pre-allocate memory for typical column counts

            std::string field;
            field.reserve(line.size());  // Reserve buffer to avoid frequent reallocations

            // Helper lambda to push the current field into the row
            auto flushField = [&]() {
                if (opt.trim_whitespace) trimInPlace(field);
                row.push_back(field);
                field.clear();
                };

            // Parse the line character by character
            for (size_t i = 0; i < line.size(); ++i) {
                char c = line[i];

                // Check if the current character matches any configured delimiter
                const bool is_delim =
                    (c == opt.delimiter) ||
                    (opt.allow_tabs && c == '\t') ||
                    (opt.allow_spaces && c == ' ');

                if (is_delim) {
                    flushField();    // Store the field accumulated so far

                    // Special handling for space delimiters:
                    // Collapse consecutive spaces into a single delimiter (common behavior for space-separated data)
                    if (opt.allow_spaces) {
                        while (i + 1 < line.size() && line[i + 1] == ' ') ++i;
                    }
                }
                else {
                    field.push_back(c);
                }
            }

            // Flush the last field after the loop ends
            flushField();

            // Optionally verify if the parsed row is effectively empty (e.g., contains only empty strings)
            if (opt.skip_empty_lines) {
                bool allEmpty = true;
                for (auto& f : row) {
                    if (!f.empty()) { allEmpty = false; break; }
                }
                if (allEmpty) continue;
            }

            // Move the constructed row into the final table
            outTable.push_back(std::move(row));
        }

        return true;
    }

} // namespace rayt::io