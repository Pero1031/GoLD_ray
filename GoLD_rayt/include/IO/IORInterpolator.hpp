#pragma once

#include <algorithm>  
#include <cmath>      
#include <complex>    
#include <fstream>    
#include <iostream>   
#include <map>        
#include <sstream>    
#include <string>     
#include <vector>     

// Class for managing Complex Refractive Index (n + ik) data.
// Handles loading from CSV and linear interpolation for arbitrary wavelengths.
class IORInterpolator {
public:
    struct DataPoint {
        double wavelength; // Wavelength in nanometers (nm)
        double n;          // Real part: Refractive index
        double k;          // Imaginary part: Extinction coefficient

        // Operator overload for sorting by wavelength
        bool operator<(const DataPoint& other) const {
            return wavelength < other.wavelength;
        }
    };

    IORInterpolator() = default;

    // Loads optical data from a CSV file.
    // Supports formats from "RefractiveIndex.info" where n and k data might be listed in separate blocks.
    bool loadCSV(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "[Error] Could not open file: " << filename << std::endl;
            return false;
        }

        // Temporary storage using maps to associate values with wavelengths.
        // Using maps allows us to handle data where 'n' and 'k' rows are separated or out of order.
        std::map<double, double> n_map;
        std::map<double, double> k_map;

        std::string line;
        bool reading_k_mode = false; // Flag to track if we are currently reading 'k' values

        while (std::getline(file, line)) {
            if (line.empty()) continue;

            // Normalize delimiters: Replace commas and tabs with spaces for easy parsing
            std::replace(line.begin(), line.end(), ',', ' ');
            std::replace(line.begin(), line.end(), '\t', ' ');

            // Header detection and mode switching.
            // If the line contains "wl" (wavelength) and "k", switch to k-mode.
            // Otherwise, assume it's an n-mode block.
            if (line.find("wl") != std::string::npos) {
                if (line.find("k") != std::string::npos) {
                    reading_k_mode = true;
                }
                else {
                    reading_k_mode = false; // Default to n-mode
                }
                continue; // Skip the header line itself
            }

            std::stringstream ss(line);
            double wl_val, val;

            // Process only if two numerical values (Wavelength, Value) are successfully extracted
            if (ss >> wl_val >> val) {
                // Unit Conversion: Micrometers (um) -> Nanometers (nm)
                // Many databases provide data in um, but rendering usually works in nm.
                double wl_nm = wl_val * 1000.0;

                // Store in the appropriate map based on the current mode
                if (reading_k_mode) {
                    k_map[wl_nm] = val;
                }
                else {
                    n_map[wl_nm] = val;
                }
            }
        }

        if (n_map.empty()) {
            std::cerr << "[Warning] No data found." << std::endl;
            return false;
        }

        // Merge data: Iterate through n_map and try to find the corresponding k value.
        data_.clear();
        for (const auto& pair : n_map) {
            double wl = pair.first;
            double n_val = pair.second;
            double k_val = 0.0;

            // Search for the corresponding k value.
            // First, check for an exact match.
            if (k_map.count(wl)) {
                k_val = k_map[wl];
            }
            else {
                // Fallback: Nearest neighbor search.
                // If wavelengths in n-data and k-data don't align perfectly, 
                // find the closest k-value within a small tolerance (1.0 nm).
                auto it = k_map.lower_bound(wl);
                if (it != k_map.end() && std::abs(it->first - wl) < 1.0) {
                    k_val = it->second;
                }
                else if (it != k_map.begin()) {
                    auto prev = std::prev(it);
                    if (std::abs(prev->first - wl) < 1.0) {
                        k_val = prev->second;
                    }
                }
            }

            DataPoint p;
            p.wavelength = wl;
            p.n = n_val;
            p.k = k_val;
            data_.push_back(p);
        }

        // Ensure data is sorted by wavelength for binary search (lower_bound)
        std::sort(data_.begin(), data_.end());
        return true;
    }

    // Evaluates the complex IOR at a specific wavelength using linear interpolation.
    std::complex<double> evaluate(double wavelength_nm) const {
        // Boundary checks
        if (data_.empty()) return { 1.0, 0.0 };
        if (wavelength_nm <= data_.front().wavelength) return { data_.front().n, data_.front().k };
        if (wavelength_nm >= data_.back().wavelength) return { data_.back().n, data_.back().k };

        // Binary search to find the insertion point
        DataPoint target = { wavelength_nm, 0.0, 0.0 };
        auto it = std::lower_bound(data_.begin(), data_.end(), target);

        // Perform linear interpolation between p1 (previous) and p2 (current)
        const DataPoint& p2 = *it;
        const DataPoint& p1 = *(it - 1);

        double t = (wavelength_nm - p1.wavelength) / (p2.wavelength - p1.wavelength);
        return {
            p1.n + (p2.n - p1.n) * t,
            p1.k + (p2.k - p1.k) * t
        };
    }

    // Debug helper: Prints information about the loaded data range
    void printInfo() const {
        if (!data_.empty()) {
            std::cout << "Loaded " << data_.size() << " combined points." << std::endl;
            // 先頭データの確認
            std::cout << "Start: " << data_.front().wavelength << "nm, n="
                << data_.front().n << ", k=" << data_.front().k << std::endl;
            // 末尾データの確認
            std::cout << "End:   " << data_.back().wavelength << "nm, n="
                << data_.back().n << ", k=" << data_.back().k << std::endl;
        }
    }

private:
    std::vector<DataPoint> data_;
};