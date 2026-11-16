// The film as a given camera sees (--camera-response); see Color/CameraResponse.hpp.
#include "pch.h"
#include "Color/CameraResponse.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "Color/Illuminants.hpp"
#include "Color/SpectralConversion.hpp"

namespace rayt::color {

    CameraResponse& cameraResponse() {
        static CameraResponse response;
        return response;
    }

    void resetCameraResponse() { cameraResponse() = CameraResponse(); }

    bool setCameraResponse(const std::string& path, std::string& error) {
        std::ifstream in(path);
        if (!in) { error = "cannot open " + path; return false; }
        std::vector<Real> nm;
        std::vector<Vector3> rgb;
        for (std::string line; std::getline(in, line);) {
            if (const size_t hash = line.find('#'); hash != std::string::npos) line.erase(hash);
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream ss(line);
            double l, r, g, b;
            if (!(ss >> l >> r >> g >> b)) continue;   // a header or an empty line
            if (!nm.empty() && !(l > nm.back())) { error = path + ": the wavelengths must increase"; return false; }
            nm.push_back(Real(l));
            rgb.push_back(Vector3(r, g, b));
        }
        std::string name = path;
        if (const size_t slash = name.find_last_of("/\\"); slash != std::string::npos) name = name.substr(slash + 1);
        return setCameraResponse(nm, rgb, name, error);
    }

    bool setCameraResponse(const std::vector<Real>& nm, const std::vector<Vector3>& rgb, const std::string& name, std::string& error) {
        if (nm.size() < 2 || nm.size() != rgb.size()) { error = name + ": needs two or more wavelengths with R, G and B"; return false; }
        CameraResponse cam;
        cam.name = name;
        const int count = CameraResponse::kLastNm - CameraResponse::kFirstNm + 1;
        cam.curves.assign(static_cast<size_t>(count), Vector3(0.0));
        for (int i = 0; i < count; ++i) {
            const Real l = Real(CameraResponse::kFirstNm + i);
            if (l < nm.front() || l > nm.back()) continue;
            const size_t k = std::min(static_cast<size_t>(std::upper_bound(nm.begin(), nm.end(), l) - nm.begin()), nm.size() - 1);
            const Real t = (l - nm[k - 1]) / (nm[k] - nm[k - 1]);
            cam.curves[static_cast<size_t>(i)] = rgb[k - 1] * (Real(1) - t) + rgb[k] * t;
        }

        // Integrals over 380..780 nm with the rule of the CIE film (1 nm, trapezoid): the raw of a spectrum, and the CIE film's RGB of it
        auto integrate = [&](auto&& spectrum, Vector3& raw, Vector3& film) {
            raw = Vector3(0.0);
            film = Vector3(0.0);
            for (int l = 380; l <= 780; ++l) {
                const Real w = (l == 380 || l == 780) ? Real(0.5) : Real(1);
                const Real e = spectrum(Real(l)) * w;
                raw += cam.curves[static_cast<size_t>(l - CameraResponse::kFirstNm)] * e;
                film += rgbSensorWeights(Real(l)) * e;
            }
        };
        Vector3 filmWhite;
        integrate([](Real l) { return evaluateD65(l); }, cam.whiteRaw, filmWhite);
        if (!(cam.whiteRaw.x > 0 && cam.whiteRaw.y > 0 && cam.whiteRaw.z > 0)) { error = name + ": a channel sees nothing of D65 between 380 and 780 nm"; return false; }

        // The colour matrix: smooth reflectances (the renderer's spectra of 6 x 6 x 6 RGB colours) under D65; each row r minimizes
        // sum (m_r . c - t_r)^2 with m_r . (1, 1, 1) = 1, so that the white stays (1, 1, 1): m = m_ls - A^-1 1 (1 . m_ls - 1) / (1 . A^-1 1).
        const Real levels[6] = { 0.04, 0.2, 0.4, 0.6, 0.8, 0.96 };
        std::vector<Vector3> balanced, target;
        for (const Real r : levels)
            for (const Real g : levels)
                for (const Real b : levels) {
                    const Vector3 colour(r, g, b);
                    Vector3 raw, film;
                    integrate([&](Real l) { return evaluateD65(l) * rgbToReflectance(colour, l); }, raw, film);
                    balanced.push_back(raw / cam.whiteRaw);
                    target.push_back(film);
                }
        Matrix3 a(Real(0));   // sum of c c^T (symmetric; glm stores columns)
        Vector3 ct[3] = { Vector3(0.0), Vector3(0.0), Vector3(0.0) };   // sum of c t_r, per output row r
        for (size_t k = 0; k < balanced.size(); ++k) {
            const Vector3& c = balanced[k];
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) a[j][i] += c[i] * c[j];
            for (int r = 0; r < 3; ++r) ct[r] += c * target[k][r];
        }
        const Matrix3 ai = glm::inverse(a);
        const Vector3 ones(1.0), aiOnes = ai * ones;
        for (int r = 0; r < 3; ++r) {
            const Vector3 ls = ai * ct[r];
            const Vector3 row = ls - aiOnes * ((glm::dot(ones, ls) - Real(1)) / glm::dot(ones, aiOnes));
            for (int j = 0; j < 3; ++j) cam.matrix[j][r] = row[j];   // row r of the matrix
        }
        double sq = 0;
        for (size_t k = 0; k < balanced.size(); ++k) sq += double(glm::dot(cam.matrix * balanced[k] - target[k], cam.matrix * balanced[k] - target[k]));
        cam.trainingCount = static_cast<int>(balanced.size());
        cam.fitRms = std::sqrt(sq / (3.0 * double(balanced.size())));

        cam.weights.resize(cam.curves.size());
        for (size_t i = 0; i < cam.curves.size(); ++i) cam.weights[i] = cam.matrix * (cam.curves[i] / cam.whiteRaw);
        cam.active = true;
        cameraResponse() = std::move(cam);
        return true;
    }

    Vector3 filmToCameraRaw(const Vector3& rgb) {
        const CameraResponse& cam = cameraResponse();
        if (!cam.active) return rgb;
        return (glm::inverse(cam.matrix) * rgb) * cam.whiteRaw / cam.whiteRaw.y;
    }

} // namespace rayt::color
