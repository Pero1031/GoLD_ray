// The film as a given camera sees (--camera-response); see Color/CameraResponse.hpp.
#include "pch.h"
#include "Color/CameraResponse.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "Color/ColorSpace.hpp"
#include "Color/Illuminants.hpp"
#include "Color/SpectralConversion.hpp"

namespace rayt::color {

    namespace {

        constexpr Real kCctA = Real(2856), kCctD65 = Real(6504);

        /// A camera's colour fit under one light: the raw values that become grey of 1 (the light's white per unit of its Y) and the
        /// matrix from white-balanced raw to the film's RGB, fitted on smooth reflectances under the light; the targets are the CIE
        /// film's colours of them, balanced from the light's white to D65 with CAT16 (`balance`; the identity for D65).
        struct Fit {
            Vector3 whiteRaw = Vector3(1.0);
            Matrix3 matrix = Matrix3(Real(1));
            Matrix3 balance = Matrix3(Real(1));
            double rms = 0;
            int count = 0;
            bool ok = false;
        };

        Fit FitUnder(const CameraResponse& cam, const std::function<Real(Real)>& light) {
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
            Fit fit;
            Vector3 rawWhite, filmWhite;
            integrate(light, rawWhite, filmWhite);
            if (!(rawWhite.x > 0 && rawWhite.y > 0 && rawWhite.z > 0)) return fit;
            fit.balance = filmBalanceMatrix(filmWhite);
            const Real y = (fit.balance * filmWhite).y;   // the balanced white is grey of the light's Y (relative to D65 of unit scale)
            if (!(y > 0)) return fit;
            fit.whiteRaw = rawWhite / y;

            // The colour matrix: smooth reflectances (the renderer's spectra of 6 x 6 x 6 RGB colours) under the light; each row r
            // minimizes sum (m_r . c - t_r)^2 with m_r . (1, 1, 1) = 1, so that the white stays grey: m = m_ls - A^-1 1 (1 . m_ls - 1) / (1 . A^-1 1).
            const Real levels[6] = { 0.04, 0.2, 0.4, 0.6, 0.8, 0.96 };
            std::vector<Vector3> balanced, target;
            for (const Real r : levels)
                for (const Real g : levels)
                    for (const Real b : levels) {
                        const Vector3 colour(r, g, b);
                        Vector3 raw, film;
                        integrate([&](Real l) { return light(l) * rgbToReflectance(colour, l); }, raw, film);
                        balanced.push_back(raw / fit.whiteRaw);
                        target.push_back(fit.balance * film);
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
                for (int j = 0; j < 3; ++j) fit.matrix[j][r] = row[j];   // row r of the matrix
            }
            double sq = 0;
            for (size_t k = 0; k < balanced.size(); ++k) sq += double(glm::dot(fit.matrix * balanced[k] - target[k], fit.matrix * balanced[k] - target[k]));
            fit.count = static_cast<int>(balanced.size());
            fit.rms = std::sqrt(sq / (3.0 * double(balanced.size())));
            fit.ok = true;
            return fit;
        }

        /// Makes a balance the film's: gains, matrix and the weights per wavelength.
        void UseBalance(CameraResponse& cam, const Vector3& whiteRaw, const Matrix3& matrix, const std::string& name) {
            cam.whiteRaw = whiteRaw;
            cam.matrix = matrix;
            cam.balancedTo = name;
            cam.weights.resize(cam.curves.size());
            for (size_t i = 0; i < cam.curves.size(); ++i) cam.weights[i] = cam.matrix * (cam.curves[i] / cam.whiteRaw);
        }

    } // namespace

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
        const Fit d65 = FitUnder(cam, [](Real l) { return evaluateD65(l); });
        const Fit a = FitUnder(cam, [](Real l) { return evaluateBlackbody(l, kCctA); });
        if (!d65.ok || !a.ok) { error = name + ": a channel sees nothing of D65 or of illuminant A between 380 and 780 nm"; return false; }
        cam.d65Raw = d65.whiteRaw;   // D65 of unit scale has Y 1: its white is its raw
        cam.matrixD65 = d65.matrix;
        cam.whiteRawA = a.whiteRaw;
        cam.matrixA = a.matrix;
        cam.catA = a.balance;
        cam.trainingCount = d65.count;
        cam.fitRms = d65.rms;
        UseBalance(cam, d65.whiteRaw, d65.matrix, "D65");
        cam.active = true;
        cameraResponse() = std::move(cam);
        return true;
    }

    bool setCameraBalance(const std::function<Real(Real)>& light, const std::string& name) {
        CameraResponse& cam = cameraResponse();
        if (!cam.active) return false;
        const Fit fit = FitUnder(cam, light);
        if (!fit.ok) return false;
        cam.trainingCount = fit.count;
        cam.fitRms = fit.rms;
        UseBalance(cam, fit.whiteRaw, fit.matrix, name);
        return true;
    }

    Real setCameraBalanceToRawWhite(const Vector3& rawWhite) {
        CameraResponse& cam = cameraResponse();
        if (!cam.active || !(rawWhite.x > 0 && rawWhite.y > 0 && rawWhite.z > 0)) return 0;
        // The white's XYZ and colour temperature depend on the matrix, the matrix on the temperature: a few rounds settle both (as a DNG
        // reader finds the matrix for its as-shot white). The matrices mix linearly in 1 / temperature between the two calibrations.
        const Matrix3 xyzA = cameraRawFromXYZ(true), xyzD65 = cameraRawFromXYZ(false);
        Real w = 0, cct = kCctD65;
        Vector3 xyz(0.0);
        for (int round = 0; round < 8; ++round) {
            xyz = glm::inverse(xyzA * w + xyzD65 * (Real(1) - w)) * rawWhite;
            cct = correlatedColourTemperature(xyz);
            w = std::clamp((Real(1) / cct - Real(1) / kCctD65) / (Real(1) / kCctA - Real(1) / kCctD65), Real(0), Real(1));
        }
        if (!(xyz.y > 0)) return 0;
        std::ostringstream name;
        name << "auto (" << std::lround(cct) << " K)";
        UseBalance(cam, rawWhite * cam.d65Raw.y / xyz.y, cam.matrixA * w + cam.matrixD65 * (Real(1) - w), name.str());
        return cct;
    }

    Vector3 filmToCameraRaw(const Vector3& rgb) {
        const CameraResponse& cam = cameraResponse();
        if (!cam.active) return rgb;
        return (glm::inverse(cam.matrix) * rgb) * cam.whiteRaw / cam.d65Raw.y;
    }

    Vector3 cameraRawToFilm(const Vector3& raw) {
        const CameraResponse& cam = cameraResponse();
        if (!cam.active) return raw;
        return cam.matrix * (raw * cam.d65Raw.y / cam.whiteRaw);
    }

    Vector3 cameraWhiteBalanceGains() {
        const CameraResponse& cam = cameraResponse();
        return cam.active ? cam.d65Raw.y / cam.whiteRaw : Vector3(1.0);
    }

    Matrix3 cameraRawFromXYZ(bool illuminantA) {
        const Matrix3 toXYZ = filmToXYZ();   // the CIE film's RGB -> XYZ (its own scale; D65 of unit scale is the film's (1, 1, 1))
        const Matrix3 toFilm = glm::inverse(toXYZ) * (toXYZ * Vector3(1.0)).y;   // XYZ with D65's Y = 1 -> the film's RGB
        const CameraResponse& cam = cameraResponse();
        if (!cam.active) return toFilm;
        // XYZ under the light -> the CIE film's RGB -> its balance to D65 (the fit's targets) -> undo the fit's matrix and gains
        const Vector3 white = illuminantA ? cam.whiteRawA : cam.d65Raw;
        const Matrix3 matrix = illuminantA ? cam.matrixA : cam.matrixD65;
        const Matrix3 balance = illuminantA ? cam.catA : Matrix3(Real(1));
        Matrix3 gains(Real(0));
        for (int i = 0; i < 3; ++i) gains[i][i] = white[i] / cam.d65Raw.y;
        return gains * glm::inverse(matrix) * balance * toFilm;
    }

    Matrix3 cameraCalibrationMatrix(bool illuminantA) {
        const CameraResponse& cam = cameraResponse();
        if (!cam.active) return Matrix3(Real(1));
        return illuminantA ? cam.matrixA : cam.matrixD65;
    }

    Matrix3 filmBalanceMatrix(const Vector3& filmWhite) {
        const Matrix3 toXYZ = filmToXYZ();
        const Vector3 src = toXYZ * filmWhite, d65 = toXYZ * Vector3(1.0);   // the film's (1, 1, 1) is D65
        if (!(src.y > Real(0)) || !cat16WhiteValid(src)) return Matrix3(Real(1));
        const Vector3 dst = d65 * (src.y / d65.y);
        return glm::inverse(toXYZ) * cat16Matrix(src, dst, Real(1)) * toXYZ;
    }

    Real correlatedColourTemperature(const Vector3& xyz) {
        const Real sum = xyz.x + xyz.y + xyz.z;
        if (!(sum > 0)) return 0;
        const Real x = xyz.x / sum, y = xyz.y / sum;
        const Real n = (x - Real(0.3320)) / (Real(0.1858) - y);
        return ((Real(449) * n + Real(3525)) * n + Real(6823.3)) * n + Real(5520.33);
    }

} // namespace rayt::color
