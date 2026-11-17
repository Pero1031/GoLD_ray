/**
 * @file include/Color/CameraResponse.hpp
 * @brief The film as a given camera sees (--camera-response): its spectral sensitivities instead of the CIE observer. docs/CAMERA.md
 *
 * A camera's raw values are the light weighted with its R, G, B sensitivities, which are not a linear mix of the CIE colour matching
 * functions (the Luther condition does not hold): spectra that match for the eye can differ for the camera, and the other way round.
 * The film then does what a camera's processing does with them, each step linear:
 *   raw = integral of L(lambda) S(lambda)  ->  white balance (D65 becomes (1, 1, 1))  ->  a 3 x 3 colour matrix to the film's RGB,
 * the matrix fitted by least squares on a set of smooth reflectances under D65 (white kept white). The three steps make one weight per
 * wavelength, used by SpectrumToRGB in place of the CIE one; filmToCameraRaw() undoes the last two for the raw file.
 */

#pragma once

#include <algorithm>
#include <string>
#include <vector>

#include "Core/Types.hpp"
#include "Color/ChromaticAdaptation.hpp"

namespace rayt::color {

    struct CameraResponse {
        static constexpr int kFirstNm = 360;
        static constexpr int kLastNm = 830;

        std::string name;
        std::vector<Vector3> curves;    ///< raw sensitivities at kFirstNm..kLastNm, 1 nm (0 outside the data)
        std::vector<Vector3> weights;   ///< the film's weights: matrix * diag(1 / whiteRaw) * curves
        Vector3 whiteRaw = Vector3(1.0);   ///< the raw values of D65 of unit scale
        Matrix3 matrix = Matrix3(Real(1));  ///< white-balanced raw -> the film's RGB
        int trainingCount = 0;          ///< reflectances the matrix was fitted on
        double fitRms = 0;              ///< RMS of the fit's error on them (the film's RGB)
        bool active = false;
    };

    /// The response the film uses (inactive by default: the CIE observer). Set it before rendering starts.
    CameraResponse& cameraResponse();

    /// Reads a CSV of wavelength (nm), R, G, B ('#' comments; a line that is not numbers is skipped) and makes it the film's response.
    bool setCameraResponse(const std::string& path, std::string& error);
    /// The same from sensitivities at the given wavelengths (increasing); linear in between, 0 outside.
    bool setCameraResponse(const std::vector<Real>& nm, const std::vector<Vector3>& rgb, const std::string& name, std::string& error);
    void resetCameraResponse();

    /// The camera's raw RGB of a film value: diag(whiteRaw) * matrix^-1 * rgb, scaled so that the G of D65 of unit scale is 1.
    /// Without a camera response the raw is the film's value itself (an ideal camera whose channels are the film's linear sRGB).
    Vector3 filmToCameraRaw(const Vector3& rgb);
    /// The inverse: the film value of a raw RGB (the camera's white balance and colour matrix).
    Vector3 cameraRawToFilm(const Vector3& raw);
    /// The white balance's gains on the raw channels: whiteRaw.y / whiteRaw, so that D65 becomes grey ((1, 1, 1) without a response).
    Vector3 cameraWhiteBalanceGains();
    /// The raw RGB of a CIE XYZ (D65 white of Y 1 under the camera's balance; a DNG's ColorMatrix before its normalization).
    Matrix3 cameraRawFromXYZ();

    /// The film's weight at a wavelength under the active response (linear between the 1 nm entries).
    inline Vector3 cameraFilmWeight(const CameraResponse& cam, Real nm) {
        const Real x = nm - Real(CameraResponse::kFirstNm);
        if (x < Real(0) || x > Real(CameraResponse::kLastNm - CameraResponse::kFirstNm)) return Vector3(0.0);
        const int i = std::min(static_cast<int>(x), CameraResponse::kLastNm - CameraResponse::kFirstNm - 1);
        const Real f = x - Real(i);
        return cam.weights[static_cast<size_t>(i)] * (Real(1) - f) + cam.weights[static_cast<size_t>(i) + 1] * f;
    }

} // namespace rayt::color
