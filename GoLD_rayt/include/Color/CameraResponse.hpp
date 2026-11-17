/**
 * @file include/Color/CameraResponse.hpp
 * @brief The film as a given camera sees (--camera-response): its spectral sensitivities instead of the CIE observer. docs/CAMERA.md
 *
 * A camera's raw values are the light weighted with its R, G, B sensitivities, which are not a linear mix of the CIE colour matching
 * functions (the Luther condition does not hold): spectra that match for the eye can differ for the camera, and the other way round.
 * The film then does what a camera's processing does with them, each step linear:
 *   raw = integral of L(lambda) S(lambda)  ->  white balance (the white becomes grey)  ->  a 3 x 3 colour matrix to the film's RGB,
 * the matrix fitted by least squares on a set of smooth reflectances under the light the camera is balanced to (white kept white). The
 * three steps make one weight per wavelength, used by SpectrumToRGB in place of the CIE one; filmToCameraRaw() undoes the last two.
 *
 * The balance is D65 by default. Balanced to another light (--wb with a scene light, or --awb), the gains make that light's white grey
 * and the matrix is fitted under it, its targets the CIE film's colours adapted from the light's white to D65 with CAT16 (what --wb does
 * for the CIE film). The raw values keep D65's scale (the G of D65 of unit scale is 1), so the light on the sensor does not depend on
 * the balance. Fits under D65 and illuminant A are kept for the raw file's two calibrations and for the auto white balance, which
 * interpolates between them by the white's correlated colour temperature, as a DNG reader does.
 */

#pragma once

#include <algorithm>
#include <functional>
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
        Vector3 d65Raw = Vector3(1.0);     ///< the raw values of D65 of unit scale (their G sets the raw scale, whatever the balance)
        Vector3 whiteRaw = Vector3(1.0);   ///< the balance: the raw values that become grey of 1 (the balanced light's white, per unit Y)
        Matrix3 matrix = Matrix3(Real(1));  ///< white-balanced raw -> the film's RGB
        std::string balancedTo = "D65";    ///< what the balance was set to (for the log)
        Vector3 whiteRawA = Vector3(1.0);  ///< whiteRaw and matrix of the fits under D65 and illuminant A
        Matrix3 matrixA = Matrix3(Real(1)), matrixD65 = Matrix3(Real(1));
        Matrix3 catA = Matrix3(Real(1));   ///< the CIE film's CAT16 balance from illuminant A to D65 (the A fit's targets)
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

    /// Balances the active camera to a light given by its spectrum (any scale): gains that make its white grey, the matrix fitted under
    /// it. False (the balance unchanged) if the camera is not active or the light has no usable white.
    bool setCameraBalance(const std::function<Real(Real)>& light, const std::string& name);
    /// Balances the active camera to a white given as raw values (filmToCameraRaw's scale), as an auto white balance does: gains from
    /// it, the matrix interpolated between the A and D65 fits by its correlated colour temperature. Returns that temperature (K).
    Real setCameraBalanceToRawWhite(const Vector3& rawWhite);

    /// The camera's raw RGB of a film value: diag(whiteRaw) * matrix^-1 * rgb, scaled so that the G of D65 of unit scale is 1.
    /// Without a camera response the raw is the film's value itself (an ideal camera whose channels are the film's linear sRGB).
    Vector3 filmToCameraRaw(const Vector3& rgb);
    /// The inverse: the film value of a raw RGB (the camera's white balance and colour matrix).
    Vector3 cameraRawToFilm(const Vector3& raw);
    /// The white balance's gains on the raw channels: d65Raw.y / whiteRaw ((1, 1, 1) without a response).
    Vector3 cameraWhiteBalanceGains();
    /// CIE XYZ under a calibration light (illuminant A, or D65) -> the camera's raw RGB: a DNG's ColorMatrix before its normalization.
    /// Without a response: XYZ -> the film's linear sRGB.
    Matrix3 cameraRawFromXYZ(bool illuminantA = false);

    /// The colour matrix of a calibration (white-balanced raw -> the film's RGB) fitted under illuminant A, or D65 (identity without a
    /// response, whose raw is the film's RGB).
    Matrix3 cameraCalibrationMatrix(bool illuminantA);

    /// The CIE film's white balance from a white (the film's RGB of it) to D65 of the same Y, with CAT16: a matrix on the film's RGB.
    /// Identity if the white has no usable cone responses.
    Matrix3 filmBalanceMatrix(const Vector3& filmWhite);

    /// The correlated colour temperature (K) of an XYZ, by McCamy's cubic (1992; within a few kelvin near the Planckian locus).
    Real correlatedColourTemperature(const Vector3& xyz);

    /// The film's weight at a wavelength under the active response (linear between the 1 nm entries).
    inline Vector3 cameraFilmWeight(const CameraResponse& cam, Real nm) {
        const Real x = nm - Real(CameraResponse::kFirstNm);
        if (x < Real(0) || x > Real(CameraResponse::kLastNm - CameraResponse::kFirstNm)) return Vector3(0.0);
        const int i = std::min(static_cast<int>(x), CameraResponse::kLastNm - CameraResponse::kFirstNm - 1);
        const Real f = x - Real(i);
        return cam.weights[static_cast<size_t>(i)] * (Real(1) - f) + cam.weights[static_cast<size_t>(i) + 1] * f;
    }

} // namespace rayt::color
