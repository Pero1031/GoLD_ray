#pragma once
/**
 * @file Sensor.hpp
 * @brief The camera's sensor (--sensor): the light on a pixel becomes electrons, the electrons become an integer. docs/CAMERA.md
 *
 * A pixel counts the photons that free an electron during the exposure. Photons arrive at random, so the count is random even under
 * steady light: Poisson distributed, its standard deviation the square root of its mean (shot noise; the signal-to-noise ratio of N
 * electrons is sqrt(N)). Reading the charge out adds read noise: one part before the amplifier that the ISO setting turns up, one part
 * after it. The analog-to-digital converter rounds the amplified signal to an integer, the digital number (DN), counted from a black
 * level so that noise below black is kept. A pixel holds at most its full well of electrons; as in a camera, the ADC's range ends there
 * (or at the white, above the base ISO), so a full pixel reads the top DN.
 *
 * The ISO setting does not change the light: it is the amplifier's gain, from electrons to DN. At the base ISO the full well is the
 * picture's white. At S times the base ISO the white is reached with 1/S of the electrons (the same light gives the same electrons, the
 * picture is S times brighter): the shot noise relative to the white grows as sqrt(S), and the read noise after the amplifier, counted
 * in electrons, shrinks as 1/S. Below the base ISO the well fills before the white.
 *
 * Colour: a single sensor sees colour through a mosaic of filters, one per pixel. The Bayer pattern RGGB repeats R G / G B over 2 x 2
 * pixels (green twice: it carries most of the luminance). Each pixel records one channel; the other two are interpolated from its
 * neighbours (demosaicing), on white-balanced values. Bilinear interpolation averages the nearest pixels of the missing colour.
 * Malvar, He and Cutler (ICASSP 2004) add to it a correction from the pixel's own colour, which follows the other colours' edges where
 * R, G and B change together, as they mostly do in natural pictures: half the false colour of bilinear at a grey edge, but more error
 * where the colours change independently.
 *
 * Values in and out are the camera's linear raw values (Color/CameraResponse::filmToCameraRaw: 1 = the G of the picture's white).
 * Not modelled: dark current, fixed-pattern noise (pixel to pixel differences of gain and offset), an optical low-pass filter.
 */

#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "Core/Types.hpp"

namespace rayt {

class Film;

enum class DemosaicMethod { Bilinear, Malvar };

struct SensorModel {
    Real fullWell = 50000;   ///< electrons a pixel holds: the picture's white at the base ISO
    Real baseIso = 100;      ///< the ISO speed whose white is the full well
    Real readNoise = 1.5;    ///< electrons RMS before the amplifier (the floor at high ISO)
    Real adcNoise = 3;       ///< electrons RMS after the amplifier, counted at the base ISO (base / ISO of it at a higher ISO)
    int bits = 14;           ///< of the analog-to-digital converter
    bool mosaic = true;      ///< a Bayer filter mosaic (RGGB): one channel per pixel, the others demosaiced; false: all three per pixel
    DemosaicMethod demosaic = DemosaicMethod::Malvar;
};

/// The sensor at one ISO speed.
struct SensorSetting {
    Real iso = 0;              ///< the ISO speed (the base ISO when the picture has no physical exposure)
    Real whiteElectrons = 0;   ///< mean electrons of a raw value of 1 (the picture's white): full well * base ISO / ISO
    Real readElectrons = 0;    ///< both read noises at this ISO (electrons RMS)
    Real dnPerElectron = 0;    ///< the conversion gain K: the DN range above black spans the white, or the full well below the base ISO
    int black = 0;             ///< DN of no light
    int maxDn = 0;             ///< 2^bits - 1
};

/// The setting of `m` at the ISO speed `iso` (0 or less: the base ISO).
SensorSetting MakeSensorSetting(const SensorModel& m, Real iso);

/// One raw value through the sensor: Poisson electrons of mean raw * whiteElectrons (a negative raw counts as no light) plus Gaussian
/// read noise, converted to DN, rounded and kept in [0, maxDn] (the top is the full well, or the white above the base ISO).
int SenseDn(Real raw, const SensorSetting& s, std::mt19937_64& rng);

/// A DN as a raw value again: (DN - black) / K / whiteElectrons.
inline Real DnToRaw(int dn, const SensorSetting& s) {
    return (Real(dn) - Real(s.black)) / (s.dnPerElectron * s.whiteElectrons);
}

/// The channel (0 R, 1 G, 2 B) that the Bayer pattern RGGB records at pixel (x, y), row 0 at the top: R G / G B.
inline int BayerChannel(int x, int y) {
    return (y & 1) == 0 ? ((x & 1) == 0 ? 0 : 1) : ((x & 1) == 0 ? 1 : 2);
}

/// The three channels of every pixel from an RGGB mosaic (w x h values, row by row; at least 3 x 3). Beyond the edges the mosaic is
/// mirrored about the edge pixel, which keeps the pattern (the pixel at -1 is the one at +1, of the same colour).
std::vector<Vector3> DemosaicBayer(const std::vector<Real>& mosaic, int w, int h, DemosaicMethod method);

/// Every pixel of the film through the sensor. The film holds radiance; `gain` turns it into the light on the sensor (raw 1 = white)
/// and back, so the film keeps its units. With the mosaic, `mosaicDn` (if given) receives the recorded DN, row by row. Returns the
/// share of pixels with a channel at the top DN. The same seed draws the same noise.
Real ApplySensor(Film& film, const SensorModel& m, const SensorSetting& s, Real gain, std::uint64_t seed, std::vector<std::uint16_t>* mosaicDn = nullptr);

/// "default", or a comma list of key=value: well, base, read, adc, bits (numbers), cfa=rggb|none, demosaic=malvar|bilinear (the others
/// keep their defaults). False with a message.
bool ParseSensorSpec(const std::string& spec, SensorModel& m, std::string& error);

} // namespace rayt
