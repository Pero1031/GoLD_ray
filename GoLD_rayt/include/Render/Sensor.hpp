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
 * Values in and out are the camera's linear raw values (Color/CameraResponse::filmToCameraRaw: 1 = the G of the picture's white).
 * Not modelled: the colour filter mosaic, dark current, fixed-pattern noise (pixel to pixel differences of gain and offset).
 */

#include <cstdint>
#include <random>
#include <string>

#include "Core/Types.hpp"

namespace rayt {

class Film;

struct SensorModel {
    Real fullWell = 50000;   ///< electrons a pixel holds: the picture's white at the base ISO
    Real baseIso = 100;      ///< the ISO speed whose white is the full well
    Real readNoise = 1.5;    ///< electrons RMS before the amplifier (the floor at high ISO)
    Real adcNoise = 3;       ///< electrons RMS after the amplifier, counted at the base ISO (base / ISO of it at a higher ISO)
    int bits = 14;           ///< of the analog-to-digital converter
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

/// Every pixel of the film through the sensor. The film holds radiance; `gain` turns it into the light on the sensor (raw 1 = white)
/// and back, so the film keeps its units. Returns the share of pixels with a channel at the top DN. The same seed draws the same noise.
Real ApplySensor(Film& film, const SensorSetting& s, Real gain, std::uint64_t seed);

/// "default", or a comma list of key=value from well, base, read, adc, bits (the others keep their defaults). False with a message.
bool ParseSensorSpec(const std::string& spec, SensorModel& m, std::string& error);

} // namespace rayt
