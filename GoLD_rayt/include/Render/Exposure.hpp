#pragma once
/**
 * @file Exposure.hpp
 * @brief Physical exposure: how bright a camera with a given f-number, shutter time and ISO speed records the scene (--iso). docs/CAMERA.md
 *
 * The film's irradiance at the centre of the picture is E = L * G, where L is the scene's luminance and G (steradians) the projected solid
 * angle of the exit pupil seen from the film centre: pi / (4 N_w^2) for a lens that meets the sine condition, N_w = N * s / (s - f) the
 * working f-number of a thin lens focused at the distance s (it grows at close range: the image moves away from the lens). A real lens
 * (Render/LensSystem) gives G from its traced rays. The exposure is H = E t (lux seconds).
 *
 * The ISO speed is the standard output sensitivity of ISO 12232: a camera of speed S records the exposure 10 / S lux seconds as mid grey,
 * 0.18 linear (118 of 255 in sRGB). So a pixel's linear value is 0.018 S H. The scene's radiance has no unit of its own; `nits` states the
 * luminance (cd/m^2) of a radiance whose Y is 1, and the gain on the radiance is 0.018 S t nits G.
 *
 * Not included: the transmittance of the glass (a real lens passes about 90 %; the T-stop) and flare. The sensor's noise and full well:
 * Render/Sensor.hpp (--sensor).
 */

#include <cmath>

#include "Core/Types.hpp"
#include "Core/Precision.hpp"

namespace rayt {

/// G of a thin lens: pi / (4 N_w^2), N_w = N s / (s - f) (s and f in the same unit; s at or below f: 0, nothing reaches a film).
inline Real ThinLensIrradianceFactor(Real fNumber, Real focalLength, Real focusDistance) {
    if (!(fNumber > 0) || !(focusDistance > focalLength)) return Real(0);
    const Real working = fNumber * focusDistance / (focusDistance - focalLength);
    return Real(3.14159265358979323846) / (Real(4) * working * working);
}

/// The gain from the scene's radiance (Y = 1 is `nits` cd/m^2) to the camera's linear value (mid grey 0.18) for the ISO speed `iso`,
/// the shutter time `seconds` and the irradiance factor G (sr) at the film centre.
inline Real PhysicalExposureGain(Real iso, Real seconds, Real nits, Real irradianceFactor) {
    return Real(0.018) * iso * seconds * nits * irradianceFactor;
}

/// The luminance (cd/m^2) that the exposure records as mid grey (0.18): the light a reflected-light meter would expose for.
inline Real MidGreyLuminance(Real iso, Real seconds, Real irradianceFactor) {
    const Real g = Real(0.018) * iso * seconds * irradianceFactor;
    return g > 0 ? Real(0.18) / g : Real(0);
}

} // namespace rayt
