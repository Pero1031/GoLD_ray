#pragma once
// The lens of the camera as the user states it: focal length, f-number, focus, shape of the aperture.
// Every field is optional; a scene (built in, or "camera" of a scene file) gives a base, and the command line
// ("--focal-length", "--fstop", ...) overrides it field by field. RenderPipeline.cpp turns the result into a Camera.

#include <optional>
#include <stdexcept>
#include <string>

#include "Core/Types.hpp"
#include "Core/Precision.hpp"

namespace rayt::cli {

struct LensSpec {
    std::optional<Real> orthoHeight;       // > 0: an orthographic view this high (scene units; replaces the field of view), 0: perspective
    std::optional<std::string> lensFile;   // a real lens (Render/LensSystem.hpp): its prescription replaces the thin lens; "" = the thin lens
    std::optional<Real> focalLengthMm;     // 35 mm equivalent: the long side of the picture is the 36 mm side of the film (replaces the field of view)
    std::optional<Real> fStop;             // f-number N: aperture diameter = focal length / N (none or 0: a pinhole, no defocus blur)
    std::optional<Real> focusDist;         // distance of the plane of focus along the view axis (scene units, metres)
    std::optional<Vector3> focusAt;        // ... or the point to focus on (its distance along the view axis)
    std::optional<int> blades;             // aperture blades: 0 or none = a round aperture, 3 or more = a regular polygon
    std::optional<Real> bladeRotationDeg;  // turns the polygon
    std::optional<Real> diffraction;       // > 0: diffraction at the aperture (Render/Diffraction.hpp), its pattern this many times its real size; 0: off
    // Aberrations (Camera::setAberrations); none = an ideal lens
    std::optional<Real> distortionK1;      // radial distortion r' = r (1 + k1 r^2 + k2 r^4); k1 > 0 barrel, k1 < 0 pincushion
    std::optional<Real> distortionK2;
    std::optional<Real> lateralCA;         // size of the picture changes by this fraction per 100 nm (spectral build)
    std::optional<Real> axialCA;           // focal length changes by this fraction per 100 nm: a plane of focus per wavelength (spectral build)

    // Vignetting (Camera::setVignetting); none = off
    std::optional<Real> vignetteFStop;     // f-number N of the opening of the lens barrel: diameter = focal length / N (smaller N = larger opening = less vignetting)
    std::optional<Real> vignetteDistMm;    // how far in front of the aperture stop that opening is (default: half the focal length)
    std::optional<Real> naturalVignetting; // the cos^4 law: 0 = off, 1 = full

    // Motion blur (Camera::setShutter, setMotion); none = off
    std::optional<Real> shutter;           // seconds the shutter stays open (1/60 s = 0.0167); objects with a velocity and a moving camera blur over it
    std::optional<Vector3> velocity;       // the camera moves with this velocity (m/s)
    std::optional<Real> panYawRate;        // the view turns to the right by this many degrees per second
    std::optional<Real> panPitchRate;      // ... and up

    // The fields that `o` states replace the ones here.
    void overrideWith(const LensSpec& o) {
        if (o.orthoHeight) orthoHeight = o.orthoHeight;
        if (o.lensFile) lensFile = o.lensFile;
        if (o.diffraction) diffraction = o.diffraction;
        if (o.shutter) shutter = o.shutter;
        if (o.velocity) velocity = o.velocity;
        if (o.panYawRate) panYawRate = o.panYawRate;
        if (o.panPitchRate) panPitchRate = o.panPitchRate;
        if (o.vignetteFStop) vignetteFStop = o.vignetteFStop;
        if (o.vignetteDistMm) vignetteDistMm = o.vignetteDistMm;
        if (o.naturalVignetting) naturalVignetting = o.naturalVignetting;
        if (o.distortionK1) distortionK1 = o.distortionK1;
        if (o.distortionK2) distortionK2 = o.distortionK2;
        if (o.lateralCA) lateralCA = o.lateralCA;
        if (o.axialCA) axialCA = o.axialCA;
        if (o.focalLengthMm) focalLengthMm = o.focalLengthMm;
        if (o.fStop) fStop = o.fStop;
        if (o.focusDist) { focusDist = o.focusDist; focusAt.reset(); }
        if (o.focusAt) { focusAt = o.focusAt; focusDist.reset(); }
        if (o.blades) blades = o.blades;
        if (o.bladeRotationDeg) bladeRotationDeg = o.bladeRotationDeg;
    }
};

// A shutter time as typed: seconds ("0.02") or the photographer's fraction ("1/50"). Throws std::invalid_argument for anything else.
inline Real ParseShutterSeconds(const std::string& text) {
    const size_t slash = text.find('/');
    if (slash == std::string::npos) return std::stod(text);
    const Real num = std::stod(text.substr(0, slash)), den = std::stod(text.substr(slash + 1));
    if (!(den > 0)) throw std::invalid_argument("a shutter time of the form 1/N needs N above 0");
    return num / den;
}

} // namespace rayt::cli
