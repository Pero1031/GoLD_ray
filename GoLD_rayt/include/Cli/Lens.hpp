#pragma once
// The lens of the camera as the user states it: focal length, f-number, focus, shape of the aperture.
// Every field is optional; a scene (built in, or "camera" of a scene file) gives a base, and the command line
// ("--focal-length", "--fstop", ...) overrides it field by field. RenderPipeline.cpp turns the result into a Camera.

#include <optional>

#include "Core/Types.hpp"
#include "Core/Precision.hpp"

namespace rayt::cli {

struct LensSpec {
    std::optional<Real> focalLengthMm;     // 35 mm equivalent: the long side of the picture is the 36 mm side of the film (replaces the field of view)
    std::optional<Real> fStop;             // f-number N: aperture diameter = focal length / N (none or 0: a pinhole, no defocus blur)
    std::optional<Real> focusDist;         // distance of the plane of focus along the view axis (scene units, metres)
    std::optional<Vector3> focusAt;        // ... or the point to focus on (its distance along the view axis)
    std::optional<int> blades;             // aperture blades: 0 or none = a round aperture, 3 or more = a regular polygon
    std::optional<Real> bladeRotationDeg;  // turns the polygon
    // Aberrations (Camera::setAberrations); none = an ideal lens
    std::optional<Real> distortionK1;      // radial distortion r' = r (1 + k1 r^2 + k2 r^4); k1 > 0 barrel, k1 < 0 pincushion
    std::optional<Real> distortionK2;
    std::optional<Real> lateralCA;         // size of the picture changes by this fraction per 100 nm (spectral build)
    std::optional<Real> axialCA;           // focal length changes by this fraction per 100 nm: a plane of focus per wavelength (spectral build)

    // The fields that `o` states replace the ones here.
    void overrideWith(const LensSpec& o) {
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

} // namespace rayt::cli
