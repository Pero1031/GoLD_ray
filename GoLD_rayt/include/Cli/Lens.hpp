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

    // The fields that `o` states replace the ones here.
    void overrideWith(const LensSpec& o) {
        if (o.focalLengthMm) focalLengthMm = o.focalLengthMm;
        if (o.fStop) fStop = o.fStop;
        if (o.focusDist) { focusDist = o.focusDist; focusAt.reset(); }
        if (o.focusAt) { focusAt = o.focusAt; focusDist.reset(); }
        if (o.blades) blades = o.blades;
        if (o.bladeRotationDeg) bladeRotationDeg = o.bladeRotationDeg;
    }
};

} // namespace rayt::cli
