/**
 * @file include/Render/Diffraction.hpp
 * @brief Diffraction at the aperture (--diffraction): the Fraunhofer point spread function of the round or polygonal aperture,
 * integrated over the pixels, convolved with the picture. docs/CAMERA.md
 *
 * The amplitude in the image plane is the Fourier transform of the aperture, so the pattern of a point is |F|^2 at the angular frequency
 * 2 pi x / (lambda N): an Airy disc for a round aperture, a star with a spike across every edge for a polygon (n spikes for an even number
 * of blades, 2n for an odd one). Its size is lambda N, a few micrometres; a pixel of the 35 mm convention (36 mm over the long side of the
 * picture) is tens of micrometres, so at the real scale the core stays inside a pixel and only the spikes of very bright points show.
 * The pattern grows with the wavelength: the kernels are those of CIE X, Y and Z under a D65-like spectrum, and the picture is convolved
 * in XYZ.
 */

#pragma once

#include <complex>
#include <vector>

#include "Core/Types.hpp"

namespace rayt {

    class Film;

    /// The aperture and the sensor of the pattern.
    struct DiffractionSetup {
        Real fNumber = 8;           ///< N = focal length / aperture diameter
        int blades = 0;             ///< 0: round (Airy); 3 or more: a regular polygon of the same area (as Camera::setApertureShape)
        Real bladeRotation = 0;     ///< radians, as Camera::setApertureShape
        Real pixelPitch = 30e-6;    ///< metres on the sensor
        Real scale = 1;             ///< 1: the real size; s > 1 makes the pattern s times larger (as if the pixels were s times smaller)
        int radius = 255;           ///< the kernel covers (2 radius + 1)^2 pixels
    };

    /// The pattern integrated over each pixel, for CIE X, Y and Z: xyz[c][(radius + y) * (2 radius + 1) + radius + x] is the share of the
    /// light of a point that lands x pixels to the right of it and y pixels below. Each kernel is renormalized to sum 1.
    struct DiffractionKernels {
        int radius = 0;
        std::vector<double> xyz[3];
        double captured = 1;        ///< the share of Y inside the kernel before the renormalization (the rest is the far tail)
    };

    DiffractionKernels ComputeDiffractionKernels(const DiffractionSetup& s);

    /// Convolves the film's picture with the kernels (in CIE XYZ; the film stays linear sRGB). Outside the picture its border pixels repeat.
    void ApplyDiffraction(Film& film, const DiffractionKernels& k);

    /// The Fourier transform  F(q) = integral over the aperture of exp(-i q.r) d^2 r  of the aperture of unit diameter (area pi / 4; a polygon of
    /// that area), at the angular frequency q (1 / unit length). Exposed for the checks.
    std::complex<double> ApertureTransform(int blades, double rotation, double qx, double qy);

} // namespace rayt
