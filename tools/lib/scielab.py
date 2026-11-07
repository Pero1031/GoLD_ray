"""S-CIELAB (Zhang & Wandell, 1996): CIELAB color differences after the image is filtered like the visual system does.

The eye resolves color much more coarsely than brightness (chromatic contrast sensitivity falls off at lower spatial frequencies), so a fine pattern of
colors blends into its surround at a normal viewing distance. S-CIELAB models this: convert the image to an opponent space (one lightness-like plane,
a red-green and a blue-yellow plane), blur each plane with the spatial filter of that channel (a sum of Gaussians, the widths in degrees of visual angle),
convert back to XYZ and take the usual CIELAB difference. The filters depend on how many pixels fall into one degree: pixels_per_degree.

The opponent matrix and the filter widths / weights are the published values of the original method, written from the paper's table (Zhang & Wandell 1996;
Zhang, Wandell, "Color image fidelity metrics evaluated using image distortion maps", Signal Processing 70, 1998). They have NOT been checked against a reference
implementation here; the study script therefore also shows how much the result depends on the viewing distance, which changes the widths in pixels the same
way an error in the table would. The kernels are normalized to unit sum per plane, so flat areas are unchanged.

Input images are XYZ with the white of the display at Y = 1 (the same convention as tools/lib/gamut.py).
"""
import numpy as np

# XYZ (D65) -> opponent planes: luminance, red-green, blue-yellow
M_XYZ_TO_OPP = np.array([[0.279, 0.720, -0.107],
                         [-0.449, 0.290, -0.077],
                         [0.086, -0.590, 0.501]])

# per plane: (weights, widths in degrees of visual angle)
FILTERS = (
    ((0.921, 0.105, -0.108), (0.0283, 0.133, 4.336)),   # luminance
    ((0.531, 0.330), (0.0392, 0.494)),                  # red-green
    ((0.488, 0.371), (0.0536, 0.386)),                  # blue-yellow
)


def _filter_plane(plane, weights, widths_px, pad):
    """Blur one plane with sum_i w_i * Gaussian(sigma_i) (sigmas in pixels), normalized to unit sum, via the Fourier transform."""
    h, w = plane.shape
    p = np.pad(plane, pad, mode="edge")
    fy = np.fft.fftfreq(p.shape[0])[:, None]
    fx = np.fft.rfftfreq(p.shape[1])[None, :]
    r2 = fx * fx + fy * fy
    gain = sum(wt * np.exp(-2.0 * np.pi ** 2 * (s ** 2) * r2) for wt, s in zip(weights, widths_px)) / sum(weights)
    out = np.fft.irfft2(np.fft.rfft2(p) * gain, s=p.shape)
    return out[pad:pad + h, pad:pad + w]


def spatial_filter_xyz(xyz, pixels_per_degree):
    """Applies the S-CIELAB spatial filtering to an (h, w, 3) XYZ image."""
    opp = xyz @ M_XYZ_TO_OPP.T
    h, w = opp.shape[:2]
    pad = int(min(max(h, w), 3 * max(s for _, ws in FILTERS for s in ws) * pixels_per_degree)) + 1
    out = np.empty_like(opp)
    for c, (weights, widths_deg) in enumerate(FILTERS):
        out[..., c] = _filter_plane(opp[..., c], weights, [s * pixels_per_degree for s in widths_deg], pad)
    return out @ np.linalg.inv(M_XYZ_TO_OPP).T


def pixels_per_degree(pixel_pitch_mm, distance_mm):
    """How many pixels fall into one degree of visual angle."""
    return distance_mm * np.tan(np.radians(1.0)) / pixel_pitch_mm
