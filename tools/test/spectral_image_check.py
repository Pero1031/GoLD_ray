"""Checks of tools/lib/hdr_io.read_envi and tools/lib/colorimetry on a render: the RGB picture computed in Python from the spectral image
must be the renderer's own RGB picture.

    out\\cli-dev-spectral\\GoLD_rayt.exe --scene spd --width 320 --spp 1024 --seed 1 --spectral-image 5,360,830 --out x.pfm
    python tools/test/spectral_image_check.py x

360..830 nm is the range the renderer samples wavelengths in, so the spectral image then holds every sample of the RGB picture.

1. With colorimetry.py's defaults (the tabulated CIE 1931 functions, which the renderer's film uses too, and sRGB) the two pictures
   must agree: what is left is the box approximation of the bands (a sample's weight is the CMF at its wavelength in the renderer, the
   mean over its band here), which grows with the band width, and float32 rounding.
2. colour-science's CIEDE2000 against tools/lib/colordiff.delta_e2000 (used by the figure scripts), the white of D65 at L* = 100 and
   at RGB (1, 1, 1).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
from colordiff import delta_e2000                                      # noqa: E402
from colorimetry import spectral_to_xyz, white_xyz, xyz_to_lab, xyz_to_rgb, delta_e   # noqa: E402
from hdr_io import read_envi, read_pfm                                 # noqa: E402


def differences(rgb, ref):
    """Image-mean difference per channel and the per-pixel difference relative to the pixel's brightness (pixels above 5 % of white)."""
    mean = rgb.reshape(-1, 3).mean(0) / ref.reshape(-1, 3).mean(0) - 1
    lum = ref.mean(-1)
    m = lum > 0.05
    rel = (np.abs(rgb - ref)[m] / lum[m][:, None]).ravel()
    return mean, np.median(rel), np.percentile(rel, 99)


def main(stem):
    ok = True

    def check(name, good, detail):
        nonlocal ok
        ok &= bool(good)
        print("  [%s] %s: %s" % ("PASS" if good else "FAIL", name, detail))

    cube, wl = read_envi(stem)
    ref = read_pfm(os.path.splitext(stem)[0] + ".pfm")
    print("%s: %d x %d, %d bands, %.1f .. %.1f nm" % (stem, cube.shape[1], cube.shape[0], cube.shape[2], wl[0], wl[-1]))
    if wl[0] - 0.5 * (wl[1] - wl[0]) > 360.5 or wl[-1] + 0.5 * (wl[1] - wl[0]) < 829.5:
        print("  note: the bands do not cover 360..830 nm, so part of the RGB picture is missing from them")

    rgb = xyz_to_rgb(spectral_to_xyz(cube, wl), "srgb")
    mean, med, p99 = differences(rgb, ref)
    step = float(np.median(np.diff(wl)))
    # The box error averages out over the picture but not in a pixel: there it scales with the band width (on the spd scene at 1024 spp
    # the 99th percentile is 0.16 % with 1 nm bands, 0.80 % with 5 nm bands). A flipped or shifted picture fails by far.
    check("Python's XYZ and sRGB reproduce the .pfm", np.abs(mean).max() < 1e-3 and p99 < 3e-3 * step,
          "image mean %s %%, pixels median %.4f %%, 99th percentile %.4f %% (limit %.2f %% for %g nm bands)"
          % (np.round(100 * mean, 4), 100 * med, 100 * p99, 0.3 * step, step))

    rng = np.random.default_rng(1)
    lab1 = np.column_stack([rng.uniform(0, 100, 10000), rng.uniform(-100, 100, 10000), rng.uniform(-100, 100, 10000)])
    lab2 = lab1 + rng.normal(0, 5, lab1.shape)
    diff = np.abs(delta_e(lab1, lab2) - delta_e2000(lab1, lab2)).max()
    check("CIEDE2000 of colour-science = colordiff.delta_e2000", diff < 1e-6, "largest difference %.2e over 10000 pairs" % diff)

    lab = xyz_to_lab(white_xyz("D65"), white_xyz("D65"))
    check("white of D65 against itself", np.allclose(lab, [100, 0, 0], atol=1e-9), "L*a*b* = %s" % np.round(lab, 9))
    rgb = xyz_to_rgb(white_xyz("D65"), "srgb")
    check("white of D65 in sRGB", np.abs(rgb - 1).max() < 1e-3, "RGB = %s" % np.round(rgb, 5))

    print("all checks passed" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1]))
