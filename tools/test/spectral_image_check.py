"""Checks of tools/lib/hdr_io.read_envi and tools/lib/colorimetry on a render: the RGB picture computed in Python from the spectral image
must be the renderer's own RGB picture.

    out\\cli-dev-spectral\\GoLD_rayt.exe --scene spd --width 320 --spp 1024 --seed 1 --spectral-image 5,360,830 --out x.pfm
    python tools/test/spectral_image_check.py x

360..830 nm is the range the renderer samples wavelengths in, so the spectral image then holds every sample of the RGB picture.

1. With the renderer's colour-matching functions (the analytic fit of Wyman et al. 2013 and its white normalisation) the two pictures
   must agree: what is left is the box approximation of the bands (a sample's weight is the CMF at its wavelength in the renderer, the
   mean over its band here), which grows with the band width, and float32 rounding.
2. With the tabulated CIE 1931 functions (colour-science) the difference is the fit's error: reported, not checked.
3. colour-science's CIEDE2000 against tools/lib/colordiff.delta_e2000 (used by the figure scripts), and the white of D65 at L* = 100.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
import colour                                                          # noqa: E402
from colordiff import delta_e2000                                      # noqa: E402
from colorimetry import spectral_to_xyz, white_xyz, xyz_to_lab, xyz_to_rgb, delta_e   # noqa: E402
from hdr_io import read_envi, read_pfm                                 # noqa: E402


def renderer_cmfs():
    """The renderer's CIE 1931 functions: the multi-lobe fit of Wyman, Sloan and Shirley (2013), Color/ColorSpace.hpp, at 1 nm."""
    def g(x, mu, kl, kr, a):
        d = (x - mu) * np.where(x < mu, kl, kr)
        return a * np.exp(-0.5 * d * d)
    x = np.arange(360.0, 831.0)
    X = g(x, 442.0, 0.0624, 0.0374, 0.362) + g(x, 599.8, 0.0264, 0.0323, 1.056) - g(x, 501.1, 0.0490, 0.0382, 0.065)
    Y = g(x, 568.8, 0.0213, 0.0247, 0.821) + g(x, 530.9, 0.0613, 0.0322, 0.286)
    Z = g(x, 437.0, 0.0845, 0.0278, 1.217) + g(x, 459.0, 0.0385, 0.0725, 0.681)
    return colour.MultiSpectralDistributions(dict(zip(x, np.stack([X, Y, Z], axis=1))), labels=("x", "y", "z"), name="renderer fit")


def renderer_white_rgb(cmfs):
    """The renderer normalises each sRGB channel by its value for D65 (1 nm trapezoid on 380..780 nm, Color/SpectralConversion.hpp)."""
    x = np.arange(380.0, 781.0)
    d65 = colour.SDS_ILLUMINANTS["D65"]
    s = np.interp(x, d65.wavelengths, d65.values)
    f = np.stack([np.interp(x, cmfs.wavelengths, cmfs.values[:, c]) for c in range(3)], axis=1)
    xyz = np.trapezoid(s[:, None] * f, x, axis=0)
    return xyz_to_rgb(xyz / xyz[1], "srgb")


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

    cmfs = renderer_cmfs()
    rgb = xyz_to_rgb(spectral_to_xyz(cube, wl, observer=cmfs), "srgb") / renderer_white_rgb(cmfs)
    mean, med, p99 = differences(rgb, ref)
    step = float(np.median(np.diff(wl)))
    # The box error averages out over the picture but not in a pixel: there it scales with the band width (on the spd scene at 1024 spp
    # the 99th percentile is 0.16 % with 1 nm bands, 0.80 % with 5 nm bands). A flipped or shifted picture fails by far.
    check("renderer CMFs reproduce the .pfm", np.abs(mean).max() < 1e-3 and p99 < 3e-3 * step,
          "image mean %s %%, pixels median %.4f %%, 99th percentile %.4f %% (limit %.2f %% for %g nm bands)"
          % (np.round(100 * mean, 4), 100 * med, 100 * p99, 0.3 * step, step))

    rgb = xyz_to_rgb(spectral_to_xyz(cube, wl), "srgb")
    mean, med, p99 = differences(rgb, ref)
    print("  [info] tabulated CIE 1931 instead of the fit: image mean %s %%, pixels median %.3f %%, 99th percentile %.3f %%"
          % (np.round(100 * mean, 3), 100 * med, 100 * p99))

    rng = np.random.default_rng(1)
    lab1 = np.column_stack([rng.uniform(0, 100, 10000), rng.uniform(-100, 100, 10000), rng.uniform(-100, 100, 10000)])
    lab2 = lab1 + rng.normal(0, 5, lab1.shape)
    diff = np.abs(delta_e(lab1, lab2) - delta_e2000(lab1, lab2)).max()
    check("CIEDE2000 of colour-science = colordiff.delta_e2000", diff < 1e-6, "largest difference %.2e over 10000 pairs" % diff)

    lab = xyz_to_lab(white_xyz("D65"), white_xyz("D65"))
    check("white of D65 against itself", np.allclose(lab, [100, 0, 0], atol=1e-9), "L*a*b* = %s" % np.round(lab, 9))

    print("all checks passed" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1]))
