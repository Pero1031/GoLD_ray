"""Colour science on spectral images (the renderer's --spectral-image, or any (H, W, bands) radiance cube), on top of colour-science
(pip install colour-science).

    from hdr_io import read_envi
    from colorimetry import spectral_to_xyz, white_xyz, xyz_to_lab, delta_e, xyz_to_rgb, to_sd

    cube, wl = read_envi("out/x")                  # out/x_spectral.hdr + .img: (H, W, bands), band centres in nm
    xyz = spectral_to_xyz(cube, wl)                # CIE 1931 2 deg; Y = 1 for a white surface under the renderer's D65 white
    lab = xyz_to_lab(xyz, white_xyz("D65"))        # or the XYZ of a white patch in the picture: xyz[mask].mean(axis=0)
    de = delta_e(lab, lab_other)                   # CIEDE2000 (other methods of colour.delta_E by name)
    rgb = xyz_to_rgb(xyz, "Display P3")            # linear RGB of a display gamut, no adaptation
    sd = to_sd(cube[120, 80], wl)                  # a colour.SpectralDistribution: CCT, CRI, CAM16, ... from colour-science itself

Bands are boxes: a band's value is the mean radiance over [centre - width/2, centre + width/2] (what the renderer writes), so each band
is weighted by the integral of the colour-matching functions over its box, not by their value at the centre. The width is the band
spacing unless given.

The renderer's own RGB uses an analytic fit of the CIE 1931 functions (Wyman et al. 2013); the tabulated functions used here differ from
it by about one percent in places (tools/test/spectral_image_check.py measures it on a render).
"""
import colour
import numpy as np

CIE1931 = "CIE 1931 2 Degree Standard Observer"
CIE1964 = "CIE 1964 10 Degree Standard Observer"

# The renderer's --gamut names and colour-science's colourspace names.
GAMUTS = {"srgb": "sRGB", "p3": "Display P3", "rec2020": "ITU-R BT.2020"}


def _table(observer):
    """(wavelengths, values (N, 3)) of colour-matching functions given by name (colour.MSDS_CMFS) or as MultiSpectralDistributions."""
    cmfs = colour.MSDS_CMFS[observer] if isinstance(observer, str) else observer
    return np.asarray(cmfs.wavelengths, dtype=np.float64), np.asarray(cmfs.values, dtype=np.float64)


def _sd_table(sd):
    """(wavelengths, values) of a spectrum given by illuminant name (colour.SDS_ILLUMINANTS), as a SpectralDistribution or as (wl, v)."""
    if isinstance(sd, str):
        sd = colour.SDS_ILLUMINANTS[sd]
    if isinstance(sd, tuple):
        return np.asarray(sd[0], dtype=np.float64), np.asarray(sd[1], dtype=np.float64)
    return np.asarray(sd.wavelengths, dtype=np.float64), np.asarray(sd.values, dtype=np.float64)


def _cumulative(x, f, at):
    """The integral of the piecewise linear table (x, f) from x[0] up to each point of `at` (constant outside the table). f: (N,) or (N, m)."""
    f = f.reshape(len(x), -1)
    seg = 0.5 * (f[1:] + f[:-1]) * np.diff(x)[:, None]
    F = np.vstack([np.zeros((1, f.shape[1])), np.cumsum(seg, axis=0)])
    at = np.clip(np.asarray(at, dtype=np.float64), x[0], x[-1])
    i = np.clip(np.searchsorted(x, at, side="right") - 1, 0, len(x) - 2)
    t = (at - x[i])[:, None]
    fa = f[i] + (f[i + 1] - f[i]) * t / (x[i + 1] - x[i])[:, None]
    return F[i] + 0.5 * (f[i] + fa) * t


def band_cmfs(wl, width=None, observer=CIE1931):
    """(bands, 3): the colour-matching functions integrated over each band's box (nm x CMF), zero outside their table."""
    wl = np.asarray(wl, dtype=np.float64)
    if width is None:
        if len(wl) < 2:
            raise ValueError("one band: give its width")
        width = np.full(len(wl), np.median(np.diff(wl)))
    width = np.broadcast_to(np.asarray(width, dtype=np.float64), wl.shape)
    x, f = _table(observer)
    return _cumulative(x, f, wl + 0.5 * width) - _cumulative(x, f, wl - 0.5 * width)


def _white_y(illuminant, observer):
    """Integral of the illuminant times y-bar over the whole table of the observer."""
    x, f = _table(observer)
    sx, sv = _sd_table(illuminant)
    s = np.interp(x, sx, sv, left=0.0, right=0.0)
    return float(np.trapezoid(s * f[:, 1], x))


def spectral_to_xyz(cube, wl, width=None, observer=CIE1931, k=None):
    """XYZ (..., 3) of spectral radiance (..., bands) in boxes centred on `wl`.

    k = None: the renderer's units, where a surface the RGB picture shows as white under D65 has the D65 spectrum (100 at 560 nm); then
    such a white has Y = 1. Otherwise XYZ = k x integral (e.g. k = 683 lm/W turns W sr^-1 m^-2 nm^-1 into cd/m^2)."""
    w = band_cmfs(wl, width, observer)
    if k is None:
        k = 1.0 / _white_y("D65", observer)
    cube = np.asarray(cube)
    xyz = cube @ w.astype(cube.dtype if cube.dtype == np.float32 else np.float64)   # no float64 copy of a large float32 cube
    return xyz.astype(np.float64) * k


def white_xyz(illuminant="D65", observer=CIE1931):
    """XYZ of a perfect white diffuser under an illuminant (name in colour.SDS_ILLUMINANTS, a SpectralDistribution or (wl, values)), Y = 1."""
    x, f = _table(observer)
    sx, sv = _sd_table(illuminant)
    s = np.interp(x, sx, sv, left=0.0, right=0.0)
    xyz = np.trapezoid(s[:, None] * f, x, axis=0)
    return xyz / xyz[1]


def xyz_to_lab(xyz, white):
    """CIELAB of XYZ against a reference white XYZ (any scale: the white itself gets L* = 100)."""
    white = np.asarray(white, dtype=np.float64)
    return colour.XYZ_to_Lab(np.asarray(xyz) / white[1], colour.XYZ_to_xy(white))


def delta_e(lab1, lab2, method="CIE 2000"):
    """Colour difference of two CIELAB arrays (colour.delta_E: "CIE 2000", "CIE 1994", "CIE 1976", "CMC")."""
    return colour.delta_E(lab1, lab2, method=method)


def xyz_to_rgb(xyz, space="sRGB"):
    """Linear RGB of a display gamut (colour-science name, or the renderer's srgb / p3 / rec2020); XYZ relative to its white, which
    goes to (1, 1, 1). No chromatic adaptation: a D65 white stays (1, 1, 1) in all three gamuts the renderer writes."""
    space = GAMUTS.get(space, space)
    return colour.XYZ_to_RGB(xyz, space, apply_cctf_encoding=False)


def to_sd(values, wl, name="pixel"):
    """A colour.SpectralDistribution of one spectrum (e.g. a pixel or the mean of a region), for the rest of colour-science."""
    return colour.SpectralDistribution(dict(zip(np.asarray(wl, dtype=np.float64), np.asarray(values, dtype=np.float64))), name=name)
