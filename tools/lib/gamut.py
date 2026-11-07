"""Output gamuts (sRGB, Display P3, Rec.2020) for the analysis scripts: matrices from the primaries, and what a display of each gamut shows.

All three share the D65 white, so a change of primaries is one 3x3 matrix on linear values (see GoLD_rayt/include/Color/Gamut.hpp, which
holds the same matrices; tools/test/gamut_check.py compares them).
"""
import numpy as np

D65_XY = (0.3127, 0.3290)
PRIMARIES = {
    "srgb":    [(0.640, 0.330), (0.300, 0.600), (0.150, 0.060)],
    "p3":      [(0.680, 0.320), (0.265, 0.690), (0.150, 0.060)],
    "rec2020": [(0.708, 0.292), (0.170, 0.797), (0.131, 0.046)],
}
NAMES = {"srgb": "sRGB", "p3": "Display P3", "rec2020": "Rec.2020"}


def rgb_to_xyz(name):
    """Linear RGB of the gamut -> XYZ (white (1,1,1) = D65 with Y = 1)."""
    p = np.array([[x / y, 1.0, (1 - x - y) / y] for x, y in PRIMARIES[name]]).T
    w = np.array([D65_XY[0] / D65_XY[1], 1.0, (1 - D65_XY[0] - D65_XY[1]) / D65_XY[1]])
    return p * np.linalg.solve(p, w)


def from_srgb(name):
    """3x3 matrix: linear sRGB (the film's values) -> linear RGB of the gamut."""
    return np.linalg.inv(rgb_to_xyz(name)) @ rgb_to_xyz("srgb")


def oetf(x):
    """The sRGB curve (also used for Display P3 and, by the renderer, for Rec.2020 stills)."""
    return np.where(x <= 0.0031308, 12.92 * x, 1.055 * np.power(np.maximum(x, 0), 1 / 2.4) - 0.055)


def eotf(v):
    return np.where(v <= 0.04045, v / 12.92, np.power((v + 0.055) / 1.055, 2.4))


def display(lin_srgb, name):
    """What Film::save writes for --gamut <name>: change of primaries, negatives to 0, exposure 0, Reinhard per channel, curve, 8 bit."""
    x = np.clip(lin_srgb @ from_srgb(name).T, 0, None)
    x = x / (1 + x)
    return np.clip(np.round(255 * oetf(x)), 0, 255).astype(np.uint8)


def displayed_xyz(img8, name):
    """XYZ (D65 white at Y = 1) of an 8-bit image shown on a display of the gamut."""
    return eotf(img8.astype(np.float64) / 255.0) @ rgb_to_xyz(name).T


def lab(xyz):
    white = rgb_to_xyz("srgb") @ np.ones(3)
    t = xyz / white
    f = np.where(t > (6 / 29) ** 3, np.cbrt(t), t / (3 * (6 / 29) ** 2) + 4 / 29)
    return np.stack([116 * f[..., 1] - 16, 500 * (f[..., 0] - f[..., 1]), 200 * (f[..., 1] - f[..., 2])], axis=-1)
