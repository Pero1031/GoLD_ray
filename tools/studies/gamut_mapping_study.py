"""
How much of the "sRGB hides the RGB-vs-spectral difference" result depends on clipping as the gamut mapping?

The recorded result (docs/METAL_ACCURACY.md 3.5) clips negative channels to 0 on the way to the display. A real pipeline may map the
out-of-gamut colors instead. This script re-measures the displayed CIEDE2000 between the RGB and the spectral renders of the gold bowl
(and of reference-model patches) with several gamut mappings, for an sRGB and a Display P3 display:

  clip        per-channel clip to 0 (what the renderer writes, and the earlier result)
  gray        desaturate toward the pixel's own luminance (linear light) until the lowest channel is 0: hue (in linear RGB) and luminance kept
  white       add the same amount to all channels until the lowest channel is 0 (raises luminance; a "veiling white")
  oklch       reduce the OKLCH chroma at constant lightness and hue until the color is inside the gamut (the idea of the CSS Color 4 gamut mapping)

Only the lower bound is mapped (negative channels). Values above 1 are handled by the Reinhard tone curve, as in the renderer's output.
Pipeline per pixel: linear sRGB -> display primaries -> mapping -> Reinhard per channel -> sRGB curve -> 8 bit -> CIEDE2000.

  python tools/studies/gamut_mapping_study.py <dir with goldbowl_<light>_<rgb|spectral>_<a|b>.pfm> <out.txt>
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
from hdr_io import read_hdr                                   # noqa: E402
from colordiff import M_RGB_TO_XYZ, delta_e2000               # noqa: E402
from gamut import NAMES, from_srgb, lab, displayed_xyz, oetf, rgb_to_xyz   # noqa: E402
from gold_compare_linear import xyz_to_lab, REGION_BOXES      # noqa: E402
from gold_compare_linear_fig import box_blur                  # noqa: E402

LIGHTS = ("d65", "tungsten", "candle")
DISPLAYS = ("srgb", "p3")
METHODS = ("clip", "gray", "white", "oklch")

# OKLab (Ottosson 2020)
M1 = np.array([[0.8189330101, 0.3618667424, -0.1288597137],
               [0.0329845436, 0.9293118715, 0.0361456387],
               [0.0482003018, 0.2643662691, 0.6338517070]])
M2 = np.array([[0.2104542553, 0.7936177850, -0.0040720468],
               [1.9779984951, -2.4285922050, 0.4505937099],
               [0.0259040371, 0.7827717662, -0.8086757660]])
M1_INV, M2_INV = np.linalg.inv(M1), np.linalg.inv(M2)


def map_gamut(x, gamut, method):
    """x: (..., 3) linear RGB in the display primaries (may have negative channels). Returns values with no channel below 0."""
    if method == "clip":
        return np.clip(x, 0.0, None)
    xmin = x.min(axis=-1)
    bad = xmin < 0
    out = x.copy()
    if not bad.any():
        return out
    xb = x[bad]
    mn = xmin[bad]
    if method == "white":
        out[bad] = xb - mn[:, None]
    elif method == "gray":
        wy = rgb_to_xyz(gamut)[1]                       # luminance weights of the display primaries
        y = xb @ wy
        t = np.where(y > 0, y / np.maximum(y - mn, 1e-12), 0.0)
        out[bad] = y[:, None] + t[:, None] * (xb - y[:, None])
    elif method == "oklch":
        to_xyz = rgb_to_xyz(gamut)
        lab_ok = (np.cbrt((xb @ to_xyz.T) @ M1.T)) @ M2.T          # (N, 3): L, a, b
        lo, hi = np.zeros(len(xb)), np.ones(len(xb))               # chroma scale: 1 = original (outside), 0 = gray
        to_rgb = np.linalg.inv(to_xyz)
        def rgb_at(s):
            c = np.stack([lab_ok[:, 0], lab_ok[:, 1] * s, lab_ok[:, 2] * s], axis=-1)
            return (((c @ M2_INV.T) ** 3) @ M1_INV.T) @ to_rgb.T
        for _ in range(28):                                         # the largest chroma scale that keeps every channel >= 0
            mid = 0.5 * (lo + hi)
            ok = rgb_at(mid).min(axis=-1) >= -1e-7
            lo = np.where(ok, mid, lo)
            hi = np.where(ok, hi, mid)
        out[bad] = np.clip(rgb_at(lo), 0.0, None)
    else:
        raise ValueError(method)
    return out


def display8(lin_srgb, gamut, method):
    x = map_gamut(lin_srgb @ from_srgb(gamut).T, gamut, method)
    x = x / (1 + x)                                                 # Reinhard per channel, exposure 0
    return np.clip(np.round(255 * oetf(x)), 0, 255).astype(np.uint8)


def de_displayed(a, b, gamut, method):
    return delta_e2000(lab(displayed_xyz(display8(a, gamut, method), gamut)), lab(displayed_xyz(display8(b, gamut, method), gamut)))


def main():
    src, out_txt = sys.argv[1], sys.argv[2]
    lines = ["Gamut-mapping sensitivity of the displayed RGB-vs-spectral difference (gold bowl, 800x450, 256 spp, mean of two renders per build, 3x3 blur; bowl box, lit pixels).",
             "CIEDE2000 between the displayed RGB and spectral images. Methods: see tools/studies/gamut_mapping_study.py.", ""]
    for light in LIGHTS:
        rep = {(m, r): read_hdr(os.path.join(src, f"goldbowl_{light}_{m}_{r}.pfm")).astype(np.float64) for m in ("rgb", "spectral") for r in "ab"}
        rgb = box_blur((rep[("rgb", "a")] + rep[("rgb", "b")]) / 2)
        spe = box_blur((rep[("spectral", "a")] + rep[("spectral", "b")]) / 2)
        blur = {k: box_blur(v) for k, v in rep.items()}
        h, w = rgb.shape[:2]

        def box(name):
            y0, y1, x0, x1 = REGION_BOXES[name]
            return slice(int(y0 / 450 * h), int(y1 / 450 * h)), slice(int(x0 / 800 * w), int(x1 / 800 * w))

        wy = (M_RGB_TO_XYZ @ rgb[box("white ball")].reshape(-1, 3).mean(0))[1]
        lit = (0.2126 * spe[..., 0] + 0.7152 * spe[..., 1] + 0.0722 * spe[..., 2]) > 0.02 * wy
        mask = np.zeros_like(lit)
        mask[box("bowl total")] = True
        mask &= lit
        lines.append(f"=== {light} ({int(mask.sum())} pixels) ===")
        lines.append(f"{'display':11s} {'mapping':8s} {'dE00 mean':>9s} {'95 pct':>7s} {'>1':>6s}   {'outside':>8s}   {'noise RGB':>9s} {'noise spec':>10s}")
        for g in DISPLAYS:
            outside = ((spe @ from_srgb(g).T).min(axis=2) < 0)[mask].mean() * 100
            for method in METHODS:
                d = de_displayed(rgb, spe, g, method)[mask]
                n_rgb = de_displayed(blur[("rgb", "a")], blur[("rgb", "b")], g, method)[mask].mean()
                n_spe = de_displayed(blur[("spectral", "a")], blur[("spectral", "b")], g, method)[mask].mean()
                lines.append(f"{NAMES[g]:11s} {method:8s} {d.mean():9.2f} {np.percentile(d, 95):7.2f} {(d > 1.0).mean() * 100:5.1f}%   {outside:7.1f}%   {n_rgb:9.2f} {n_spe:10.2f}")
        lines.append("")
    text = "\n".join(lines)
    print(text)
    os.makedirs(os.path.dirname(os.path.abspath(out_txt)), exist_ok=True)
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()
