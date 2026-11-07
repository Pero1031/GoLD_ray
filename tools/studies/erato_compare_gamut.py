"""RGB build vs spectral build on the gold statue (Erato): how much of the difference does each display gamut show?

The statue counterpart of gold_compare_gamut.py: per-pixel CIEDE2000 between the RGB render and the spectral render, measured on what a display of
each gamut shows (Film::save for --gamut: change of primaries, negatives to 0, Reinhard per channel, sRGB curve, 8 bit). Pixels and box as in
tools/figures/compare_figure.py: the middle of the picture, lit pixels only (a displayed sRGB channel > 12 in either image; the same pixels for
every gamut). Two parts: A the denoised renders (single, no blur), B the raw renders (mean of a and b, 3x3 box blur). The "noise" columns (the
same dE00 between the repeats a and b of one build) are the floor of the measurement.

THE DENOISER CLIPS NEGATIVE VALUES to 0, so part A can not show a color outside sRGB (the candle's negative blue is gone); part B is the one to use.

Input:  <dir>/erato_<light>_<rgb|spectral>_<a|b>.pfm and ..._denoised.pfm   (`GoLD_rayt scenes/erato_gold.json --denoise --out <dir>/erato_<light>_<build>_<rep>.pfm`,
        the light changed in the scene file: out/gamut_study/erato_<light>.json were made by replacing "tungsten")

    python tools/studies/erato_compare_gamut.py <dir> <out.txt> d65 tungsten candle
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
from hdr_io import read_hdr                                    # noqa: E402
from gamut import NAMES, display, from_srgb                    # noqa: E402
from gold_compare_gamut import GAMUTS, displayed_de, stats     # noqa: E402

BOX = (0.15, 0.85, 0.05, 0.98)   # x0 x1 y0 y1 as fractions, like compare_figure.py


def box_blur(im):
    p = np.pad(im, ((1, 1), (1, 1), (0, 0)), mode="edge")
    return sum(p[dy:dy + im.shape[0], dx:dx + im.shape[1]] for dy in range(3) for dx in range(3)) / 9.0


def analyze(src, light, mode, lines):
    """mode 'denoised': single denoised renders, no blur.  mode 'raw': the raw renders, mean of a and b, 3x3 box blur (noise floor from a vs b blurred)."""
    head = (f"{'light':9s} {'display':11s} {'dE00 mean':>9s} {'median':>7s} {'95 pct':>7s} {'>1':>6s} {'>2.3':>6s}   {'outside':>8s}   {'noise RGB':>9s} {'noise spec':>10s}")
    suffix = "_denoised" if mode == "denoised" else ""
    im = {(m, r): read_hdr(os.path.join(src, f"erato_{light}_{m}_{r}{suffix}.pfm")).astype(np.float64) for m in ("rgb", "spectral") for r in "ab"}
    if mode == "raw":
        a = {k: box_blur(v) for k, v in im.items()}                                   # for the noise floor
        rgb, spe = box_blur((im[("rgb", "a")] + im[("rgb", "b")]) / 2), box_blur((im[("spectral", "a")] + im[("spectral", "b")]) / 2)
    else:
        a = im
        rgb, spe = im[("rgb", "a")], im[("spectral", "a")]
    h, w = rgb.shape[:2]
    lit = (display(rgb, "srgb").max(axis=2) > 12) | (display(spe, "srgb").max(axis=2) > 12)
    mask = np.zeros_like(lit)
    mask[int(BOX[2] * h):int(BOX[3] * h), int(BOX[0] * w):int(BOX[1] * w)] = True
    mask &= lit
    lines += [f"=== {light} ({w}x{h}, {int(mask.sum())} pixels) ===", head]
    for g in GAMUTS:
        d = displayed_de(rgb, spe, g)[mask]
        mean, p95, over1, over23 = stats(d)
        outside = ((spe @ from_srgb(g).T).min(axis=2) < 0)[mask].mean() * 100
        n_rgb = displayed_de(a[("rgb", "a")], a[("rgb", "b")], g)[mask].mean()
        n_spe = displayed_de(a[("spectral", "a")], a[("spectral", "b")], g)[mask].mean()
        lines.append(f"{light:9s} {NAMES[g]:11s} {mean:9.2f} {np.median(d):7.2f} {p95:7.2f} {over1:5.1f}% {over23:5.1f}%   {outside:7.1f}%   {n_rgb:9.2f} {n_spe:10.2f}")
    lines.append("")


def main():
    src, out_txt = sys.argv[1], sys.argv[2]
    lights = sys.argv[3:]
    lines = ["RGB build vs spectral build, gold statue (Erato, 600x800, 256 spp; middle box, lit pixels).",
             "CIEDE2000 measured on what a display of each gamut shows (see the header of tools/studies/erato_compare_gamut.py).",
             "",
             "##### A. denoised renders (single renders, no blur) #####",
             "NOTE: the denoiser (OIDN) clips negative values to 0. A color outside sRGB (the candle: blue < 0) is therefore already clipped in",
             "these images, and a wider display can not show it. Use part B for anything that depends on out-of-sRGB colors.", ""]
    for light in lights:
        analyze(src, light, "denoised", lines)
    lines += ["##### B. raw renders (no denoiser; mean of the repeats a and b, 3x3 box blur; noise = a vs b blurred) #####", ""]
    for light in lights:
        analyze(src, light, "raw", lines)
    text = "\n".join(lines)
    print(text)
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()
