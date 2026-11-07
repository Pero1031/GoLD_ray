"""RGB build vs spectral build on the gold bowl: how much of the difference does each display gamut show?

Same renders and the same statistics as gold_compare_linear_fig.py (mean of the repeats a and b, 3 x 3 box blur, bowl box, lit pixels),
but the images are "displayed" in sRGB, Display P3 and Rec.2020 (Film::save for --gamut: change of primaries, negatives to 0, Reinhard
per channel, sRGB curve, 8 bit) and the CIEDE2000 is measured on what that display shows (decoded with the gamut's own primaries, D65).

Input:  <dir>/goldbowl_<light>_<rgb|spectral>_<a|b>.pfm  (linear sRGB, written by `GoLD_rayt --out x.pfm`; a and b are independent renders)
Output: the table below, printed and written to <out.txt>.

    python tools/studies/gold_compare_gamut.py <dir> <out.txt> d65 tungsten candle

Columns: dE00 mean (and 95th percentile) between the RGB and the spectral image over the bowl; the share of pixels above 1 and above 2.3 (about
the just-noticeable difference); "outside" = share of the spectral pixels whose color lies outside the gamut (a channel < 0 before the clip);
"noise" = the same dE00 between the two repeats a and b of one build (each is a single render, so the noise of the averaged images that are
compared is about 1/sqrt(2) of it; the value is for the displayed, blurred images). The mean dE00 of the linear images (no clip) is gamut independent.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
from hdr_io import read_hdr                        # noqa: E402
from colordiff import M_RGB_TO_XYZ, delta_e2000    # noqa: E402
from gamut import NAMES, display, displayed_xyz, from_srgb, lab   # noqa: E402
from gold_compare_linear import xyz_to_lab, REGION_BOXES          # noqa: E402
from gold_compare_linear_fig import box_blur                      # noqa: E402

GAMUTS = ("srgb", "p3", "rec2020")


def displayed_de(lin_a, lin_b, gamut):
    """Per-pixel CIEDE2000 between two linear-sRGB images as shown on a display of the gamut."""
    return delta_e2000(lab(displayed_xyz(display(lin_a, gamut), gamut)), lab(displayed_xyz(display(lin_b, gamut), gamut)))


def stats(d):
    return d.mean(), np.percentile(d, 95), (d > 1.0).mean() * 100, (d > 2.3).mean() * 100


def main():
    src, out_txt = sys.argv[1], sys.argv[2]
    lights = sys.argv[3:]
    lines = ["RGB build vs spectral build, gold bowl (800x450, 256 spp, mean of two renders per build, 3x3 box blur; bowl box, lit pixels).",
             "CIEDE2000 measured on what a display of each gamut shows (see the header of tools/studies/gold_compare_gamut.py).", ""]
    head = (f"{'light':9s} {'display':11s} {'dE00 mean':>9s} {'95 pct':>7s} {'>1':>6s} {'>2.3':>6s}   {'outside':>8s}   {'noise RGB':>9s} {'noise spec':>10s}")
    for light in lights:
        rep = {(m, r): read_hdr(os.path.join(src, f"goldbowl_{light}_{m}_{r}.pfm")).astype(np.float64) for m in ("rgb", "spectral") for r in "ab"}
        rgb = box_blur((rep[("rgb", "a")] + rep[("rgb", "b")]) / 2)
        spe = box_blur((rep[("spectral", "a")] + rep[("spectral", "b")]) / 2)
        blur = {k: box_blur(v) for k, v in rep.items()}
        h, w = rgb.shape[:2]

        def box(name):
            y0, y1, x0, x1 = REGION_BOXES[name]
            return slice(int(y0 / 450 * h), int(y1 / 450 * h)), slice(int(x0 / 800 * w), int(x1 / 800 * w))

        wy = (M_RGB_TO_XYZ @ rgb[box("white ball")].reshape(-1, 3).mean(0))[1]     # luminance of the white reference ball
        lit = (0.2126 * spe[..., 0] + 0.7152 * spe[..., 1] + 0.0722 * spe[..., 2]) > 0.02 * wy
        mask = np.zeros_like(lit)
        mask[box("bowl total")] = True
        mask &= lit

        white = (M_RGB_TO_XYZ @ np.ones(3)) * wy
        de_lin = delta_e2000(xyz_to_lab(rgb @ M_RGB_TO_XYZ.T, white), xyz_to_lab(spe @ M_RGB_TO_XYZ.T, white))
        lines += [f"=== {light} ({w}x{h}, {int(mask.sum())} pixels) ===", head]
        lines.append(f"{light:9s} {'linear':11s} {de_lin[mask].mean():9.2f} {np.percentile(de_lin[mask], 95):7.2f} {'':>6s} {'':>6s}   {'':>8s}   (no clip, white = D65 at the luminance of the white ball)")
        for g in GAMUTS:
            d = displayed_de(rgb, spe, g)[mask]
            mean, p95, over1, over23 = stats(d)
            outside = ((spe @ from_srgb(g).T).min(axis=2) < 0)[mask].mean() * 100
            n_rgb = displayed_de(blur[("rgb", "a")], blur[("rgb", "b")], g)[mask].mean()
            n_spe = displayed_de(blur[("spectral", "a")], blur[("spectral", "b")], g)[mask].mean()
            lines.append(f"{light:9s} {NAMES[g]:11s} {mean:9.2f} {p95:7.2f} {over1:5.1f}% {over23:5.1f}%   {outside:7.1f}%   {n_rgb:9.2f} {n_spe:10.2f}")
        lines.append("")
    text = "\n".join(lines)
    print(text)
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()
