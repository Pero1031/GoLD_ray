"""RGB build vs spectral build on the gold statue, with S-CIELAB: does the difference survive the eye's spatial filtering?

The per-pixel CIEDE2000 map of the statue shows large differences in the deep folds (many reflections) and small ones on the flat parts, yet the two
renders look different as a whole (the spectral one more vivid). S-CIELAB (tools/lib/scielab.py) blurs the color planes of both images as the visual
system does before taking the difference, so a fine pattern of colors blends into its surround. This script compares the plain per-pixel dE00 with the
S-CIELAB dE00 (same images, same mask), for several viewing distances, and with the difference of the object's mean color.

Input:  <dir>/erato_<light>_<rgb|spectral>_<a|b>.pfm  (raw renders, not denoised: see tools/studies/erato_compare_gamut.py; mean of a and b, no blur;
        the plain dE00 is shown on the 3x3-blurred means like the other results, and also unblurred)
Output: the table (printed and written to <out.txt>) and, with --fig DIR, a figure per light and gamut.

    python tools/studies/erato_scielab.py <dir> <out.txt> d65 tungsten candle [--fig DIR]

"Display" is as in Film::save (change of primaries, negatives to 0, Reinhard per channel, sRGB curve, 8 bit); the gamuts shown are sRGB and Display P3.
Viewing distance: pixels per degree 20 / 40 / 80 (40 = a 0.25 mm pixel pitch seen from 60 cm; 20 = the image seen from half the distance, or on a coarser screen).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
from hdr_io import read_hdr                                          # noqa: E402
from colordiff import delta_e2000                                    # noqa: E402
from gamut import NAMES, display, displayed_xyz, lab                 # noqa: E402
from scielab import spatial_filter_xyz                               # noqa: E402
from erato_compare_gamut import BOX, box_blur                        # noqa: E402

PPDS = (20, 40, 80)
GAMUTS = ("srgb", "p3")
VMAX = 10.0


def stats(d):
    return d.mean(), np.median(d), np.percentile(d, 95), (d > 2.3).mean() * 100


def de_maps(lin_a, lin_b, gamut, ppd=None, blur=False):
    """dE00 map between two linear-sRGB images as displayed in `gamut`; with ppd: after the S-CIELAB filtering."""
    if blur:
        lin_a, lin_b = box_blur(lin_a), box_blur(lin_b)
    xa, xb = displayed_xyz(display(lin_a, gamut), gamut), displayed_xyz(display(lin_b, gamut), gamut)
    if ppd:
        xa, xb = spatial_filter_xyz(xa, ppd), spatial_filter_xyz(xb, ppd)
    return delta_e2000(lab(xa), lab(xb))


def main():
    args = sys.argv[1:]
    fig_dir = None
    if "--fig" in args:
        i = args.index("--fig")
        fig_dir = args[i + 1]
        del args[i:i + 2]
    src, out_txt = args[:2]
    lights = args[2:]
    lines = ["RGB build vs spectral build, gold statue (Erato, 600x800, 256 spp, raw renders, mean of two; middle box, lit pixels).",
             "dE00 between the two images as displayed; plain = per pixel, S-CIELAB = after the spatial filtering of the visual system at the given pixels per degree.",
             "'noise' = the same measure between the repeats a and b of one build (a single render each). 'object mean' = dE00 between the mean colors of the masked pixels.",
             "The S-CIELAB filter widths are written from the published table and not checked against a reference implementation (see tools/lib/scielab.py).", ""]
    head = f"{'measure':22s} {'mean':>6s} {'median':>7s} {'95 pct':>7s} {'>2.3':>6s}   {'noise RGB':>9s} {'noise spec':>10s}"
    for light in lights:
        im = {(m, r): read_hdr(os.path.join(src, f"erato_{light}_{m}_{r}.pfm")).astype(np.float64) for m in ("rgb", "spectral") for r in "ab"}
        rgb, spe = (im[("rgb", "a")] + im[("rgb", "b")]) / 2, (im[("spectral", "a")] + im[("spectral", "b")]) / 2
        h, w = rgb.shape[:2]
        lit = (display(box_blur(rgb), "srgb").max(axis=2) > 12) | (display(box_blur(spe), "srgb").max(axis=2) > 12)
        mask = np.zeros_like(lit)
        mask[int(BOX[2] * h):int(BOX[3] * h), int(BOX[0] * w):int(BOX[1] * w)] = True
        mask &= lit
        for g in GAMUTS:
            lines += [f"=== {light}, shown in {NAMES[g]} ({int(mask.sum())} pixels) ===", head]
            rows = [("plain, 3x3 blur", dict(blur=True)), ("plain, no blur", dict())] + [(f"S-CIELAB {p} px/deg", dict(ppd=p)) for p in PPDS]
            maps = {}
            for name, kw in rows:
                d = de_maps(rgb, spe, g, **kw)
                maps[name] = d
                n_rgb = de_maps(im[("rgb", "a")], im[("rgb", "b")], g, **kw)[mask].mean()
                n_spe = de_maps(im[("spectral", "a")], im[("spectral", "b")], g, **kw)[mask].mean()
                m, med, p95, over = stats(d[mask])
                lines.append(f"{name:22s} {m:6.2f} {med:7.2f} {p95:7.2f} {over:5.1f}%   {n_rgb:9.2f} {n_spe:10.2f}")
            mean_xyz = lambda x: displayed_xyz(display(x, g), g)[mask].mean(axis=0)   # noqa: E731
            d_obj = float(delta_e2000(lab(mean_xyz(rgb)), lab(mean_xyz(spe))))
            lines += [f"{'object mean':22s} {d_obj:6.2f}", ""]
            if fig_dir:
                save_figure(os.path.join(fig_dir, f"erato_scielab_{light}_{g}.png"), rgb, spe, g, maps["plain, 3x3 blur"], maps["S-CIELAB 40 px/deg"], mask)
    text = "\n".join(lines)
    print(text)
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(text + "\n")


def save_figure(path, rgb, spe, gamut, plain, scie, mask):
    from PIL import Image
    from matplotlib import colormaps
    os.makedirs(os.path.dirname(path), exist_ok=True)
    cmap = colormaps["magma"]

    def heat(d):
        out = (cmap(np.clip(d / VMAX, 0, 1))[..., :3] * 255).astype(np.uint8)
        out[~mask] = (out[~mask] * 0.35).astype(np.uint8)
        return out
    # the images as the display shows them (8 bit; for P3 the file would be tagged, here only the numbers matter)
    panels = [display(box_blur(rgb), gamut), display(box_blur(spe), gamut), heat(plain), heat(scie)]
    Image.fromarray(np.concatenate(panels, axis=1)).save(path)


if __name__ == "__main__":
    main()
