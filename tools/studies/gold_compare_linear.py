"""RGB build vs spectral build on the goldbowl scene, measured on LINEAR .hdr renders (no tone mapping, no clipping, no denoiser).

Input: <dir>/goldbowl_<light>_<rgb|spectral>_<a|b>.pfm  (32-bit float; written by `GoLD_rayt --out x.pfm`. .hdr cannot be used: RGBE drops
       negative values and quantizes the small channels of strongly colored pixels).
       a and b are independent renders of the same build: the difference a-b is the noise floor of that build.

    python tools/studies/gold_compare_linear.py <dir> <out.txt> d65 tungsten candle

Colour difference: CIEDE2000 between the region-mean colours (noise averages out over thousands of pixels), under two white points:
  "D65"      : white = the D65 white of the display, with the luminance of the white reference ball (the neutral Lambertian ball)
  "adapted"  : white = the light's own white (what the white ball reflects), i.e. the eye / camera is adapted to the light
Also reported: the ratio spectral/RGB of the linear R, G, B means, and how much of the spectral image is out of the sRGB gamut (negative channel).
"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
from hdr_io import read_hdr                       # noqa: E402
from colordiff import M_RGB_TO_XYZ, delta_e2000   # noqa: E402


def xyz_to_lab(xyz, white):
    t = xyz / white
    f = np.where(t > (6 / 29) ** 3, np.cbrt(t), t / (3 * (6 / 29) ** 2) + 4 / 29)   # defined for negative t as well
    return np.stack([116 * f[..., 1] - 16, 500 * (f[..., 0] - f[..., 1]), 200 * (f[..., 1] - f[..., 2])], axis=-1)


def regions(h, w):
    """Boxes (y0, y1, x0, x1) for an 800 x 450 image, scaled."""
    sy, sx = h / 450.0, w / 800.0
    r = {}
    for name, (y0, y1, x0, x1) in REGION_BOXES.items():
        r[name] = (int(y0 * sy), int(y1 * sy), int(x0 * sx), int(x1 * sx))
    return r


# name -> (y0, y1, x0, x1) in an 800 x 450 image (read off the picture: the white reference ball is the right-hand ball)
REGION_BOXES = {
    "white ball": (255, 300, 575, 630),    # neutral Lambertian reference: sets the exposure and the adapted white
    "floor": (330, 430, 20, 120),          # neutral gray floor
    "bowl inside": (120, 200, 340, 480),   # inner wall of the gold bowl (light bounces between the walls)
    "bowl outside": (295, 318, 270, 380),  # outer bottom of the bowl (reflects the floor)
    "gold ball": (200, 240, 285, 325),     # the gold ball in the bowl
    "bowl total": (100, 325, 150, 495),    # everything of the bowl and its ball
}


def main():
    src, out_txt = sys.argv[1], sys.argv[2]
    lights = sys.argv[3:]
    lines = []
    for light in lights:
        im = {(m, rep): read_hdr(os.path.join(src, f"goldbowl_{light}_{m}_{rep}.pfm")).astype(np.float64)
              for m in ("rgb", "spectral") for rep in ("a", "b")}
        h, w = im[("rgb", "a")].shape[:2]
        boxes = regions(h, w)
        lines.append(f"=== {light} ({w}x{h}) ===")
        # the white reference ball in the RGB render sets the exposure and the adapted white
        y0, y1, x0, x1 = boxes["white ball"]
        white_rgb = np.mean([im[("rgb", r)][y0:y1, x0:x1].reshape(-1, 3).mean(0) for r in "ab"], axis=0)
        white_adapt = M_RGB_TO_XYZ @ white_rgb
        white_d65 = (M_RGB_TO_XYZ @ np.ones(3)) * white_adapt[1]
        lines.append(f"white ball (RGB render) linear = ({white_rgb[0]:.4f}, {white_rgb[1]:.4f}, {white_rgb[2]:.4f})")
        for name, (y0, y1, x0, x1) in boxes.items():
            def mean_rgb(m, rep):
                px = im[(m, rep)][y0:y1, x0:x1].reshape(-1, 3)
                return px.mean(0)
            rgb_a, rgb_b, spe_a, spe_b = (mean_rgb("rgb", "a"), mean_rgb("rgb", "b"), mean_rgb("spectral", "a"), mean_rgb("spectral", "b"))
            rgb, spe = (rgb_a + rgb_b) / 2, (spe_a + spe_b) / 2
            out = [f"  {name:12s} mean linear RGB: RGB ({rgb[0]:.4f} {rgb[1]:.4f} {rgb[2]:.4f})  spectral ({spe[0]:.4f} {spe[1]:.4f} {spe[2]:.4f})"
                   f"  spectral/RGB ({spe[0] / rgb[0]:.3f} {spe[1] / rgb[1]:.3f} {spe[2] / rgb[2]:.3f})"]
            for wname, white in (("D65", white_d65), ("adapted", white_adapt)):
                lab = lambda c: xyz_to_lab(M_RGB_TO_XYZ @ c, white)            # noqa: E731
                de = float(delta_e2000(lab(rgb), lab(spe)))
                floor_rgb = float(delta_e2000(lab(rgb_a), lab(rgb_b)))
                floor_spe = float(delta_e2000(lab(spe_a), lab(spe_b)))
                # the same, but with the spectral colour clipped into the sRGB gamut (what a display shows)
                de_clip = float(delta_e2000(lab(np.clip(rgb, 0, None)), lab(np.clip(spe, 0, None))))
                out.append(f"      {wname:8s} dE00 RGB vs spectral {de:6.2f}   (clipped at 0: {de_clip:6.2f})   noise floor: RGB a-b {floor_rgb:5.2f}  spectral a-b {floor_spe:5.2f}")
            px = np.concatenate([im[("spectral", r)][y0:y1, x0:x1].reshape(-1, 3) for r in "ab"])
            neg = px.min(axis=1) < 0
            out.append(f"      spectral pixels with a negative channel: {100 * neg.mean():.1f} %"
                       f"   (mean of the negative values per channel: R {px[:, 0][px[:, 0] < 0].mean() if (px[:, 0] < 0).any() else 0:.4f}"
                       f" G {px[:, 1][px[:, 1] < 0].mean() if (px[:, 1] < 0).any() else 0:.4f} B {px[:, 2][px[:, 2] < 0].mean() if (px[:, 2] < 0).any() else 0:.4f})")
            lines += out
    text = "\n".join(lines)
    print(text)
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()
