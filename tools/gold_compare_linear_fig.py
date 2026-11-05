"""Figure for the linear RGB-vs-spectral gold comparison: does a person see the difference?

For every light, a row with
  1. the RGB render as a display shows it (exposure -> Reinhard per channel -> sRGB -> clip, like Film::save for PNG),
  2. the spectral render, the same way,
  3. CIEDE2000 between the two LINEAR renders (no clipping; white = D65 white at the luminance of the white reference ball) = what is measured,
  4. CIEDE2000 between the two DISPLAYED images = what a viewer of the screen can see.
The renders are the mean of the two independent repeats (a, b) and a 3 x 3 box blur, so that the maps show the difference and not the noise.

    python tools/gold_compare_linear_fig.py <dir with the .pfm files> <out.png> d65 tungsten candle
"""
import os
import sys

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hdr_io import read_hdr                                              # noqa: E402
from gold_compare import M_RGB_TO_XYZ, delta_e2000, font, to_lab, LABEL   # noqa: E402
from gold_compare_linear import xyz_to_lab, REGION_BOXES                  # noqa: E402
from matplotlib import colormaps                                          # noqa: E402

VMAX = 10.0


def box_blur(im):
    p = np.pad(im, ((1, 1), (1, 1), (0, 0)), mode="edge")
    return sum(p[dy:dy + im.shape[0], dx:dx + im.shape[1]] for dy in range(3) for dx in range(3)) / 9.0


def display(lin):
    x = np.clip(lin, 0, None)
    x = x / (1 + x)                                       # Reinhard, per channel
    x = np.where(x <= 0.0031308, 12.92 * x, 1.055 * np.power(x, 1 / 2.4) - 0.055)
    return np.clip(np.round(255 * x), 0, 255).astype(np.uint8)


def main():
    src, out_png = sys.argv[1], sys.argv[2]
    lights = sys.argv[3:]
    cmap = colormaps["magma"]
    rows, stat_lines = [], []
    for light in lights:
        im = {m: (read_hdr(os.path.join(src, f"goldbowl_{light}_{m}_a.pfm")).astype(np.float64)
                  + read_hdr(os.path.join(src, f"goldbowl_{light}_{m}_b.pfm")).astype(np.float64)) / 2 for m in ("rgb", "spectral")}
        rgb, spe = box_blur(im["rgb"]), box_blur(im["spectral"])
        h, w = rgb.shape[:2]
        y0, y1, x0, x1 = REGION_BOXES["white ball"]
        wb = rgb[int(y0 / 450 * h):int(y1 / 450 * h), int(x0 / 800 * w):int(x1 / 800 * w)].reshape(-1, 3).mean(0)
        wy = (M_RGB_TO_XYZ @ wb)[1]
        white = (M_RGB_TO_XYZ @ np.ones(3)) * wy
        de_lin = delta_e2000(xyz_to_lab(rgb @ M_RGB_TO_XYZ.T, white), xyz_to_lab(spe @ M_RGB_TO_XYZ.T, white))
        d_rgb, d_spe = display(rgb), display(spe)
        de_disp = delta_e2000(to_lab(d_rgb), to_lab(d_spe))
        # statistics over the bowl (box) only
        by0, by1, bx0, bx1 = REGION_BOXES["bowl total"]
        sl = (slice(int(by0 / 450 * h), int(by1 / 450 * h)), slice(int(bx0 / 800 * w), int(bx1 / 800 * w)))
        lit = (0.2126 * spe[..., 0] + 0.7152 * spe[..., 1] + 0.0722 * spe[..., 2]) > 0.02 * wy   # skip the black / deep shadow, as gold_compare.py does
        m = np.zeros_like(lit)
        m[sl] = True
        m &= lit
        st = (float(de_lin[m].mean()), float(de_disp[m].mean()), float((de_disp[m] > 1.0).mean()) * 100, float((de_disp[m] > 2.3).mean()) * 100)
        p95 = (float(np.percentile(de_lin[m], 95)), float(np.percentile(de_disp[m], 95)))
        line = (f"{light:9s} bowl box, lit pixels: dE00 mean linear {st[0]:5.2f} (95th pct {p95[0]:5.2f})   displayed {st[1]:5.2f} (95th pct {p95[1]:5.2f})"
                f"   displayed pixels > 1: {st[2]:4.1f} %   > 2.3: {st[3]:4.1f} %")
        print(line)
        stat_lines.append(line)
        maps = [(cmap(np.clip(d / VMAX, 0, 1))[..., :3] * 255).astype(np.uint8) for d in (de_lin, de_disp)]
        rows.append((light, d_rgb, d_spe, maps[0], maps[1], st))

    h, w = rows[0][1].shape[:2]
    s = 0.6                                              # panel scale
    pw, ph = int(w * s), int(h * s)
    gap, head = 6, 40
    canvas = Image.new("RGB", (4 * pw + 5 * gap, head + len(rows) * (ph + gap) + gap), (18, 18, 18))
    d = ImageDraw.Draw(canvas)
    titles = ["RGB renderer (as displayed)", "Spectral renderer (as displayed)",
              f"dE00 of the linear renders (0 - {VMAX:g})", f"dE00 of the displayed images (0 - {VMAX:g})"]
    for c, t in enumerate(titles):
        d.text((gap + c * (pw + gap) + 6, 10), t, font=font(17, True), fill=(235, 235, 235))
    for r, (light, a, b, m1, m2, st) in enumerate(rows):
        y = head + r * (ph + gap)
        for c, im in enumerate((a, b, m1, m2)):
            canvas.paste(Image.fromarray(im).resize((pw, ph), Image.LANCZOS), (gap + c * (pw + gap), y))
        for c, text in ((0, LABEL.get(light, light)), (2, f"bowl mean {st[0]:.1f}"), (3, f"bowl mean {st[1]:.1f}   > 2.3: {st[3]:.0f} % of px")):
            over = Image.new("RGBA", (pw, 26), (0, 0, 0, 150))
            canvas.paste(over, (gap + c * (pw + gap), y + ph - 26), over)
            d.text((gap + c * (pw + gap) + 8, y + ph - 23), text, font=font(16), fill=(255, 255, 255))
    canvas.save(out_png, optimize=True)
    print(out_png, canvas.size)
    with open(os.path.splitext(out_png)[0] + "_stats.txt", "w", encoding="utf-8") as f:
        f.write("Per-pixel CIEDE2000 statistics of the figure (mean of two repeats, 3x3 box blur; bowl box, lit pixels)\n" + "\n".join(stat_lines) + "\n")


if __name__ == "__main__":
    main()
