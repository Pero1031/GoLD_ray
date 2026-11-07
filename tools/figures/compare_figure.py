"""Figure for any scene: the RGB build against the spectral build under several lights, with a CIEDE2000 color-difference map.

Input: for every light L two denoised linear renders (.pfm), from
    GoLD_rayt scene.json --denoise --out <dir>/<prefix>_<L>_rgb.pfm        (RGB build)
    GoLD_rayt scene.json --denoise --out <dir>/<prefix>_<L>_spectral.pfm   (spectral build)
which write <prefix>_<L>_rgb_denoised.pfm and <prefix>_<L>_spectral_denoised.pfm. Same size, same camera.

    python tools/figures/compare_figure.py <dir> <prefix> <out.png> <out.txt> d65 tungsten [--box x0 x1 y0 y1]

--box is the part of the picture that is measured (fractions of the width and height; default the middle of the picture).
The color difference is measured on the displayed images (exposure -> Reinhard per channel -> sRGB, 8 bit; white = D65 display white),
as in tools/studies/gold_compare_linear_fig.py: what a person looking at the two screens could tell apart. It is not a comparison with a photograph.
"""
import os
import sys

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
from hdr_io import read_hdr                                    # noqa: E402
from colordiff import delta_e2000, font, to_lab, LABEL      # noqa: E402
from matplotlib import colormaps                               # noqa: E402

VMAX = 10.0


def display(lin):
    x = np.clip(lin, 0, None)
    x = x / (1 + x)                                            # Reinhard, per channel
    x = np.where(x <= 0.0031308, 12.92 * x, 1.055 * np.power(x, 1 / 2.4) - 0.055)
    return np.clip(np.round(255 * x), 0, 255).astype(np.uint8)


def main():
    args = sys.argv[1:]
    box = (0.15, 0.85, 0.05, 0.98)
    if "--box" in args:
        i = args.index("--box")
        box = tuple(float(v) for v in args[i + 1:i + 5])
        del args[i:i + 5]
    src, prefix, out_png, out_txt = args[:4]
    lights = args[4:]
    cmap = colormaps["magma"]
    rows, lines = [], []
    for light in lights:
        rgb = display(read_hdr(os.path.join(src, f"{prefix}_{light}_rgb_denoised.pfm")).astype(np.float64))
        spe = display(read_hdr(os.path.join(src, f"{prefix}_{light}_spectral_denoised.pfm")).astype(np.float64))
        de = delta_e2000(to_lab(rgb), to_lab(spe))
        h, w = de.shape
        lit = (rgb.max(axis=2) > 12) | (spe.max(axis=2) > 12)
        region = np.zeros_like(lit)
        region[int(box[2] * h):int(box[3] * h), int(box[0] * w):int(box[1] * w)] = True
        m = lit & region
        mean, p50, p95 = float(de[m].mean()), float(np.percentile(de[m], 50)), float(np.percentile(de[m], 95))
        lines.append(f"{light:9s} measured box: dE00 mean {mean:5.2f}   median {p50:5.2f}   95th percentile {p95:5.2f}   max {float(de[m].max()):5.2f}"
                     f"   | whole image: mean {float(de[lit].mean()):5.2f}   95th percentile {float(np.percentile(de[lit], 95)):5.2f}")
        heat = (cmap(np.clip(de / VMAX, 0, 1))[..., :3] * 255).astype(np.uint8)
        heat[~lit] = (heat[~lit] * 0.35).astype(np.uint8)
        rows.append((light, rgb, spe, heat, mean, p95))

    h, w = rows[0][1].shape[:2]
    gap, head = 6, 34
    canvas = Image.new("RGB", (3 * w + 4 * gap, head + len(rows) * (h + gap) + gap), (18, 18, 18))
    d = ImageDraw.Draw(canvas)
    for c, t in enumerate(["RGB renderer", "Spectral renderer", f"CIEDE2000 difference (0 - {VMAX:g})"]):
        d.text((gap + c * (w + gap) + 6, 6), t, font=font(20, True), fill=(235, 235, 235))
    for r, (light, rgb, spe, heat, mean, p95) in enumerate(rows):
        y = head + r * (h + gap)
        for c, im in enumerate((rgb, spe, heat)):
            canvas.paste(Image.fromarray(im), (gap + c * (w + gap), y))
        over = Image.new("RGBA", (w, 30), (0, 0, 0, 150))
        canvas.paste(over, (gap, y + h - 30), over)
        d.text((gap + 8, y + h - 27), LABEL.get(light, light), font=font(18), fill=(255, 255, 255))
        d.text((gap + 2 * (w + gap) + 8, y + h - 27), f"mean {mean:.1f}   95th pct {p95:.1f}", font=font(18), fill=(255, 255, 255))
    canvas.save(out_png, optimize=True)
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(f"RGB build vs spectral build, {prefix}: CIEDE2000 on the displayed sRGB images (lit pixels only)\n" + "\n".join(lines) + "\n")
    print("\n".join(lines))
    print(out_png, canvas.size)


if __name__ == "__main__":
    main()
