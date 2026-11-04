"""Side-by-side figure: the same gold scene rendered by the RGB build and by the spectral build, under several lights, with a CIEDE2000
color-difference map.

Input: for every light L a pair of renders  <dir>/goldbowl_<L>_rgb_denoised.png  and  <dir>/goldbowl_<L>_spectral_denoised.png  (both from
`GoLD_rayt --scene goldbowl --light L --denoise --out goldbowl_<L>_<mode>.png`, the RGB build and the spectral build; same size, same camera).

Output: the figure (PNG) and a text file with the numbers.

    python tools/gold_compare.py <dir> <out.png> <out.txt> d65 tungsten candle

The color difference is measured on the displayed (sRGB, 8 bit) images, white = D65 display white. It is not a comparison with a photograph.
"""
import os
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont
import matplotlib

matplotlib.use("Agg")
from matplotlib import colormaps

M_RGB_TO_XYZ = np.array([[0.4124564, 0.3575761, 0.1804375],
                         [0.2126729, 0.7151522, 0.0721750],
                         [0.0193339, 0.1191920, 0.9503041]])
WHITE = M_RGB_TO_XYZ @ np.ones(3)


def srgb_to_linear(v):
    return np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4)


def to_lab(img8):
    lin = srgb_to_linear(img8.astype(np.float64) / 255.0)
    xyz = lin @ M_RGB_TO_XYZ.T / WHITE
    f = np.where(xyz > (6 / 29) ** 3, np.cbrt(xyz), xyz / (3 * (6 / 29) ** 2) + 4 / 29)
    return np.stack([116 * f[..., 1] - 16, 500 * (f[..., 0] - f[..., 1]), 200 * (f[..., 1] - f[..., 2])], axis=-1)


def delta_e2000(lab1, lab2):
    L1, a1, b1 = lab1[..., 0], lab1[..., 1], lab1[..., 2]
    L2, a2, b2 = lab2[..., 0], lab2[..., 1], lab2[..., 2]
    C1, C2 = np.hypot(a1, b1), np.hypot(a2, b2)
    Cm = (C1 + C2) / 2
    G = 0.5 * (1 - np.sqrt(Cm ** 7 / (Cm ** 7 + 25.0 ** 7)))
    a1p, a2p = (1 + G) * a1, (1 + G) * a2
    C1p, C2p = np.hypot(a1p, b1), np.hypot(a2p, b2)
    h1p = np.degrees(np.arctan2(b1, a1p)) % 360
    h2p = np.degrees(np.arctan2(b2, a2p)) % 360
    dLp, dCp = L2 - L1, C2p - C1p
    dh = h2p - h1p
    dh = np.where(dh > 180, dh - 360, np.where(dh < -180, dh + 360, dh))
    dh = np.where(C1p * C2p == 0, 0, dh)
    dHp = 2 * np.sqrt(C1p * C2p) * np.sin(np.radians(dh) / 2)
    Lpm, Cpm = (L1 + L2) / 2, (C1p + C2p) / 2
    hs = h1p + h2p
    hpm = np.where(np.abs(h1p - h2p) <= 180, hs / 2, np.where(hs < 360, (hs + 360) / 2, (hs - 360) / 2))
    hpm = np.where(C1p * C2p == 0, hs, hpm)
    T = (1 - 0.17 * np.cos(np.radians(hpm - 30)) + 0.24 * np.cos(np.radians(2 * hpm))
         + 0.32 * np.cos(np.radians(3 * hpm + 6)) - 0.20 * np.cos(np.radians(4 * hpm - 63)))
    dTheta = 30 * np.exp(-(((hpm - 275) / 25) ** 2))
    Rc = 2 * np.sqrt(Cpm ** 7 / (Cpm ** 7 + 25.0 ** 7))
    Sl = 1 + 0.015 * (Lpm - 50) ** 2 / np.sqrt(20 + (Lpm - 50) ** 2)
    Sc, Sh = 1 + 0.045 * Cpm, 1 + 0.015 * Cpm * T
    Rt = -np.sin(np.radians(2 * dTheta)) * Rc
    return np.sqrt((dLp / Sl) ** 2 + (dCp / Sc) ** 2 + (dHp / Sh) ** 2 + Rt * (dCp / Sc) * (dHp / Sh))


def font(size, bold=False):
    for name in (("segoeuib.ttf" if bold else "segoeui.ttf"), "arial.ttf"):
        p = os.path.join(r"C:\Windows\Fonts", name)
        if os.path.exists(p):
            return ImageFont.truetype(p, size)
    return ImageFont.load_default()


LABEL = {"d65": "D65 daylight", "tungsten": "Tungsten (illuminant A, 2856 K)", "candle": "Candle (1850 K)", "led": "Cool-white LED"}


def main():
    src, out_png, out_txt = sys.argv[1], sys.argv[2], sys.argv[3]
    lights = sys.argv[4:]
    cmap = colormaps["magma"]
    VMAX = 8.0
    rows, lines = [], []
    for light in lights:
        rgb = np.array(Image.open(os.path.join(src, f"goldbowl_{light}_rgb_denoised.png")).convert("RGB"))
        spe = np.array(Image.open(os.path.join(src, f"goldbowl_{light}_spectral_denoised.png")).convert("RGB"))
        de = delta_e2000(to_lab(rgb), to_lab(spe))
        lum = 0.2126 * srgb_to_linear(spe[..., 0] / 255.0) + 0.7152 * srgb_to_linear(spe[..., 1] / 255.0) + 0.0722 * srgb_to_linear(spe[..., 2] / 255.0)
        lit = lum > 0.02   # skip the black background / deep shadow
        # the bowl and the ball in it: a box around them (the floor and the white ball are the same in both renderers and would dilute the numbers);
        # the box is given for an 800 x 450 image and scaled
        hh, ww = de.shape
        box = np.zeros_like(lit)
        box[int(125 / 450 * hh):int(365 / 450 * hh), int(150 / 800 * ww):int(510 / 800 * ww)] = True
        gold = lit & box
        mean, p50, p95 = float(de[gold].mean()), float(np.percentile(de[gold], 50)), float(np.percentile(de[gold], 95))
        lines.append(f"{light:9s} bowl box: dE00 mean {mean:5.2f}   median {p50:5.2f}   95th percentile {p95:5.2f}   max {float(de[gold].max()):5.2f}"
                     f"   | whole image: mean {float(de[lit].mean()):5.2f}   95th percentile {float(np.percentile(de[lit], 95)):5.2f}")
        heat = (cmap(np.clip(de / VMAX, 0, 1))[..., :3] * 255).astype(np.uint8)
        heat[~lit] = (heat[~lit] * 0.35).astype(np.uint8)
        Image.fromarray(heat).save(os.path.join(src, f"goldbowl_{light}_dE.png"))   # kept separately (the README hero uses it)
        rows.append((light, rgb, spe, heat, mean, p95))

    h, w = rows[0][1].shape[:2]
    gap, head, side = 6, 34, 0
    W = 3 * w + 4 * gap
    H = head + len(rows) * (h + gap) + gap
    canvas = Image.new("RGB", (W, H), (18, 18, 18))
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
        d.text((gap + 2 * (w + gap) + 8, y + h - 27), f"bowl: mean {mean:.1f}   95th pct {p95:.1f}", font=font(18), fill=(255, 255, 255))
    canvas.save(out_png, optimize=True)
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("RGB build vs spectral build, scene goldbowl: CIEDE2000 on the displayed sRGB images (lit pixels only)\n")
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(out_png, canvas.size)


if __name__ == "__main__":
    main()
