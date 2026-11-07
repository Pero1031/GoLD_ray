"""Thin-film figure: the RGB build against the spectral build of a scene, side by side, with labels.

Input: two denoised linear renders (.pfm) of the same scene,
    GoLD_rayt scene.json --denoise --out <dir>/<name>_rgb.pfm        (RGB build)  -> <name>_rgb_denoised.pfm
    GoLD_rayt scene.json --denoise --out <dir>/<name>_spectral.pfm   (spectral build)

    python tools/figures/thinfilm_figure.py <dir> <out.png> NAME[,NAME...] ["title one|title two|..."]

Each scene is one row: the RGB build on the left (the color of the film under white light, from the film sensor), the spectral build on the right.
The tone mapping is the one of Film::save for a PNG (Reinhard per channel, sRGB curve).
"""
import os
import sys

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
from hdr_io import read_hdr          # noqa: E402
from colordiff import font           # noqa: E402


def display(lin):
    x = np.clip(lin, 0, None)
    x = x / (1 + x)
    x = np.where(x <= 0.0031308, 12.92 * x, 1.055 * np.power(x, 1 / 2.4) - 0.055)
    return np.clip(np.round(255 * x), 0, 255).astype(np.uint8)


def load(src, name, build):
    path = os.path.join(src, f"{name}_{build}_denoised.pfm")
    if not os.path.exists(path):
        path = os.path.join(src, f"{name}_{build}.pfm")
    return Image.fromarray(display(read_hdr(path).astype(np.float64)))


def main():
    src, out = sys.argv[1], sys.argv[2]
    names = sys.argv[3].split(",")
    titles = sys.argv[4].split("|") if len(sys.argv) > 4 else names
    rows = []
    for n in names:
        rows.append((load(src, n, "rgb"), load(src, n, "spectral")))
    width = 900
    gap, head = 6, 42
    panels = []
    for a, b in rows:
        h = int(a.size[1] * (width / 2) / a.size[0])
        panels.append((a.resize((width // 2, h), Image.LANCZOS), b.resize((width // 2, h), Image.LANCZOS)))
    total_h = head + sum(p[0].size[1] + gap for p in panels)
    canvas = Image.new("RGB", (width + gap, total_h), (18, 18, 18))
    d = ImageDraw.Draw(canvas)
    d.text((10, 10), "RGB build (color under white light)", font=font(20, True), fill=(235, 235, 235))
    d.text((width // 2 + gap + 10, 10), "Spectral build (interference per wavelength)", font=font(20, True), fill=(235, 235, 235))
    y = head
    for (a, b), t in zip(panels, titles):
        canvas.paste(a, (0, y))
        canvas.paste(b, (width // 2 + gap, y))
        over = Image.new("RGBA", (width // 2, 28), (0, 0, 0, 150))
        canvas.paste(over, (0, y + a.size[1] - 28), over)
        d.text((8, y + a.size[1] - 25), t, font=font(18), fill=(255, 255, 255))
        y += a.size[1] + gap
    canvas.save(out, optimize=True)
    print(out, canvas.size)


if __name__ == "__main__":
    main()
