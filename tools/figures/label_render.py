"""Tone-maps a linear render (.pfm) the way the renderer does for a PNG and writes labels on it.

    python tools/figures/label_render.py render_denoised.pfm out.png "Wax|220 /m, albedo 1 / 0.997 / 0.972, g 0.4@0.20" "Marble|...@0.50" ...

Each label is  title|second line@x  with x the horizontal position of its centre as a fraction of the width.
The tone mapping is the one of gold_compare_linear_fig.py: Reinhard per channel, then sRGB.
"""
import os
import sys

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
from hdr_io import read_hdr                  # noqa: E402
from colordiff import font                # noqa: E402


def display(lin):
    x = np.clip(lin, 0, None)
    x = x / (1 + x)
    x = np.where(x <= 0.0031308, 12.92 * x, 1.055 * np.power(x, 1 / 2.4) - 0.055)
    return np.clip(np.round(255 * x), 0, 255).astype(np.uint8)


def main():
    pfm, out = sys.argv[1], sys.argv[2]
    labels = sys.argv[3:]
    img = Image.fromarray(display(read_hdr(pfm).astype(np.float64)))
    w, h = img.size
    over = Image.new("RGBA", (w, 74), (0, 0, 0, 140))
    img = img.convert("RGBA")
    img.alpha_composite(over, (0, h - 74))
    d = ImageDraw.Draw(img)
    for lab in labels:
        text, x = lab.rsplit("@", 1)
        title, sub = (text.split("|", 1) + [""])[:2]
        cx = float(x) * w
        t1, t2 = font(26, True), font(16)
        d.text((cx - d.textlength(title, font=t1) / 2, h - 68), title, font=t1, fill=(255, 255, 255, 255))
        d.text((cx - d.textlength(sub, font=t2) / 2, h - 32), sub, font=t2, fill=(215, 215, 215, 255))
    img.convert("RGB").save(out, optimize=True)
    print(out, img.size)


if __name__ == "__main__":
    main()
