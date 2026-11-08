"""Figures of docs/HAIR.md: a strip of the hair lobes, and captions on the color comparison.

    python tools/figures/hair_figures.py lobes lobe_R.png lobe_TT.png lobe_TRT.png lobe_all.png docs/img/hair_lobes.png
    python tools/figures/hair_figures.py caption colors.png docs/img/hair_colors.png "Black@0.19" "Brown@0.35" ...

lobes: the four renders of scenes/hair_lock.json with "lobes" set to R, TT, TRT and (nothing) are cropped to the lock and set side by side.
caption: writes captions on a picture; each is  text@x  with x the horizontal position of its centre as a fraction of the width.
"""
import os
import sys

from PIL import Image, ImageDraw

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
from colordiff import font  # noqa: E402

LOBE_TITLES = [
    ("R", "reflection off the surface"),
    ("TT", "through the fiber"),
    ("TRT", "one internal reflection (lights 6x)"),
    ("All lobes", "R + TT + TRT + rest"),
]


def lobes(inputs, out, crop=(270, 0, 720, 506)):
    tiles = [Image.open(p).convert("RGB").crop(crop) for p in inputs]
    w, h = tiles[0].size
    gap = 6
    sheet = Image.new("RGB", (len(tiles) * w + (len(tiles) - 1) * gap, h), (20, 20, 20))
    for i, t in enumerate(tiles):
        sheet.paste(t, (i * (w + gap), 0))
    d = ImageDraw.Draw(sheet, "RGBA")
    d.rectangle((0, h - 62, sheet.size[0], h), fill=(0, 0, 0, 150))
    for i, (title, sub) in enumerate(LOBE_TITLES[: len(tiles)]):
        cx = i * (w + gap) + w // 2
        for text, f, y, col in ((title, font(26, True), h - 56, (255, 255, 255)), (sub, font(15), h - 24, (210, 210, 210))):
            tw = d.textlength(text, font=f)
            d.text((cx - tw / 2, y), text, font=f, fill=col)
    sheet.save(out)


def caption(src, out, labels):
    img = Image.open(src).convert("RGBA")
    w, h = img.size
    img.alpha_composite(Image.new("RGBA", (w, 54), (0, 0, 0, 150)), (0, h - 54))
    d = ImageDraw.Draw(img)
    f = font(24, True)
    for lab in labels:
        text, x = lab.rsplit("@", 1)
        tw = d.textlength(text, font=f)
        d.text((float(x) * w - tw / 2, h - 42), text, font=f, fill=(255, 255, 255))
    img.convert("RGB").save(out)


if __name__ == "__main__":
    if sys.argv[1] == "lobes":
        lobes(sys.argv[2:-1], sys.argv[-1])
    elif sys.argv[1] == "caption":
        caption(sys.argv[2], sys.argv[3], sys.argv[4:])
    else:
        sys.exit(__doc__)
