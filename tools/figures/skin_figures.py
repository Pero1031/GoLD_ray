"""Figures of docs/SKIN.md.

    python tools/figures/skin_figures.py compare left.png right.png out.png "Left caption" "Right caption"

compare: two renders of the same patch side by side, a caption under each.
"""
import os
import sys

from PIL import Image, ImageDraw

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
from colordiff import font  # noqa: E402


def compare(left, right, out, cap_left, cap_right):
    a, b = Image.open(left).convert("RGB"), Image.open(right).convert("RGB")
    h = min(a.size[1], b.size[1])
    a = a.resize((round(a.size[0] * h / a.size[1]), h), Image.LANCZOS)
    b = b.resize((round(b.size[0] * h / b.size[1]), h), Image.LANCZOS)
    gap = 8
    sheet = Image.new("RGB", (a.size[0] + b.size[0] + gap, h), (20, 20, 20))
    sheet.paste(a, (0, 0))
    sheet.paste(b, (a.size[0] + gap, 0))
    d = ImageDraw.Draw(sheet, "RGBA")
    d.rectangle((0, h - 52, sheet.size[0], h), fill=(0, 0, 0, 150))
    f = font(24, True)
    for text, cx in ((cap_left, a.size[0] / 2), (cap_right, a.size[0] + gap + b.size[0] / 2)):
        tw = d.textlength(text, font=f)
        d.text((cx - tw / 2, h - 42), text, font=f, fill=(255, 255, 255))
    sheet.save(out)


if __name__ == "__main__":
    if sys.argv[1] == "compare":
        compare(*sys.argv[2:7])
    else:
        sys.exit(__doc__)
