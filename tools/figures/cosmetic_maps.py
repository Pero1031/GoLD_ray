"""The reflectance of thick foundation layers as coloured patches, the way Lanza et al. (2024) plot their measurements and fits
(their figures 4 and 6): the angle of observation across (-40 ... +70 degrees, + is the mirror side), the angle of incidence down
(20 ... 70 degrees). A patch shows pi f(theta_i, theta_o) in sRGB: white = as bright as a white diffuser, or brighter.

    GoLD_rayt.exe --scene cosmetic > cosmetic.log
    python tools/figures/cosmetic_maps.py cosmetic.log OUT.png

Reads the lines "[Cosmetic]   map NAME | THETA_I | r,g,b r,g,b ..." of the self-check.
"""
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

maps = {}
for line in open(sys.argv[1], encoding="utf-8", errors="replace"):
    if "[Cosmetic]   map " not in line:
        continue
    name, theta, values = (part.strip() for part in line.split("map ", 1)[1].split("|"))
    maps.setdefault(name, {})[int(theta)] = np.array([[float(c) for c in v.split(",")] for v in values.split()])
if not maps:
    sys.exit("no map lines in " + sys.argv[1])


def srgb(lin):
    lin = np.clip(lin, 0, 1)
    return np.where(lin <= 0.0031308, 12.92 * lin, 1.055 * lin ** (1 / 2.4) - 0.055)


cell_w, cell_h, pad, head = 22, 26, 34, 30
names = list(maps)
thetas = sorted(next(iter(maps.values())))
n_o = len(next(iter(maps.values()))[thetas[0]])
panel_w, panel_h = n_o * cell_w, len(thetas) * cell_h
try:
    font = ImageFont.truetype("arial.ttf", 15)
    small = ImageFont.truetype("arial.ttf", 12)
except OSError:
    font = small = ImageFont.load_default()
sheet = Image.new("RGB", (pad + len(names) * (panel_w + pad), head + panel_h + 44), (255, 255, 255))
draw = ImageDraw.Draw(sheet)
for k, name in enumerate(names):
    x0 = pad + k * (panel_w + pad)
    draw.text((x0, 6), name, fill=(0, 0, 0), font=font)
    for r, th in enumerate(thetas):
        rgb = (srgb(maps[name][th]) * 255 + 0.5).astype(np.uint8)
        for c in range(n_o):
            draw.rectangle([x0 + c * cell_w, head + r * cell_h, x0 + (c + 1) * cell_w - 1, head + (r + 1) * cell_h - 1], fill=tuple(int(v) for v in rgb[c]))
        if k == 0:
            draw.text((6, head + r * cell_h + 6), str(th), fill=(0, 0, 0), font=small)
    for c, label in ((0, "-40"), (4, "0"), (8, "+40"), (11, "+70")):
        draw.text((x0 + c * cell_w + 2, head + panel_h + 4), label, fill=(0, 0, 0), font=small)
draw.text((pad, head + panel_h + 24), "across: angle of observation (+ = mirror side); down: angle of incidence; white = a white diffuser or brighter", fill=(0, 0, 0), font=small)
sheet.save(sys.argv[2])
print("wrote", sys.argv[2], sheet.size)
