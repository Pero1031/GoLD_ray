"""The small assets of scenes/gpu_textures.json (docs/GPU.md), made reproducibly:

  assets/Textures/gpu_pattern.png      64 x 48 sRGB image: a coloured grid over a diagonal ramp (bilinear filtering and the repeat show)
  assets/Textures/gpu_cloth.sptex      32 x 32 spectral texture (Core/SpectralTexture.hpp): stripes of two smooth dye spectra
  assets/Meshes/gpu_uv_tile.obj        a wavy 12 x 12 tile with texture coordinates from 0 to 2 (the texture repeats twice)

    python tools/test/make_gpu_texture_assets.py
"""
import math
import os
import sys

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools", "data"))
from hsm_texture import write_sptex  # noqa: E402

ASSETS = os.path.join(ROOT, "GoLD_rayt", "assets")


def pattern():
    w, h = 64, 48
    img = np.zeros((h, w, 3), np.uint8)
    for y in range(h):
        for x in range(w):
            ramp = (x + y) / (w + h)
            cell = ((x // 8) + (y // 8)) % 4
            base = [(220, 60, 50), (60, 170, 80), (50, 90, 210), (230, 200, 60)][cell]
            img[y, x] = [int(c * (0.45 + 0.55 * ramp)) for c in base]
            if x % 8 == 0 or y % 8 == 0:
                img[y, x] = (245, 245, 240)
    path = os.path.join(ASSETS, "Textures", "gpu_pattern.png")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    Image.fromarray(img).save(path)
    return path


def cloth():
    w, h = 32, 32
    lam = np.arange(380, 781, 10, dtype=np.float64)
    dye_a = 0.08 + 0.75 / (1 + np.exp(-(lam - 590) / 18))          # a red dye
    dye_b = 0.06 + 0.55 * np.exp(-0.5 * ((lam - 470) / 35) ** 2)     # a blue dye
    mean = 0.5 * (dye_a + dye_b)
    basis = np.stack([dye_a - mean, (dye_b - mean) * 0.5])
    coef = np.zeros((h * w, 2))
    for y in range(h):
        for x in range(w):
            t = 0.5 + 0.5 * math.sin(2 * math.pi * (x + 0.5 * y) / 8)   # slanted stripes
            coef[y * w + x] = [t - 0.5 + 0.5, (1 - t) - 0.5 + 0.5]
    # reflectance = mean + coef @ basis; keep it in [0, 1]
    spectra = mean + coef @ basis
    assert spectra.min() >= 0 and spectra.max() <= 1, (spectra.min(), spectra.max())
    path = os.path.join(ASSETS, "Textures", "gpu_cloth.sptex")
    write_sptex(path, w, h, 380.0, 10.0, mean.astype(np.float32), basis.astype(np.float32), coef.astype(np.float32))
    return path


def tile():
    n = 12
    path = os.path.join(ASSETS, "Meshes", "gpu_uv_tile.obj")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="\n") as f:
        f.write("# a wavy tile with texture coordinates from 0 to 2 (tools/test/make_gpu_texture_assets.py)\n")
        for j in range(n + 1):
            for i in range(n + 1):
                x, z = i / n - 0.5, j / n - 0.5
                f.write("v %.6f %.6f %.6f\n" % (x, 0.03 * math.sin(2 * math.pi * x * 1.5) * math.cos(2 * math.pi * z), z))
        for j in range(n + 1):
            for i in range(n + 1):
                f.write("vt %.6f %.6f\n" % (2 * i / n, 2 * j / n))
        for j in range(n):
            for i in range(n):
                a = j * (n + 1) + i + 1
                b, c, d = a + 1, a + n + 2, a + n + 1
                f.write("f %d/%d %d/%d %d/%d\n" % (a, a, d, d, c, c))
                f.write("f %d/%d %d/%d %d/%d\n" % (a, a, c, c, b, b))
    return path


if __name__ == "__main__":
    for p in (pattern(), cloth(), tile()):
        print(p, os.path.getsize(p), "bytes")
