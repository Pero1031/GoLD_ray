"""
The fluorescence demo (docs/FLUORESCENCE.md): six samples of the RIT bispectral database as flat patches under three lights.

  1. daylight: D65 with its ultraviolet (light preset "d65uv")
  2. the same daylight with the fluorescence of every sample switched off (each sample keeps its reflectance, the matrix diagonal)
  3. a black light: assets/Illuminants/blacklight_365nm.csv (a 365 nm band and a weak violet leak)

Samples, left to right: a white ceramic tile, white printer paper (optical brightener), a highlighter on an index card, a fluorescent
orange pigment, an orange ceramic tile, a yellow textile. The matrices must be on this PC (tools/data/rit_bispectral.py fetch).
Needs the spectral build. Writes the scene files and renders to <out>, and the figure to --fig.

  python tools/studies/fluorescence_figure.py --exe out/cli-dev-spectral/GoLD_rayt.exe --out out/renders/fluorescence [--spp 512]
                                              [--fig docs/img/fluorescence.png]
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools", "lib"))
from hdr_io import read_pfm                                  # noqa: E402

DATA = "assets/Fluorescence/rit_bispectral/Labsphere/"       # relative to GoLD_rayt/, where the renderer runs
SAMPLES = [("BCRA white tile", "bcraday01/BCRWH1WO.BFC"), ("white paper", "PHP8HPWH.BFC"), ("highlighter", "IXCSSTHL.BFC"),
           ("fluorescent orange", "HERPIORA.BFC"), ("orange tile", "bcraday01/BCROR1WR.BFC"), ("yellow textile", "TEXTYELL.BFC")]
LIGHTS = {"daylight": {"preset": "d65uv"}, "blacklight": {"file": "assets/Illuminants/blacklight_365nm.csv"}}


def scene(light, fluorescence, radiance):
    mats = {"floor": {"type": "diffuse", "albedo": 0.06},
            "lamp": {"type": "emission", "spectrum": LIGHTS[light], "radiance": radiance}}
    objects = [{"shape": "quad", "corner": [-1, 0, -1], "edge1": [0, 0, 2], "edge2": [2, 0, 0], "material": "floor"},
               {"shape": "quad", "corner": [-0.5, 0.9, -0.4], "edge1": [1.0, 0, 0], "edge2": [0, 0, 0.8], "material": "lamp"}]
    size, gap = 0.1, 0.13
    x0 = -gap * (len(SAMPLES) - 1) / 2 - size / 2
    for i, (label, path) in enumerate(SAMPLES):
        name = f"s{i}"
        mats[name] = {"type": "fluorescent", "file": DATA + path, "fluorescence": fluorescence}
        x = x0 + i * gap
        objects.append({"shape": "quad", "corner": [x, 0.001, -size / 2], "edge1": [0, 0, size], "edge2": [size, 0, 0], "material": name})
    return {"options": {"width": 960, "height": 300},
            "camera": {"from": [0, 0.42, 0.52], "at": [0, 0, 0], "fov": 24},
            "materials": mats, "objects": objects}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--spp", type=int, default=512)
    ap.add_argument("--fig", default="")
    ap.add_argument("--day-radiance", type=float, default=3.0)
    ap.add_argument("--black-radiance", type=float, default=0.05)
    a = ap.parse_args()
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)
    runs = [("daylight", True, a.day_radiance, "Daylight (D65 with its ultraviolet)"),
            ("daylight", False, a.day_radiance, "The same daylight, fluorescence switched off"),
            ("blacklight", True, a.black_radiance, "Black light (365 nm)")]
    images = []
    for light, fluo, rad, title in runs:
        tag = f"{light}_{'on' if fluo else 'off'}"
        sf = os.path.join(out, f"fluorescence_{tag}.json")
        with open(sf, "w", newline="\n") as f:
            json.dump(scene(light, fluo, rad), f, indent=1)
        pfm = os.path.join(out, f"fluorescence_{tag}.pfm")
        cmd = [os.path.abspath(a.exe), sf, "--spp", str(a.spp), "--seed", "1", "--out", pfm]
        r = subprocess.run(cmd, cwd=os.path.join(ROOT, "GoLD_rayt"), capture_output=True, text=True, errors="replace")
        if r.returncode != 0:
            sys.exit(r.stdout[-3000:] + r.stderr[-2000:])
        images.append((title, read_pfm(pfm).astype(np.float64)))
        print("rendered", tag)
    if a.fig:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        def show(x):                                         # linear sRGB -> display: Reinhard on luminance keeps the hues, sRGB curve
            lum = np.maximum(x @ np.array([0.2126, 0.7152, 0.0722]), 1e-9)[..., None]
            y = np.clip(x * (1 / (1 + lum)) * 1.6, 0, 1)
            return np.where(y <= 0.0031308, 12.92 * y, 1.055 * np.power(y, 1 / 2.4) - 0.055)

        fig, axes = plt.subplots(len(images), 1, figsize=(9.6, 3.0 * len(images) + 0.5))
        for ax, (title, img) in zip(axes, images):
            ax.imshow(show(img))
            ax.set_title(title, fontsize=11)
            ax.axis("off")
        fig.text(0.5, 0.01, "Left to right: " + ", ".join(s for s, _ in SAMPLES) + "  (RIT bispectral database, Gonzalez & Fairchild 2000)",
                 ha="center", fontsize=8.5)
        fig.tight_layout(rect=(0, 0.03, 1, 1))
        fig.savefig(a.fig, dpi=100)
        print("wrote", a.fig)


if __name__ == "__main__":
    main()
