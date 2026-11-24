"""
Noise and time of the multiple-bounce models in a rendering (docs/MULTIBOUNCE.md).

scenes/multibounce_gold.json with all four spheres switched to one model, rendered twice with different seeds; the noise of a
render is half the mean squared difference of the two (per pixel, on the gold only), so it does not depend on any reference or on
the small bias of a model. Efficiency = 1 / (noise variance x time): higher is better.

Usage (from the repository root, after building):  python tools/studies/multibounce_noise.py [--exe out/cli-dev-spectral/GoLD_rayt.exe]
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools", "lib"))
from hdr_io import read_pfm  # noqa: E402

OUT = os.path.join(ROOT, "out", "studies", "multibounce")


def strip_comments(text):
    return "\n".join(re.sub(r"^\s*//.*$", "", line) for line in text.splitlines())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default=os.path.join(ROOT, "out", "cli-dev-spectral", "GoLD_rayt.exe"))
    ap.add_argument("--spp", type=int, default=64)
    ap.add_argument("--width", type=int, default=400)
    args = ap.parse_args()
    args.exe = os.path.abspath(args.exe)
    os.makedirs(OUT, exist_ok=True)
    cwd = os.path.join(ROOT, "GoLD_rayt")
    with open(os.path.join(cwd, "scenes", "multibounce_gold.json"), encoding="utf-8") as f:
        scene = json.loads(strip_comments(f.read()))

    results = {}
    for model in ("off", "kulla-conty", "zhang", "heitz"):
        s = json.loads(json.dumps(scene))
        for name, m in s["materials"].items():
            if m.get("type") == "conductor":
                m["multipleScattering"] = model
        path = os.path.join(OUT, f"noise_{model}.json")
        with open(path, "w", encoding="utf-8") as f:
            json.dump(s, f, indent=1)
        imgs, secs = [], []
        for seed in (1, 2):
            out = os.path.join(OUT, f"noise_{model}_{seed}.pfm")
            t0 = time.perf_counter()
            subprocess.run([args.exe, path, "--width", str(args.width), "--spp", str(args.spp), "--seed", str(seed), "--out", out],
                           cwd=cwd, check=True, stdout=subprocess.DEVNULL)
            secs.append(time.perf_counter() - t0)
            imgs.append(read_pfm(out))
        a, b = imgs
        lum = lambda x: 0.2126 * x[..., 0] + 0.7152 * x[..., 1] + 0.0722 * x[..., 2]
        mean = 0.5 * (lum(a) + lum(b))
        h, w = mean.shape
        # the gold: the band of rows through the spheres' centres, where the spheres cover most of the width
        rows = slice(int(h * 0.30), int(h * 0.62))
        mask = np.zeros_like(mean, dtype=bool)
        mask[rows, :] = True
        var = 0.5 * np.mean(((lum(a) - lum(b)) ** 2)[mask])
        level = np.mean(mean[mask])
        t = float(np.mean(secs))
        results[model] = (level, var, t)
        print(f"{model:12s} mean {level:.4f}  noise variance {var:.3e}  time {t:6.2f} s  efficiency {1 / (var * t):.3e}")
    base = results["heitz"]
    print("\nrelative to Heitz: variance x time = " + ", ".join(f"{m} {r[1] * r[2] / (base[1] * base[2]):.2f}" for m, r in results.items()))


if __name__ == "__main__":
    main()
