"""Makes the dye concentration table c(t) of the engine's "fading from RGB" (a "fading" material with "rgb" or "image") from a measured series.

Model (Kubelka-Munk, one parameter): the dye adds D(lambda) = max(0, K/S_color - K/S_paper) to the K/S of the bare paper, and fading lowers it:
    K/S(t) = K/S_paper + c(t) * D(lambda),     K/S = (1 - R)^2 / (2 R)
For every color of the series, c is fitted by least squares to its measured spectra (440..680 nm) in every round; the table is the median
over the colors. It is meant for dyed paper (the dye series); the pigment and origami series do not follow this model.

Input:  <series folder>/hours.csv and r<N>/<color>.csv (written by tools/data/fading_refangle.py).
Output: a CSV "hours,c" with the per-color fits in the comments.

Usage:
  python tools/data/fading_dye_model.py GoLD_rayt/assets/Fading_measured/dye_paper --out GoLD_rayt/assets/Fading_measured/dye_rgb_model.csv [--paper 0.9]
"""
import argparse
import glob
import os

import numpy as np


def read_curve(path):
    wl, v = [], []
    for line in open(path, encoding="utf-8"):
        p = line.strip().split(",")
        try:
            wl.append(float(p[0])); v.append(float(p[1]))
        except (ValueError, IndexError):
            continue
    return np.array(wl), np.array(v)


def ks(r):
    r = np.clip(r, 0.02, 0.98)
    return (1 - r) ** 2 / (2 * r)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("series", help="series folder with hours.csv and r1, r2, ...")
    ap.add_argument("--out", required=True)
    ap.add_argument("--paper", type=float, default=0.9, help="reflectance of the bare paper under the dye (flat)")
    ap.add_argument("--exclude", nargs="*", default=["white"])
    a = ap.parse_args(argv)
    hours = {}
    for line in open(os.path.join(a.series, "hours.csv"), encoding="utf-8"):
        p = line.strip().split(",")
        try:
            hours[int(p[0])] = float(p[1])
        except (ValueError, IndexError):
            continue
    rounds = sorted(hours)
    colors = sorted(os.path.splitext(os.path.basename(f))[0] for f in glob.glob(os.path.join(a.series, f"r{rounds[0]}", "*.csv")))
    colors = [c for c in colors if c not in a.exclude]
    sub = ks(np.array(a.paper))
    fits = {}
    for c in colors:
        wl, r0 = read_curve(os.path.join(a.series, f"r{rounds[0]}", c + ".csv"))
        band = (wl >= 440) & (wl <= 680)
        d = np.maximum(ks(r0) - sub, 0)[band]
        row = []
        for k in rounds:
            _, rk = read_curve(os.path.join(a.series, f"r{k}", c + ".csv"))
            row.append(1.0 if k == rounds[0] else float(np.clip((((ks(rk) - sub)[band]) * d).sum() / (d * d).sum(), 0, 1)))
        fits[c] = row
    median = np.median(np.array([fits[c] for c in colors]), axis=0)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8", newline="\n") as f:
        f.write("# Dye concentration c(t) of the series " + os.path.basename(os.path.normpath(a.series)) + ": K/S = K/S_paper + c(t) max(0, K/S_color - K/S_paper)\n")
        f.write(f"# bare paper reflectance {a.paper:g} (flat); median over the colors of the least-squares fits to the measured spectra, 440..680 nm\n")
        f.write("# hours of exposure: " + " ".join(f"{hours[k]:g}" for k in rounds) + "\n")
        for c in colors:
            f.write(f"# {c:9s} " + " ".join(f"{x:.3f}" for x in fits[c]) + "\n")
        f.write("hours,c\n")
        for k, m in zip(rounds, median):
            f.write(f"{hours[k]:g},{m:.4f}\n")
    print(f"{a.out}: " + ", ".join(f"{hours[k]:g} h: {m:.3f}" for k, m in zip(rounds, median)))


if __name__ == "__main__":
    main()
