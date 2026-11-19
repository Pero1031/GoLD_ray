"""
Condition map: when does an RGB renderer differ visibly from the spectral result, and what does the output gamut / gamut mapping do to it?

Reference model of tools/studies/metal_accuracy.py (spectral = 1 nm integral of F^n L against the observer; RGB = the renderer's
eta/k reduction in a working basis), run over  light (blackbody 1500..6500 K + D65)  x  metal  x  number of reflections  x
working basis  x  output (display gamut + mapping).  Normal incidence, the same metal at every bounce.

Display pipeline of the renderer's PNG output (exposure 0): XYZ -> display primaries -> gamut mapping (clip / oklch) -> Reinhard per
channel -> CIEDE2000 on the displayed color (no 8 bit rounding). 'linear' = CIEDE2000 of the two XYZ, no display at all.

  python tools/studies/condition_map.py <out_dir>        (writes results.txt and condition_map.png)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
import metal_accuracy as ma                                   # noqa: E402
import rgb_basis_study as rb                                  # noqa: E402
import gamut as gm                                            # noqa: E402
from gamut_mapping_study import map_gamut                     # noqa: E402
from colordiff import delta_e2000                             # noqa: E402

LAM, W_FULL = ma.LAM, ma.W_FULL
BASES = {"srgb": "sRGB (Rec.709)", "p3": "Display P3", "rec2020": "Rec.2020", "sharp": "Sharp"}
CCTS = [1500, 1850, 2200, 2856, 3500, 4500, 5500, 6500]
BOUNCES = (1, 2, 3, 4, 6)
METALS = ("Au", "Cu", "Ag")


def light_spd(cct):
    spd = ma.blackbody(LAM, float(cct))
    return spd / ((rb.CMF[:, 1] * spd * ma.W_VIS).sum() / rb.YN)          # Y = 1 relative to D65


class Model:
    """Spectral reference and RGB result (XYZ) of one metal in all working bases."""
    def __init__(self):
        self.data = {m: ma.load_metal(m) for m in METALS}
        self.T = {b: rb.basis_matrix(b) for b in BASES}
        self.S = {b: rb.sensors(self.T[b]) for b in BASES}
        self.red = {(m, b): rb.reduce_to_rgb(self.S[b], *self.data[m]) for m in METALS for b in BASES}

    def xyz_pair(self, metal, spd, n, basis):
        eta, k = self.data[metal]
        ref = rb.xyz_of(ma.fresnel_conductor(1.0, eta, k) ** n * spd)
        eta_c, k_c = self.red[(metal, basis)]
        light_rgb = (self.S[basis] * (spd * ma.W_VIS)[:, None]).sum(0)
        rgb = light_rgb * ma.fresnel_conductor(1.0, eta_c, k_c) ** n
        return ref, np.linalg.inv(self.T[basis]) @ rgb


def shown(xyz, out, mapping):
    """What a display of gamut `out` shows for an XYZ (D65 white, Y = 1): returns the displayed XYZ."""
    x = (np.linalg.inv(gm.rgb_to_xyz(out)) @ xyz)[None, :]
    x = map_gamut(x, out, mapping)
    x = x / (1 + x)
    return (x @ gm.rgb_to_xyz(out).T)[0]


def de_lab(xyz_a, xyz_b):
    return float(delta_e2000(gm.lab(xyz_a), gm.lab(xyz_b)))


def outside(xyz, name):
    return bool((np.linalg.inv(gm.rgb_to_xyz(name)) @ xyz).min() < 0)


def main():
    out_dir = sys.argv[1]
    os.makedirs(out_dir, exist_ok=True)
    M = Model()
    L = []

    def emit(s=""):
        print(s)
        L.append(s)

    emit("Condition map (reference model, normal incidence, same metal at every bounce). CIEDE2000 between the RGB renderer and the spectral result.")
    emit("'linear' = no display. Display = gamut primaries, mapping, Reinhard per channel (exposure 0). 'out' = the spectral color is outside that gamut (a channel < 0).")
    emit("Working basis of the RGB renderer: see each table. Written by tools/studies/condition_map.py.")
    emit()

    grid = {}       # (cct, n) -> dict of values, for Au in the sRGB working space (the figure)
    emit("=== A. Au, RGB renderer in sRGB (Rec.709) primaries: light x reflections ===")
    emit(f"{'light':>8s} {'n':>2s} {'linear':>7s} {'sRGB clip':>10s} {'sRGB oklch':>11s} {'P3':>6s} {'Rec2020':>8s} {'hidden by clip':>15s}  out sRGB / P3")
    lights = [(f"{c} K", light_spd(c)) for c in CCTS] + [("D65", ma.d65(LAM) / ((rb.CMF[:, 1] * ma.d65(LAM) * ma.W_VIS).sum() / rb.YN))]
    for name, spd in lights:
        for n in BOUNCES:
            ref, rgb = M.xyz_pair("Au", spd, n, "srgb")
            lin = de_lab(ref, rgb)
            vals = {o: de_lab(shown(ref, out, mp), shown(rgb, out, mp)) for o, (out, mp) in
                    {"clip": ("srgb", "clip"), "oklch": ("srgb", "oklch"), "p3": ("p3", "clip"), "rec2020": ("rec2020", "clip")}.items()}
            hidden = 1 - vals["clip"] / vals["p3"] if vals["p3"] > 0.05 else float("nan")
            grid[(name, n)] = dict(lin=lin, **vals, out_s=outside(ref, "srgb"), out_p=outside(ref, "p3"))
            emit(f"{name:>8s} {n:>2d} {lin:7.2f} {vals['clip']:10.2f} {vals['oklch']:11.2f} {vals['p3']:6.2f} {vals['rec2020']:8.2f} {hidden * 100:14.0f}%  {'yes' if outside(ref, 'srgb') else 'no':>3s} / {'yes' if outside(ref, 'p3') else 'no'}")
        emit()

    emit("=== B. Working space x output (Au; the RGB renderer works in the row's basis, the display is the column) ===")
    emit("Columns: output gamut and mapping. 'linear' needs no display. Same pairs as above.")
    for lname, cct in (("tungsten 2856 K", 2856), ("candle 1850 K", 1850)):
        spd = light_spd(cct)
        for n in (1, 3):
            emit(f"-- {lname}, {n} reflection(s) --")
            emit(f"{'working basis':>16s} {'linear':>7s} {'sRGB clip':>10s} {'sRGB oklch':>11s} {'P3':>6s} {'Rec2020':>8s}")
            for b, label in BASES.items():
                ref, rgb = M.xyz_pair("Au", spd, n, b)
                row = [de_lab(ref, rgb)] + [de_lab(shown(ref, o, mp), shown(rgb, o, mp)) for o, mp in (("srgb", "clip"), ("srgb", "oklch"), ("p3", "clip"), ("rec2020", "clip"))]
                emit(f"{label:>16s} " + " ".join(f"{v:{w}.2f}" for v, w in zip(row, (7, 10, 11, 6, 8))))
        emit()

    emit("=== C. Metals (RGB renderer in sRGB primaries) ===")
    emit(f"{'metal':>5s} {'light':>10s} {'n':>2s} {'linear':>7s} {'sRGB clip':>10s} {'sRGB oklch':>11s} {'P3':>6s}")
    for m in METALS:
        for lname, cct in (("2856 K", 2856), ("1850 K", 1850)):
            spd = light_spd(cct)
            for n in (1, 3):
                ref, rgb = M.xyz_pair(m, spd, n, "srgb")
                row = [de_lab(ref, rgb)] + [de_lab(shown(ref, o, mp), shown(rgb, o, mp)) for o, mp in (("srgb", "clip"), ("srgb", "oklch"), ("p3", "clip"))]
                emit(f"{m:>5s} {lname:>10s} {n:>2d} " + " ".join(f"{v:{w}.2f}" for v, w in zip(row, (7, 10, 11, 6))))
    text = "\n".join(L)
    with open(os.path.join(out_dir, "results.txt"), "w", encoding="utf-8") as f:
        f.write(text + "\n")

    # figure: Au, sRGB working space; panels = linear error, sRGB clip, sRGB oklch, P3
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = [f"{c} K" for c in CCTS] + ["D65"]
    panels = [("lin", "rendering error (no display)"), ("clip", "sRGB display, clip"), ("oklch", "sRGB display, chroma reduction"), ("p3", "Display P3 display")]
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.6), sharey=True)
    for ax, (key, title) in zip(axes, panels):
        img = np.array([[grid[(r, n)][key] for n in BOUNCES] for r in rows])
        im = ax.imshow(img, cmap="magma_r", vmin=0, vmax=10, aspect="auto")
        for i in range(len(rows)):
            for j in range(len(BOUNCES)):
                ax.text(j, i, f"{img[i, j]:.1f}", ha="center", va="center", fontsize=8, color="white" if img[i, j] > 5 else "black")
        ax.set_xticks(range(len(BOUNCES))); ax.set_xticklabels([str(n) for n in BOUNCES])
        ax.set_xlabel("reflections"); ax.set_title(title, fontsize=10)
    axes[0].set_yticks(range(len(rows))); axes[0].set_yticklabels(rows)
    fig.suptitle("Gold: CIEDE2000 between the RGB renderer (sRGB primaries) and the spectral result, by light and reflections (about 1 = just visible)", fontsize=11)
    fig.colorbar(im, ax=axes, shrink=0.8, label="ΔE00")
    fig.savefig(os.path.join(out_dir, "condition_map.png"), dpi=150, bbox_inches="tight")
    print("wrote", out_dir)


if __name__ == "__main__":
    main()
