"""
Does a wider RGB working space (or a smarter RGB renderer) remove the RGB-vs-spectral color error of metals?

The question behind docs/data/metal_accuracy_results.txt: that study fixes the RGB renderer to sRGB primaries. A reader can
fairly ask whether the error is a property of "RGB rendering" or of "sRGB primaries". This script keeps the model of
tools/studies/metal_accuracy.py (spectral reference = integral of F(lambda)^n L(lambda) against the observer; Johnson & Christy
tables; the renderer's eta/k reduction) and changes only the RGB basis the "RGB renderer" works in:

  sRGB, Display P3, Rec.2020     working RGB = linear RGB of that gamut (D65 white = (1,1,1));
                                 eta/k are reduced per channel with that basis (F0 matched under D65, as the renderer does).
  Sharp                          the "Sharp" sensors of Ward & Eydelberg-Vileshin (2002): non-negative, narrow channels
                                 (the matrix is written from memory of that paper: check before quoting).
  + prefilter                    the same basis, but the material's per-channel reflectance is computed under the actual
                                 light (spectral prefiltering): rho_c = Integral(s_c F L) / Integral(s_c L), normal incidence;
                                 exact for one reflection, the best an RGB renderer can do when it knows the light.

The color difference is CIEDE2000 between the spectral XYZ and the XYZ of the RGB result, white = D65 (Y = 1), linear, no clipping
(a pure measure of the rendering error, independent of any display).

Usage:  python tools/studies/rgb_basis_study.py [--out docs/data/rgb_basis_study_results.txt]
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "tools", "lib"))
import metal_accuracy as ma          # noqa: E402
import gamut as gm                   # noqa: E402
import colour                        # noqa: E402

LAM, VIS, W_VIS, W_FULL = ma.LAM, ma.VIS, ma.W_VIS, ma.W_FULL
CMF = ma.cmf_wyman(LAM)                                       # (N, 3) XYZ observer (Wyman fit, like the renderer's film)
D65_SPD = ma.d65(LAM)
XYZ_D65 = (CMF * (D65_SPD * W_VIS)[:, None]).sum(0)
YN = XYZ_D65[1]                                               # D65 -> Y = 1 after dividing by YN
WHITE_XYZ = XYZ_D65 / YN

# XYZ -> sharpened sensors (Ward & Eydelberg-Vileshin 2002, "Sharp"); white is renormalized to D65 below
M_XYZ_TO_SHARP = np.array([[1.2694, -0.0988, -0.1706],
                           [-0.8364, 1.8006, 0.0357],
                           [0.0297, -0.0315, 1.0018]])


def basis_matrix(name):
    """T: XYZ (D65 white, Y = 1) -> working RGB, with D65 -> (1, 1, 1)."""
    if name == "sharp":
        m = M_XYZ_TO_SHARP
        return np.diag(1.0 / (m @ WHITE_XYZ)) @ m
    return np.linalg.inv(gm.rgb_to_xyz(name))


BASES = [("srgb", "sRGB"), ("p3", "Display P3"), ("rec2020", "Rec.2020"), ("sharp", "Sharp")]


def sensors(T):
    """(N, 3): the working-RGB channel responses s_c(lambda) for a radiance sampled on LAM and normalized so that D65 (Y = 1) -> (1, 1, 1)."""
    return (CMF @ T.T) / YN


def d65_average(S, values):
    """The renderer's spectrumToRGB_D65: sensor response of `values` under D65, divided by the response of 1."""
    w = (D65_SPD * W_VIS)[:, None]
    return (S * w * values[:, None]).sum(0) / (S * w).sum(0)


def reduce_to_rgb(S, eta, k):
    """integrateIorToRgb_CMF in an arbitrary basis: eta_rgb = D65 average of eta, k_rgb so that F0 matches the D65 average of F0."""
    r = ma.f0(eta, k)
    r_max = min(r[VIS].max(), 0.9999)
    f0_rgb = d65_average(S, r)
    eta_c = np.maximum(d65_average(S, eta), 0.01)
    f0_c = np.clip(f0_rgb, 0.0, r_max)
    k2 = (f0_c * (eta_c + 1) ** 2 - (eta_c - 1) ** 2) / (1 - f0_c)
    return eta_c, np.where(k2 > 0, np.sqrt(np.maximum(k2, 0)), 0.0)


def light_set():
    raw = {"D65": D65_SPD, "A (2856K)": ma.blackbody(LAM, 2856.0), "candle (1850K)": ma.blackbody(LAM, 1850.0), "LED cool-white": ma.load_led()}
    out = {}
    for name, spd in raw.items():
        y = (CMF[:, 1] * spd * W_VIS).sum() / YN       # Y of the light relative to D65's (D65 -> 1)
        out[name] = spd / y
    return out


def xyz_of(spec, weights=W_FULL):
    return (CMF * (spec * weights)[:, None]).sum(0) / YN


def de2000(xyz_a, xyz_b):
    lab_a = colour.XYZ_to_Lab(xyz_a, ma._XY_D65)
    lab_b = colour.XYZ_to_Lab(xyz_b, ma._XY_D65)
    return float(colour.delta_E(lab_a, lab_b, method="CIE 2000"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="")
    ap.add_argument("--fig", default="", help="write the bar figure here")
    ap.add_argument("--lang", default="en", choices=["en", "ja"])
    a = ap.parse_args()
    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit("RGB working space and the metal color error (CIEDE2000, linear, white = D65, normal incidence, same metal at every bounce).")
    emit("Reference: spectral (1 nm integral). Models: the renderer's RGB pipeline in other bases; 'prefilter' = reflectance integrated under the light.")
    emit("Written by tools/studies/rgb_basis_study.py. The Sharp matrix is from memory of Ward & Eydelberg-Vileshin 2002: verify before quoting.")
    emit()
    lights = light_set()
    data = {m: ma.load_metal(m) for m in ("Au", "Cu", "Ag")}
    bounce_list = (1, 2, 3, 4)
    table = {}      # (metal, light, n, model) -> dE

    for metal, (eta, k) in data.items():
        emit(f"=== {metal} ===")
        header = f"{'light':16s} {'n':>2s}  " + "  ".join(f"{label:>11s}" for _, label in BASES) + "  |  " + "  ".join(f"{label + ' +pf':>15s}" for _, label in BASES)
        emit(header)
        F0 = ma.fresnel_conductor(1.0, eta, k)
        for ln, spd in lights.items():
            for n in bounce_list:
                ref = xyz_of(F0 ** n * spd)
                row_plain, row_pf = [], []
                for bname, _ in BASES:
                    T = basis_matrix(bname)
                    S = sensors(T)
                    eta_c, k_c = reduce_to_rgb(S, eta, k)
                    light_rgb = (S * (spd * W_VIS)[:, None]).sum(0)
                    rgb = light_rgb * ma.fresnel_conductor(1.0, eta_c, k_c) ** n
                    plain = de2000(np.linalg.inv(T) @ rgb, ref)
                    # spectral prefiltering: per-channel reflectance of the actual light
                    rho = (S * (F0 * spd * W_FULL)[:, None]).sum(0) / (S * (spd * W_FULL)[:, None]).sum(0)
                    rgb_pf = (S * (spd * W_FULL)[:, None]).sum(0) * rho ** n
                    pf = de2000(np.linalg.inv(T) @ rgb_pf, ref)
                    row_plain.append(plain); row_pf.append(pf)
                    table[(metal, ln, n, bname, "plain")] = plain
                    table[(metal, ln, n, bname, "pf")] = pf
                emit(f"{ln:16s} {n:>2d}  " + "  ".join(f"{v:11.2f}" for v in row_plain) + "  |  " + "  ".join(f"{v:15.2f}" for v in row_pf))
        emit()

    # a short summary for the paper figure
    emit("Summary (Au): mean over the 3 lights without the LED, n = 1 and n = 3")
    for n in (1, 3):
        for kind, lab in (("plain", "reduced eta/k"), ("pf", "prefiltered")):
            vals = {b: np.mean([table[("Au", ln, n, b, kind)] for ln in ("D65", "A (2856K)", "candle (1850K)")]) for b, _ in BASES}
            emit(f"  n = {n}, {lab:14s}: " + ", ".join(f"{label} {vals[b]:.2f}" for b, label in BASES))

    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        with open(a.out, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(lines) + "\n")
        print("wrote", a.out)
    if a.fig:
        figure(table, a.fig, a.lang)


TXT = {"en": dict(lights=["D65 (daylight)", "Tungsten (2856 K)", "Candle (1850 K)"], plain="reduced eta/k (the renderer's RGB mode)", pf="with spectral prefiltering",
                  title="Gold: color error of an RGB renderer by working space (CIEDE2000 to the spectral result)", ylabel="CIEDE2000", nlab="{} reflection(s)",
                  jnd="visible-difference level (about 1)", off="off scale"),
       "ja": dict(lights=["D65（昼光）", "タングステン（2856 K）", "ろうそく（1850 K）"], plain="eta/k を縮約（今の RGB モード）", pf="分光プレフィルタあり",
                  title="金：RGB 計算の色誤差と作業色空間（分光の結果との CIEDE2000）", ylabel="CIEDE2000", nlab="{} 回反射",
                  jnd="差が見える目安（約 1）", off="範囲外")}


def figure(table, path, lang):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    if lang == "ja":
        for name in ("Yu Gothic", "Meiryo", "BIZ UDGothic", "MS Gothic"):
            if any(f.name == name for f in font_manager.fontManager.ttflist):
                plt.rcParams["font.family"] = name
                break
    t = TXT[lang]
    cols = {"srgb": "#6c7a89", "p3": "#d9822b", "rec2020": "#2b7bb9", "sharp": "#3f9b5b"}
    lights = ["D65", "A (2856K)", "candle (1850K)"]
    ymax = 14.0
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for ax, n in zip(axes, (1, 3)):
        w = 0.1
        for li, ln in enumerate(lights):
            for bi, (b, label) in enumerate(BASES):
                for ki, kind in enumerate(("plain", "pf")):
                    x = li + (bi * 2 + ki - 3.5) * w
                    v = table[("Au", ln, n, b, kind)]
                    ax.bar(x, min(v, ymax), w * 0.92, color=cols[b], hatch="" if kind == "plain" else "///", edgecolor="white" if kind == "plain" else "#222222", lw=0.6,
                           label=(label if (li == 0 and ki == 0) else None))
                    if v > ymax:
                        ax.text(x, ymax + 0.15, f"{v:.0f}", ha="center", fontsize=7, rotation=90, va="bottom")
        ax.axhline(1.0, color="#b03030", ls="--", lw=1, label=t["jnd"])
        ax.set_xticks(range(3)); ax.set_xticklabels(t["lights"])
        ax.set_title(t["nlab"].format(n)); ax.set_ylim(0, ymax + 2.2)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(t["ylabel"])
    from matplotlib.patches import Patch
    h, lab = axes[0].get_legend_handles_labels()
    h += [Patch(fc="#cccccc", ec="white", label=t["plain"]), Patch(fc="#cccccc", ec="#222222", hatch="///", label=t["pf"])]
    lab += [t["plain"], t["pf"]]
    fig.legend(h, lab, loc="lower center", ncol=4, frameon=False, fontsize=8.5)
    fig.suptitle(t["title"], fontsize=11.5)
    fig.tight_layout(rect=(0, 0.1, 1, 0.95))
    fig.savefig(path, dpi=170)
    print("wrote", path)


if __name__ == "__main__":
    main()
