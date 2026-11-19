"""
Figures and a table for the RGB-vs-spectral / display-gamut study, made from the recorded results in docs/data
(no rendering): gold_compare_gamut_results.txt (bowl), erato_compare_gamut_results.txt (statue, raw renders = part B),
metal_accuracy_results.txt (single and multiple reflections of gold).

  python tools/studies/gamut_figures.py <out_dir> [--lang en|ja]

Writes in <out_dir>:
  gamut_de00_bars.png    CIEDE2000 between the RGB and the spectral image, as shown on an sRGB / Display P3 / Rec.2020 display
  gamut_de00_table.png   the same numbers as a table (mean, share above 1, noise floor)
  gamut_chromaticity.png chromaticity of gold (spectral vs RGB) under D65 / tungsten / candle, with the three gamuts

The CIE 1931 observer for the spectral locus is the multi-lobe fit of Wyman, Sloan & Shirley (2013), the same fit the renderer uses.
"""
import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "docs", "data")
LIGHTS = ["d65", "tungsten", "candle"]
DISPLAYS = ["sRGB", "Display P3", "Rec.2020"]

TXT = {
    "en": dict(light={"d65": "D65 (daylight)", "tungsten": "Tungsten (2856 K)", "candle": "Candle (1850 K)"},
               obj={"bowl": "Gold bowl", "statue": "Gold statue"},
               bars_title="Color difference between RGB and spectral rendering, as shown on each display",
               ylabel="CIEDE2000 (mean over lit pixels)", jnd="visible-difference level (about 1)", noise="noise floor (same build, two renders)",
               tab_cols=["Light", "Object", "sRGB", "Display P3", "Rec.2020", "noise"],
               tab_note="mean CIEDE2000 (share of pixels above 1). Raw renders, 256 spp, mean of two renders per build, 3x3 blur.",
               chroma_title="Gold under three lights: spectral vs RGB renderer", spec="spectral", rgbm="RGB", x="x", y="y",
               one="1 reflection", three="3 reflections", locus="spectral locus"),
    "ja": dict(light={"d65": "D65（昼光）", "tungsten": "タングステン（2856 K）", "candle": "ろうそく（1850 K）"},
               obj={"bowl": "金の器", "statue": "金の像"},
               bars_title="RGB と分光の描画の色差（各ディスプレイで表示した画像）",
               ylabel="CIEDE2000（明るい画素の平均）", jnd="差が見える目安（約 1）", noise="ノイズの下限（同じ版を2回描いた差）",
               tab_cols=["光源", "対象", "sRGB", "Display P3", "Rec.2020", "ノイズ"],
               tab_note="ΔE00 の平均（1 を超える画素の割合）。デノイズなしの描画、256 spp、各版 2 回の平均、3×3 ぼかし。",
               chroma_title="3 つの光源の下の金：分光（正しい色）と RGB", spec="分光", rgbm="RGB", x="x", y="y",
               one="1 回反射", three="3 回反射", locus="スペクトル軌跡"),
}


def setup_font(lang):
    if lang == "ja":
        for name in ("Yu Gothic", "Meiryo", "BIZ UDGothic", "MS Gothic"):
            if any(f.name == name for f in font_manager.fontManager.ttflist):
                plt.rcParams["font.family"] = name
                return
        print("no Japanese font found", file=sys.stderr)
    plt.rcParams["axes.unicode_minus"] = False


def parse_compare(path, raw_section_marker=None):
    """{(light, display): dict(mean, gt23, noise)} from a *_compare_gamut_results.txt; the section after `raw_section_marker` if given."""
    text = open(path, encoding="utf-8").read()
    if raw_section_marker:
        text = text[text.index(raw_section_marker):]
    out = {}
    for line in text.splitlines():
        m = re.match(r"^(d65|tungsten|candle)\s+(sRGB|Display P3|Rec\.2020)\s+(.*)$", line)
        if not m:
            continue
        nums = [float(v.strip("%")) for v in re.findall(r"[-+]?\d+\.?\d*%?", m.group(3))]
        if len(nums) == 7:      # mean, p95, >1, >2.3, outside, noise rgb, noise spec
            mean, gt23, nr, ns = nums[0], nums[2], nums[5], nums[6]
        elif len(nums) == 8:    # mean, median, p95, >1, >2.3, outside, noise rgb, noise spec
            mean, gt23, nr, ns = nums[0], nums[3], nums[6], nums[7]
        else:
            continue
        out[(m.group(1), m.group(2))] = dict(mean=mean, gt23=gt23, noise=0.5 * (nr + ns))
    return out


def figure_bars(data, t, path):
    colors = {"sRGB": "#6c7a89", "Display P3": "#d9822b", "Rec.2020": "#2b7bb9"}
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), sharey=True)
    objs = ["bowl", "statue"]
    w = 0.25
    for ax, light in zip(axes, LIGHTS):
        for i, d in enumerate(DISPLAYS):
            vals = [data[o][(light, d)]["mean"] for o in objs]
            ax.bar(np.arange(2) + (i - 1) * w, vals, w, color=colors[d], label=d)
            for x, v in zip(np.arange(2) + (i - 1) * w, vals):
                ax.text(x, v + 0.08, f"{v:.1f}", ha="center", va="bottom", fontsize=8)
        for k, o in enumerate(objs):
            n = data[o][(light, "sRGB")]["noise"]
            ax.plot([k - 1.5 * w, k + 1.5 * w], [n, n], color="k", lw=1.2, ls=":", label=t["noise"] if (light == LIGHTS[0] and k == 0) else None)
        ax.axhline(1.0, color="#b03030", lw=1, ls="--", label=t["jnd"] if light == LIGHTS[0] else None)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([t["obj"][o] for o in objs])
        ax.set_title(t["light"][light], fontsize=11)
        ax.set_ylim(0, 6.3)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(t["ylabel"])
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=5, frameon=False, fontsize=9)
    fig.suptitle(t["bars_title"], fontsize=12)
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    fig.savefig(path, dpi=170)
    plt.close(fig)


def figure_table(data, t, path):
    rows = []
    for light in LIGHTS:
        for o in ("bowl", "statue"):
            cells = [t["light"][light] if o == "bowl" else "", t["obj"][o]]
            for d in DISPLAYS:
                r = data[o][(light, d)]
                cells.append(f"{r['mean']:.2f} ({r['gt23']:.0f} %)")
            cells.append(f"{data[o][(light, 'sRGB')]['noise']:.2f}")
            rows.append(cells)
    fig, ax = plt.subplots(figsize=(10.5, 3.3))
    ax.axis("off")
    tb = ax.table(cellText=rows, colLabels=t["tab_cols"], loc="upper center", cellLoc="center",
                  colWidths=[0.22, 0.14, 0.16, 0.16, 0.16, 0.1])
    tb.auto_set_font_size(False)
    tb.set_fontsize(10)
    tb.scale(1, 1.55)
    for (r, c), cell in tb.get_celld().items():
        cell.set_edgecolor("#bbbbbb")
        if r == 0:
            cell.set_facecolor("#2f3b4a"); cell.get_text().set_color("white"); cell.get_text().set_weight("bold")
        elif c in (2, 3, 4):
            v = float(rows[r - 1][c].split()[0])
            cell.set_facecolor("#fbe3c8" if v >= 1.0 else "white")
        if r > 0 and (r - 1) % 2 == 1:
            cell.set_edgecolor("#888888")
    ax.text(0.5, 0.02, t["tab_note"], transform=ax.transAxes, ha="center", fontsize=8.5, color="#444444")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


# --- chromaticity diagram -----------------------------------------------------------------------------------------
def wyman_xyz(lam):
    def g(l, mu, s1, s2):
        s = np.where(l < mu, s1, s2)
        return np.exp(-0.5 * ((l - mu) / s) ** 2)
    x = 1.056 * g(lam, 599.8, 37.9, 31.0) + 0.362 * g(lam, 442.0, 16.0, 26.7) - 0.065 * g(lam, 501.1, 20.4, 26.2)
    y = 0.821 * g(lam, 568.8, 46.9, 40.5) + 0.286 * g(lam, 530.9, 16.3, 31.1)
    z = 1.217 * g(lam, 437.0, 11.8, 36.0) + 0.681 * g(lam, 459.0, 26.0, 13.8)
    return x, y, z


M_SRGB = np.array([[0.4124564, 0.3575761, 0.1804375], [0.2126729, 0.7151522, 0.0721750], [0.0193339, 0.1191920, 0.9503041]])


def lin_to_xy(rgb):
    X, Y, Z = M_SRGB @ np.asarray(rgb)
    s = X + Y + Z
    return X / s, Y / s


def gold_points(path):
    """(light, n) -> (spectral linear sRGB, RGB linear sRGB) for gold at 0 degrees, from E6 of metal_accuracy_results.txt."""
    text = open(path, encoding="utf-8").read()
    text = text[text.index("== E3."):text.index("== E4.")] + text[text.index("== E6."):text.index("== E7.")]   # n = 1 (all lights), n = 2.. (D65, tungsten)
    out = {}
    for line in text.splitlines():
        m = re.match(r"^Au\s+(D65|A \(2856K\)|candle \(1850K\))\s+0\s+(\d+)\s+\(([^)]*)\)\s+\(([^)]*)\)", line)
        if m:
            key = ({"D65": "d65", "A (2856K)": "tungsten", "candle (1850K)": "candle"}[m.group(1)], int(m.group(2)))
            out[key] = (np.array([float(v) for v in m.group(3).split(",")]), np.array([float(v) for v in m.group(4).split(",")]))
    return out


def figure_chromaticity(t, path):
    gp = gold_points(os.path.join(DATA, "metal_accuracy_results.txt"))
    lam = np.linspace(400, 700, 301)
    X, Y, Z = wyman_xyz(lam)
    lx, ly = X / (X + Y + Z), Y / (X + Y + Z)
    gamuts = {"sRGB": ([0.64, 0.30, 0.15], [0.33, 0.60, 0.06], "#6c7a89"),
              "Display P3": ([0.680, 0.265, 0.150], [0.320, 0.690, 0.060], "#d9822b"),
              "Rec.2020": ([0.708, 0.170, 0.131], [0.292, 0.797, 0.046], "#2b7bb9")}
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2), gridspec_kw=dict(width_ratios=[1, 1.15]))
    for ai, ax in enumerate(axes):
        ax.plot(np.r_[lx, lx[0]], np.r_[ly, ly[0]], color="#999999", lw=1, label=t["locus"])
        for name, (px, py, col) in gamuts.items():
            ax.plot(px + [px[0]], py + [py[0]], color=col, lw=1.6, label=name)
        for light, col in zip(LIGHTS, ("#555555", "#c97a10", "#b03030")):
            for n, mk in ((1, "o"), (3, "s")):
                if (light, n) not in gp:
                    continue
                s, r = gp[(light, n)]
                xs, ys = lin_to_xy(s)
                xr, yr = lin_to_xy(r)
                ax.annotate("", xy=(xr, yr), xytext=(xs, ys), arrowprops=dict(arrowstyle="->", color=col, lw=1))
                ax.plot(xs, ys, mk, color=col, ms=6, label=(f"{t['light'][light]}, {t['one'] if n == 1 else t['three']}: {t['spec']}" if ai == 1 else None))
                ax.plot(xr, yr, mk, mfc="white", mec=col, ms=6)
        ax.set_xlabel(t["x"]); ax.set_ylabel(t["y"])
        ax.set_aspect("equal")
        ax.grid(alpha=0.25)
    axes[0].set_xlim(0, 0.8); axes[0].set_ylim(0, 0.9)
    axes[1].set_xlim(0.46, 0.70); axes[1].set_ylim(0.34, 0.46)
    axes[1].set_aspect("auto")
    axes[0].add_patch(plt.Rectangle((0.46, 0.34), 0.24, 0.12, fill=False, ec="k", lw=0.8, ls=":"))
    axes[1].legend(fontsize=7.5, loc="upper right", framealpha=0.9)
    fig.suptitle(t["chroma_title"] + "  (" + ("filled: " if t["spec"] == "spectral" else "塗り: ") + t["spec"] + ("; open: " if t["spec"] == "spectral" else "、白抜き: ") + t["rgbm"] + ")", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--lang", default="en", choices=["en", "ja"])
    a = ap.parse_args()
    setup_font(a.lang)
    t = TXT[a.lang]
    os.makedirs(a.out_dir, exist_ok=True)
    data = {"bowl": parse_compare(os.path.join(DATA, "gold_compare_gamut_results.txt")),
            "statue": parse_compare(os.path.join(DATA, "erato_compare_gamut_results.txt"), "##### B.")}
    for o in data:
        assert len(data[o]) == 9, (o, len(data[o]))
    figure_bars(data, t, os.path.join(a.out_dir, "gamut_de00_bars.png"))
    figure_table(data, t, os.path.join(a.out_dir, "gamut_de00_table.png"))
    figure_chromaticity(t, os.path.join(a.out_dir, "gamut_chromaticity.png"))
    print("wrote 3 figures to", a.out_dir)


if __name__ == "__main__":
    main()
