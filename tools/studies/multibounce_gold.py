"""
Colour of rough gold with each multiple-bounce model (docs/MULTIBOUNCE.md).

Input: the CSVs of `GoLD_rayt.exe --scene multibounce-gold` from both builds (run from GoLD_rayt/):
    ../out/studies/multibounce/gold_albedo_spectral.csv   directional albedo per wavelength (380..810 nm, 10 nm)
    ../out/studies/multibounce/gold_albedo_rgb.csv        the same per channel of the RGB build
The albedo is the light a gold surface sends back into the whole hemisphere for light arriving at theta_i (0 or 60 degrees),
as estimated by RoughConductor's own sample(): single bounce only, Kulla-Conty lobe, Zhang et al. 2026, Heitz et al. 2016.

Colour of an albedo under an illuminant: the XYZ of R(l) S(l) (CIE 1931 2 deg), white = the illuminant itself (Y = 1), CIELAB
relative to that white. The RGB build: linear sRGB albedo times the illuminant's linear sRGB (what the RGB renderer does), to XYZ
with the sRGB matrix, same white. The reference is Heitz's model in the spectral build.

Usage:  python tools/studies/multibounce_gold.py [--figs DIR]
"""

import argparse
import csv
import os

import colour
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "out", "studies", "multibounce")
CMFS = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"]
ILLUM = {"D65": colour.SDS_ILLUMINANTS["D65"], "A": colour.SDS_ILLUMINANTS["A"]}
SRGB = colour.RGB_COLOURSPACES["sRGB"]
M_RGB2XYZ = SRGB.matrix_RGB_to_XYZ
M_XYZ2RGB = SRGB.matrix_XYZ_to_RGB
MODELS = ["single", "kulla-conty", "zhang", "heitz"]


def read(path):
    rows = {}
    with open(path, newline="") as f:
        r = csv.reader(f)
        head = next(r)
        for row in r:
            rows[(row[0], float(row[1]), float(row[2]))] = np.array([float(v) for v in row[3:]])
    return head[3:], rows


def spectral_xyz(wl, refl, illum):
    """XYZ of a reflectance on the grid wl (nm) under an illuminant, white Y = 1 (box bands of 10 nm, 1 nm integration)."""
    lam = np.arange(380, 781)
    r = np.interp(lam, wl, refl)
    s = illum.copy().align(colour.SpectralShape(380, 780, 1)).values
    cmf = CMFS.copy().align(colour.SpectralShape(380, 780, 1)).values
    k = 1.0 / np.sum(s * cmf[:, 1])
    return k * np.sum((r * s)[:, None] * cmf, axis=0)


def illum_rgb(illum):
    """The illuminant's linear sRGB as the RGB renderer would use it (XYZ with Y = 1, then the sRGB matrix)."""
    return M_XYZ2RGB @ spectral_xyz(np.array([380.0, 780.0]), np.array([1.0, 1.0]), illum)


def lab(xyz, white):
    return colour.XYZ_to_Lab(xyz, colour.XYZ_to_xy(white))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figs", default=DATA)
    args = ap.parse_args()

    head, spec = read(os.path.join(DATA, "gold_albedo_spectral.csv"))
    wl = np.array([float(h) for h in head])
    _, rgb = read(os.path.join(DATA, "gold_albedo_rgb.csv"))
    alphas = sorted({k[1] for k in spec})
    thetas = sorted({k[2] for k in spec})

    colours = {}   # (source, model, alpha, theta, illum) -> (XYZ, Lab)
    for name, illum in ILLUM.items():
        white = spectral_xyz(wl, np.ones_like(wl), illum)
        lrgb = illum_rgb(illum)
        white_rgb = M_RGB2XYZ @ lrgb
        for key, refl in spec.items():
            xyz = spectral_xyz(wl, refl, illum)
            colours[("spectral",) + key + (name,)] = (xyz, lab(xyz, white))
        for key, a in rgb.items():
            xyz = M_RGB2XYZ @ (a * lrgb)
            colours[("rgb",) + key + (name,)] = (xyz, lab(xyz, white_rgb))

    def de(a, b):
        return float(colour.delta_E(a[1], b[1], method="CIE 2000"))

    # ---- tables -------------------------------------------------------------------------------------------------------------
    for name in ILLUM:
        for th in thetas:
            print(f"\nCIEDE2000 against Heitz (spectral), illuminant {name}, theta_i = {th:g} deg")
            print("alpha  | single  kulla-conty  zhang  | RGB build: heitz  kulla-conty | L* a* b* of the reference")
            for a in alphas:
                ref = colours[("spectral", "heitz", a, th, name)]
                cells = [de(colours[("spectral", m, a, th, name)], ref) for m in ("single", "kulla-conty", "zhang")]
                cells += [de(colours[("rgb", m, a, th, name)], ref) for m in ("heitz", "kulla-conty")]
                L, A_, B = ref[1]
                print(f"{a:4.1f}   | {cells[0]:6.2f}  {cells[1]:10.2f}  {cells[2]:6.2f} | {cells[3]:15.2f}  {cells[4]:10.2f}  | {L:5.1f} {A_:5.1f} {B:5.1f}")

    # saturation (chroma) and lightness of the reference against Kulla-Conty
    print("\nChroma C*ab and lightness L* at theta_i = 0 under A: Heitz / Kulla-Conty (spectral)")
    for a in alphas:
        h = colours[("spectral", "heitz", a, 0.0, "A")][1]
        k = colours[("spectral", "kulla-conty", a, 0.0, "A")][1]
        print(f"{a:4.1f}  C* {np.hypot(h[1], h[2]):5.1f} / {np.hypot(k[1], k[2]):5.1f}   L* {h[0]:5.1f} / {k[0]:5.1f}")

    # ---- figure -------------------------------------------------------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(args.figs, exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    styles = {"single": ("#999999", "--"), "kulla-conty": ("#d9822b", "-"), "zhang": ("#2b6cd9", "-"), "heitz": ("#111111", ":")}
    for m in MODELS:
        c, ls = styles[m]
        ax[0].plot(wl, spec[(m, 1.0, 0.0)], color=c, ls=ls, lw=2, label=m)
    ax[0].set_title("albedo of gold, GGX alpha = 1, normal incidence")
    ax[0].set_xlabel("wavelength (nm)")
    ax[0].set_xlim(380, 780)
    ax[0].legend()
    for th, mk in zip(thetas, ("o", "s")):
        for m in ("single", "kulla-conty", "zhang"):
            c, _ = styles[m]
            ax[1].plot(alphas, [de(colours[("spectral", m, a, th, "A")], colours[("spectral", "heitz", a, th, "A")]) for a in alphas],
                       color=c, marker=mk, ms=4, lw=1.5, label=f"{m}, theta {th:g}")
        ax[1].plot(alphas, [de(colours[("rgb", "heitz", a, th, "A")], colours[("spectral", "heitz", a, th, "A")]) for a in alphas],
                   color="#b03030", marker=mk, ms=4, lw=1.5, label=f"Heitz in the RGB build, theta {th:g}")
    ax[1].set_title("CIEDE2000 against Heitz (spectral), illuminant A")
    ax[1].set_xlabel("GGX alpha")
    ax[1].legend(fontsize=7)
    # swatches: rows = models (spectral) + Heitz RGB build, columns = alpha, colour under A shown white-balanced to D65 (Bradford)
    rows = [("spectral", m) for m in MODELS] + [("rgb", "heitz")]
    img = np.zeros((len(rows), len(alphas), 3))
    white_a = spectral_xyz(wl, np.ones_like(wl), ILLUM["A"])
    for i, (src, m) in enumerate(rows):
        for j, a in enumerate(alphas):
            xyz = colours[(src, m, a, 0.0, "A")][0]
            if src == "rgb":
                xyz = xyz * (white_a[1] / (M_RGB2XYZ @ illum_rgb(ILLUM["A"]))[1])
            xyz_d65 = colour.chromatic_adaptation(xyz, colour.xy_to_XYZ(colour.XYZ_to_xy(white_a)),
                                                  colour.xy_to_XYZ(SRGB.whitepoint), method="Von Kries", transform="Bradford")
            img[i, j] = np.clip(colour.cctf_encoding(np.clip(M_XYZ2RGB @ xyz_d65, 0, None)), 0, 1)
    ax[2].imshow(img, aspect="auto", interpolation="nearest")
    ax[2].set_yticks(range(len(rows)), [m if s == "spectral" else "heitz (RGB build)" for s, m in rows])
    ax[2].set_xticks(range(len(alphas)), [f"{a:g}" for a in alphas])
    ax[2].set_xlabel("GGX alpha")
    ax[2].set_title("albedo colour under A (adapted to D65, sRGB)")
    fig.tight_layout()
    out = os.path.join(args.figs, "multibounce_gold.png")
    fig.savefig(out, dpi=130)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
