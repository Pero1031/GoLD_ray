"""
Reference solution for metal color: spectral (1 nm integral) vs the RGB approximation of the renderer.

What this models (a faithful port of the renderer's own pipeline):
  * Film sensor   : s_c(l) = (M * cmf(l))_c / W_c, cmf = Wyman/Sloan/Shirley 2013 multi-lobe fit,
                    M = the renderer's XYZ->linear sRGB matrix, W = sensor response to D65 (so D65 -> (1,1,1)).
  * Light         : SPD scaled to luminance Y = 1 (EmissionSpectrum::FromCurve); RGB mode keeps only its
                    sensor RGB  L_c = Integral(L s_c).
  * Conductor     : unpolarized Fresnel F(l, theta) from the Johnson table (union grid, linear interpolation,
                    clamped at the ends, like ComplexIorTable).
  * Spectral mode (what the path tracer converges to):  RGB_ref_c = Integral(F^n L s_c dl) over 360..830 nm.
  * RGB mode      : eta/k reduced per channel by integrateIorToRgb_CMF (F0 matched under D65), then
                    RGB_c = L_c * F_c(theta; eta_c, k_c)^n.

Everything uses the film's own CMF, so the differences measured here are the RGB approximation only,
not the (separate) accuracy of the CMF fit.  `--true-cmf` additionally compares to tabulated CIE 1931.

Usage:  python tools/studies/metal_accuracy.py [--true-cmf] [--figs DIR]
        python tools/studies/metal_accuracy.py --compare probe_rgb.log probe_spectral.log   (engine check)
"""

import argparse
import csv
import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ASSETS = os.path.join(ROOT, "GoLD_rayt", "assets")
HEADER_ILLUM = os.path.join(ROOT, "GoLD_rayt", "include", "Color", "Illuminants.hpp")

LAM = np.arange(360, 831, dtype=np.float64)          # film / spectral renderer range, 1 nm
VIS = (LAM >= 380) & (LAM <= 780)                      # range used by the RGB reduction and light RGB

# ---------------------------------------------------------------------------------------------
# Colorimetry (copied from Color/ColorSpace.hpp)
# ---------------------------------------------------------------------------------------------


def _g(l, mu, kl, kr, a):
    d = (l - mu) * np.where(l < mu, kl, kr)
    return a * np.exp(-0.5 * d * d)


def cmf_wyman(l):
    x = _g(l, 442.0, 0.0624, 0.0374, 0.362) + _g(l, 599.8, 0.0264, 0.0323, 1.056) - _g(l, 501.1, 0.0490, 0.0382, 0.065)
    y = _g(l, 568.8, 0.0213, 0.0247, 0.821) + _g(l, 530.9, 0.0613, 0.0322, 0.286)
    z = _g(l, 437.0, 0.0845, 0.0278, 1.217) + _g(l, 459.0, 0.0385, 0.0725, 0.681)
    return np.stack([x, y, z], axis=-1)


M_XYZ_TO_RGB = np.array([[3.2406, -1.5372, -0.4986],
                         [-0.9689, 1.8758, 0.0415],
                         [0.0557, -0.2040, 1.0570]])
M_RGB_TO_XYZ = np.linalg.inv(M_XYZ_TO_RGB)


def load_d65_table():
    text = open(HEADER_ILLUM, encoding="utf-8").read()
    body = text[text.index("kD65_SPD"):]
    body = body[body.index("{"):body.index("};")]
    body = re.sub(r"//[^\n]*", "", body)
    vals = [float(v) for v in re.findall(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?", body)]
    assert len(vals) == 81, len(vals)
    return np.array(vals)


D65_TABLE = load_d65_table()


def d65(l):
    return np.interp(l, np.arange(380, 781, 5), D65_TABLE, left=0.0, right=0.0)


def trapz_weights(mask):
    w = mask.astype(np.float64)
    idx = np.nonzero(mask)[0]
    w[idx[0]] *= 0.5
    w[idx[-1]] *= 0.5
    return w


W_VIS = trapz_weights(VIS)           # 1 nm trapezoid over 380..780
W_FULL = np.ones_like(LAM)           # 360..830 (spectral film)
W_FULL[0] = W_FULL[-1] = 0.5

CMF = cmf_wyman(LAM)                                   # (N, 3)
RGB_CMF = CMF @ M_XYZ_TO_RGB.T                         # (N, 3)  M * cmf
WHITE = (RGB_CMF * (d65(LAM) * W_VIS)[:, None]).sum(0)  # sensor response to D65
SENSOR = RGB_CMF / WHITE                               # s_c(l), D65 -> (1,1,1)


def sensor_rgb(spec, weights=W_FULL):
    """Film RGB of a radiance spectrum sampled on LAM."""
    return (SENSOR * (spec * weights)[:, None]).sum(0)


# ---------------------------------------------------------------------------------------------
# Light sources
# ---------------------------------------------------------------------------------------------


def blackbody(l, kelvin):
    c2 = 1.438777e7
    return 1e15 / (l ** 5 * np.expm1(c2 / (l * kelvin)))


def load_led():
    nm, v = [], []
    with open(os.path.join(ASSETS, "Illuminants", "led_synthetic_coolwhite.csv"), encoding="utf-8") as f:
        for row in csv.reader(f):
            if not row or row[0].startswith("#") or not row[0][0].isdigit():
                continue
            nm.append(float(row[0]))
            v.append(float(row[1]))
    return np.interp(LAM, nm, v, left=0.0, right=0.0)


def make_light(spec):
    """Scale an SPD to luminance Y = 1 (EmissionSpectrum::FromCurve); returns (spd, rgb_for_RGB_mode)."""
    rgb = sensor_rgb(spec, W_VIS)                       # FromCurve integrates the sensor over 380..780
    y = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return spec / y, rgb / y


def lights():
    out = {"D65": make_light(d65(LAM)),
           "A (2856K)": make_light(blackbody(LAM, 2856.0)),
           "candle (1850K)": make_light(blackbody(LAM, 1850.0)),
           "LED cool-white": make_light(load_led())}
    return out


# ---------------------------------------------------------------------------------------------
# Metal tables and Fresnel
# ---------------------------------------------------------------------------------------------


def load_metal(name):
    path = os.path.join(ASSETS, "Metal_data", f"Johnson_{name}.csv")
    sections, cur = {}, None
    for row in csv.reader(open(path, encoding="utf-8")):
        if not row:
            continue
        if row[0] == "wl":
            cur = row[1]
            sections[cur] = []
        else:
            sections[cur].append((float(row[0]) * 1000.0, float(row[1])))   # um -> nm
    n, k = np.array(sections["n"]), np.array(sections["k"])
    grid = np.union1d(n[:, 0], k[:, 0])                  # ComplexIorTable::fromCurves
    eta_g = np.interp(grid, n[:, 0], n[:, 1])
    k_g = np.interp(grid, k[:, 0], k[:, 1])
    return np.interp(LAM, grid, eta_g), np.interp(LAM, grid, k_g)   # clamped at the ends


def fresnel_conductor(cos_i, eta, k):
    """Unpolarized conductor Fresnel (port of fresnel::fresnelConductor)."""
    cos_i = np.clip(np.abs(cos_i), 0.0, 1.0)
    c2 = cos_i ** 2
    s2 = 1.0 - c2
    eta2, k2 = eta ** 2, k ** 2
    t0 = eta2 - k2 - s2
    a2b2 = np.sqrt(t0 * t0 + 4.0 * eta2 * k2)
    t1 = a2b2 + c2
    a = np.sqrt(np.maximum(0.5 * (a2b2 + t0), 0.0))
    t2 = 2.0 * cos_i * a
    rs = (t1 - t2) / (t1 + t2)
    t3 = c2 * a2b2 + s2 * s2
    t4 = t2 * s2
    rp = rs * (t3 - t4) / (t3 + t4)
    return 0.5 * (rs + rp)


def f0(eta, k):
    return ((eta - 1) ** 2 + k ** 2) / ((eta + 1) ** 2 + k ** 2)


def d65_rgb_of(spec_fn_values):
    """spectrumToRGB_D65: sensor RGB of S under D65, divided by the white (S == 1)."""
    w = d65(LAM) * W_VIS
    num = (RGB_CMF * (spec_fn_values * w)[:, None]).sum(0)
    den = (RGB_CMF * w[:, None]).sum(0)
    return num / den


def reduce_to_rgb(eta, k):
    """integrateIorToRgb_CMF: eta_rgb by D65-CMF averaging, k_rgb solved so that F0 matches R's color."""
    r = f0(eta, k)
    r_max = min(r[VIS].max(), 0.9999)
    f0_rgb = d65_rgb_of(r)
    eta_rgb = d65_rgb_of(eta)
    eta_c = np.maximum(eta_rgb, 0.01)
    f0_c = np.clip(f0_rgb, 0.0, r_max)
    k2 = (f0_c * (eta_c + 1) ** 2 - (eta_c - 1) ** 2) / (1 - f0_c)
    k_rgb = np.where(k2 > 0, np.sqrt(np.maximum(k2, 0)), 0.0)
    return eta_c, k_rgb, f0_rgb


def reduce_3wavelength(eta, k, lam=(650.0, 550.0, 450.0)):
    """approxIorToRgb_3Wavelengths: point samples of the table (R, G, B order)."""
    idx = [int(round(l - LAM[0])) for l in lam]
    return eta[idx], k[idx]


# ---------------------------------------------------------------------------------------------
# Rendering models
# ---------------------------------------------------------------------------------------------


def spectral_rgb(eta, k, cos_i, light_spd, bounces=1):
    """What the spectral path tracer converges to for n specular bounces off the same metal."""
    f = fresnel_conductor(cos_i, eta, k) ** bounces
    return sensor_rgb(f * light_spd, W_FULL)


def rgb_model(eta_rgb, k_rgb, cos_i, light_rgb, bounces=1):
    return light_rgb * fresnel_conductor(cos_i, eta_rgb, k_rgb) ** bounces


# ---------------------------------------------------------------------------------------------
# Color difference (CIEDE2000, white = D65 of the output sRGB)
# ---------------------------------------------------------------------------------------------

try:
    import colour
    _XY_D65 = np.array([0.3127, 0.3290])
except ImportError:
    colour = None


def delta_e2000(rgb_a, rgb_b):
    """CIEDE2000 between two linear-sRGB colors (Y=1 is the white of luminance-1 D65). Negative
    channels are allowed (the conversion goes through XYZ, no clipping)."""
    xyz_a = M_RGB_TO_XYZ @ np.asarray(rgb_a)
    xyz_b = M_RGB_TO_XYZ @ np.asarray(rgb_b)
    lab_a = colour.XYZ_to_Lab(xyz_a, _XY_D65)
    lab_b = colour.XYZ_to_Lab(xyz_b, _XY_D65)
    return float(colour.delta_E(lab_a, lab_b, method="CIE 2000"))


def delta_e2000_displayed(rgb_a, rgb_b):
    """CIEDE2000 after clipping both colors to the sRGB gamut [0, 1] (what a PNG / monitor shows)."""
    return delta_e2000(np.clip(rgb_a, 0, 1), np.clip(rgb_b, 0, 1))


def to_srgb8(rgb):
    v = np.clip(rgb, 0, 1)
    v = np.where(v <= 0.0031308, 12.92 * v, 1.055 * v ** (1 / 2.4) - 0.055)
    return np.round(255 * v).astype(int)


def white_balanced(rgb, light_rgb):
    """Film white balance of the renderer (--wb): divide by the light's color (luminance preserved)."""
    return rgb / light_rgb


# ---------------------------------------------------------------------------------------------
# Experiments
# ---------------------------------------------------------------------------------------------

METALS = ["Au", "Ag", "Cu"]


def fmt(v, p=3):
    return "(" + ", ".join(f"{x:.{p}f}" for x in v) + ")"


def exp_reduction(data):
    print("== E1. RGB reduction of the tables (integrateIorToRgb_CMF) ==")
    for m in METALS:
        eta, k = data[m]
        eta_rgb, k_rgb, f0_rgb = reduce_to_rgb(eta, k)
        f0_model = f0(eta_rgb, k_rgb)
        clamp = f0_rgb - np.clip(f0_rgb, 0, min(f0(eta, k)[VIS].max(), 0.9999))
        e3, k3 = reduce_3wavelength(eta, k)
        print(f"{m}: F0 target (spectrum under D65) {fmt(f0_rgb)}   F0 of reduced eta/k {fmt(f0_model)}"
              f"   eta_rgb {fmt(eta_rgb)}   k_rgb {fmt(k_rgb)}   clamped by {fmt(clamp)}")
        print(f"    3-wavelength approx: eta {fmt(e3)} k {fmt(k3)} -> F0 {fmt(f0(e3, k3))}")


def exp_single(data, lts, thetas, bounces=1, quiet=False):
    """Returns rows: metal, light, theta, bounces, ref rgb, model rgb, dE (raw), dE (white balanced)."""
    rows = []
    for m in METALS:
        eta, k = data[m]
        eta_rgb, k_rgb, _ = reduce_to_rgb(eta, k)
        for ln, (spd, lrgb) in lts.items():
            for th in thetas:
                c = np.cos(np.radians(th))
                ref = spectral_rgb(eta, k, c, spd, bounces)
                mod = rgb_model(eta_rgb, k_rgb, c, lrgb, bounces)
                rows.append(dict(metal=m, light=ln, theta=th, bounces=bounces, ref=ref, rgb=mod,
                                 dE=delta_e2000(ref, mod),
                                 dE_clip=delta_e2000_displayed(ref, mod),
                                 dE_wb=delta_e2000(white_balanced(ref, lrgb), white_balanced(mod, lrgb)),
                                 rel=(mod - ref) / ref))
    return rows


def print_rows(rows, title):
    print(f"== {title} ==")
    print(f"{'metal':5s} {'light':15s} {'th':>3s} {'n':>2s} {'spectral (lin)':>26s} {'RGB (lin)':>26s}"
          f" {'dE00':>6s} {'clip':>6s} {'dE00 wb':>8s} {'rel err % (R,G,B)':>24s}")
    for r in rows:
        rel = ", ".join(f"{100 * x:+5.1f}" for x in r["rel"])
        print(f"{r['metal']:5s} {r['light']:15s} {r['theta']:3.0f} {r['bounces']:2d} {fmt(r['ref']):>26s} "
              f"{fmt(r['rgb']):>26s} {r['dE']:6.2f} {r['dE_clip']:6.2f} {r['dE_wb']:8.2f} {rel:>24s}")


def exp_decompose(data, lts, thetas):
    """Per channel: ideal reflectance rho_c = ref_c / L_c  vs  the RGB model's F_c."""
    print("== E4. Decomposition: ideal per-channel reflectance rho_c(theta; light) vs RGB-model F_c(theta) ==")
    for m in ["Au", "Cu"]:
        eta, k = data[m]
        eta_rgb, k_rgb, _ = reduce_to_rgb(eta, k)
        for ln in ["D65", "A (2856K)"]:
            spd, lrgb = lts[ln]
            print(f"{m} / {ln}")
            print("   theta | rho_c (spectral)            | F_c (RGB model)             | RGB - spectral (abs)")
            for th in thetas:
                c = np.cos(np.radians(th))
                rho = spectral_rgb(eta, k, c, spd) / lrgb
                fc = fresnel_conductor(c, eta_rgb, k_rgb)
                print(f"   {th:5.0f} | {fmt(rho, 4):27s} | {fmt(fc, 4):27s} | {fmt(fc - rho, 4)}")


def exp_progression(data, lts):
    """Effective per-bounce reflectance r_c(n) = (ref_c(n) / L_c)^(1/n): constant in n for a flat spectrum.
    The RGB model uses one fixed number per channel (F_c) for every bounce."""
    print("== E8. Effective per-bounce reflectance r_c(n) under D65 (theta = 0): spectral vs the RGB model's fixed F_c ==")
    spd, lrgb = lts["D65"]
    for m in METALS:
        eta, k = data[m]
        eta_rgb, k_rgb, _ = reduce_to_rgb(eta, k)
        fc = fresnel_conductor(1.0, eta_rgb, k_rgb)
        print(f"{m}: RGB model F_c = {fmt(fc, 4)}")
        for n in (1, 2, 3, 4, 6, 8):
            ref = spectral_rgb(eta, k, 1.0, spd, n) / lrgb
            eff = np.sign(ref) * np.abs(ref) ** (1.0 / n)
            print(f"    n={n}: spectral r_c = {fmt(eff, 4)}   spectral value {fmt(ref, 4)}   RGB value {fmt(fc ** n, 4)}")


def exp_blackbody(data, kelvins, thetas=(0.0, 60.0)):
    print("== E5. Blackbody sweep (single reflection) ==")
    print(f"{'metal':5s} {'K':>6s} {'th':>3s} {'dE00':>6s} {'dE00 wb':>8s} {'rel err % (R,G,B)':>24s} {'spectral (lin)':>26s}")
    for m in METALS:
        eta, k = data[m]
        eta_rgb, k_rgb, _ = reduce_to_rgb(eta, k)
        for kel in kelvins:
            spd, lrgb = make_light(blackbody(LAM, float(kel)))
            for th in thetas:
                c = np.cos(np.radians(th))
                ref = spectral_rgb(eta, k, c, spd)
                mod = rgb_model(eta_rgb, k_rgb, c, lrgb)
                rel = ", ".join(f"{100 * x:+5.1f}" for x in (mod - ref) / ref)
                print(f"{m:5s} {kel:6.0f} {th:3.0f} {delta_e2000(ref, mod):6.2f} "
                      f"{delta_e2000(white_balanced(ref, lrgb), white_balanced(mod, lrgb)):8.2f} {rel:>24s} {fmt(ref):>26s}")


def exp_true_cmf(data):
    """How far is the film's CMF fit from tabulated CIE 1931 2deg (D65 from colour, exact sRGB matrix)?"""
    print("== E7. Absolute check: film CMF (Wyman fit) vs tabulated CIE 1931 / D65 ==")
    cmfs = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"]
    d65t = colour.SDS_ILLUMINANTS["D65"]
    lam = np.arange(380, 781, 1.0)
    xyz_true = cmfs.copy().align(colour.SpectralShape(380, 780, 1)).values
    d65v = d65t.copy().align(colour.SpectralShape(380, 780, 1)).values
    wts = trapz_weights(np.ones_like(lam, dtype=bool))
    m_exact = colour.models.RGB_COLOURSPACE_sRGB.matrix_XYZ_to_RGB
    wh = (xyz_true * (d65v * wts)[:, None]).sum(0)
    wh_rgb = m_exact @ wh
    print(f"(white point of the table-based chain in linear sRGB: {fmt(wh_rgb / wh_rgb[1], 4)}; exact would be (1,1,1))")
    for m in METALS:
        eta, k = data[m]
        r = f0(eta, k)[(LAM >= 380) & (LAM <= 780)]
        # film-convention result (Wyman fit, renormalized by its own white)
        film = d65_rgb_of(f0(eta, k))
        # tabulated chain, normalized by its own white (so D65 -> (1,1,1)) like the film
        true = (m_exact @ (xyz_true * (r * d65v * wts)[:, None]).sum(0)) / wh_rgb
        print(f"{m}: F0 color  film {fmt(film, 4)}   tabulated CIE {fmt(true, 4)}   dE00 {delta_e2000(film, true):.2f}"
              f"   sRGB8 film {to_srgb8(film)} tab {to_srgb8(true)}")


LIGHT_KEYS = {"d65": "D65", "tungsten": "A (2856K)", "candle": "candle (1850K)", "led": "LED cool-white"}


def parse_probe_log(path):
    """Reads the output of `main_cli --scene metalprobe` -> (metal, light, roughness, mode, [(mean cos, rgb)])."""
    metal = light = mode = None
    rough = 0.01
    rows = []
    for line in open(path, encoding="utf-8", errors="replace"):
        m = re.search(r"\[MetalProbe\] (\w+) under (\w+)(?: r=([\d.]+))? \[(\w+)", line)
        if m:
            metal, light, mode = m.group(1), m.group(2), m.group(4).lower()
            rough = float(m.group(3)) if m.group(3) else 0.01
        m = re.search(r"mean cos ([\d.]+)\s+RGB \(([-\d.]+), ([-\d.]+), ([-\d.]+)\)", line)
        if m:
            rows.append((float(m.group(1)), np.array([float(m.group(i)) for i in (2, 3, 4)])))
    return metal, light, rough, mode, rows


def compare_probe(paths, data, lts):
    print("== Engine (metalprobe render) vs this script's prediction ==")
    logs = {}
    summary = []
    for path in paths:
        metal, light, rough, mode, rows = parse_probe_log(path)
        if not rows:
            print(f"{path}: no probe rows found")
            continue
        logs[(metal, light, rough, mode)] = rows
        if rough > 0.05:    # the prediction below is the smooth (single bounce) model
            continue
        eta, k = data[metal]
        eta_rgb, k_rgb, _ = reduce_to_rgb(eta, k)
        spd, lrgb = lts[LIGHT_KEYS[light]]
        print(f"{metal} under {light}, engine mode {mode}")
        print("   theta   engine RGB                  prediction                  rel diff % (R,G,B)        dE00")
        worst = 0.0
        for cos_t, rgb in rows:
            pred = rgb_model(eta_rgb, k_rgb, cos_t, lrgb) if mode == "rgb" else spectral_rgb(eta, k, cos_t, spd)
            rel = ", ".join(f"{100 * x:+6.2f}" for x in (rgb - pred) / pred)
            theta = np.degrees(np.arccos(min(cos_t, 1.0)))
            de = delta_e2000(rgb, pred)
            if theta < 80:      # the 85 deg bin averages a steep F(theta) over +-1.5 deg: not a point value
                worst = max(worst, de)
            print(f"   {theta:5.1f}   {fmt(rgb, 4):27s} {fmt(pred, 4):27s} {rel:26s} {de:5.2f}")
        summary.append((metal, light, mode, worst))

    print()
    print("== Summary: engine vs prediction (worst dE00 over theta < 80 deg) ==")
    for metal, light, mode, worst in summary:
        print(f"   {metal} {light:9s} {mode:9s} {worst:5.2f}")

    print()
    print("== Engine RGB build vs engine spectral build (dE00 per angle of incidence) ==")
    for (metal, light, rough, mode), rows in sorted(logs.items()):
        other = logs.get((metal, light, rough, "spectral"))
        if mode != "rgb" or other is None:
            continue
        des = [delta_e2000(a[1], b[1]) for a, b in zip(rows, other) if abs(a[0] - b[0]) < 1e-3]
        thetas = [np.degrees(np.arccos(min(a[0], 1.0))) for a, b in zip(rows, other) if abs(a[0] - b[0]) < 1e-3]
        cells = "  ".join(f"{t:4.0f}:{d:5.2f}" for t, d in zip(thetas, des))
        print(f"   {metal} {light:9s} r={rough:4.2f} {cells}")


def make_figures(outdir, data, lts):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(outdir, exist_ok=True)
    cols = {"Au": "#c9a227", "Ag": "#8a8f98", "Cu": "#b8693d"}
    chan = ("#d62728", "#2ca02c", "#1f77b4")

    # Fig 1: reflectance spectra, the film sensor, and the flat per-channel value the RGB model uses
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    for ax, m in zip(axes, METALS):
        eta, k = data[m]
        eta_rgb, k_rgb, _ = reduce_to_rgb(eta, k)
        ax.plot(LAM, f0(eta, k), color=cols[m], lw=2.2, label="R(l) from the table")
        fc = fresnel_conductor(1.0, eta_rgb, k_rgb)
        spans = ((560, 780), (480, 620), (380, 520))
        for c, name in enumerate("RGB"):
            ax.hlines(fc[c], *spans[c], colors=chan[c], lw=2.5, linestyles="--", label=f"RGB model {name} = {fc[c]:.3f}")
            s = SENSOR[:, c] / np.abs(SENSOR).max() * 0.35
            ax.fill_between(LAM, 0.0, np.clip(s, 0, None), color=chan[c], alpha=0.12)
        ax.set_title(m)
        ax.set_xlabel("wavelength (nm)")
        ax.set_xlim(360, 830)
        ax.set_ylim(0, 1.05)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="lower right")
    axes[0].set_ylabel("normal-incidence reflectance")
    fig.suptitle("Metal reflectance vs the flat per-channel value of the RGB model (shaded: positive lobes of the film sensor)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "metal_reflectance_spectra.png"), dpi=130)
    plt.close(fig)

    # Fig 2: color error vs number of bounces (D65) and vs light color temperature (one bounce)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.8))
    for m in METALS:
        eta, k = data[m]
        eta_rgb, k_rgb, _ = reduce_to_rgb(eta, k)
        spd, lrgb = lts["D65"]
        ns = np.arange(1, 9)
        de = [delta_e2000(spectral_rgb(eta, k, 1.0, spd, n), rgb_model(eta_rgb, k_rgb, 1.0, lrgb, n)) for n in ns]
        a1.plot(ns, de, "o-", color=cols[m], label=m)
        kel = np.array([1500, 1850, 2200, 2856, 3500, 4000, 5000, 6500, 8000, 10000], dtype=float)
        de2 = []
        for kk in kel:
            sp, lr = make_light(blackbody(LAM, kk))
            de2.append(delta_e2000(spectral_rgb(eta, k, 1.0, sp), rgb_model(eta_rgb, k_rgb, 1.0, lr)))
        a2.plot(kel, de2, "o-", color=cols[m], label=m)
    a1.set_xlabel("number of specular bounces (D65, normal incidence)")
    a1.set_ylabel("CIEDE2000, RGB model vs spectral")
    a2.set_xscale("log")
    a2.set_xlabel("blackbody temperature (K), one bounce")
    for a in (a1, a2):
        a.axhline(1, color="k", lw=0.6, ls=":")
        a.grid(alpha=0.3)
        a.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "metal_rgb_error.png"), dpi=130)
    plt.close(fig)

    # Fig 3: swatches as displayed (clipped, sRGB encoded): spectral on top, RGB model below
    cases = [("D65", 1), ("D65", 3), ("A (2856K)", 1), ("candle (1850K)", 1)]
    fig, ax = plt.subplots(figsize=(9, 3.2))
    x = 0.0
    for m in METALS:
        eta, k = data[m]
        eta_rgb, k_rgb, _ = reduce_to_rgb(eta, k)
        for ln, n in cases:
            spd, lrgb = lts[ln]
            ref = spectral_rgb(eta, k, 1.0, spd, n)
            mod = rgb_model(eta_rgb, k_rgb, 1.0, lrgb, n)
            for row, v in ((1, ref), (0, mod)):
                ax.add_patch(plt.Rectangle((x, row), 0.9, 0.9, color=to_srgb8(v) / 255.0))
            ax.text(x + 0.45, -0.25, f"{m}\n{ln.split(' ')[0]}\nn={n}", ha="center", va="top", fontsize=6.5)
            ax.text(x + 0.45, 2.0, f"dE {delta_e2000_displayed(ref, mod):.1f}", ha="center", fontsize=6.5)
            x += 1
        x += 0.4
    ax.set_xlim(-0.2, x)
    ax.set_ylim(-1.3, 2.4)
    ax.text(-0.25, 1.45, "spectral", ha="right", fontsize=7)
    ax.text(-0.25, 0.45, "RGB", ha="right", fontsize=7)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "metal_swatches.png"), dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--true-cmf", action="store_true", help="also compare with tabulated CIE 1931 (needs `colour`)")
    ap.add_argument("--figs", default=None, help="directory for PNG figures")
    ap.add_argument("--compare", nargs="+", default=None, metavar="LOG",
                    help="only compare logs of `main_cli --scene metalprobe` with the prediction")
    args = ap.parse_args()
    if colour is None:
        sys.exit("this script needs the `colour-science` package for CIEDE2000")

    data = {m: load_metal(m) for m in METALS}
    lts = lights()
    if args.compare:
        compare_probe(args.compare, data, lts)
        return
    for n, (spd, lrgb) in lts.items():
        print(f"light {n:15s} sensor RGB (Y=1): {fmt(lrgb, 4)}")
    print()

    exp_reduction(data)
    print()
    angles = [0, 30, 45, 60, 70, 80, 85]
    rows = exp_single(data, {"D65": lts["D65"]}, angles)
    print_rows(rows, "E2. Single reflection under D65 (the case the reduction is built for)")
    print()
    rows = exp_single(data, lts, [0, 60])
    print_rows(rows, "E3. Single reflection, all lights, theta 0 / 60 deg")
    print()
    exp_decompose(data, lts, [0, 30, 60, 80])
    print()
    exp_blackbody(data, [1500, 1850, 2856, 4000, 5000, 6500, 10000])
    print()
    exp_progression(data, lts)
    print()
    rows = []
    for n in (1, 2, 3, 4, 6):
        rows += exp_single(data, {k: lts[k] for k in ("D65", "A (2856K)")}, [0, 45], bounces=n)
    print_rows(rows, "E6. n specular bounces (same metal each time)")
    if args.true_cmf:
        print()
        exp_true_cmf(data)
    if args.figs:
        make_figures(args.figs, data, lts)
        print(f"\nfigures written to {args.figs}")


if __name__ == "__main__":
    main()
