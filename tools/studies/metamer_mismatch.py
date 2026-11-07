"""How far can the color under another light be from the truth when only the color under D65 (or m spectral moments) is known?

An RGB renderer, and every RGB-to-spectrum conversion, knows a surface only through its RGB under D65 (the film sensor of the engine). The true
reflectance is one of all bounded reflectances 0 <= R(lambda) <= 1 with that RGB: a convex set (the "metamer set"). Under another light the members of
the set have different colors (metamer mismatch, Wyszecki & Stiles; Finlayson & Morovic 2005; Logvinenko, Funt & Godau 2014). The spread of that
color is the error that NO conversion can avoid: whichever spectrum it picks from the set, the truth may be another member.

This script measures the spread with linear programs (scipy HiGHS): for many directions d in the film's RGB space under the light L, maximize
d . (color of R under L) over R in the set. The extreme colors found are the vertices of the mismatch body. Reported per (spectrum, light):

  radius   the smallest worst-case dE00 of a single guess: min over the candidate guesses of max over the body of dE00(guess, body color);
           an approximation from above (candidates: the vertices and their centroid). No conversion can beat it without knowing more.
  to truth the largest dE00 between the TRUE color and a point of the body (what an unlucky pick costs).
  the error of the engine's current conversion (rgb_to_spectrum.py) of the same RGB, for comparison.

With --moments the set is instead the bounded reflectances that share the first m real trigonometric moments (cosine moments of the phase of the
wavelength, the mirrored-signal variant of Peters et al. 2019, optionally with their wavelength warp from the reference code) with the truth:
m = 3 is as many numbers as an RGB. The radius as a function of m says how many moments guarantee a color under the other light.

    python tools/studies/metamer_mismatch.py [--moments 3,4,5,6,8,10] [--warp PATH_TO_ReferenceImplementation.py] [--out FILE]

The spectra: the measured cloth of the lab (assets/Metal_measured, not in the repository; skipped when absent) and the Fresnel reflectance of gold
(normal incidence and 60 degrees, one and three bounces). Cloth spectra are cleaned (430-700 nm kept, clipped to [0, 1], 5 nm moving average,
the ends held constant): the instrument's ratio is unreliable outside that range.
"""
import argparse
import glob
import importlib.util
import os
import sys

import numpy as np
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "tools", "data"))
sys.path.insert(0, os.path.join(ROOT, "tools", "lib"))
import metal_accuracy as ma          # noqa: E402  (film sensor, lights, Fresnel of the engine's reference model)
import rgb_to_spectrum as r2s        # noqa: E402  (the engine's RGB -> spectrum conversion)
from colordiff import delta_e2000    # noqa: E402  (vectorized CIEDE2000)
from gamut import lab as xyz_to_lab  # noqa: E402

GRID = np.arange(380.0, 781.0, 5.0)                       # 81 wavelengths: the LP variables
SENS = np.array([ma.SENSOR[int(l - ma.LAM[0])] for l in GRID])        # (81, 3) film sensor weights s_c(lambda)
WEIGHT = np.where((GRID == 380) | (GRID == 780), 0.5, 1.0) * 5.0   # trapezoid, 5 nm steps
D65 = ma.d65(GRID)


def color_rows(light):
    """(3, 81): film RGB of a reflectance R under `light` (a spectrum on ma.LAM, luminance 1) is  rows @ R."""
    L = np.interp(GRID, ma.LAM, light)
    return (SENS * (L * WEIGHT)[:, None]).T


def to_lab(rgb):
    """Film RGB (linear sRGB, light of luminance 1) -> CIELAB with the D65 white at Y = 1; works on (..., 3) arrays."""
    return xyz_to_lab(np.asarray(rgb, dtype=float) @ ma.M_RGB_TO_XYZ.T)


def lab_de(rgb_a, rgb_b):
    return delta_e2000(to_lab(rgb_a), to_lab(rgb_b))


def directions(n, seed=1):
    """Roughly even directions on the sphere (Fibonacci), plus the axes."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    d = np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], axis=1)
    return np.concatenate([d, np.eye(3), -np.eye(3)])


def second_difference():
    """(79, 81): R[i-1] - 2 R[i] + R[i+1], the curvature of a spectrum on the 5 nm grid."""
    n = len(GRID)
    D = np.zeros((n - 2, n))
    for i in range(n - 2):
        D[i, i:i + 3] = (1.0, -2.0, 1.0)
    return D


def mismatch_vertices(A_eq, b_eq, B, n_dir=200, kappa=None):
    """Extreme colors B R of the set {0 <= R <= 1, A_eq R = b_eq, (|curvature of R| <= kappa)}: one LP per direction. Returns (k, 3)."""
    pts = []
    A_ub = b_ub = None
    if kappa is not None:
        D = second_difference()
        A_ub, b_ub = np.vstack([D, -D]), np.full(2 * D.shape[0], kappa)
    for d in directions(n_dir):
        res = linprog(-(d @ B), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=(0.0, 1.0), method="highs")
        if res.status == 0:
            pts.append(B @ res.x)
    return np.unique(np.round(np.array(pts), 9), axis=0)


def guess_radius(vertices, truth):
    """(radius, worst dE00 to the truth, the best guess): candidates are the vertices and their centroid."""
    cands = np.vstack([vertices, vertices.mean(axis=0)])
    lc, lv = to_lab(cands), to_lab(vertices)
    de = delta_e2000(lc[:, None, :], lv[None, :, :])                      # (candidates, vertices)
    worst = de.max(axis=1)
    j = int(np.argmin(worst))
    to_truth = float(delta_e2000(to_lab(truth)[None, :], lv).max())
    return worst[j], to_truth, cands[j]


def cloth_spectra():
    import sr5100
    base = os.path.join(ROOT, "GoLD_rayt", "assets", "Metal_measured", "Cloth_measured", "布_反射率")
    out = {}
    for f in sorted(glob.glob(os.path.join(base, "CSV", "*_flame.csv")) + glob.glob(os.path.join(base, "CSV_2回目", "*_trans.csv"))):
        d = sr5100.read(f)
        for i, s in enumerate(d.spectra):
            v = np.interp(GRID, d.wavelengths, np.nan_to_num(s, nan=0.0))
            keep = (GRID >= 430) & (GRID <= 700)
            v = np.clip(v, 0.0, 1.0)
            v[GRID < 430] = v[keep][0]
            v[GRID > 700] = v[keep][-1]
            v = np.convolve(np.pad(v, 2, mode="edge"), np.ones(5) / 5, mode="valid")
            name = os.path.basename(f).replace("nuno_", "").replace(".csv", "") + (f"#{i + 1}" if d.spectra.shape[0] > 1 else "")
            out[name] = v
    return out


def gold_spectra():
    eta, k = ma.load_metal("Au")
    out = {}
    for theta in (0, 60):
        F = ma.fresnel_conductor(np.cos(np.radians(theta)), eta, k)
        for n in (1, 3):
            out[f"Au theta={theta} n={n}"] = np.clip(np.interp(GRID, ma.LAM, F ** n), 0.0, 1.0)
    return out


def phase_of(wl, warp_path):
    if warp_path:
        spec = importlib.util.spec_from_file_location("peters_ref", warp_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod.WavelengthToPhase(wl, 360.0, 830.0, True, True)           # mirrored signal, with the paper's warp
    return -np.pi * (wl - 380.0) / 400.0                                      # mirrored signal, linear phase in [-pi, 0]


def moment_rows(m, warp_path):
    """(m, 81): real trigonometric moments  integral R(phi) cos(k phi) dphi, k = 0..m-1, of the reflectance as a function of the phase.

    The measure is d(phi), not d(lambda): with the paper's warp the phase runs unevenly over the wavelengths, and the moments of the bounded
    signal are taken over the phase (trapezoid rule on the phase values of the 5 nm grid)."""
    ph = phase_of(GRID, warp_path)
    dphi = np.zeros_like(ph)
    dphi[1:-1] = (ph[2:] - ph[:-2]) / 2.0
    dphi[0], dphi[-1] = (ph[1] - ph[0]) / 2.0, (ph[-1] - ph[-2]) / 2.0
    return np.array([np.cos(k * ph) * dphi for k in range(m)])


def current_conversion_error(spec, rows_L, truth_L, A_rgb):
    """dE00 under the light of the engine's conversion of the spectrum's own RGB (the 'estimate' of the current renderer)."""
    rgb = A_rgb @ spec
    est = np.clip(r2s.reflectance_spectrum(np.clip(rgb, 0, 1), GRID), 0, 1)
    return lab_de(truth_L, rows_L @ est)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--moments", default="", help="comma-separated moment counts m (cosine moments); empty: none")
    ap.add_argument("--smooth", default="", help="comma-separated factors f: also the sets with |curvature| <= f x (the spectrum's own largest curvature)")
    ap.add_argument("--warp", default=None, help="path to Peters et al. ReferenceImplementation.py: use their wavelength warp for the moments")
    ap.add_argument("--dirs", type=int, default=120, help="directions per LP set (more: a finer body)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    spectra = {}
    try:
        spectra.update(cloth_spectra())
    except Exception as e:          # the lab data is not in the repository
        print("[cloth spectra skipped: %s]" % e)
    spectra.update(gold_spectra())
    lts = {k: v[0] for k, v in ma.lights().items() if not k.startswith("D65")}      # the lights other than D65
    A_rgb = color_rows(ma.d65(ma.LAM))                                                # the color under D65 (film sensor)
    A_rgb = A_rgb / A_rgb.sum(axis=1, keepdims=True)                                  # a flat reflectance of 1 gives (1, 1, 1)
    D2 = second_difference()

    # constraint sets: (label, equality rows, curvature factor or None)
    modes = [("D65 color", A_rgb, None)]
    modes += [(f"D65 color + smooth x{f}", A_rgb, float(f)) for f in args.smooth.split(",") if f]
    modes += [(f"{m} moments" + (" (warp)" if args.warp else ""), moment_rows(m, args.warp), None) for m in (int(x) for x in args.moments.split(",") if x)]

    rows = []
    lines = [f"{'spectrum':22s} {'light':16s} {'constraints':26s} {'radius':>7s} {'to truth':>9s} {'engine conv.':>13s}"]
    for name, spec in spectra.items():
        kappa0 = float(np.abs(D2 @ spec).max())
        for lname, light in lts.items():
            B = color_rows(light)
            truth = B @ spec
            for cname, A, f in modes:
                verts = mismatch_vertices(A, A @ spec, B, args.dirs, None if f is None else max(f * kappa0, 1e-6))
                if len(verts) < 2:
                    continue
                rad, to_truth, _ = guess_radius(verts, truth)
                eng = current_conversion_error(spec, B, truth, A_rgb) if cname == "D65 color" else float("nan")
                rows.append((name, lname, cname, rad, to_truth, eng))
                lines.append(f"{name:22s} {lname:16s} {cname:26s} {rad:7.2f} {to_truth:9.2f} {eng:13.2f}")

    # summary: mean radius over the spectra of each class, per constraint set and light
    summary = ["", "Mean radius (dE00) over the spectra of a class; 'engine' = mean dE00 of the engine's current conversion of the same color.",
               f"{'class':8s} {'constraints':26s} " + " ".join(f"{k:>16s}" for k in lts)]
    for cls, pred in (("cloth", lambda n: not n.startswith("Au")), ("gold", lambda n: n.startswith("Au"))):
        for cname, _, _ in modes:
            vals = [np.mean([r[3] for r in rows if pred(r[0]) and r[1] == l and r[2] == cname] or [np.nan]) for l in lts]
            summary.append(f"{cls:8s} {cname:26s} " + " ".join(f"{v:16.2f}" for v in vals))
        eng = [np.mean([r[5] for r in rows if pred(r[0]) and r[1] == l and r[2] == "D65 color"] or [np.nan]) for l in lts]
        summary.append(f"{cls:8s} {'engine conversion':26s} " + " ".join(f"{v:16.2f}" for v in eng))
    text = "\n".join(summary + [""] + lines)
    print("\n".join(summary))
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(text + "\n")


if __name__ == "__main__":
    main()
