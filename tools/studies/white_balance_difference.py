"""
The RGB-vs-spectral difference after white balance: does balancing to the light hide it, show it, or enlarge it?

Sections 3.5-3.8 of docs/METAL_ACCURACY.md compare the two renders with the light's color cast left in (no white balance). A photograph
or a film is usually balanced to the light. The renderer's --wb (CIE film) is CAT16 from the light's white to D65 of the same luminance,
a 3x3 matrix on the linear film values (Color/EmissionSpectrum.hpp setWhiteBalance); D < 1 adapts only partly (Color/ChromaticAdaptation.hpp).
This script applies that matrix to both results and measures the displayed CIEDE2000, through the renderer's own output (Reinhard,
clip) and through the view transforms of view_transform_hiding.py (Blender 4.5 and ACES 1.3, with OpenColorIO).

Part A  reference model of working_space_hiding.py (gold, normal incidence, 9 lights x 1-6 reflections), RGB renderer in sRGB primaries
        (and ACEScg primaries), white balance with D = 1, 0.8, 0.6 and none.
Part B  why: the two results in CAT16's cone-like space (L, M, S) before the balance, the balance's gains, and the CIEDE2000 after the
        balance when only S (or only L and M) is taken from the RGB result.
Part C  gold bowl renders of the current builds (goldbowl_<light>_<rgb|spectral>_<a|b>.pfm without --wb, their .log next to them, from
        which the light's film RGB is read); the matrix is applied afterwards (equal to rendering with --wb within 0.1 %, checked on
        small renders with the same seed), so one set of renders serves every D.

  python tools/studies/white_balance_difference.py --ocio-python <python with PyOpenColorIO> [--renders <folder>]
         [--out docs/data/white_balance_difference_results.txt]
"""
import argparse
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
import gamut as gm                                            # noqa: E402
from colordiff import delta_e2000                             # noqa: E402
import metal_accuracy as ma                                   # noqa: E402
import working_space_hiding as wsh                            # noqa: E402
import view_transform_hiding as vth                           # noqa: E402

M16 = ma.M16
DEGREES = (None, 0.6, 0.8, 1.0)          # None = no white balance
SHOW = [("Reinhard + clip (renderer)", "srgb"), ("Reinhard + clip (renderer)", "p3"), ("Blender Standard", "srgb"),
        ("Blender AgX", "srgb"), ("Blender AgX", "p3"), ("ACES 1.3 SDR", "srgb"), ("ACES 1.3 SDR", "p3")]


def cat16_xyz(src_white, degree):
    """XYZ -> XYZ: CAT16 from src_white to D65 of the same Y (the renderer's --wb), degree D."""
    dst = wsh.rb.WHITE_XYZ * (src_white[1] / wsh.rb.WHITE_XYZ[1])
    scale = degree * (M16 @ dst) / (M16 @ src_white) + (1.0 - degree)
    return np.linalg.inv(M16) @ np.diag(scale) @ M16


def dlabel(d):
    return "no WB" if d is None else f"D = {d:g}"


def part_a(ocio, emit):
    M = wsh.Model("Au")
    lights = [(f"{c} K", wsh.light_spd(c)) for c in wsh.CCTS] + [("D65", wsh.light_spd("D65"))]
    conds = [(ln, spd, n) for ln, spd in lights for n in wsh.BOUNCES]
    idx = {(ln, n): i for i, (ln, _, n) in enumerate(conds)}
    picks = [("1850 K", 1, "candle x1"), ("1500 K", 1, "1500 K x1"), ("2856 K", 1, "tungsten x1"), ("1850 K", 3, "candle x3"), ("D65", 3, "D65 x3")]
    for basis in ("srgb", "acescg"):
        base = [(M.xyz_pair(spd, n, basis), wsh.rb.xyz_of(spd)) for _, spd, n in conds]
        emit(f"=== A ({wsh.BASES[basis]} working basis). Gold, reference model, 45 conditions, exposure 0 ===")
        emit("Columns: CIEDE2000 between the RGB and the spectral result after the same white balance; 'out' = share of the 45 balanced spectral colors outside sRGB.")
        for d in DEGREES:
            ref, rgb, lin = [], [], []
            for (r, g), lw in base:
                A = np.eye(3) if d is None else cat16_xyz(lw, d)
                ref.append(wsh.to_srgb(A @ r)); rgb.append(wsh.to_srgb(A @ g)); lin.append(wsh.de(A @ r, A @ g))
            ref, rgb = np.array(ref), np.array(rgb)
            codes = vth.all_codes(ocio, np.concatenate([ref, rgb]), 0.0)
            dd = {k: delta_e2000(vth.displayed_lab(v[:len(ref)], k[1]), vth.displayed_lab(v[len(ref):], k[1])) for k, v in codes.items()}
            emit(f"--- {dlabel(d)} (out of sRGB: {np.mean(ref.min(axis=1) < 0) * 100:.0f} %) ---")
            emit(f"{'display':>10s} {'view':>27s} " + " ".join(f"{p[2]:>11s}" for p in picks) + f" {'mean':>6s}")
            emit(f"{'':>10s} {'(linear, no display)':>27s} " + " ".join(f"{lin[idx[(ln, n)]]:11.2f}" for ln, n, _ in picks) + f" {np.mean(lin):6.2f}")
            for label, disp in SHOW:
                v = dd[(label, disp)]
                emit(f"{vth.DISPLAY_NAME[disp]:>10s} {label:>27s} " + " ".join(f"{v[idx[(ln, n)]]:11.2f}" for ln, n, _ in picks) + f" {v.mean():6.2f}")
            if d == 1.0 and basis == "srgb":
                i = idx[("1850 K", 1)]
                emit(f"  candle x1 after the balance, linear sRGB: spectral ({ref[i][0]:.4f}, {ref[i][1]:.4f}, {ref[i][2]:.4f}), RGB ({rgb[i][0]:.4f}, {rgb[i][1]:.4f}, {rgb[i][2]:.4f})")
                i = idx[("2856 K", 1)]
                emit(f"  tungsten x1 after the balance, linear sRGB: spectral ({ref[i][0]:.4f}, {ref[i][1]:.4f}, {ref[i][2]:.4f}), RGB ({rgb[i][0]:.4f}, {rgb[i][1]:.4f}, {rgb[i][2]:.4f})")
        emit()


def part_b(emit):
    M = wsh.Model("Au")
    emit("=== B. Why: the two results in CAT16's cone-like space (sRGB working basis, one reflection, normal incidence) ===")
    emit("LMS of the result divided by the light's LMS (= the effective reflectance each cone sees); gain = the balance's scale (D = 1);")
    emit("'S only' = dE after the balance when only S is taken from the RGB result (L, M from the spectral one); 'L,M only' = the reverse.")
    emit(f"{'light':>8s} {'spectral L/M/S':>24s} {'RGB L/M/S':>24s} {'gain L/M/S':>21s} {'dE balanced':>12s} {'S only':>7s} {'L,M only':>9s} {'unbalanced':>11s}")
    for ln, cct in (("1500 K", 1500), ("1850 K", 1850), ("2200 K", 2200), ("2856 K", 2856), ("4500 K", 4500), ("D65", "D65")):
        spd = wsh.light_spd(cct)
        r, g = M.xyz_pair(spd, 1, "srgb")
        lw = wsh.rb.xyz_of(spd)
        lms_w, lms_r, lms_g = M16 @ lw, M16 @ r, M16 @ g
        dst = wsh.rb.WHITE_XYZ * (lw[1] / wsh.rb.WHITE_XYZ[1])
        gain = (M16 @ dst) / lms_w
        A = cat16_xyz(lw, 1.0)
        inv = np.linalg.inv(M16)
        s_only = inv @ np.array([lms_r[0], lms_r[1], lms_g[2]])
        lm_only = inv @ np.array([lms_g[0], lms_g[1], lms_r[2]])
        f = lambda v: "/".join(f"{x:.3f}" for x in v)
        emit(f"{ln:>8s} {f(lms_r / lms_w):>24s} {f(lms_g / lms_w):>24s} {f(gain):>21s} {wsh.de(A @ r, A @ g):12.2f} "
             f"{wsh.de(A @ r, A @ s_only):7.2f} {wsh.de(A @ r, A @ lm_only):9.2f} {wsh.de(r, g):11.2f}")
    emit()
    # gold's reflectance as each cone sees it under the light: the spectral result weights F by cone x light, the RGB renderer by
    # the sRGB channel sensors under D65 (its eta/k reduction) and then mixes the channels through the light's RGB
    cmf = wsh.rb.CMF
    lms_sens = cmf @ M16.T                                      # (N, 3) cone-like sensors
    F = ma.fresnel_conductor(1.0, M.eta, M.k)
    for ln, cct in (("1850 K", 1850), ("D65", "D65")):
        spd = wsh.light_spd(cct)
        w = spd * wsh.W_VIS
        cone_peak = wsh.LAM[np.argmax(lms_sens[:, 2])]
        s_weight = lms_sens[:, 2] * w
        mean_nm = (wsh.LAM * s_weight).sum() / s_weight.sum()
        emit(f"  {ln}: the S sensor x light has its weighted mean at {mean_nm:.0f} nm (S peak {cone_peak:.0f} nm); gold's F there {np.interp(mean_nm, wsh.LAM, F):.3f}")
    emit("Which working-RGB channels carry the candle's (1850 K) S stimulus, and gold's S reflectance in each working basis (one reflection):")
    spd = wsh.light_spd(1850)
    lw = wsh.rb.xyz_of(spd)
    for basis in ("srgb", "rec2020", "acescg"):
        light_rgb = (M.S[basis] * (spd * wsh.W_VIS)[:, None]).sum(0)
        contrib = (M16 @ np.linalg.inv(M.T[basis]))[2] * light_rgb
        r, g = M.xyz_pair(spd, 1, basis)
        f0 = ma.fresnel_conductor(1.0, *M.red[basis])
        emit(f"  {wsh.BASES[basis]:>16s}: light ({', '.join(f'{v:.4f}' for v in light_rgb)}), S share R/G/B {'/'.join(f'{v * 100:.0f}%' for v in contrib / contrib.sum())}, "
             f"gold F0 ({', '.join(f'{v:.3f}' for v in f0)}), S reflectance spectral {(M16 @ r)[2] / (M16 @ lw)[2]:.3f} / RGB {(M16 @ g)[2] / (M16 @ lw)[2]:.3f}")
    emit()


def read_light_rgb(log_path):
    with open(log_path, encoding="utf-8", errors="replace") as f:
        m = re.search(r"sensor RGB of the light = \(([-0-9.e]+), ([-0-9.e]+), ([-0-9.e]+)\)", f.read())
    return np.array([float(m.group(i)) for i in (1, 2, 3)])


def part_c(ocio, renders, emit):
    from hdr_io import read_hdr
    from colordiff import M_RGB_TO_XYZ
    from gold_compare_linear import REGION_BOXES
    from gold_compare_linear_fig import box_blur
    emit("=== C. Gold bowl renders, current builds (800x450, 256 spp, two renders per build with seeds 1 and 2, raw .pfm, 3x3 blur, bowl box, lit pixels; 8-bit) ===")
    emit("The balance is applied to the renders afterwards (the same matrix as --wb, computed from the light's film RGB in the render's log).")
    for light in ("candle", "tungsten"):
        rep = {(m, r): box_blur(read_hdr(os.path.join(renders, f"goldbowl_{light}_{m}_{r}.pfm")).astype(np.float64)) for m in ("rgb", "spectral") for r in "ab"}
        light_rgb = read_light_rgb(os.path.join(renders, f"goldbowl_{light}_rgb_a.log"))
        h, w = rep[("rgb", "a")].shape[:2]

        def box(name):
            y0, y1, x0, x1 = REGION_BOXES[name]
            return slice(int(y0 / 450 * h), int(y1 / 450 * h)), slice(int(x0 / 800 * w), int(x1 / 800 * w))

        spe0 = (rep[("spectral", "a")] + rep[("spectral", "b")]) / 2
        rgb0 = (rep[("rgb", "a")] + rep[("rgb", "b")]) / 2
        wy = (M_RGB_TO_XYZ @ rgb0[box("white ball")].reshape(-1, 3).mean(0))[1]
        mask = np.zeros((h, w), bool)
        mask[box("bowl total")] = True
        mask &= (0.2126 * spe0[..., 0] + 0.7152 * spe0[..., 1] + 0.0722 * spe0[..., 2]) > 0.02 * wy
        emit(f"--- {light} (film RGB of the light {light_rgb.round(4).tolist()}; {int(mask.sum())} pixels) ---")
        emit(f"{'balance':>8s} {'display':>10s} {'view':>27s} {'dE00 mean':>9s} {'95 pct':>7s} {'>1':>6s} {'noise RGB':>10s} {'noise spec':>10s} {'out sRGB':>9s}")
        for d in DEGREES:
            Mf = np.eye(3) if d is None else ma.cat16_matrix(light_rgb, d)
            imgs = [x @ Mf.T for x in (rgb0, spe0, rep[("rgb", "a")], rep[("rgb", "b")], rep[("spectral", "a")], rep[("spectral", "b")])]
            out_s = ((imgs[1] @ gm.from_srgb("srgb").T).min(axis=2) < 0)[mask].mean() * 100
            codes = vth.all_codes(ocio, np.stack(imgs), 0.0)
            for label, disp in SHOW:
                c = codes[(label, disp)]
                lab = [vth.displayed_lab(c[i], disp, round8=True) for i in range(6)]
                dE = delta_e2000(lab[0], lab[1])[mask]
                emit(f"{dlabel(d):>8s} {vth.DISPLAY_NAME[disp]:>10s} {label:>27s} {dE.mean():9.2f} {np.percentile(dE, 95):7.2f} {(dE > 1).mean() * 100:5.1f}% "
                     f"{delta_e2000(lab[2], lab[3])[mask].mean():10.2f} {delta_e2000(lab[4], lab[5])[mask].mean():10.2f} {out_s:8.0f}%")
            # the white ball after the balance (should be neutral in both builds)
            wb = [imgs[i][box("white ball")].reshape(-1, 3).mean(0) for i in (0, 1)]
            emit(f"{'':>8s} white ball, linear sRGB: RGB build ({wb[0][0]:.3f}, {wb[0][1]:.3f}, {wb[0][2]:.3f}), spectral ({wb[1][0]:.3f}, {wb[1][1]:.3f}, {wb[1][2]:.3f})")
        emit()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ocio-python", required=True)
    ap.add_argument("--blender-config", default="")
    ap.add_argument("--renders", default="")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    config = a.blender_config or os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(a.ocio_python)))),
                                              "datafiles", "colormanagement", "config.ocio")
    ocio = vth.Ocio(a.ocio_python, config)
    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit("The RGB-vs-spectral difference after white balance (CAT16 from the light's white to D65, as the renderer's --wb). Written by")
    emit("tools/studies/white_balance_difference.py. Views as in view_transform_hiding.py; CIEDE2000 with the display's white (D65).")
    emit()
    part_a(ocio, emit)
    part_b(emit)
    if a.renders:
        part_c(ocio, a.renders, emit)
    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        with open(a.out, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(lines) + "\n")
        print("wrote", a.out)


if __name__ == "__main__":
    main()
