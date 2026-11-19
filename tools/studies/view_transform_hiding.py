"""
Does an sRGB display still hide the RGB-vs-spectral difference when the image goes through a real view transform?

The earlier results (docs/METAL_ACCURACY.md 3.5 and 3.7) show the renderer's own output: Reinhard per channel, negative channels clipped.
Production images go through a view transform instead. This script shows the same pairs of colors (the RGB renderer's and the spectral
result) through the view transforms of Blender 4.5 (Standard, Filmic, AgX, Khronos PBR Neutral; on Display P3: Standard, AgX) and the
ACES 1.3 SDR output transform (OCIO's built-in CG config; sRGB and Display P3), with OpenColorIO itself (tools/studies/ocio_apply.py, run
in a Python that has PyOpenColorIO, e.g. Blender's), and measures CIEDE2000 between the two displayed colors (display-encoded values
decoded with the sRGB piecewise curve and the display's primaries; white = the display's white, D65).

Part A  the reference model of working_space_hiding.py (gold, normal incidence, 9 lights x 1-6 reflections), RGB renderer in sRGB
        primaries (and in ACEScg primaries for comparison); no 8-bit rounding.
Part B  the gold bowl renders of 3.5 (RGB and spectral builds, two renders each, raw .pfm), bowl pixels, 8-bit rounding, with the
        noise floor (two renders of the same build).

  python tools/studies/view_transform_hiding.py --ocio-python <python with PyOpenColorIO> [--blender-config <config.ocio>]
         [--renders <folder of goldbowl_<light>_<rgb|spectral>_<a|b>.pfm>] [--out docs/data/view_transform_hiding_results.txt]
"""
import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
import gamut as gm                                            # noqa: E402
from colordiff import delta_e2000                             # noqa: E402
import working_space_hiding as wsh                            # noqa: E402

# (label, OCIO key or "reinhard", display gamut)
PIPELINES = [
    ("Reinhard + clip (renderer)", "reinhard", "srgb"),
    ("Blender Standard", "blender|sRGB|Standard", "srgb"),
    ("Blender Filmic", "blender|sRGB|Filmic", "srgb"),
    ("Blender AgX", "blender|sRGB|AgX", "srgb"),
    ("Khronos PBR Neutral", "blender|sRGB|Khronos PBR Neutral", "srgb"),
    ("ACES 1.3 SDR", "aces|sRGB - Display|ACES 1.0 - SDR Video", "srgb"),
    ("Reinhard + clip (renderer)", "reinhard", "p3"),
    ("Blender Standard", "blender|Display P3|Standard", "p3"),
    ("Blender AgX", "blender|Display P3|AgX", "p3"),
    ("ACES 1.3 SDR", "aces|Display P3 - Display|ACES 1.0 - SDR Video", "p3"),
]
FAMILIES = ["Reinhard + clip (renderer)", "Blender Standard", "Blender AgX", "ACES 1.3 SDR"]     # those with an sRGB and a P3 version
DISPLAY_NAME = {"srgb": "sRGB", "p3": "Display P3"}


class Ocio:
    def __init__(self, python, config):
        self.python, self.config = python, config
        self.tmp = tempfile.mkdtemp(prefix="ocio_")

    def apply(self, lin709, ev=0.0):
        """lin709: (..., 3) scene-linear Rec.709. Returns {key: display-encoded values clipped to [0, 1]} for every OCIO pipeline."""
        inp, out = os.path.join(self.tmp, "in.npy"), os.path.join(self.tmp, "out.npz")
        np.save(inp, lin709.astype(np.float32))
        keys = ";".join(k for _, k, _ in PIPELINES if k != "reinhard")
        subprocess.run([self.python, "-I", os.path.join(HERE, "ocio_apply.py"), inp, out, "--blender-config", self.config,
                        "--ev", str(ev), "--keys", keys], check=True, stdout=subprocess.DEVNULL)
        with np.load(out) as d:
            return {k: d[k].astype(np.float64) for k in d.files}


def reinhard(lin709, display, ev=0.0):
    """The renderer's PNG: display primaries, negatives to 0, exposure, Reinhard per channel, sRGB curve (code values in [0, 1])."""
    x = np.clip(lin709 @ gm.from_srgb(display).T, 0.0, None) * 2.0 ** ev
    return gm.oetf(x / (1 + x))


def displayed_lab(code, display, round8=False):
    if round8:
        code = np.round(code * 255.0) / 255.0
    return gm.lab(gm.eotf(code) @ gm.rgb_to_xyz(display).T)


def all_codes(ocio, lin709, ev=0.0):
    """{(label, display): code values} for every pipeline."""
    oc = ocio.apply(lin709, ev)
    return {(label, disp): (reinhard(lin709, disp, ev) if key == "reinhard" else oc[key]) for label, key, disp in PIPELINES}


def part_a(ocio, emit):
    M = wsh.Model("Au")
    lights = [(f"{c} K", wsh.light_spd(c)) for c in wsh.CCTS] + [("D65", wsh.light_spd("D65"))]
    conds = [(ln, spd, n) for ln, spd in lights for n in wsh.BOUNCES]
    picks = [("1850 K", 1, "candle x1"), ("1500 K", 1, "1500 K x1"), ("2856 K", 1, "tungsten x1"), ("1850 K", 3, "candle x3"), ("D65", 1, "D65 x1")]

    for basis in ("srgb", "acescg"):
        ref, rgb, lin_de = [], [], []
        for ln, spd, n in conds:
            r, g = M.xyz_pair(spd, n, basis)
            ref.append(wsh.to_srgb(r)); rgb.append(wsh.to_srgb(g)); lin_de.append(wsh.de(r, g))
        ref, rgb = np.array(ref), np.array(rgb)
        out_srgb = ref.min(axis=1) < 0
        evs = (-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0)
        res = {}
        for ev in evs:
            codes = all_codes(ocio, np.concatenate([ref, rgb]), ev)
            res[ev] = {k: delta_e2000(displayed_lab(v[:len(ref)], k[1]), displayed_lab(v[len(ref):], k[1])) for k, v in codes.items()}
        d0 = res[0.0]
        idx = {(ln, n): i for i, (ln, _, n) in enumerate(conds)}
        vis = d0[("Reinhard + clip (renderer)", "p3")] >= wsh.VISIBLE

        emit(f"=== A ({wsh.BASES[basis]} working basis). Gold, reference model, 45 conditions, exposure 0 ===")
        emit(f"'below 1' counts the {int(vis.sum())} conditions that are visible (>= {wsh.VISIBLE:g}) on P3 with the renderer's own output (section 3.7).")
        emit(f"{'display':>10s} {'view':>27s} " + " ".join(f"{p[2]:>11s}" for p in picks) + f" {'mean':>6s} {'below 1':>8s}")
        for label, _, disp in PIPELINES:
            d = d0[(label, disp)]
            emit(f"{DISPLAY_NAME[disp]:>10s} {label:>27s} " + " ".join(f"{d[idx[(ln, n)]]:11.2f}" for ln, n, _ in picks)
                 + f" {d.mean():6.2f} {int((d[vis] < 1).sum()):>4d}/{int(vis.sum())}")
        emit(f"{'(linear, no display)':>38s} " + " ".join(f"{lin_de[idx[(ln, n)]]:11.2f}" for ln, n, _ in picks) + f" {np.mean(lin_de):6.2f}")
        emit()
        emit("The view transforms are not linear, so the result depends on the exposure (scene-linear multiplier 2^EV before the view;")
        emit("at EV 0 the light has Y = 1, so a white diffuser reads 1 and candle x1 about 0.86; half the lit bowl pixels sit 1.3 to 2.2 stops below that).")
        emit("Candle x1 at EV " + " / ".join(f"{e:+g}" for e in evs) + ":")
        for label, _, disp in PIPELINES:
            emit(f"  {DISPLAY_NAME[disp]:>10s} {label:27s} " + " ".join(f"{res[ev][(label, disp)][idx[('1850 K', 1)]]:5.2f}" for ev in evs))
        emit()
        emit("Share of the P3 difference that the sRGB version of the same view removes (mean of 1 - dE sRGB / dE P3 over the conditions whose")
        emit("spectral color is outside sRGB), at EV " + " / ".join(f"{e:+g}" for e in evs) + ":")
        for fam in FAMILIES:
            vals = []
            for ev in evs:
                s_, p_ = res[ev][(fam, "srgb")], res[ev][(fam, "p3")]
                m = out_srgb & (p_ > 0.05)
                vals.append(np.mean(1 - s_[m] / p_[m]) * 100)
            emit(f"  {fam:27s} " + " ".join(f"{v:5.0f}%" for v in vals))
        emit()
        if basis == "srgb":
            emit("Display-encoded values of candle x1 (spectral | RGB renderer), R G B in [0, 1]:")
            codes = all_codes(ocio, np.stack([ref[idx[("1850 K", 1)]], rgb[idx[("1850 K", 1)]]]), 0.0)
            for label, _, disp in PIPELINES:
                c = codes[(label, disp)]
                emit(f"  {DISPLAY_NAME[disp]:>10s} {label:27s} ({c[0, 0]:.3f}, {c[0, 1]:.3f}, {c[0, 2]:.3f}) | ({c[1, 0]:.3f}, {c[1, 1]:.3f}, {c[1, 2]:.3f})")
            emit()


def part_b(ocio, renders, emit):
    from hdr_io import read_hdr
    from colordiff import M_RGB_TO_XYZ
    from gold_compare_linear import REGION_BOXES
    from gold_compare_linear_fig import box_blur
    emit("=== B. Gold bowl renders (800x450, 256 spp, two renders per build, 3x3 blur, bowl box, lit pixels; 8-bit rounding) ===")
    emit("dE00 between the displayed RGB and spectral images; noise = the same build rendered twice (RGB / spectral).")
    for light in ("d65", "tungsten", "candle"):
        rep = {(m, r): box_blur(read_hdr(os.path.join(renders, f"goldbowl_{light}_{m}_{r}.pfm")).astype(np.float64)) for m in ("rgb", "spectral") for r in "ab"}
        rgb = (rep[("rgb", "a")] + rep[("rgb", "b")]) / 2
        spe = (rep[("spectral", "a")] + rep[("spectral", "b")]) / 2
        h, w = rgb.shape[:2]

        def box(name):
            y0, y1, x0, x1 = REGION_BOXES[name]
            return slice(int(y0 / 450 * h), int(y1 / 450 * h)), slice(int(x0 / 800 * w), int(x1 / 800 * w))

        wy = (M_RGB_TO_XYZ @ rgb[box("white ball")].reshape(-1, 3).mean(0))[1]
        mask = np.zeros((h, w), bool)
        mask[box("bowl total")] = True
        mask &= (0.2126 * spe[..., 0] + 0.7152 * spe[..., 1] + 0.0722 * spe[..., 2]) > 0.02 * wy
        stack = np.stack([rgb, spe, rep[("rgb", "a")], rep[("rgb", "b")], rep[("spectral", "a")], rep[("spectral", "b")]])
        codes = all_codes(ocio, stack, 0.0)
        emit(f"--- {light} ({int(mask.sum())} pixels; spectral outside sRGB: {((spe @ gm.from_srgb('srgb').T).min(axis=2) < 0)[mask].mean() * 100:.0f} %) ---")
        emit(f"{'display':>10s} {'view':>27s} {'dE00 mean':>9s} {'95 pct':>7s} {'>1':>6s} {'noise RGB':>10s} {'noise spec':>10s}")
        for label, _, disp in PIPELINES:
            c = codes[(label, disp)]
            lab = [displayed_lab(c[i], disp, round8=True) for i in range(6)]
            d = delta_e2000(lab[0], lab[1])[mask]
            n_rgb = delta_e2000(lab[2], lab[3])[mask].mean()
            n_spe = delta_e2000(lab[4], lab[5])[mask].mean()
            emit(f"{DISPLAY_NAME[disp]:>10s} {label:>27s} {d.mean():9.2f} {np.percentile(d, 95):7.2f} {(d > 1).mean() * 100:5.1f}% {n_rgb:10.2f} {n_spe:10.2f}")
        emit()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ocio-python", required=True, help="a Python with PyOpenColorIO, e.g. <blender>/4.5/python/bin/python.exe")
    ap.add_argument("--blender-config", default="", help="Blender's config.ocio (default: next to the Python of --ocio-python)")
    ap.add_argument("--renders", default="", help="folder with goldbowl_<light>_<rgb|spectral>_<a|b>.pfm (part B is skipped without it)")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    config = a.blender_config or os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(a.ocio_python)))),
                                              "datafiles", "colormanagement", "config.ocio")
    ocio = Ocio(a.ocio_python, config)
    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit("View transforms and the sRGB 'hiding' of the RGB-vs-spectral difference. Written by tools/studies/view_transform_hiding.py.")
    emit("Views: Blender 4.5's config.ocio (Standard, Filmic, AgX, Khronos PBR Neutral) and OCIO's built-in ACES 1.3 CG config (SDR Video), applied with")
    emit("OpenColorIO (tools/studies/ocio_apply.py). Input = the renderer's film, scene-linear Rec.709. Display values decoded with the sRGB curve and the")
    emit("display's primaries; CIEDE2000 with the display's white (D65, Y = 1). 'Reinhard + clip' is the renderer's own PNG output (sections 3.5 and 3.7).")
    emit()
    part_a(ocio, emit)
    if a.renders:
        part_b(ocio, a.renders, emit)
    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        with open(a.out, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(lines) + "\n")
        print("wrote", a.out)


if __name__ == "__main__":
    main()
