"""
Is "an sRGB display hides the RGB-vs-spectral difference" a property of the RGB renderer's working space, and why does it happen?

Reference model of tools/studies/metal_accuracy.py and condition_map.py: the spectral result is the 1 nm integral of F(lambda)^n L(lambda)
against the observer; the RGB renderer reduces eta/k and the light to three channels of its working basis and multiplies per channel.
Normal incidence, the same metal at every bounce. The display pipeline is the renderer's PNG output (exposure 0): XYZ -> display primaries ->
negative channels clipped to 0 (or OKLCH chroma reduction) -> Reinhard per channel -> CIEDE2000 on the displayed color.

Part A  counts, for each working basis, the conditions (light x reflections) in which the difference is visible on a Display P3 display
        (CIEDE2000 >= 2) but not on an sRGB display with clipping (< 1), and gives the share of the P3 difference that the sRGB clip removes.
Part B  splits the error of the sRGB-basis renderer into its linear-sRGB channels: the CIEDE2000 (linear, no display) of the spectral color
        against the spectral color with only the blue channel (or only red and green) taken from the RGB result. If the error sits in the
        blue channel and blue is the channel that is negative (clipped), per-channel clipping removes the error together with the blue.
Part C  is Part B for the Rec.2020 and ACEScg bases, whose results are shown on the same sRGB display.

ACEScg here means the AP1 primaries with the white scaled to D65 (as a D65 renderer that adopts the AP1 primaries would use them); the
ACES D60 white and its adaptation are not modeled. The Sharp matrix is from memory of Ward & Eydelberg-Vileshin 2002 (see rgb_basis_study.py).

  python tools/studies/working_space_hiding.py [--out docs/data/working_space_hiding_results.txt]
"""
import argparse
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

LAM, W_VIS = ma.LAM, ma.W_VIS
AP1 = [(0.713, 0.293), (0.165, 0.830), (0.128, 0.044)]       # ACES AP1 (ACEScg) primaries
BASES = {"srgb": "sRGB (Rec.709)", "p3": "Display P3", "rec2020": "Rec.2020", "acescg": "ACEScg (AP1)", "sharp": "Sharp"}
CCTS = [1500, 1850, 2200, 2856, 3500, 4500, 5500, 6500]
BOUNCES = (1, 2, 3, 4, 6)
VISIBLE, HIDDEN = 2.0, 1.0      # "visible on P3" and "hidden on sRGB" thresholds (CIEDE2000; about 1 = just visible)


def basis_matrix(name):
    """T: XYZ (D65 white, Y = 1) -> working RGB with D65 -> (1, 1, 1)."""
    if name == "acescg":
        p = np.array([[x / y, 1.0, (1 - x - y) / y] for x, y in AP1]).T      # columns = primaries' XYZ (Y = 1)
        m = np.linalg.inv(p)
        return np.diag(1.0 / (m @ rb.WHITE_XYZ)) @ m
    return rb.basis_matrix(name)


def light_spd(cct):
    spd = ma.blackbody(LAM, float(cct)) if cct != "D65" else ma.d65(LAM)
    return spd / ((rb.CMF[:, 1] * spd * W_VIS).sum() / rb.YN)          # Y = 1 relative to D65


class Model:
    def __init__(self, metal="Au"):
        self.eta, self.k = ma.load_metal(metal)
        self.T = {b: basis_matrix(b) for b in BASES}
        self.S = {b: rb.sensors(self.T[b]) for b in BASES}
        self.red = {b: rb.reduce_to_rgb(self.S[b], self.eta, self.k) for b in BASES}

    def xyz_pair(self, spd, n, basis):
        """(spectral XYZ, RGB renderer's XYZ)."""
        ref = rb.xyz_of(ma.fresnel_conductor(1.0, self.eta, self.k) ** n * spd)
        eta_c, k_c = self.red[basis]
        light_rgb = (self.S[basis] * (spd * W_VIS)[:, None]).sum(0)
        rgb = light_rgb * ma.fresnel_conductor(1.0, eta_c, k_c) ** n
        return ref, np.linalg.inv(self.T[basis]) @ rgb


def shown(xyz, out, mapping="clip"):
    x = (np.linalg.inv(gm.rgb_to_xyz(out)) @ xyz)[None, :]
    x = map_gamut(x, out, mapping)
    x = x / (1 + x)
    return (x @ gm.rgb_to_xyz(out).T)[0]


def de(xyz_a, xyz_b):
    return float(delta_e2000(gm.lab(xyz_a), gm.lab(xyz_b)))


def to_srgb(xyz):
    return np.linalg.inv(gm.rgb_to_xyz("srgb")) @ xyz


def from_srgb(rgb):
    return gm.rgb_to_xyz("srgb") @ rgb


def channel_split(ref, rgb):
    """Linear CIEDE2000 of the spectral color against itself with only blue / only red+green taken from the RGB result (linear sRGB)."""
    s, r = to_srgb(ref), to_srgb(rgb)
    only_b = np.array([s[0], s[1], r[2]])
    only_rg = np.array([r[0], r[1], s[2]])
    return de(ref, from_srgb(only_b)), de(ref, from_srgb(only_rg))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    M = Model("Au")
    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    lights = [(f"{c} K", light_spd(c)) for c in CCTS] + [("D65", light_spd("D65"))]
    emit("Working space and the sRGB 'hiding' of the RGB-vs-spectral difference (gold, reference model, normal incidence).")
    emit("CIEDE2000 between the RGB renderer and the spectral result. Display: primaries, negatives clipped to 0, Reinhard per channel, exposure 0.")
    emit("Written by tools/studies/working_space_hiding.py; see its docstring for the model.")
    emit()

    # ---------------- A ----------------
    per_basis = {}
    for metal in ("Au", "Cu"):
        Mm = M if metal == "Au" else Model(metal)
        emit(f"=== A ({metal}). Hidden on sRGB (clip, dE < {HIDDEN:g}) but visible on Display P3 (dE >= {VISIBLE:g}); {len(lights)} lights x {len(BOUNCES)} reflection counts ===")
        emit("'removed by sRGB clip' = mean of 1 - dE(sRGB clip) / dE(P3) over the conditions whose spectral color is outside sRGB. Means over all conditions.")
        emit(f"{'working basis':>16s} {'hidden / visible on P3':>23s} {'removed by sRGB clip':>21s} {'mean dE linear':>15s} {'P3':>6s} {'sRGB clip':>10s} {'sRGB oklch':>11s}")
        for b, label in BASES.items():
            rows = []
            for lname, spd in lights:
                for n in BOUNCES:
                    ref, rgb = Mm.xyz_pair(spd, n, b)
                    rows.append(dict(light=lname, n=n, lin=de(ref, rgb), clip=de(shown(ref, "srgb"), shown(rgb, "srgb")),
                                     oklch=de(shown(ref, "srgb", "oklch"), shown(rgb, "srgb", "oklch")), p3=de(shown(ref, "p3"), shown(rgb, "p3")),
                                     out=to_srgb(ref).min() < 0))
            if metal == "Au":
                per_basis[b] = rows
            hid = [r for r in rows if r["clip"] < HIDDEN and r["p3"] >= VISIBLE]
            vis = [r for r in rows if r["p3"] >= VISIBLE]
            outs = [r for r in rows if r["out"] and r["p3"] > 0.05]
            removed = np.mean([1 - r["clip"] / r["p3"] for r in outs]) if outs else float("nan")
            mean = lambda key: np.mean([r[key] for r in rows])
            emit(f"{label:>16s} {len(hid):>11d} / {len(vis):<10d} {removed * 100:20.0f}% {mean('lin'):15.2f} {mean('p3'):6.2f} {mean('clip'):10.2f} {mean('oklch'):11.2f}")
            if hid:
                emit(f"{'':>16s}   hidden: " + "; ".join(f"{r['light']} x{r['n']}: {r['clip']:.2f} -> {r['p3']:.2f}" for r in hid))
        emit()
    emit("--- A1. Au, one reflection only (the case the bowl mostly shows): dE sRGB clip / dE P3 ---")
    emit(f"{'light':>8s} " + " ".join(f"{BASES[b]:>18s}" for b in BASES))
    for i, (lname, _) in enumerate(lights):
        idx = i * len(BOUNCES)
        emit(f"{lname:>8s} " + " ".join(f"{per_basis[b][idx]['clip']:8.2f} / {per_basis[b][idx]['p3']:6.2f}  " for b in BASES))
    emit()
    emit("--- A2. Au, every condition: dE sRGB clip / dE P3 ---")
    emit(f"{'light':>8s} {'n':>2s} " + " ".join(f"{BASES[b]:>18s}" for b in BASES))
    for i, (lname, _) in enumerate(lights):
        for j, n in enumerate(BOUNCES):
            idx = i * len(BOUNCES) + j
            emit(f"{lname:>8s} {n:>2d} " + " ".join(f"{per_basis[b][idx]['clip']:8.2f} / {per_basis[b][idx]['p3']:6.2f}  " for b in BASES))
    emit()

    # ---------------- B ----------------
    emit("=== B. Where the error of the sRGB-basis renderer sits (linear sRGB channels; dE linear, no display) ===")
    emit("s = spectral, r = RGB renderer (linear sRGB, D65 white = 1). 'blue only' = dE(s, s with r's blue); 'R,G only' = dE(s, s with r's red and green).")
    emit(f"{'light':>8s} {'n':>2s} {'s (R, G, B)':>26s} {'r (R, G, B)':>26s} {'dE lin':>7s} {'blue only':>10s} {'R,G only':>9s} {'sRGB clip':>10s} {'P3':>6s}")
    for lname, spd in lights:
        for n in (1, 3):
            ref, rgb = M.xyz_pair(spd, n, "srgb")
            s, r = to_srgb(ref), to_srgb(rgb)
            d_b, d_rg = channel_split(ref, rgb)
            emit(f"{lname:>8s} {n:>2d} ({s[0]:7.4f},{s[1]:7.4f},{s[2]:8.4f}) ({r[0]:7.4f},{r[1]:7.4f},{r[2]:8.4f}) {de(ref, rgb):7.2f} {d_b:10.2f} {d_rg:9.2f} "
                 f"{de(shown(ref, 'srgb'), shown(rgb, 'srgb')):10.2f} {de(shown(ref, 'p3'), shown(rgb, 'p3')):6.2f}")
    emit()
    # the blue channel of the light and of the observer's sRGB blue sensor
    sb = M.S["srgb"][:, 2]
    neg = LAM[(sb < 0) & (LAM >= 380) & (LAM <= 780)]
    emit(f"sRGB blue sensor (CMF -> linear sRGB): negative from {neg.min():.0f} to {neg.max():.0f} nm (minimum {sb.min() / sb.max():.3f} of its peak at {LAM[np.argmax(sb)]:.0f} nm).")
    for lname, spd in lights:
        lb = (M.S["srgb"] * (spd * W_VIS)[:, None]).sum(0)
        if lname in ("1500 K", "1850 K", "2856 K", "D65"):
            emit(f"  light {lname:>7s}: linear sRGB ({lb[0]:.4f}, {lb[1]:.4f}, {lb[2]:.4f})")
    eta_c, k_c = M.red["srgb"]
    emit(f"  gold, reduced to sRGB channels: F0 = ({', '.join(f'{v:.4f}' for v in ma.fresnel_conductor(1.0, eta_c, k_c))})")
    # the colour temperature below which a blackbody light itself is outside sRGB (its blue < 0): the RGB renderer's blue is the light's
    # blue times gold's, so it can be <= 0 only below this temperature
    temps = np.arange(1000.0, 4000.0, 5.0)
    blue = np.array([(M.S["srgb"][:, 2] * light_spd(t) * W_VIS).sum() for t in temps])
    t0 = temps[np.argmax(blue > 0)]
    emit(f"  a blackbody's sRGB blue is negative below about {t0:.0f} K (the light itself is outside sRGB).")
    emit()

    # ---------------- C ----------------
    emit("=== C. The same split for the wide working bases, shown on the sRGB display ===")
    for b in ("rec2020", "acescg"):
        emit(f"-- {BASES[b]} --")
        emit(f"{'light':>8s} {'n':>2s} {'r in linear sRGB (R, G, B)':>28s} {'dE lin':>7s} {'blue only':>10s} {'R,G only':>9s} {'sRGB clip':>10s} {'P3':>6s}")
        for lname, spd in lights:
            if lname not in ("1500 K", "1850 K", "2856 K", "4500 K", "D65"):
                continue
            for n in (1, 3):
                ref, rgb = M.xyz_pair(spd, n, b)
                r = to_srgb(rgb)
                d_b, d_rg = channel_split(ref, rgb)
                emit(f"{lname:>8s} {n:>2d} ({r[0]:7.4f},{r[1]:7.4f},{r[2]:8.4f})   {de(ref, rgb):7.2f} {d_b:10.2f} {d_rg:9.2f} "
                     f"{de(shown(ref, 'srgb'), shown(rgb, 'srgb')):10.2f} {de(shown(ref, 'p3'), shown(rgb, 'p3')):6.2f}")
        emit()

    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        with open(a.out, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(lines) + "\n")
        print("wrote", a.out)


if __name__ == "__main__":
    main()
