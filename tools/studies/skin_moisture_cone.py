"""Which mechanism of a lotion changes how much light a thin part of skin lets through? (docs/SKIN_MOISTURE.md, step D)

A cone of skin (20 mm at the base, 3 mm near the tip) is lit from behind (a thin strip: the glow) or from the front (a softbox: the veil), one state of the skin at a time:
  dry     flakes scatter white light (rho0 = 0.2, not measured)
  m1      the flakes are gone (nothing else changes); also the reference with no surface layers
  m2      a water film changes the specular reflection (the flakes are kept)
  m3      the scatter distance is 1.3 x (a hydrated stratum corneum scatters less inside)
  peak    a lotion 200 s after it was applied (everything together; scatter distance 1.24 x)
  hour    a lotion 3600 s after it was applied (scatter distance 1.038 x)

The BSSRDF cannot show how the transmission depends on thickness (docs/SKIN_MOISTURE.md 5.4), so the bulk scattering is taken from a random walk in a medium (walk_*), and
the effect of the surface layers from the BSSRDF as a ratio (with / without the layers; they act at the entrance and the exit, where the BSSRDF is right):
    glow(state) = walk(bulk of the state) x BSSRDF_back(state) / BSSRDF_back(the same bulk, no surface layers).

    python tools/studies/skin_moisture_cone.py make  OUTDIR      writes the scenes (about 20 files)
    (render every OUTDIR/cone_*_{back,front}.json and OUTDIR/walk_s*.json from GoLD_rayt\\ with   --out FILE.pfm   about 40 s each)
    python tools/studies/skin_moisture_cone.py analyze OUTDIR     the tables
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

CONE = [[0, 0], [0.010, 0], [0.0015, 0.06], [0, 0.0615]]
CAMERA = {"from": [0, 0.03, 0.2], "at": [0, 0.03, 0], "up": [1, 0, 0], "fov": 11.5}
W_PX, H_PX = 480, 270
OPTIONS = {"width": W_PX, "spp": 256, "depth": 12, "denoise": False}
BASE_SKIN = {"type": "skin", "albedo": [0.80, 0.52, 0.42], "roughness": 0.3, "specular": 0.5, "scatterDistanceMm": [3.67, 1.37, 0.68], "relief": {"wrinkles": {"depthScale": 1.6}}}
FLAKES = {"rho0": 0.2, "threshold": 0.5, "softness": 0.06}
STATES = {
    "dry": {"scatter": FLAKES},
    "m1": {"scatter": {"rho0": 0.0}},
    "m2": {"scatter": dict(FLAKES, gapUm=1e6), "film": {"thicknessUm": 8, "ior": 1.33, "roughness": 0.02}},
    "m3": {"scatter": FLAKES, "sssScale": 1.3},
    "peak": {"scatter": dict(FLAKES, kappaH=0.3), "lotion": {"timeS": 200, "bulk": 0.3}},
    "hour": {"scatter": dict(FLAKES, kappaH=0.3), "lotion": {"timeS": 3600, "bulk": 0.3}},
    # a rich emulsion (20 um of liquid, 15 % oil, 1 um droplets) just after it was applied and an hour later; the same hydration as the water lotion
    "emul0": {"scatter": dict(FLAKES, kappaH=0.3), "lotion": {"timeS": 0, "bulk": 0.3, "doseUm": 20, "oil": 0.15}},
    "emul3600": {"scatter": dict(FLAKES, kappaH=0.3), "lotion": {"timeS": 3600, "bulk": 0.3, "doseUm": 20, "oil": 0.15}},
    "plain124": {"scatter": {"rho0": 0.0}, "sssScale": 1.24},
    "plain1038": {"scatter": {"rho0": 0.0}, "sssScale": 1.038},
}
BACK = {"shape": "quad", "material": "back", "corner": [-0.0008, 0.0, -0.03], "edge1": [0.0016, 0, 0], "edge2": [0, 0.0615, 0]}
FRONT = {"shape": "quad", "material": "front", "corner": [-0.2, 0.2, 0.25], "edge1": [0.2, 0, 0], "edge2": [0, -0.2, 0.1]}
BACK_RADIANCE, FRONT_RADIANCE = 100.0, 1.2
BULK_OF = {"dry": "1000", "m1": "1000", "m2": "1000", "m3": "1300", "peak": "1240", "hour": "1038", "emul0": "1000", "emul3600": "1038"}
REFERENCE_OF = {"dry": "m1", "m1": "m1", "m2": "m1", "peak": "plain124", "hour": "plain1038", "emul0": "m1", "emul3600": "plain1038"}
EDGES = (3, 5, 7, 9, 11, 13, 15, 17.5)   # diameter bins (mm)
Y_W = np.array([0.2126, 0.7152, 0.0722])
M = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])
WHITE = np.array([0.9505, 1.0, 1.089])


def make(out):
    out.mkdir(parents=True, exist_ok=True)
    for name, st in STATES.items():
        for light, quad in (("back", BACK), ("front", FRONT)):
            skin = dict(BASE_SKIN)
            skin.update({k: v for k, v in st.items() if k in ("scatter", "film", "lotion", "sssScale")})
            scene = {"options": OPTIONS, "camera": CAMERA,
                     "materials": {"skin": skin, "back": {"type": "emission", "spectrum": "d65", "radiance": BACK_RADIANCE},
                                   "front": {"type": "emission", "spectrum": "d65", "radiance": FRONT_RADIANCE}},
                     "objects": [{"shape": "lathe", "material": "skin", "segments": 96, "translate": [0, 0, 0], "profile": CONE}, quad]}
            (out / f"cone_{name}_{light}.json").write_text(json.dumps(scene, indent=1) + "\n", encoding="utf-8")
    for scale in (1.0, 1.038, 1.24, 1.3):   # the same optical properties as a medium (sigma_t / scale)
        mat = {"type": "dielectric", "ior": 1.4, "roughness": 0.45,
               "medium": {"sigmaT": [round(s / scale, 3) for s in (272, 730, 1470)], "albedo": [0.9906, 0.921, 0.8652], "g": 0}}
        scene = {"options": OPTIONS, "camera": CAMERA, "materials": {"skin": mat, "back": {"type": "emission", "spectrum": "d65", "radiance": BACK_RADIANCE}},
                 "objects": [{"shape": "lathe", "material": "skin", "segments": 96, "translate": [0, 0, 0], "profile": CONE}, BACK]}
        (out / f"walk_s{int(round(scale * 1000))}.json").write_text(json.dumps(scene, indent=1) + "\n", encoding="utf-8")
    print("wrote", len(list(out.glob("*.json"))), "scenes to", out)


def read_pfm(path):
    with open(path, "rb") as f:
        kind = f.readline().strip()
        w, h = map(int, f.readline().split())
        scale = float(f.readline())
        data = np.frombuffer(f.read(), dtype="<f4" if scale < 0 else ">f4")
    return data.reshape(h, w, 3 if kind == b"PF" else 1)[::-1].astype(float)


def profile(path, per_radiance):
    """Mean linear RGB along the axis of the cone, in bins of its diameter, per unit radiance of the light."""
    strip = read_pfm(path)[H_PX // 2 - 4: H_PX // 2 + 5].mean(0)
    half_w = np.tan(np.radians(CAMERA["fov"] / 2)) * 0.2 * (W_PX / H_PX)
    ys = 0.03 + half_w - (np.arange(W_PX) + 0.5) / (W_PX / (2 * half_w))   # +y runs to the left
    dm = 2 * (0.010 + (0.0015 - 0.010) * ys / 0.06) * 1000
    ok = (ys > 0.001) & (ys < 0.0605)
    return np.array([strip[ok & (dm >= lo) & (dm < hi)].mean(0) for lo, hi in zip(EDGES[:-1], EDGES[1:])]) / per_radiance


def lab(rgb, white):
    xyz = (rgb / white) @ M.T / WHITE
    f = np.where(xyz > 216 / 24389, np.cbrt(xyz), (24389 / 27 * xyz + 16) / 116)
    return 116 * f[..., 1] - 16, 500 * (f[..., 0] - f[..., 1]), 200 * (f[..., 1] - f[..., 2])


def analyze(d):
    back = {s: profile(d / f"cone_{s}_back.pfm", BACK_RADIANCE) for s in ("dry", "m1", "m2", "peak", "hour", "emul0", "emul3600", "plain124", "plain1038")}
    front = {s: profile(d / f"cone_{s}_front.pfm", FRONT_RADIANCE) for s in ("dry", "m1", "m2", "m3", "peak", "hour", "emul0", "emul3600")}
    walk = {b: profile(d / f"walk_s{b}.pfm", BACK_RADIANCE) for b in ("1000", "1038", "1240", "1300")}
    mids = [0.5 * (a + b) for a, b in zip(EDGES[:-1], EDGES[1:])]
    glow = {s: walk[BULK_OF[s]] * (back[s] / back[REFERENCE_OF[s]] if s in REFERENCE_OF else 1.0) for s in BULK_OF}
    Y = lambda a: a @ Y_W
    names = ("m1", "m2", "m3", "peak", "hour", "emul0", "emul3600")
    print("glow (luminance per unit radiance of the back light): dry, and the ratio to it")
    for i, m in enumerate(mids):
        print(f"  {m:5.1f} mm  {Y(glow['dry'][i]):8.4f} | " + " | ".join(f"{s} {Y(glow[s][i]) / Y(glow['dry'][i]):5.2f}" for s in names))
    print("veil (front-lit reflection): dry, and the ratio to it")
    for i, m in enumerate(mids):
        print(f"  {m:5.1f} mm  {Y(front['dry'][i]):8.5f} | " + " | ".join(f"{s} {Y(front[s][i]) / Y(front['dry'][i]):5.2f}" for s in names))
    i8 = mids.index(8.0)
    k = Y(front["dry"][i8]) / Y(glow["dry"][i8])
    total = {s: k * glow[s] + front[s] for s in glow}
    white = 2.0 * Y(total["dry"][i8])
    print(f"both lights (the glow weighted by {k:.4g} so that at 8 mm the dry glow and the dry veil are equally bright): CIELAB, dE*ab to the dry skin")
    for i, m in enumerate(mids):
        b = lab(total["dry"][i], white)
        line = f"  {m:5.1f} mm dry L* {b[0]:5.1f} C* {np.hypot(b[1], b[2]):5.1f} |"
        for s in names:
            L, a, bb = lab(total[s][i], white)
            line += f" {s} L* {L:5.1f} C* {np.hypot(a, bb):5.1f} dE {np.sqrt((L - b[0]) ** 2 + (a - b[1]) ** 2 + (bb - b[2]) ** 2):4.1f} |"
        print(line)
    with open(d / "cone_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["diameter_mm", "state", "glow_Y", "veil_Y"])
        for i, m in enumerate(mids):
            for s in glow:
                if s in front:
                    w.writerow([m, s, Y(glow[s][i]), Y(front[s][i])])


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in ("make", "analyze"):
        print(__doc__)
        sys.exit(1)
    (make if sys.argv[1] == "make" else analyze)(Path(sys.argv[2]))
