"""A lotion, then an emulsion: what do the second product and the wait before it change? (docs/SKIN_MOISTURE.md 5.6)

Four ways of caring for the skin of the cone of skin_moisture_cone.py (the same lights, the same dry skin):
  L       a water lotion (10 um, 5 % humectant)
  E       a rich emulsion (20 um, 15 % oil)
  LE60    the lotion, and the emulsion 60 s later (the lotion is still wet)
  LE600   the lotion, and the emulsion 600 s later (the lotion has dried)
each seen 0 / 200 / 1800 / 3600 s after its last product, and LE60 an hour later with an occlusion of 0 and 0.6 instead of 0.3.

As in skin_moisture_cone.py the glow (light from behind) is a random walk for the bulk times the effect of the surface layers from the BSSRDF:
    glow(state) = walk(scatter distance x (1 + BULK x hydration)) x BSSRDF_back(state, bulk off) / BSSRDF_back(no surface layers).
The walk is interpolated between the scatter distances that were rendered (a power law per diameter and channel); the hydration of a state is read from
the "[Lotion]" line of its render log. The veil (light from the front) is the BSSRDF with the bulk on.
Also a close-up series (scenes/skin_macro.json, 450 px, 96 spp): the lotion alone, and with the emulsion at 60 s.

    python tools/studies/skin_moisture_layers.py make OUTDIR             writes the scenes (about 50 files)
    (render every OUTDIR/cone_*.json and walk_*.json from GoLD_rayt with  --out FILE.pfm,  every macro_*.json with  --out FILE.png,
     and keep what each run prints in FILE.log; about 40 s each)
    python tools/studies/skin_moisture_layers.py analyze OUTDIR CONEDIR   the tables (CONEDIR: the renders of skin_moisture_cone.py)
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import skin_moisture_cone as cone  # noqa: E402

BULK = 0.3
LOTION = {"doseUm": 10}
EMULSION = {"doseUm": 20, "oil": 0.15}
SEQUENCES = {"L": [(0, LOTION)], "E": [(0, EMULSION)], "LE60": [(0, LOTION), (60, EMULSION)], "LE600": [(0, LOTION), (600, EMULSION)]}
AFTER = (0, 200, 1800, 3600)   # seconds after the last product
OCCLUSIONS = (0.0, 0.6)        # besides the default 0.3, for LE60 an hour later
WALK_SCALES = ("1000", "1038", "1120", "1240", "1300")
MACRO_TIMES = {"L": (0, 120, 360, 960, 1860, 3660), "LE60": (60, 120, 360, 960, 1860, 3660)}   # on the clock of the first product
MACRO_FLAKES = dict(cone.FLAKES, kappaH=0.3)   # the flakes of the close-up series of docs/SKIN_MOISTURE.md 5.5


def lotion(seq, time_s, bulk, occlusion=None):
    out = {"timeS": time_s, "layers": [dict(p, atS=at) for at, p in SEQUENCES[seq]]}
    if bulk:
        out["bulk"] = bulk
    if occlusion is not None:
        out["occlusion"] = occlusion
    return out


def states():
    """name -> (sequence, time on the clock of the first product, occlusion or None)"""
    out = {}
    for seq, layers in SEQUENCES.items():
        for after in AFTER:
            out[f"{seq}_{after}"] = (seq, layers[-1][0] + after, None)
    for occ in OCCLUSIONS:
        out[f"LE60_3600_o{int(occ * 100)}"] = ("LE60", 60 + 3600, occ)
    return out


def make(out):
    out.mkdir(parents=True, exist_ok=True)
    for name, (seq, time_s, occ) in states().items():
        for light, quad, bulk in (("back", cone.BACK, 0), ("front", cone.FRONT, BULK)):
            skin = dict(cone.BASE_SKIN, scatter=dict(cone.FLAKES, kappaH=0.3), lotion=lotion(seq, time_s, bulk, occ))
            scene = {"options": cone.OPTIONS, "camera": cone.CAMERA,
                     "materials": {"skin": skin, "back": {"type": "emission", "spectrum": "d65", "radiance": cone.BACK_RADIANCE},
                                   "front": {"type": "emission", "spectrum": "d65", "radiance": cone.FRONT_RADIANCE}},
                     "objects": [{"shape": "lathe", "material": "skin", "segments": 96, "translate": [0, 0, 0], "profile": cone.CONE}, quad]}
            (out / f"cone_{name}_{light}.json").write_text(json.dumps(scene, indent=1) + "\n", encoding="utf-8")
    scale = 1.12   # one more scatter distance of the random walk, between the ones of skin_moisture_cone.py
    mat = {"type": "dielectric", "ior": 1.4, "roughness": 0.45,
           "medium": {"sigmaT": [round(s / scale, 3) for s in (272, 730, 1470)], "albedo": [0.9906, 0.921, 0.8652], "g": 0}}
    scene = {"options": cone.OPTIONS, "camera": cone.CAMERA, "materials": {"skin": mat, "back": {"type": "emission", "spectrum": "d65", "radiance": cone.BACK_RADIANCE}},
             "objects": [{"shape": "lathe", "material": "skin", "segments": 96, "translate": [0, 0, 0], "profile": cone.CONE}, cone.BACK]}
    (out / "walk_s1120.json").write_text(json.dumps(scene, indent=1) + "\n", encoding="utf-8")

    macro_path = Path(__file__).resolve().parents[2] / "GoLD_rayt" / "scenes" / "skin_macro.json"
    macro = json.loads("".join(line for line in macro_path.read_text(encoding="utf-8").splitlines(True) if not line.lstrip().startswith("//")))
    macro["options"].update({"width": 450, "spp": 96, "denoise": False})
    macro["materials"]["skin"]["scatter"] = MACRO_FLAKES
    (out / "macro_dry.json").write_text(json.dumps(macro, indent=1) + "\n", encoding="utf-8")
    for seq, times in MACRO_TIMES.items():
        for t in times:
            macro["materials"]["skin"]["lotion"] = lotion(seq, t, 0)
            (out / f"macro_{seq}_{t}.json").write_text(json.dumps(macro, indent=1) + "\n", encoding="utf-8")
    print("wrote", len(list(out.glob("*.json"))), "scenes to", out)


def hydration(log):
    m = re.search(r"\[Lotion\].*hydration ([0-9.eE+-]+)", log.read_text(encoding="ascii", errors="ignore"))
    if not m:
        raise SystemExit(f"no [Lotion] line with a hydration in {log}")
    return float(m.group(1).rstrip("."))


def image_lab(path):
    """Mean L* and mean chroma C* of a picture (display sRGB)."""
    from PIL import Image
    v = np.asarray(Image.open(path).convert("RGB"), dtype=float) / 255.0
    lin = np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4)
    L, a, b = cone.lab(lin, 1.0)
    return L.mean(), np.hypot(a, b).mean()


def analyze(d, cone_dir):
    Y = lambda a: a @ cone.Y_W
    mids = [0.5 * (a + b) for a, b in zip(cone.EDGES[:-1], cone.EDGES[1:])]
    i8 = mids.index(8.0)
    thick = [i for i, m in enumerate(mids) if 4 <= m <= 16]
    walk = {s: cone.profile((d if s == "1120" else cone_dir) / f"walk_s{s}.pfm", cone.BACK_RADIANCE) for s in WALK_SCALES}

    def walk_at(scale, nodes=WALK_SCALES):
        xs = np.log([int(s) / 1000 for s in nodes])
        logs = np.log(np.array([walk[s] for s in nodes]))   # node, diameter, channel
        out = np.empty_like(walk["1000"])
        for i in range(out.shape[0]):
            for c in range(3):
                out[i, c] = np.exp(np.interp(np.log(scale), xs, logs[:, i, c]))
        return out

    plain_back = cone.profile(cone_dir / "cone_m1_back.pfm", cone.BACK_RADIANCE)
    dry_back = cone.profile(cone_dir / "cone_dry_back.pfm", cone.BACK_RADIANCE)
    dry_front = cone.profile(cone_dir / "cone_dry_front.pfm", cone.FRONT_RADIANCE)
    dry_glow = walk["1000"] * dry_back / plain_back

    # How good are the two shortcuts? (1) the walk at 1.12 x from its neighbours; (2) the surface factor with the bulk off against the one of skin_moisture_cone.py (bulk on, own reference)
    guess = walk_at(1.12, ("1000", "1038", "1240", "1300"))
    r = [Y(guess[i]) / Y(walk["1120"][i]) for i in thick]
    print(f"check 1: random walk at 1.12 x, interpolated without it / rendered: {min(r):.3f} .. {max(r):.3f} (4-16 mm)")
    for new, old, ref in (("L_200", "peak", "plain124"), ("L_3600", "hour", "plain1038")):
        a = cone.profile(d / f"cone_{new}_back.pfm", cone.BACK_RADIANCE) / plain_back
        b = cone.profile(cone_dir / f"cone_{old}_back.pfm", cone.BACK_RADIANCE) / cone.profile(cone_dir / f"cone_{ref}_back.pfm", cone.BACK_RADIANCE)
        r = [Y(a[i]) / Y(b[i]) for i in thick]
        print(f"check 2: surface factor of {new} (bulk off) / of '{old}' (bulk on): {min(r):.3f} .. {max(r):.3f} (4-16 mm)")

    k = Y(dry_front[i8]) / Y(dry_glow[i8])
    white = 2.0 * Y((k * dry_glow + dry_front)[i8])
    base = cone.lab((k * dry_glow + dry_front)[i8], white)
    print("state: hydration, scatter distance, glow / dry (4-16 mm; at 8 mm), veil / dry (4-16 mm), both lights at 8 mm: L*, C*, dE*ab to dry")
    print(f"  dry: L* {base[0]:.1f}, C* {np.hypot(base[1], base[2]):.1f}")
    rows = []
    for name in states():
        h = hydration(d / f"cone_{name}_front.log")
        scale = 1 + BULK * h
        glow = walk_at(scale) * cone.profile(d / f"cone_{name}_back.pfm", cone.BACK_RADIANCE) / plain_back
        front = cone.profile(d / f"cone_{name}_front.pfm", cone.FRONT_RADIANCE)
        g = [Y(glow[i]) / Y(dry_glow[i]) for i in thick]
        v = [Y(front[i]) / Y(dry_front[i]) for i in thick]
        L, a, b = cone.lab((k * glow + front)[i8], white)
        de = np.sqrt((L - base[0]) ** 2 + (a - base[1]) ** 2 + (b - base[2]) ** 2)
        rows.append((name, h, scale, min(g), max(g), Y(glow[i8]) / Y(dry_glow[i8]), min(v), max(v), L, np.hypot(a, b), de))
        print(f"  {name:14s} H {h:5.3f}  x{scale:5.3f}  glow {min(g):5.2f}..{max(g):5.2f} ({rows[-1][5]:5.2f})  veil {min(v):5.2f}..{max(v):5.2f}  L* {L:5.1f} C* {np.hypot(a, b):5.1f} dE {de:4.1f}")
    with open(d / "layers_summary.csv", "w", encoding="ascii") as fh:
        fh.write("state,hydration,scatter_scale,glow_min,glow_max,glow_8mm,veil_min,veil_max,L,C,dE\n")
        for row in rows:
            fh.write(",".join([row[0]] + [f"{x:.4f}" for x in row[1:]]) + "\n")

    if (d / "macro_dry.png").exists():
        L0, C0 = image_lab(d / "macro_dry.png")
        print(f"close-up: dry L* {L0:.2f}, C* {C0:.2f}; the difference to it (dL* / dC*) on the clock of the first product")
        for seq, times in MACRO_TIMES.items():
            line = f"  {seq:5s}"
            for t in times:
                L, C = image_lab(d / f"macro_{seq}_{t}.png")
                line += f" | {t} s {L - L0:+.2f} / {C - C0:+.2f}"
            print(line)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "make":
        make(Path(sys.argv[2]))
    elif len(sys.argv) == 4 and sys.argv[1] == "analyze":
        analyze(Path(sys.argv[2]), Path(sys.argv[3]))
    else:
        print(__doc__)
        sys.exit(1)
