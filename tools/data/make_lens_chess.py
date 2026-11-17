"""Writes the scenes of a lens comparison: one chess position at night, shot through four lenses.

    python tools/data/make_lens_chess.py OUT_DIR [--width W --height H --spp N]

Writes OUT_DIR/lens_chess_<frame>.json for the frames

  fisheye   assets/Lenses/fisheye_10mm.txt  wide open (f/4)       the straight lines of the board and of a row of lights bend
  dgauss2   assets/Lenses/dgauss_50mm.txt   wide open (f/2)       shallow focus, the lights become discs, cat's eyes at the edges
  dgauss16  assets/Lenses/dgauss_50mm.txt   f/16 + diffraction    nearly everything sharp, the lights become stars (enlarged)
  singlet   assets/Lenses/singlet_bk7_100mm.txt (f/5.6)           an uncorrected lens: soft, colored fringes (spectral build)

The camera keeps the gold queen the same size in every frame (it moves closer for the shorter lenses). The scene:
a board of 8 x 8 squares (5 cm) on a dark glossy table, gold and silver pieces (surfaces of revolution), a straight row
of small lights behind the board, two strings of lights farther back and scattered lights far away. The small lights do
not light the scene (they are only seen); a key and a rim light do.
"""
import argparse
import json
import math
import os
import random

ap = argparse.ArgumentParser()
ap.add_argument("out_dir")
ap.add_argument("--width", type=int, default=1200)
ap.add_argument("--height", type=int, default=800)
ap.add_argument("--spp", type=int, default=256)
args = ap.parse_args()
rng = random.Random(11)

SQ = 0.05            # square size (m)
TOP = 0.015          # height of the board's top above the table

# ---------------------------------------------------------------- piece profiles (radius, height), counter-clockwise


def chaikin(pts, iters=3):
    """Corner cutting of an open polyline; the end points (on the axis) stay."""
    for _ in range(iters):
        out = [pts[0]]
        for a, b in zip(pts[:-1], pts[1:]):
            out.append((0.75 * a[0] + 0.25 * b[0], 0.75 * a[1] + 0.25 * b[1]))
            out.append((0.25 * a[0] + 0.75 * b[0], 0.25 * a[1] + 0.75 * b[1]))
        out.append(pts[-1])
        pts = out
    return pts


def arc(cy, R, deg0, deg1, n=10):
    """Points of a circle of radius R centred on the axis at height cy, from angle deg0 to deg1 (90 = top)."""
    return [(R * math.cos(math.radians(d)), cy + R * math.sin(math.radians(d)))
            for d in [deg0 + (deg1 - deg0) * i / n for i in range(n + 1)]]


def base(r):
    return [(0, 0), (r, 0), (r, 0.004), (r * 0.82, 0.007), (r * 0.66, 0.009), (r * 0.7, 0.011)]


def profile(name):
    if name == "pawn":
        p = base(0.017) + [(0.0075, 0.014), (0.0055, 0.027), (0.0095, 0.029), (0.0095, 0.031), (0.005, 0.033)]
        p += arc(0.041, 0.0085, -55, 90)
    elif name == "rook":
        p = base(0.019) + [(0.0115, 0.015), (0.0105, 0.040), (0.0145, 0.045), (0.0145, 0.057), (0.0115, 0.057),
                           (0.0110, 0.052), (0, 0.052)]
    elif name == "bishop":
        p = base(0.018) + [(0.008, 0.016), (0.0055, 0.042), (0.0105, 0.045), (0.0105, 0.047), (0.006, 0.049)]
        p += [(0.0085, 0.053), (0.0095, 0.058), (0.0085, 0.064), (0.005, 0.069), (0.002, 0.071)]
        p += arc(0.0745, 0.0035, -40, 90, 6)
    elif name == "queen":
        p = base(0.019) + [(0.0085, 0.017), (0.006, 0.055), (0.0115, 0.059), (0.0115, 0.061), (0.007, 0.063)]
        p += [(0.0085, 0.070), (0.0125, 0.078), (0.0105, 0.079), (0.006, 0.079)]
        p += arc(0.0835, 0.0048, -60, 90, 8)
    elif name == "king":
        p = base(0.020) + [(0.009, 0.018), (0.0065, 0.062), (0.012, 0.066), (0.012, 0.068), (0.0075, 0.070)]
        p += [(0.0095, 0.078), (0.0115, 0.085), (0.008, 0.087), (0.003, 0.088), (0.003, 0.094), (0.0055, 0.095),
              (0.0055, 0.098), (0.003, 0.099), (0.003, 0.104), (0, 0.104)]
    else:
        raise ValueError(name)
    if p[-1][0] != 0:
        p.append((0, p[-1][1]))
    return [[round(r, 5), round(h, 5)] for r, h in chaikin(p)]


def square_center(file_idx, rank):
    """file_idx 0..7 = a..h along +x, rank 1..8: rank 1 nearest the camera (+z)."""
    return (-3.5 * SQ + file_idx * SQ, 3.5 * SQ - (rank - 1) * SQ)


# ---------------------------------------------------------------- materials

materials = {
    "gold":   {"type": "conductor", "metal": "Au", "roughness": 0.07},
    "silver": {"type": "conductor", "metal": "Ag", "roughness": 0.06},
    "light_sq": {"type": "layered", "coat": {"ior": 1.5, "roughness": 0.06}, "base": {"albedo": [0.72, 0.66, 0.55]}},
    "dark_sq":  {"type": "layered", "coat": {"ior": 1.5, "roughness": 0.06}, "base": {"albedo": [0.10, 0.05, 0.03]}},
    "frame":    {"type": "layered", "coat": {"ior": 1.5, "roughness": 0.08}, "base": {"albedo": [0.20, 0.09, 0.04]}},
    "table":    {"type": "layered", "coat": {"ior": 1.5, "roughness": 0.12}, "base": {"albedo": [0.05, 0.03, 0.02]}},
    "key":  {"type": "emission", "spectrum": "tungsten", "radiance": 22},
    "rim":  {"type": "emission", "spectrum": {"blackbody": 5500}, "radiance": 40},
}
lamp_temps = [1900, 2300, 2700, 3200, 4500, 6500]
for T in lamp_temps:
    materials["lamp%d" % T] = {"type": "emission", "spectrum": {"blackbody": T}, "radiance": 300}

objects = []


def quad(corner, e1, e2, m):
    objects.append({"shape": "quad", "corner": [round(x, 5) for x in corner], "edge1": [round(x, 5) for x in e1],
                    "edge2": [round(x, 5) for x in e2], "material": m})


def lamp(c, r, T):
    objects.append({"shape": "sphere", "center": [round(x, 4) for x in c], "radius": r, "material": "lamp%d" % T,
                    "light": False})


# ---------------------------------------------------------------- table, board

quad([-6, 0, 3], [12, 0, 0], [0, 0, -14], "table")          # faces up: cross(+x, -z) = +y
half = 4 * SQ
for f in range(8):
    for rnk in range(1, 9):
        x, z = square_center(f, rnk)
        m = "dark_sq" if (f + rnk) % 2 == 0 else "light_sq"   # a1 dark
        quad([x - SQ / 2, TOP, z + SQ / 2], [SQ, 0, 0], [0, 0, -SQ], m)
B = 0.025  # border width
o = half + B
quad([-o, TOP, o], [2 * o, 0, 0], [0, 0, -B], "frame")                    # near border
quad([-o, TOP, -half], [2 * o, 0, 0], [0, 0, -B], "frame")                # far border
quad([-o, TOP, half], [B, 0, 0], [0, 0, -2 * half], "frame")              # left
quad([half, TOP, half], [B, 0, 0], [0, 0, -2 * half], "frame")            # right
quad([-o, 0, o], [2 * o, 0, 0], [0, TOP, 0], "frame")                     # sides (facing out)
quad([o, 0, -o], [-2 * o, 0, 0], [0, TOP, 0], "frame")
quad([-o, 0, -o], [0, 0, 2 * o], [0, TOP, 0], "frame")
quad([o, 0, o], [0, 0, -2 * o], [0, TOP, 0], "frame")

# ---------------------------------------------------------------- pieces (a middle game)

position = [
    ("queen", "gold", "d3"), ("king", "gold", "g1"), ("rook", "gold", "a1"), ("rook", "gold", "f1"),
    ("bishop", "gold", "c4"), ("pawn", "gold", "a2"), ("pawn", "gold", "b3"), ("pawn", "gold", "e4"),
    ("pawn", "gold", "f2"), ("pawn", "gold", "g2"), ("pawn", "gold", "h3"),
    ("king", "silver", "g8"), ("queen", "silver", "d8"), ("rook", "silver", "a8"), ("rook", "silver", "e8"),
    ("bishop", "silver", "b7"), ("bishop", "silver", "e6"), ("pawn", "silver", "a7"), ("pawn", "silver", "c6"),
    ("pawn", "silver", "d5"), ("pawn", "silver", "f7"), ("pawn", "silver", "g6"), ("pawn", "silver", "h7"),
]
for piece, mat, sq in position:
    x, z = square_center("abcdefgh".index(sq[0]), int(sq[1]))
    objects.append({"shape": "lathe", "profile": profile(piece), "segments": 96, "translate": [round(x, 4), TOP, round(z, 4)],
                    "material": mat})
queen = square_center(3, 3)
subject = (queen[0], TOP + 0.05, queen[1])   # the queen's waist: the point in focus

# ---------------------------------------------------------------- lights

objects.append({"shape": "sphere", "center": [-0.9, 0.9, 0.7], "radius": 0.2, "material": "key"})
objects.append({"shape": "sphere", "center": [0.3, 5.0, -0.7], "radius": 0.5, "material": "rim"})

# a straight row behind the board (shows the distortion)
for i in range(81):
    lamp((-4 + 0.1 * i, 0.32, -1.4), 0.006, 2700)
# two strings with a sag, farther back
for zs, ylo, yhi, n in ((-2.8, 0.30, 0.75, 70), (-4.2, 0.45, 1.10, 70)):
    for i in range(n):
        t = i / (n - 1)
        x = -4.5 + 9 * t
        y = ylo + (yhi - ylo) * (2 * t - 1) ** 2
        lamp((x, y, zs), 0.007, rng.choice([1900, 2300, 2700]))
# scattered lights far away
for _ in range(45):
    lamp((rng.uniform(-5, 5), rng.uniform(0.05, 1.2), rng.uniform(-9, -5.5)), 0.012, rng.choice(lamp_temps))

# ---------------------------------------------------------------- frames

frames = {
    # name: (lens, f-number or None, subject distance (m), vertical field (deg, for aiming), extra camera keys)
    "fisheye":  ("assets/Lenses/fisheye_10mm.txt", None, 0.12, 0.0, {}),
    "dgauss2":  ("assets/Lenses/dgauss_50mm.txt", None, 0.6, 27.0, {"blades": 7}),
    "dgauss16": ("assets/Lenses/dgauss_50mm.txt", 16, 0.6, 27.0, {"blades": 7, "diffraction": 2}),
    "singlet":  ("assets/Lenses/singlet_bk7_100mm.txt", 5.6, 1.2, 13.7, {}),
}
ELEV = math.radians(9)     # camera above the subject
AZ = math.radians(-12)     # and a little to the left
os.makedirs(args.out_dir, exist_ok=True)
for name, (lens, fstop, dist, vfov, extra) in frames.items():
    d = (math.sin(AZ) * math.cos(ELEV), math.sin(ELEV), math.cos(AZ) * math.cos(ELEV))
    cam = [subject[i] + dist * d[i] for i in range(3)]
    # aim above the subject so that it sits in the lower part of the picture
    lift = dist * math.tan(math.radians(0.16 * vfov))
    at = [subject[0], subject[1] + lift, subject[2]]
    camera = {"from": [round(x, 4) for x in cam], "at": [round(x, 4) for x in at], "up": [0, 1, 0], "lens": lens,
              "focusAt": [round(x, 4) for x in subject]}
    if fstop:
        camera["fStop"] = fstop
    camera.update(extra)
    scene = {
        "comment": "Written by tools/data/make_lens_chess.py (frame %s)." % name,
        "options": {"width": args.width, "height": args.height, "spp": args.spp, "depth": 12, "denoise": True},
        "camera": camera,
        "materials": materials,
        "objects": objects,
    }
    path = os.path.join(args.out_dir, "lens_chess_%s.json" % name)
    # one key, material and object per line (the profiles of the pieces would make an indented file long)
    lines = ["{"]
    for key in ("comment", "options", "camera"):
        lines.append(' "%s": %s,' % (key, json.dumps(scene[key])))
    lines.append(' "materials": {')
    lines.append(",\n".join('  "%s": %s' % (k, json.dumps(v)) for k, v in materials.items()))
    lines.append(' },')
    lines.append(' "objects": [')
    lines.append(",\n".join("  " + json.dumps(o) for o in objects))
    lines.append(' ]')
    lines.append("}")
    text = "\n".join(lines) + "\n"
    assert json.loads(text) == scene
    with open(path, "w") as fp:
        fp.write(text)
    print(path)
