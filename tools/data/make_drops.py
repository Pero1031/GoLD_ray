"""Water drops on a glass plate: writes the meshes and the mask a scene needs (docs/WATERDROPS.md).

    python tools/data/make_drops.py GoLD_rayt/assets/Drops --plate 0.20 0.15 --coverage 0.25 --tilt 60 --seed 1

Output (metres, the plate in the x-y plane at z = 0, drops on the +z side, "down" = -y):
    caps.ply        the drops' free surfaces (water against air):  material  { "type": "dielectric", "ior": 1.333 }
    contacts.ply    the contact faces (water against glass), normals -z:  { "type": "dielectric", "ior": 1.333, "exteriorIor": 1.5 }
    plate_top.ply   the top face of the glass, with uv, to be used with "cutout": plate_mask.png (holes under the drops)
    plate_mask.png  RGBA, alpha 0 under every drop
    drops.csv       x, y, contact radius, tail length, height, apex shift (mm) and the kind of every drop

Model (phenomenological, after Kaneda et al. 1993's sphere sections; physics only as guidance):
  - contact radius a: log-normal (median, sigma), clipped; drops are placed without overlap (dart throwing) until the coverage is reached;
    a second population of mist (tiny droplets between the drops, --mist) with radii from --mist-min to --mist-max;
  - a drop is a spherical cap with the contact angle theta (<= 90 deg): sphere radius a / sin(theta), height a tan(theta / 2);
    larger drops are flattened by gravity, height x 1 / sqrt(1 + (a / l_c)^2) with the capillary length l_c = 2.7 mm;
  - on a tilted plate (tilt from horizontal, s = sin(tilt)) a drop that is big enough to feel gravity becomes a teardrop: round head
    downhill, a tapering tail uphill (tail length a (1 + 0.8 s b), pointedness 0.6 s b, b = min(1, (a / l_c)^2)), and its apex moves
    downhill by 0.3 s a min(1, a / l_c): the downhill side is steeper (the advancing angle), the uphill side flatter;
  - a fraction of the big drops on a steep plate become streaks: a long thin teardrop, apex near the lower end.
Hydrophobic surfaces (contact angle above 90 deg, overhanging rims) and running rivulets are not modelled.
"""
import argparse
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mesh_clean import smooth_normals, write_ply   # noqa: E402

L_CAP = 2.7e-3   # capillary length of water (m)


def cap_profile(rho, theta):
    """Height of a spherical cap of contact radius 1 and contact angle theta at the normalized radius rho (0 = center, 1 = rim)."""
    R = 1.0 / math.sin(theta)
    return np.sqrt(np.maximum(R * R - rho * rho, 0.0)) - R * math.cos(theta)


def outline(phi, a, tail, point):
    """The contact line of a teardrop, star-shaped from the origin: a circle of radius a below (sin phi <= 0, the head, downhill = -y),
    above it a tapering tail of length `tail` whose width shrinks with `point` (0 = a half-ellipse, 1 = a sharp tip)."""
    s = np.sin(phi)
    up = s > 0
    x = a * np.cos(phi) * np.where(up, 1.0 - point * s, 1.0)
    y = np.where(up, tail * s, a * s)
    return x, y


def make_drop(a, tail, point, h, shift, theta, segments, rings):
    """One drop: vertices (n, 3) of the cap (rim at z = 0), its triangles, and the rim polygon (counter-clockwise seen from +z)."""
    rho = np.sin(np.linspace(0, 1, rings + 1)[1:] * math.pi / 2)   # finer toward the rim
    phi = np.linspace(0, 2 * math.pi, segments, endpoint=False)
    ox, oy = outline(phi, a, tail, point)
    z0 = cap_profile(np.array([0.0]), theta)[0]
    verts = [(0.0, -shift, h)]                                        # apex (shifted downhill)
    for r in rho:
        z = h * cap_profile(np.array([r]), theta)[0] / z0
        for k in range(segments):
            verts.append((ox[k] * r, oy[k] * r - shift * (1 - r * r), z))
    verts = np.array(verts)
    faces = []
    for k in range(segments):                                         # the fan around the apex
        faces.append((0, 1 + k, 1 + (k + 1) % segments))
    for i in range(rings - 1):                                        # the quads between rings
        base, nxt = 1 + i * segments, 1 + (i + 1) * segments
        for k in range(segments):
            k1 = (k + 1) % segments
            faces.append((base + k, nxt + k, nxt + k1))
            faces.append((base + k, nxt + k1, base + k1))
    rim = list(range(1 + (rings - 1) * segments, 1 + rings * segments))
    return verts, np.array(faces, dtype=np.int64), rim


class Grid:
    """Occupancy of discs for the overlap test (cell = the largest disc)."""

    def __init__(self, W, H, cell):
        self.cell = cell
        self.nx, self.ny = max(1, int(W / cell) + 1), max(1, int(H / cell) + 1)
        self.cells = {}
        self.W, self.H = W, H

    def key(self, x, y):
        return (int((x + self.W / 2) / self.cell), int((y + self.H / 2) / self.cell))

    def free(self, x, y, r):
        kx, ky = self.key(x, y)
        for i in range(kx - 1, kx + 2):
            for j in range(ky - 1, ky + 2):
                for (px, py, pr) in self.cells.get((i, j), ()):
                    if (x - px) ** 2 + (y - py) ** 2 < (r + pr) ** 2:
                        return False
        return True

    def add(self, x, y, r):
        self.cells.setdefault(self.key(x, y), []).append((x, y, r))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--plate", type=float, nargs=2, default=[0.20, 0.15], metavar=("W", "H"), help="plate size in metres")
    ap.add_argument("--coverage", type=float, default=0.25, help="fraction of the plate covered by drops")
    ap.add_argument("--median-mm", type=float, default=1.2, help="median contact radius (mm)")
    ap.add_argument("--sigma", type=float, default=0.6, help="log-normal spread of the radius")
    ap.add_argument("--min-mm", type=float, default=0.25)
    ap.add_argument("--max-mm", type=float, default=4.0)
    ap.add_argument("--mist", type=float, default=0.02, help="fraction of the plate covered by mist (tiny droplets); 0 = none")
    ap.add_argument("--mist-min", type=float, default=0.05, help="mist radius range (mm), log-uniform")
    ap.add_argument("--mist-max", type=float, default=0.25)
    ap.add_argument("--contact-angle", type=float, default=65.0, help="degrees, up to 90 (clean glass 20-40, a used windshield 50-70)")
    ap.add_argument("--tilt", type=float, default=0.0, help="inclination of the plate from horizontal, degrees (teardrops and streaks)")
    ap.add_argument("--streaks", type=float, default=0.15, help="fraction of the big drops that become streaks on a steep plate")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--mask-res", type=int, default=4096, help="texels along the longer side of the mask")
    ap.add_argument("--segments", type=int, default=40)
    ap.add_argument("--rings", type=int, default=10)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    W, H = args.plate
    theta = math.radians(min(args.contact_angle, 90.0))
    s = math.sin(math.radians(args.tilt))

    # ---- the drops: sizes, shapes, positions (no overlaps) ----
    drops = []        # dicts: x, y, a, tail, point, h, shift, kind
    grid = Grid(W, H, args.max_mm * 1e-3 * 2.5)
    covered, tries = 0.0, 0
    target = args.coverage * W * H
    while covered < target and tries < 400000:
        tries += 1
        a = float(np.clip(rng.lognormal(math.log(args.median_mm * 1e-3), args.sigma), args.min_mm * 1e-3, args.max_mm * 1e-3))
        big = min(1.0, (a / L_CAP) ** 2)
        h = a * math.tan(theta / 2) / math.sqrt(1 + (a / L_CAP) ** 2)
        tail, point = a * (1 + 0.8 * s * big), 0.6 * s * big
        shift = 0.3 * s * a * min(1.0, a / L_CAP)
        kind = "drop"
        if s > 0.5 and a > 1e-3 and rng.random() < args.streaks:   # a streak: a long thin teardrop, the apex near the lower end
            a *= 0.7
            tail, point = a * rng.uniform(4.0, 8.0), 0.85
            shift, h, kind = 0.3 * a, h * 0.6, "streak"
        r = max(a, tail) * 1.05
        x = rng.uniform(-W / 2 + r, W / 2 - r)
        y = rng.uniform(-H / 2 + r, H / 2 - r)
        if not grid.free(x, y, r):
            continue
        grid.add(x, y, r)
        drops.append(dict(x=x, y=y, a=a, tail=tail, point=point, h=h, shift=shift, kind=kind))
        covered += math.pi * a * (a + tail) / 2
    n_drops = len(drops)
    # mist between the drops
    mist_target = args.mist * W * H
    mist_area, tries = 0.0, 0
    while args.mist > 0 and mist_area < mist_target and tries < 2000000:
        tries += 1
        a = math.exp(rng.uniform(math.log(args.mist_min * 1e-3), math.log(args.mist_max * 1e-3)))
        x, y = rng.uniform(-W / 2 + a, W / 2 - a), rng.uniform(-H / 2 + a, H / 2 - a)
        if not grid.free(x, y, a * 1.3):
            continue
        grid.add(x, y, a * 1.3)
        drops.append(dict(x=x, y=y, a=a, tail=a, point=0.0, h=a * math.tan(theta / 2), shift=0.0, kind="mist"))
        mist_area += math.pi * a * a
    print("drops: %d (coverage %.1f %% of %.0f x %.0f mm, radius %.2f-%.2f mm), mist: %d (%.1f %%)" % (
        n_drops, 100 * covered / (W * H), W * 1e3, H * 1e3, min(d["a"] for d in drops[:n_drops]) * 1e3, max(d["a"] for d in drops[:n_drops]) * 1e3,
        len(drops) - n_drops, 100 * mist_area / (W * H)))

    # ---- meshes ----
    cap_v, cap_f, con_v, con_f, rims = [], [], [], [], []
    nv_cap = nv_con = 0
    for d in drops:
        if d["kind"] == "mist":
            seg, rings = 10, 3
        else:
            seg, rings = int(np.clip(12 + 20 * d["a"] * 1e3, 16, args.segments)), args.rings
        v, f, rim = make_drop(d["a"], d["tail"], d["point"], d["h"], d["shift"], theta, seg, rings)
        v = v + np.array([d["x"], d["y"], 0.0])
        cap_v.append(v); cap_f.append(f + nv_cap); nv_cap += len(v)
        rim_pts = v[rim]
        rims.append(rim_pts[:, :2])
        # the contact face: the rim polygon fanned from its center (star-shaped), wound clockwise seen from +z so the normal points -z
        n = len(rim_pts)
        con_v.append(np.vstack([[[d["x"], d["y"], 0.0]], rim_pts]))
        con_f.append(np.array([(nv_con, nv_con + 1 + (k + 1) % n, nv_con + 1 + k) for k in range(n)], dtype=np.int64))
        nv_con += n + 1
    cap_v = np.vstack(cap_v); cap_f = np.vstack(cap_f)
    con_v = np.vstack(con_v); con_f = np.vstack(con_f)
    os.makedirs(args.outdir, exist_ok=True)
    write_ply(os.path.join(args.outdir, "caps.ply"), cap_v, smooth_normals(cap_v, cap_f), cap_f, None, "water drops, free surfaces (tools/data/make_drops.py)")
    write_ply(os.path.join(args.outdir, "contacts.ply"), con_v, np.tile(np.array([[0.0, 0.0, -1.0]]), (len(con_v), 1)), con_f, None,
              "water drops, contact faces against the glass (normals -z)")

    # the top face of the plate, with uv (s = x / W + 1/2, t = y / H + 1/2): the renderer flips t, so texture row 0 is y = +H/2
    pv = np.array([[-W / 2, -H / 2, 0], [W / 2, -H / 2, 0], [W / 2, H / 2, 0], [-W / 2, H / 2, 0]], dtype=np.float64)
    puv = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=np.float64)
    pf = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    write_ply(os.path.join(args.outdir, "plate_top.ply"), pv, np.tile(np.array([[0.0, 0.0, 1.0]]), (4, 1)), pf, puv, "glass plate, top face with uv for the cutout mask")

    # ---- the mask: alpha 0 inside every contact polygon ----
    from PIL import Image, ImageDraw
    if W >= H:
        rx, ry = args.mask_res, max(1, int(round(args.mask_res * H / W)))
    else:
        ry, rx = args.mask_res, max(1, int(round(args.mask_res * W / H)))
    alpha = Image.new("L", (rx, ry), 255)
    draw = ImageDraw.Draw(alpha)
    for pts in rims:
        px = (pts[:, 0] + W / 2) / W * rx
        py = (H / 2 - pts[:, 1]) / H * ry                 # row 0 = y = +H/2
        draw.polygon(list(zip(px.tolist(), py.tolist())), fill=0)
    rgba = np.zeros((ry, rx, 4), dtype=np.uint8)
    rgba[..., :3] = 255
    rgba[..., 3] = np.asarray(alpha)
    Image.fromarray(rgba, "RGBA").save(os.path.join(args.outdir, "plate_mask.png"))

    with open(os.path.join(args.outdir, "drops.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_mm", "y_mm", "a_mm", "tail_mm", "height_mm", "shift_mm", "kind"])
        for d in drops:
            w.writerow(["%.3f" % (d[k] * 1e3) for k in ("x", "y", "a", "tail", "h", "shift")] + [d["kind"]])
    print("wrote caps.ply (%d triangles), contacts.ply (%d), plate_top.ply, plate_mask.png (%d x %d), drops.csv in %s" % (
        len(cap_f), len(con_f), rx, ry, args.outdir))


if __name__ == "__main__":
    main()
