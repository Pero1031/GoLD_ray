"""Water drops on a glass plate: writes the meshes and the mask a scene needs (docs/WATERDROPS.md).

    python tools/data/make_drops.py GoLD_rayt/assets/Drops --plate 0.20 0.15 --coverage 0.25 --tilt 60 --seed 1

Output (metres, the plate in the x-y plane at z = 0, drops on the +z side, "down" = -y):
    caps.ply        the drops' free surfaces (water against air):  material  { "type": "dielectric", "ior": 1.333 }
    contacts.ply    the contact discs (water against glass), normals -z:  { "type": "dielectric", "ior": 1.333, "exteriorIor": 1.5 }
    plate_top.ply   the top face of the glass, with uv, to be used with "cutout": plate_mask.png (holes under the drops)
    plate_mask.png  RGBA, alpha 0 under every drop
    drops.csv       x, y, half axes a_x a_y, height, apex shift (mm) of every drop

Model (phenomenological, after Kaneda et al. 1993's sphere sections; physics only as guidance):
  - contact radius a: log-normal (median, sigma), clipped; drops are placed without overlap (dart throwing) until the coverage is reached;
  - a drop is a spherical cap with the contact angle theta (<= 90 deg): sphere radius a / sin(theta), height a tan(theta / 2);
    larger drops are flattened by gravity, height x 1 / sqrt(1 + (a / l_c)^2) with the capillary length l_c = 2.7 mm;
  - on a tilted plate (tilt from horizontal, s = sin(tilt)) a drop is stretched downhill (half axis a_y = a (1 + 0.5 s min(1, (a / l_c)^2)))
    and its apex moves downhill by 0.3 s a min(1, a / l_c): the downhill side is steeper (the advancing angle), the uphill side flatter;
  - a fraction of the big drops on a steep plate become streaks: 3 to 6 times longer than wide, thin, apex near the lower end.
Hydrophobic surfaces (contact angle above 90 deg, overhanging rims) are not modelled.
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


def make_drop(a, ax, ay, h, shift, theta, segments, rings):
    """One drop: vertices (n, 3) of the cap (rim at z = 0), its triangles, and the rim polygon (counter-clockwise seen from +z)."""
    rho = np.sin(np.linspace(0, 1, rings + 1)[1:] * math.pi / 2)   # finer toward the rim
    phi = np.linspace(0, 2 * math.pi, segments, endpoint=False)
    z0 = cap_profile(np.array([0.0]), theta)[0]
    verts = [(0.0, -shift, h)]                                        # apex (shifted downhill)
    for r in rho:
        z = h * cap_profile(np.array([r]), theta)[0] / z0
        for p in phi:
            verts.append((ax * r * math.cos(p), ay * r * math.sin(p) - shift * (1 - r * r), z))
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--plate", type=float, nargs=2, default=[0.20, 0.15], metavar=("W", "H"), help="plate size in metres")
    ap.add_argument("--coverage", type=float, default=0.25, help="fraction of the plate covered by drops")
    ap.add_argument("--median-mm", type=float, default=1.2, help="median contact radius (mm)")
    ap.add_argument("--sigma", type=float, default=0.6, help="log-normal spread of the radius")
    ap.add_argument("--min-mm", type=float, default=0.25)
    ap.add_argument("--max-mm", type=float, default=4.0)
    ap.add_argument("--contact-angle", type=float, default=65.0, help="degrees, up to 90 (clean glass 20-40, a used windshield 50-70)")
    ap.add_argument("--tilt", type=float, default=0.0, help="inclination of the plate from horizontal, degrees (asymmetry and streaks)")
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
    target = args.coverage * W * H

    # ---- the drops: sizes, shapes, positions (no overlaps) ----
    drops = []        # (x, y, a, ax, ay, h, shift)
    covered = 0.0
    tries = 0
    while covered < target and tries < 200000:
        tries += 1
        a = float(np.clip(rng.lognormal(math.log(args.median_mm * 1e-3), args.sigma), args.min_mm * 1e-3, args.max_mm * 1e-3))
        big = min(1.0, (a / L_CAP) ** 2)
        h = a * math.tan(theta / 2) / math.sqrt(1 + (a / L_CAP) ** 2)
        ax, ay = a, a * (1 + 0.5 * s * big)
        shift = 0.3 * s * a * min(1.0, a / L_CAP)
        if s > 0.5 and a > 1e-3 and rng.random() < args.streaks:   # a streak: long, thin, the apex near the lower end
            e = rng.uniform(3.0, 6.0)
            ax, ay = a * 0.8, a * e
            shift = 0.35 * ay
            h *= 0.6
        x = rng.uniform(-W / 2 + ay, W / 2 - ay)
        y = rng.uniform(-H / 2 + ay, H / 2 - ay)
        r = max(ax, ay) * 1.05
        if any((x - d[0]) ** 2 + (y - d[1]) ** 2 < (r + max(d[3], d[4]) * 1.05) ** 2 for d in drops):
            continue
        drops.append((x, y, a, ax, ay, h, shift))
        covered += math.pi * ax * ay
    print("drops: %d, coverage %.1f %% of %.0f x %.0f mm (asked %.1f %%), contact radius %.2f-%.2f mm, %d tries" % (
        len(drops), 100 * covered / (W * H), W * 1e3, H * 1e3, 100 * args.coverage, min(d[2] for d in drops) * 1e3, max(d[2] for d in drops) * 1e3, tries))

    # ---- meshes ----
    cap_v, cap_f, con_v, con_f = [], [], [], []
    for (x, y, a, ax, ay, h, shift) in drops:
        seg = int(np.clip(12 + 20 * a * 1e3, 12, args.segments))
        v, f, rim = make_drop(a, ax, ay, h, shift, theta, seg, args.rings)
        v = v + np.array([x, y, 0.0])
        base = sum(len(c) for c in cap_v)
        cap_v.append(v)
        cap_f.append(f + base)
        # the contact disc: the rim polygon, fanned from its center, wound clockwise seen from +z so that the normal points -z (into the glass)
        rim_pts = v[rim]
        cb = sum(len(c) for c in con_v)
        con_v.append(np.vstack([[[x, y, 0.0]], rim_pts]))
        n = len(rim_pts)
        con_f.append(np.array([(cb, cb + 1 + (k + 1) % n, cb + 1 + k) for k in range(n)], dtype=np.int64))
    cap_v = np.vstack(cap_v); cap_f = np.vstack(cap_f)
    con_v = np.vstack(con_v); con_f = np.vstack(con_f)
    os.makedirs(args.outdir, exist_ok=True)
    cap_n = smooth_normals(cap_v, cap_f)
    write_ply(os.path.join(args.outdir, "caps.ply"), cap_v, cap_n, cap_f, None, "water drops, free surfaces (tools/data/make_drops.py)")
    con_n = np.tile(np.array([[0.0, 0.0, -1.0]]), (len(con_v), 1))
    write_ply(os.path.join(args.outdir, "contacts.ply"), con_v, con_n, con_f, None, "water drops, contact discs against the glass (normals -z)")

    # the top face of the plate, with uv (s = x / W + 1/2, t = y / H + 1/2): the renderer flips t, so texture row 0 is y = +H/2
    pv = np.array([[-W / 2, -H / 2, 0], [W / 2, -H / 2, 0], [W / 2, H / 2, 0], [-W / 2, H / 2, 0]], dtype=np.float64)
    puv = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=np.float64)
    pf = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)
    write_ply(os.path.join(args.outdir, "plate_top.ply"), pv, np.tile(np.array([[0.0, 0.0, 1.0]]), (4, 1)), pf, puv, "glass plate, top face with uv for the cutout mask")

    # ---- the mask: alpha 0 inside every drop's contact ellipse ----
    from PIL import Image
    if W >= H:
        rx, ry = args.mask_res, max(1, int(round(args.mask_res * H / W)))
    else:
        ry, rx = args.mask_res, max(1, int(round(args.mask_res * W / H)))
    xs = (np.arange(rx) + 0.5) / rx * W - W / 2
    ys = H / 2 - (np.arange(ry) + 0.5) / ry * H          # row 0 = y = +H/2
    alpha = np.full((ry, rx), 255, dtype=np.uint8)
    for (x, y, a, ax, ay, h, shift) in drops:
        i0, i1 = np.searchsorted(-ys, [-(y + ay), -(y - ay)])
        j0, j1 = np.searchsorted(xs, [x - ax, x + ax])
        i0, i1, j0, j1 = max(i0 - 1, 0), min(i1 + 1, ry), max(j0 - 1, 0), min(j1 + 1, rx)
        yy, xx = np.meshgrid(ys[i0:i1], xs[j0:j1], indexing="ij")
        inside = ((xx - x) / ax) ** 2 + ((yy - y) / ay) ** 2 <= 1.0
        alpha[i0:i1, j0:j1][inside] = 0
    rgba = np.zeros((ry, rx, 4), dtype=np.uint8)
    rgba[..., :3] = 255
    rgba[..., 3] = alpha
    Image.fromarray(rgba, "RGBA").save(os.path.join(args.outdir, "plate_mask.png"))

    with open(os.path.join(args.outdir, "drops.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_mm", "y_mm", "a_mm", "ax_mm", "ay_mm", "height_mm", "shift_mm"])
        for d in drops:
            w.writerow(["%.3f" % (v * 1e3) for v in d])
    print("wrote caps.ply (%d triangles), contacts.ply (%d), plate_top.ply, plate_mask.png (%d x %d), drops.csv in %s" % (
        len(cap_f), len(con_f), rx, ry, args.outdir))


if __name__ == "__main__":
    main()
