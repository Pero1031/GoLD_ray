"""
Prepares a scanned leaf (photogrammetry) for the leaf material (docs/LEAF.md).

  python tools/data/leaf_scan.py sheet in.obj out.obj [--axis z]
      A scan closes each leaf into a thin shell (an upper and a lower surface). The leaf material is a single sheet that both
      reflects and transmits, so only the faces whose normal points along +axis are kept (the upper surface).

  python tools/data/leaf_scan.py pigments texture.jpg out.png [--size 2048] [--N 1.6]
      Estimates the pigments of every texel from the scan's colour: PROSPECT's diffuse reflectance (without the first face's
      reflection: a scan's texture is taken without gloss) for a grid of chlorophyll, carotenoid,
      anthocyanin and brown-pigment contents, as linear sRGB under D65, and for each texel the weighted mean of the 8 grid points
      nearest in CIELAB (brown and carotenoids far from 8 ug/cm^2 cost a little, to choose among colours that several mixtures
      give). The scan's colour is the reflectance times an unknown factor k (exposure, white level); k is the one that fits the
      leaf texels best. Writes an RGBA PNG: chlorophyll / 80, carotenoids / 20, anthocyanins / 40, brown pigments (0..1), in
      ug/cm^2, which the leaf material reads as "pigmentMap".

  python tools/data/leaf_scan.py split scan.obj sheet.obj leaves.obj twig.obj [--width 0.25] [--pixel 0.01] [--axis z]
      Separates the twig and the petioles from the leaves, so that they can have a material of their own (bark is opaque).
      The sheet is drawn from above into a mask (--pixel per pixel, in the mesh's units). The parts of it narrower than
      --width are the thin parts (a morphological opening); of these, the large connected pieces are the twig and the
      petioles, the small ones are tips of lobes and teeth, which stay leaf. leaves.obj gets the sheet's faces outside the
      twig; twig.obj gets the scan's faces (both surfaces: a closed tube) over it. With --each, every leaf is also written to a file
      of its own (leaves_0.obj, leaves_1.obj ..., from left to right), so that each can have its own colour.

  python tools/data/leaf_scan.py extend twig.obj branch.obj [--end -x] [--toward -0.6 -0.8 0] [--length 14] [--thicken 3]
      A scanned sprig ends in a cut. This continues it: a tube that starts inside the twig's end with the twig's radius and
      direction (fitted to its last --fit units), thickens and turns towards --toward, with a slow wobble and slightly uneven
      bark, so that the branch can run out of the picture. Smooth normals and uv (around, along).
"""
import argparse
import os
import sys

import colour
import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools", "lib"))
from prospect import load_table, prospect, tav  # noqa: E402

Image.MAX_IMAGE_PIXELS = None
SCALE = np.array([80.0, 20.0, 40.0, 1.0])   # chl, car, ant, brown -> 0..1 in the PNG


def read_obj(path):
    v, vt, f = [], [], []
    for line in open(path):
        if line.startswith("v "):
            v.append([float(x) for x in line.split()[1:4]])
        elif line.startswith("vt "):
            vt.append([float(x) for x in line.split()[1:3]])
        elif line.startswith("f "):
            f.append([[int(p) if p else 0 for p in t.split("/")[:2]] for t in line.split()[1:4]])
    return np.array(v), np.array(vt) if vt else None, np.array(f) - 1


def write_obj(path, v, vt, f, comment):
    used_v, inv_v = np.unique(f[:, :, 0], return_inverse=True)
    with open(path, "w") as out:
        out.write(f"# {comment}\n")
        np.savetxt(out, v[used_v], fmt="v %.6f %.6f %.6f")
        if vt is not None:
            used_t, inv_t = np.unique(f[:, :, 1], return_inverse=True)
            np.savetxt(out, vt[used_t], fmt="vt %.6f %.6f")
            idx = np.stack([inv_v.reshape(-1, 3), inv_t.reshape(-1, 3)], axis=2).reshape(-1, 6) + 1
            np.savetxt(out, idx, fmt="f %d/%d %d/%d %d/%d")
        else:
            np.savetxt(out, inv_v.reshape(-1, 3) + 1, fmt="f %d %d %d")


def split(args):
    from PIL import ImageDraw
    v, vt, f = read_obj(args.scan)
    sv, svt, sf = read_obj(args.sheet)
    a, b = [i for i in range(3) if i != {"x": 0, "y": 1, "z": 2}[args.axis]]
    lo = np.minimum(v[:, [a, b]].min(0), sv[:, [a, b]].min(0)) - 4 * args.pixel
    hi = np.maximum(v[:, [a, b]].max(0), sv[:, [a, b]].max(0)) + 4 * args.pixel
    size = np.ceil((hi - lo) / args.pixel).astype(int)

    def px(p):
        return (p[..., [a, b]] - lo) / args.pixel

    img = Image.new("L", (int(size[0]), int(size[1])), 0)
    draw = ImageDraw.Draw(img)
    for tri in px(sv[sf[:, :, 0]]):
        draw.polygon([tuple(q) for q in tri], fill=1)
    mask = np.array(img, dtype=bool)
    r = 0.5 * args.width / args.pixel
    inside = ndimage.distance_transform_edt(mask)
    opened = ndimage.distance_transform_edt(~(inside > r)) <= r
    thin = mask & ~opened
    lab, n = ndimage.label(thin, structure=np.ones((3, 3)))
    area = ndimage.sum(thin, lab, index=np.arange(1, n + 1))
    min_area = args.min_length * args.width / args.pixel ** 2
    keep = np.zeros(n + 1, dtype=bool)
    keep[1:] = area >= min_area
    twig = ndimage.binary_dilation(keep[lab], iterations=2)

    def over_twig(pts, fs):
        c = px(pts[fs[:, :, 0]].mean(axis=1)).astype(int)
        c = np.clip(c, 0, [twig.shape[1] - 1, twig.shape[0] - 1])
        return twig[c[:, 1], c[:, 0]]

    t_sheet = over_twig(sv, sf)
    t_scan = over_twig(v, f)
    write_obj(args.leaves, sv, svt, sf[~t_sheet], f"leaves of {os.path.basename(args.sheet)}: {(~t_sheet).sum()} faces")
    if args.each:
        # one file per leaf: the pieces of the sheet that hang together once the twig is gone (vertices at the same place are
        # joined); pieces smaller than 2 % of the largest join the piece whose centre is nearest
        fs = sf[~t_sheet]
        _, vid = np.unique(np.round(sv, 5), axis=0, return_inverse=True)
        vid = vid.reshape(-1)
        parent = np.arange(vid.max() + 1)

        def root(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for tri in vid[fs[:, :, 0]]:
            r0 = root(tri[0])
            for q in tri[1:]:
                parent[root(q)] = r0
        face_lab = np.array([root(q) for q in vid[fs[:, 0, 0]]])
        labs, counts = np.unique(face_lab, return_counts=True)
        big = labs[counts >= 0.02 * counts.max()]
        centres = {l: sv[fs[face_lab == l][:, :, 0]].mean(axis=(0, 1)) for l in labs}
        for l in labs:
            if l not in big:
                face_lab[face_lab == l] = min(big, key=lambda b: np.linalg.norm(centres[b] - centres[l]))
        cx = {l: centres[l][a] for l in big}
        stem, ext = os.path.splitext(args.leaves)
        for k, l in enumerate(sorted(big, key=lambda q: cx[q])):
            sel = fs[face_lab == l]
            centre = sv[sel[:, :, 0]].mean(axis=(0, 1))
            write_obj(f"{stem}_{k}{ext}", sv, svt, sel, f"leaf {k} of {os.path.basename(args.sheet)}: {len(sel)} faces")
            print(f"  leaf {k}: {len(sel)} faces, centre {np.round(centre, 2).tolist()}")
    write_obj(args.twig, v, vt, f[t_scan], f"twig and petioles of {os.path.basename(args.scan)}: {t_scan.sum()} faces")
    print(f"{int(keep.sum())} thin pieces of {n} kept as twig; leaves {(~t_sheet).sum()} faces, twig {t_scan.sum()} faces")
    if args.preview:
        rgb = np.zeros(mask.shape + (3,), np.uint8)
        rgb[mask] = (200, 60, 40)
        rgb[thin] = (240, 200, 60)
        rgb[twig] = (60, 200, 255)
        Image.fromarray(rgb[::-1]).save(args.preview)


def extend(args):
    """A tapered, gently curved tube that continues the twig from its cut end (so the branch can leave the picture)."""
    v, _, _ = read_obj(args.twig)
    axis = {"x": 0, "y": 1, "z": 2}[args.end[1]]
    sign = -1.0 if args.end[0] == "-" else 1.0
    tip = v[np.argmax(sign * v[:, axis])]
    near = v[np.linalg.norm(v - tip, axis=1) < args.fit]
    c = near.mean(axis=0)
    d0 = np.linalg.svd(near - c)[2][0]
    if np.dot(d0, tip - c) < 0:
        d0 = -d0
    # the end of the twig's axis, and its radius: distance of the last points from the axis
    last = near[np.dot(near - c, d0) > np.dot(tip - c, d0) - 0.3]
    rel = last - c
    r0 = float(np.percentile(np.linalg.norm(rel - np.outer(rel @ d0, d0), axis=1), 75))
    p0 = c + d0 * float(np.max((near - c) @ d0))
    toward = np.array(args.toward, float)
    toward /= np.linalg.norm(toward)

    rng = np.random.default_rng(args.seed)
    L = args.length
    rings = int(np.ceil((L + args.overlap) / args.step)) + 1
    s_all = np.linspace(-args.overlap, L, rings)
    # direction: from the twig's own towards `toward`, smoothly; plus a slow wobble
    wob = [rng.normal(size=3) for _ in range(4)]
    pts, dirs = [p0 - d0 * args.overlap], []
    for i, s in enumerate(s_all):
        t = max(0.0, s) / L
        w = t * t * (3 - 2 * t)
        d = (1 - w) * d0 + w * toward
        for k, a in enumerate(wob):
            d = d + args.wobble * a * np.sin((k + 1) * np.pi * t + k)
        d /= np.linalg.norm(d)
        dirs.append(d)
        if i > 0:
            pts.append(pts[-1] + d * (s_all[i] - s_all[i - 1]))
    pts, dirs = np.array(pts), np.array(dirs)
    # frames by parallel transport
    n = np.cross(dirs[0], [0.0, 0.0, 1.0] if abs(dirs[0][2]) < 0.9 else [1.0, 0.0, 0.0])
    n /= np.linalg.norm(n)
    frames = []
    for d in dirs:
        n = n - d * np.dot(n, d)
        n /= np.linalg.norm(n)
        frames.append((n, np.cross(d, n)))
    sides = args.sides
    ang = 2 * np.pi * np.arange(sides) / sides
    # bark: ridges along the branch (higher harmonics around it, drifting slowly along it) and lumps
    ridges = [(k, rng.normal(), rng.uniform(0, 2 * np.pi), rng.normal() * 0.6) for k in range(3, 14)]
    lumps = [(k, rng.normal(), rng.uniform(0, 2 * np.pi), rng.normal() * 3.0) for k in range(1, 4)]
    # nodes: swellings where buds and side shoots were, every --node units on average
    nodes = []
    s_node = args.node * rng.uniform(0.3, 1.0)
    while s_node < L:
        nodes.append((s_node, rng.uniform(0.6, 1.0) * args.node_swell, rng.uniform(0, 2 * np.pi)))
        s_node += args.node * rng.uniform(0.6, 1.4)
    verts, uvs = [], []
    for i, s in enumerate(s_all):
        t = max(0.0, s) / L
        r = r0 * (1 + (args.thicken - 1) * t ** args.grow)
        n1, n2 = frames[i]
        for j, a in enumerate(ang):
            ridge = sum(c * np.sin(k * a + p + drift * s / max(r, 1e-6) * 0.05) for k, c, p, drift in ridges) / np.sqrt(len(ridges))
            lump = sum(c * np.sin(k * a + p + drift * s) for k, c, p, drift in lumps) / np.sqrt(len(lumps))
            node = sum(sw * np.exp(-((s - sn) / (2.0 * r)) ** 2) * (0.6 + 0.4 * np.cos(a - pa)) for sn, sw, pa in nodes)
            rr = r * (1 + args.rough * ridge * min(1.0, 4 * t) + 0.5 * args.rough * lump + node)
            verts.append(pts[i] + rr * (np.cos(a) * n1 + np.sin(a) * n2))
            uvs.append((j / sides, s / L))
    verts = np.array(verts)
    # smooth normals from the faces (the bark's relief included)
    faces = []
    for i in range(rings - 1):
        for j in range(sides):
            a_, b_ = i * sides + j, i * sides + (j + 1) % sides
            faces += [(a_, b_, b_ + sides), (a_, b_ + sides, a_ + sides)]
    faces = np.array(faces)
    fn = np.cross(verts[faces[:, 1]] - verts[faces[:, 0]], verts[faces[:, 2]] - verts[faces[:, 0]])
    norms = np.zeros_like(verts)
    for c in range(3):
        np.add.at(norms, faces[:, c], fn)
    norms /= np.maximum(np.linalg.norm(norms, axis=1, keepdims=True), 1e-12)
    centre = np.repeat(pts, sides, axis=0)
    if np.mean(np.sum(norms * (verts - centre), axis=1)) < 0:
        norms = -norms
    with open(args.output, "w") as out:
        out.write(f"# branch continuing {os.path.basename(args.twig)} from its end {np.round(p0, 3).tolist()}, radius {r0:.4f} -> {r0 * args.thicken:.4f}, length {L}\n")
        np.savetxt(out, np.array(verts), fmt="v %.6f %.6f %.6f")
        np.savetxt(out, np.array(norms), fmt="vn %.6f %.6f %.6f")
        np.savetxt(out, np.array(uvs), fmt="vt %.6f %.6f")
        f = faces + 1
        np.savetxt(out, np.stack([f, f, f], axis=2).reshape(-1, 9), fmt="f %d/%d/%d %d/%d/%d %d/%d/%d")
    print(f"end {np.round(p0, 3).tolist()}, direction {np.round(d0, 3).tolist()}, radius {r0:.4f} -> {args.output} ({rings} rings)")


def sheet(args):
    v, vt, f = read_obj(args.input)
    axis = {"x": 0, "y": 1, "z": 2}[args.axis]
    p = v[f[:, :, 0]]
    n = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    keep = n[:, axis] > 0
    f = f[keep]
    used_v, inv_v = np.unique(f[:, :, 0], return_inverse=True)
    with open(args.output, "w") as out:
        out.write(f"# upper surface (+{args.axis}) of {os.path.basename(args.input)}: {keep.sum()} of {len(keep)} faces\n")
        np.savetxt(out, v[used_v], fmt="v %.6f %.6f %.6f")
        if vt is not None:
            used_t, inv_t = np.unique(f[:, :, 1], return_inverse=True)
            np.savetxt(out, vt[used_t], fmt="vt %.6f %.6f")
            idx = np.stack([inv_v.reshape(-1, 3), inv_t.reshape(-1, 3)], axis=2).reshape(-1, 6) + 1
            np.savetxt(out, idx, fmt="f %d/%d %d/%d %d/%d")
        else:
            np.savetxt(out, inv_v.reshape(-1, 3) + 1, fmt="f %d %d %d")
    print(f"kept {keep.sum()} of {len(keep)} faces -> {args.output}")


def grid_colours(N):
    t = load_table()
    lam = t["lambda"]
    vis = (lam >= 380) & (lam <= 780)
    cmf = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"]
    d65 = colour.SDS_ILLUMINANTS["D65"]
    W = np.array([cmf[x] for x in lam[vis]]) * np.array([d65[x] for x in lam[vis]])[:, None]
    W /= W[:, 1].sum()
    rface = 1 - tav(90, t["nrefrac"])   # the first face's reflection: a scan's texture has no gloss, it shows the light from inside
    params, xyz = [], []
    for chl in (0, 1, 2, 4, 7, 10, 15, 20, 30, 40, 55, 70):
        for car in (1, 3, 5, 8, 12, 16):
            for ant in (0, 0.5, 1, 2, 4, 6, 9, 12, 16, 20, 26, 32, 40):
                for brown in (0, 0.1, 0.25, 0.5, 0.8):
                    R, _ = prospect(t, N=N, chl=chl, car=car, ant=ant, brown=brown, ewt=0.012, lma=0.006, alpha=90)
                    params.append((chl, car, ant, brown))
                    xyz.append((R - rface)[vis] @ W)
    return np.array(params, dtype=float), np.array(xyz), colour.XYZ_to_xy(W.sum(0))


def pigments(args):
    im = Image.open(args.input).convert("RGB")
    if args.size and im.size[0] > args.size:
        im = im.resize((args.size, args.size * im.size[1] // im.size[0]), Image.LANCZOS)
    srgb = np.asarray(im).reshape(-1, 3)
    uniq, inv = np.unique(srgb, axis=0, return_inverse=True)
    inv = inv.ravel()
    print(f"{len(srgb)} texels, {len(uniq)} distinct colours")
    lin = colour.cctf_decoding(uniq / 255.0)
    xyz_tex = colour.RGB_to_XYZ(lin, "sRGB", apply_cctf_decoding=False)

    params, xyz_grid, white = grid_colours(args.N)
    print(f"{len(params)} PROSPECT mixtures (N = {args.N})")
    # extra coordinates: what a mixture costs besides its colour (brown pigments, carotenoids far from 8)
    penalty = np.column_stack([6.0 * params[:, 3], 0.4 * (params[:, 1] - 8.0)])

    def fit(k):
        lab_t = colour.XYZ_to_Lab(xyz_tex / k, white)
        lab_g = colour.XYZ_to_Lab(xyz_grid, white)
        tree = cKDTree(np.column_stack([lab_g, penalty]))
        d, i = tree.query(np.column_stack([lab_t, np.zeros((len(lab_t), 2))]), k=8)
        return d, i

    # the leaf texels: reddish or greenish, not the dark gaps of the atlas
    weight = np.bincount(inv, minlength=len(uniq))
    leafy = (lin.max(1) > 0.02)
    best = None
    for k in (0.6, 0.8, 1.0, 1.25, 1.6, 2.0, 2.5, 3.2, 4.0):
        d, _ = fit(k)
        e = np.average(d[leafy, 0], weights=weight[leafy])
        print(f"  k = {k:4.2f}: mean distance to the nearest mixture {e:.2f}")
        if best is None or e < best[0]:
            best = (e, k)
    k = best[1]
    d, i = fit(k)
    w = 1.0 / np.maximum(d, 0.3) ** 2
    est = (params[i] * w[..., None]).sum(1) / w.sum(1)[:, None]
    out = np.clip(est / SCALE, 0, 1)[inv].reshape(im.size[1], im.size[0], 4)
    Image.fromarray((out * 255 + 0.5).astype(np.uint8)).save(args.output)
    q = est[inv][np.repeat(leafy, 1)[inv]]
    print(f"k = {k}; leaf texels: chlorophyll {np.median(q[:, 0]):.1f}, carotenoids {np.median(q[:, 1]):.1f}, "
          f"anthocyanins {np.median(q[:, 2]):.1f}, brown {np.median(q[:, 3]):.2f} (medians) -> {args.output}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("sheet")
    s.add_argument("input")
    s.add_argument("output")
    s.add_argument("--axis", default="z", choices=["x", "y", "z"])
    p = sub.add_parser("pigments")
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--size", type=int, default=2048)
    p.add_argument("--N", type=float, default=1.6)
    t = sub.add_parser("split")
    t.add_argument("scan")
    t.add_argument("sheet")
    t.add_argument("leaves")
    t.add_argument("twig")
    t.add_argument("--width", type=float, default=0.25, help="parts narrower than this are thin (mesh units)")
    t.add_argument("--min-length", type=float, default=1.5, help="thin pieces shorter than this stay leaf (mesh units)")
    t.add_argument("--pixel", type=float, default=0.01)
    t.add_argument("--axis", default="z", choices=["x", "y", "z"])
    t.add_argument("--each", action="store_true", help="also write each leaf on its own: leaves_0.obj, leaves_1.obj ... (left to right)")
    t.add_argument("--preview", help="PNG of the mask: leaf red, thin yellow, twig blue")
    e = sub.add_parser("extend")
    e.add_argument("twig")
    e.add_argument("output")
    e.add_argument("--end", default="-x", choices=["-x", "+x", "-y", "+y", "-z", "+z"], help="the cut end is the twig's extreme point this way")
    e.add_argument("--toward", type=float, nargs=3, default=[-0.6, -0.8, 0.0], help="direction the branch turns to (mesh coordinates)")
    e.add_argument("--length", type=float, default=14.0)
    e.add_argument("--thicken", type=float, default=3.0, help="radius at the far end / radius at the cut")
    e.add_argument("--fit", type=float, default=1.0, help="the twig's direction is fitted to its points this near the end")
    e.add_argument("--overlap", type=float, default=0.15, help="the tube starts this far inside the twig")
    e.add_argument("--wobble", type=float, default=0.08)
    e.add_argument("--grow", type=float, default=0.8, help="radius grows as (s / length)^grow: below 1 it thickens quickly near the cut")
    e.add_argument("--rough", type=float, default=0.06, help="relief of the bark (fraction of the radius)")
    e.add_argument("--node", type=float, default=3.0, help="mean distance between nodes (swellings)")
    e.add_argument("--node-swell", type=float, default=0.15, help="how much a node swells (fraction of the radius)")
    e.add_argument("--step", type=float, default=0.05)
    e.add_argument("--sides", type=int, default=20)
    e.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    {"sheet": sheet, "pigments": pigments, "split": split, "extend": extend}[args.cmd](args)


if __name__ == "__main__":
    main()
