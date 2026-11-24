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
    args = ap.parse_args()
    {"sheet": sheet, "pigments": pigments, "split": split}[args.cmd](args)


if __name__ == "__main__":
    main()
