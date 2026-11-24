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
"""
import argparse
import os
import sys

import colour
import numpy as np
from PIL import Image
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
    args = ap.parse_args()
    sheet(args) if args.cmd == "sheet" else pigments(args)


if __name__ == "__main__":
    main()
