"""
Cross-check against Mitsuba 3: the same scenes rendered by GoLD_rayt and by Mitsuba, compared within the Monte Carlo noise.

Each scene is described once (camera, area lights, a uniform environment, spheres and quads, materials) and written both as a
GoLD_rayt scene file and as a Mitsuba scene dictionary. Both renderers draw it twice with independent seeds; the noise of each pixel
is estimated from the two renders, and the images are compared over blocks of BxB pixels:

  z = (Mitsuba - GoLD_rayt) / sqrt(var_M + var_G)     per block, from the means of the two renders
  If the two renderers converge to the same image, |z| > 3 happens in about 0.3 % of the blocks, and the relative differences shrink
  as the spp grows. A systematic difference shows as many blocks with |z| > 3 of one sign, and as a mean ratio away from 1.

Conventions matched here (see docs of each renderer): vertical field of view; a box pixel filter; GGX with alpha = roughness^2 (the
GoLD_rayt "roughness" is sqrt(alpha)); the multiple-scattering compensation of the GoLD_rayt conductor off (Mitsuba has none);
white lights only (rgb (1,1,1) x radiance); a uniform environment of radiance `env`; no limit on the path length in either renderer
(Russian roulette only); a medium boundary with index 1 is a Mitsuba "null" BSDF.

  python tools/test/mitsuba_compare.py --exe out/cli-dev/GoLD_rayt.exe --spp 256 [--scenes a,b] [--variant scalar_rgb]
                                       [--work <folder>] [--out docs/data/mitsuba_compare_results.txt]
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools", "lib"))
from hdr_io import read_pfm                                  # noqa: E402

# ------------------------------------------------------------------------------------------------ scene descriptions
FLOOR = dict(shape="quad", corner=[-2, 0, -2], edge1=[0, 0, 4], edge2=[4, 0, 0], material="floor")
LIGHT = dict(corner=[-0.25, 1.5, -0.25], edge1=[0.5, 0, 0], edge2=[0, 0, 0.5], radiance=12.0)   # faces down
SPHERE = dict(shape="sphere", center=[0, 0.3, 0], radius=0.3, material="obj")
CAMERA = dict(origin=[0, 0.8, 2.2], target=[0, 0.3, 0], up=[0, 1, 0], fov=35.0)


def scene(name, obj, env=0.0, light=True, floor=True, note=""):
    return dict(name=name, camera=CAMERA, env=env, lights=[LIGHT] if light else [],
                materials={"floor": dict(type="diffuse", albedo=[0.5, 0.5, 0.5]), "obj": obj},
                objects=[FLOOR, SPHERE] if floor else [SPHERE], note=note)


MEDIUM = dict(sigmaT=[20.0, 20.0, 20.0], albedo=[0.99, 0.9, 0.7], g=0.5)
DENSE = dict(sigmaT=[200.0, 200.0, 200.0], albedo=[0.995, 0.98, 0.9], g=0.3)
DENSE_GRAY = dict(sigmaT=[200.0, 200.0, 200.0], albedo=[0.98, 0.98, 0.98], g=0.3)
SCENES = [
    scene("diffuse", dict(type="diffuse", albedo=[0.8, 0.5, 0.3]), note="diffuse sphere and floor, one area light"),
    scene("env", dict(type="diffuse", albedo=[0.8, 0.5, 0.3]), env=1.0, light=False, note="the same under a uniform environment of radiance 1"),
    scene("conductor", dict(type="conductor", eta=0.2, k=3.5, roughness=0.5), env=0.2,
          note="rough conductor (GGX, alpha 0.25, n 0.2, k 3.5), area light + environment 0.2"),
    scene("glass", dict(type="dielectric", ior=1.5, roughness=0.0), env=0.2, note="smooth glass sphere (n 1.5), area light + environment 0.2"),
    scene("roughglass", dict(type="dielectric", ior=1.5, roughness=0.3), env=0.2, note="rough glass sphere (alpha 0.09)"),
    scene("medium", dict(type="dielectric", ior=1.5, roughness=0.0, medium=MEDIUM), env=0.2,
          note="smooth glass boundary (n 1.5) filled with a medium: sigma_t 20 /m, albedo (0.99, 0.9, 0.7), HG g 0.5"),
    scene("roughmedium", dict(type="dielectric", ior=1.5, roughness=0.3, medium=MEDIUM), env=0.2, note="the same medium, rough boundary (alpha 0.09)"),
    scene("matchedmedium", dict(type="dielectric", ior=1.0, roughness=0.0, medium=MEDIUM), env=0.2, note="the same medium, index-matched boundary (n 1)"),
    scene("glassfurnace", dict(type="dielectric", ior=1.5, roughness=0.0), env=1.0, light=False, floor=False,
          note="white furnace: smooth glass sphere alone in a uniform environment of radiance 1 (exactly 1 everywhere)"),
    scene("roughglassfurnace", dict(type="dielectric", ior=1.5, roughness=0.3), env=1.0, light=False, floor=False,
          note="white furnace: rough glass sphere (alpha 0.09) alone (below 1 by the energy a single-scattering microfacet model loses)"),
    scene("densemedium", dict(type="dielectric", ior=1.5, roughness=0.0, medium=DENSE), env=0.2,
          note="smooth boundary (n 1.5), dense medium: sigma_t 200 /m, albedo (0.995, 0.98, 0.9), HG g 0.3 (wax-like)"),
    scene("densegray", dict(type="dielectric", ior=1.5, roughness=0.0, medium=DENSE_GRAY), env=0.2,
          note="smooth boundary (n 1.5), dense gray medium: sigma_t 200 /m, albedo 0.98, HG g 0.3 (no weights between channels)"),
    scene("densegraymatched", dict(type="dielectric", ior=1.0, roughness=0.0, medium=DENSE_GRAY), env=0.2,
          note="the same dense gray medium, index-matched boundary (n 1)"),
    scene("densewhitefurnace", dict(type="dielectric", ior=1.0, roughness=0.0, medium=dict(DENSE_GRAY, albedo=[1.0, 1.0, 1.0])), env=1.0,
          light=False, floor=False, note="white furnace: index-matched sphere of a dense non-absorbing medium (exactly 1 everywhere)"),
    scene("densegrayfurnace", dict(type="dielectric", ior=1.0, roughness=0.0, medium=DENSE_GRAY), env=1.0, light=False, floor=False,
          note="furnace: index-matched sphere of the dense gray medium (albedo 0.98): a pixel = the chance that its light escapes"),
    scene("roughdensemedium", dict(type="dielectric", ior=1.5, roughness=0.3, medium=DENSE), env=0.2, note="the same dense medium, rough boundary (alpha 0.09)"),
]


# ------------------------------------------------------------------------------------------------ GoLD_rayt scene file
def write_gold(sc, folder):
    mats = {}
    for name, m in sc["materials"].items():
        if m["type"] == "diffuse":
            mats[name] = {"type": "diffuse", "albedo": m["albedo"]}
        elif m["type"] == "conductor":
            csv = os.path.join(folder, f"nk_{name}.csv")
            with open(csv, "w", newline="\n") as f:     # a flat n, k table (RefractiveIndex.info layout)
                f.write("wl,n\n0.30,{0}\n0.90,{0}\n\nwl,k\n0.30,{1}\n0.90,{1}\n".format(m["eta"], m["k"]))
            mats[name] = {"type": "conductor", "table": csv.replace("\\", "/"), "roughness": m["roughness"], "multipleScattering": False}
        elif m["type"] == "dielectric":
            d = {"type": "dielectric", "ior": m["ior"], "roughness": m["roughness"]}
            if "medium" in m:
                d["medium"] = {"sigmaT": m["medium"]["sigmaT"], "albedo": m["medium"]["albedo"], "g": m["medium"]["g"]}
            mats[name] = d
    objects = []
    for i, L in enumerate(sc["lights"]):
        mats[f"light{i}"] = {"type": "emission", "spectrum": {"rgb": [1, 1, 1]}, "radiance": L["radiance"]}
        objects.append({"shape": "quad", "corner": L["corner"], "edge1": L["edge1"], "edge2": L["edge2"], "material": f"light{i}"})
    for o in sc["objects"]:
        objects.append({k: v for k, v in o.items()})
    cam = sc["camera"]
    doc = {"camera": {"from": cam["origin"], "at": cam["target"], "up": cam["up"], "fov": cam["fov"]},
           "environment": {"type": "uniform", "scale": sc["env"]} if sc["env"] > 0 else "none",
           "materials": mats, "objects": objects}
    path = os.path.join(folder, f"{sc['name']}.json")
    with open(path, "w", newline="\n") as f:
        json.dump(doc, f, indent=1)
    return path


def render_gold(exe, scene_file, out, w, h, spp, seed):
    cmd = [exe, scene_file, "--width", str(w), "--height", str(h), "--spp", str(spp), "--seed", str(seed), "--depth", "100000", "--out", out]
    t = time.time()
    r = subprocess.run(cmd, cwd=os.path.join(ROOT, "GoLD_rayt"), capture_output=True, text=True, errors="replace")
    if r.returncode != 0 or not os.path.exists(out):
        raise RuntimeError(f"GoLD_rayt failed ({r.returncode}):\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    return read_pfm(out).astype(np.float64), time.time() - t, r.stdout


# ------------------------------------------------------------------------------------------------ Mitsuba scene dictionary
def quad_to_world(mi, corner, e1, e2):
    c, e1, e2 = np.array(corner, float), np.array(e1, float), np.array(e2, float)
    n = np.cross(e1, e2)
    n /= np.linalg.norm(n)
    m = np.eye(4)
    m[:3, 0], m[:3, 1], m[:3, 2], m[:3, 3] = e1 / 2, e2 / 2, n, c + e1 / 2 + e2 / 2
    return mi.ScalarTransform4f(m.tolist())


# "volpathmis" (spectral MIS, Miller et al. 2019) rather than "volpath": with a chromatic medium "volpath" weights the channels by the
# density of the one channel it sampled, which gives fireflies in the channel of highest albedo (red pixels up to 117 in densemedium).
MITSUBA_VOLUME_INTEGRATOR = "volpathmis"
# Mitsuba's Russian roulette keeps a path with probability min(throughput, 0.95) at every interaction from rr_depth on, medium
# scatterings included. In a dense medium of high albedo a path scatters thousands of times with a throughput near 1, so it survives
# with about 0.95^n and its weight grows as 1/0.95^n: unbiased in expectation, but the survivors are too rare to be sampled, and a
# finite render comes out dark (9 % on the sphere of densewhitefurnace, which must be 1). The roulette is therefore left out.
MITSUBA_RR_DEPTH = 1 << 30


def mitsuba_dict(mi, sc, w, h, spp):
    def rgb(v):
        return {"type": "rgb", "value": [float(x) for x in v]}

    has_medium = any("medium" in m for m in sc["materials"].values())
    cam = sc["camera"]
    d = {"type": "scene",
         "integrator": {"type": MITSUBA_VOLUME_INTEGRATOR if has_medium else "path", "max_depth": -1, "rr_depth": MITSUBA_RR_DEPTH},
         "sensor": {"type": "perspective", "fov": cam["fov"], "fov_axis": "y",
                    "to_world": mi.ScalarTransform4f().look_at(origin=cam["origin"], target=cam["target"], up=cam["up"]),
                    "sampler": {"type": "independent", "sample_count": spp},
                    "film": {"type": "hdrfilm", "width": w, "height": h, "pixel_format": "rgb", "rfilter": {"type": "box"}}}}
    if sc["env"] > 0:
        d["env"] = {"type": "constant", "radiance": rgb([sc["env"]] * 3)}
    bsdfs = {}
    for name, m in sc["materials"].items():
        if m["type"] == "diffuse":
            bsdfs[name] = {"type": "diffuse", "reflectance": rgb(m["albedo"])}
        elif m["type"] == "conductor":
            alpha = m["roughness"] ** 2
            base = {"eta": rgb([m["eta"]] * 3), "k": rgb([m["k"]] * 3)}
            bsdfs[name] = dict(type="roughconductor", distribution="ggx", alpha=alpha, **base) if alpha > 0 else dict(type="conductor", **base)
        elif m["type"] == "dielectric":
            alpha = m["roughness"] ** 2
            if m["ior"] == 1.0:
                bsdfs[name] = {"type": "null"}
            elif alpha > 0:
                bsdfs[name] = {"type": "roughdielectric", "distribution": "ggx", "alpha": alpha, "int_ior": m["ior"], "ext_ior": 1.0}
            else:
                bsdfs[name] = {"type": "dielectric", "int_ior": m["ior"], "ext_ior": 1.0}
    for i, L in enumerate(sc["lights"]):
        d[f"light{i}"] = {"type": "rectangle", "to_world": quad_to_world(mi, L["corner"], L["edge1"], L["edge2"]),
                          "bsdf": {"type": "diffuse", "reflectance": rgb([0, 0, 0])},
                          "emitter": {"type": "area", "radiance": rgb([L["radiance"]] * 3)}}
    for i, o in enumerate(sc["objects"]):
        m = sc["materials"][o["material"]]
        if o["shape"] == "sphere":
            s = {"type": "sphere", "center": [float(x) for x in o["center"]], "radius": float(o["radius"])}
        else:
            s = {"type": "rectangle", "to_world": quad_to_world(mi, o["corner"], o["edge1"], o["edge2"])}
        s["bsdf"] = bsdfs[o["material"]]
        if "medium" in m:
            md = m["medium"]
            st = max(md["sigmaT"])                       # an RGB value must lie in [0, 1]: normalize it and move the size to "scale"
            s["interior"] = {"type": "homogeneous", "sigma_t": rgb([v / st for v in md["sigmaT"]]), "scale": float(st),
                             "albedo": rgb(md["albedo"]), "phase": {"type": "hg", "g": float(md["g"])}}
        d[f"obj{i}"] = s
    return d


def render_mitsuba(mi, sc, w, h, spp, seed):
    t = time.time()
    s = mi.load_dict(mitsuba_dict(mi, sc, w, h, spp))
    img = np.array(mi.render(s, spp=spp, seed=seed), dtype=np.float64)[..., :3]
    return img, time.time() - t


# ------------------------------------------------------------------------------------------------ comparison
def block_means(img, b):
    h, w = img.shape[0] // b * b, img.shape[1] // b * b
    return img[:h, :w].reshape(h // b, b, w // b, b, -1).mean(axis=(1, 3))


def compare(g, m, b):
    """g, m: (2, H, W, 3) two renders each. Returns a dict of statistics."""
    G, M = g.mean(0), m.mean(0)
    # per-pixel variance of one render from the pair, then the variance of the block mean of the two renders
    var_g = block_means((g[0] - g[1]) ** 2 / 2, b) / (2 * b * b)
    var_m = block_means((m[0] - m[1]) ** 2 / 2, b) / (2 * b * b)
    bg, bm = block_means(G, b), block_means(M, b)
    lum = lambda x: x @ np.array([0.2126, 0.7152, 0.0722])
    var = var_g + var_m
    # blocks with no noise in either renderer (background that sees only the environment) cannot give a z: they are compared
    # directly (largest relative difference), and the z statistics use the noisy, non-black blocks only
    exact = (var.max(axis=-1) <= (1e-7 * np.maximum(bg, bm).max(axis=-1)) ** 2)
    keep = (lum(bg) > 0.02 * lum(bg).max()) & ~exact
    z = (bm - bg) / np.sqrt(var + 1e-30)
    rel = (bm - bg) / np.maximum(bg, 1e-12)
    zk, rk = z[keep], rel[keep]
    mean_g, mean_m = G.reshape(-1, 3).mean(0), M.reshape(-1, 3).mean(0)
    # standard error of the image mean, from the pixel variances
    se = np.sqrt((((g[0] - g[1]) ** 2 / 2).reshape(-1, 3).mean(0) + ((m[0] - m[1]) ** 2 / 2).reshape(-1, 3).mean(0)) / (2 * G.shape[0] * G.shape[1]))
    return dict(mean_g=mean_g, mean_m=mean_m, ratio=mean_m / mean_g, z_mean=(mean_m - mean_g) / se,
                blocks=int(keep.sum()), frac_z3=float((np.abs(zk) > 3).mean()), mean_z=float(zk.mean()),
                med_rel=float(np.median(np.abs(rk))), p95_rel=float(np.percentile(np.abs(rk), 95)),
                noise_rel=float(np.median(np.sqrt(var_g + var_m)[keep] / np.maximum(bg[keep], 1e-12))), z=z, keep=keep,
                exact=int(exact.sum()), exact_rel=float(np.abs(rel[exact]).max()) if exact.any() else 0.0)


def sphere_mask(sc, w, h):
    """Pixels whose center ray (pinhole, vertical fov) hits the scene's sphere."""
    cam = sc["camera"]
    o, t, up = (np.array(cam[k], float) for k in ("origin", "target", "up"))
    fwd = (t - o) / np.linalg.norm(t - o)
    right = np.cross(fwd, up); right /= np.linalg.norm(right)
    upv = np.cross(right, fwd)
    th = np.tan(np.radians(cam["fov"]) / 2)
    ys, xs = np.mgrid[0:h, 0:w]
    u = (2 * (xs + 0.5) / w - 1) * th * w / h
    v = (1 - 2 * (ys + 0.5) / h) * th
    d = fwd[None, None] + u[..., None] * right + v[..., None] * upv
    d /= np.linalg.norm(d, axis=-1, keepdims=True)
    sp = next(x for x in sc["objects"] if x["shape"] == "sphere")
    oc = o - np.array(sp["center"], float)
    b = d @ oc
    return b * b - (oc @ oc - sp["radius"] ** 2) > 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True, help="GoLD_rayt.exe (the RGB build for the RGB variants)")
    ap.add_argument("--variant", default="scalar_rgb")
    ap.add_argument("--spp", type=int, default=256)
    ap.add_argument("--width", type=int, default=320)
    ap.add_argument("--height", type=int, default=240)
    ap.add_argument("--block", type=int, default=16)
    ap.add_argument("--scenes", default="")
    ap.add_argument("--work", default=os.path.join(ROOT, "out", "studies", "mitsuba_compare"))
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    import mitsuba as mi
    mi.set_variant(a.variant)
    exe = os.path.abspath(a.exe)
    a.work = os.path.abspath(a.work)          # GoLD_rayt runs in its project folder, so its paths must be absolute
    os.makedirs(a.work, exist_ok=True)
    chosen = [s for s in SCENES if not a.scenes or s["name"] in a.scenes.split(",")]
    lines = []

    def emit(s=""):
        print(s, flush=True)
        lines.append(s)

    emit(f"GoLD_rayt against Mitsuba {mi.__version__} ({a.variant}). {a.width}x{a.height}, {a.spp} spp, two renders per renderer (seeds 1, 2), blocks of {a.block}x{a.block} pixels.")
    emit("Written by tools/test/mitsuba_compare.py. ratio = Mitsuba / GoLD_rayt (image mean, R G B) with its z; |z|>3 = share of blocks (expected about 0.3 %);")
    emit("mean z = average block z (sign = which is brighter); |rel| = relative block difference (median, 95th percentile); noise = median relative noise of a block difference.")
    emit()
    for sc in chosen:
        folder = os.path.join(a.work, sc["name"])
        os.makedirs(folder, exist_ok=True)
        sf = write_gold(sc, folder)
        g, m, tg, tm = [], [], 0.0, 0.0
        for seed in (1, 2):
            img, t, _ = render_gold(exe, sf, os.path.join(folder, f"gold_{seed}.pfm"), a.width, a.height, a.spp, seed)
            g.append(img); tg += t
            img, t = render_mitsuba(mi, sc, a.width, a.height, a.spp, seed)
            m.append(img); tm += t
            np.save(os.path.join(folder, f"mitsuba_{seed}.npy"), img.astype(np.float32))
        g, m = np.stack(g), np.stack(m)
        r = compare(g, m, a.block)
        mask = sphere_mask(sc, a.width, a.height)
        sg, sm = g.mean(0)[mask].mean(0), m.mean(0)[mask].mean(0)
        sse = np.sqrt((((g[0] - g[1]) ** 2 / 2)[mask].mean(0) + ((m[0] - m[1]) ** 2 / 2)[mask].mean(0)) / (2 * mask.sum()))
        emit(f"--- {sc['name']}: {sc['note']} (GoLD_rayt {tg:.0f} s, Mitsuba {tm:.0f} s) ---")
        emit(f"  image mean  GoLD_rayt ({', '.join(f'{v:.4f}' for v in r['mean_g'])})  Mitsuba ({', '.join(f'{v:.4f}' for v in r['mean_m'])})")
        emit(f"  ratio ({', '.join(f'{v:.4f}' for v in r['ratio'])})  z ({', '.join(f'{v:+.1f}' for v in r['z_mean'])})")
        emit(f"  sphere ({int(mask.sum())} px): ratio ({', '.join(f'{v:.4f}' for v in sm / sg)}) +- ({', '.join(f'{v:.4f}' for v in sse / sg)})  z ({', '.join(f'{v:+.1f}' for v in (sm - sg) / sse)})")
        if r["exact"]:
            emit(f"  {r['exact']} blocks without noise (environment only): largest relative difference {r['exact_rel'] * 100:.3f} %")
        emit(f"  blocks {r['blocks']}: |z|>3 {r['frac_z3'] * 100:.1f} %, mean z {r['mean_z']:+.2f}, |rel| median {r['med_rel'] * 100:.2f} % / 95th {r['p95_rel'] * 100:.2f} %, noise {r['noise_rel'] * 100:.2f} %")
        np.save(os.path.join(folder, "block_z.npy"), r["z"].astype(np.float32))
    if a.out:
        with open(a.out, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(lines) + "\n")
        print("wrote", a.out)


if __name__ == "__main__":
    main()
