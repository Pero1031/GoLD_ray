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
(no Russian roulette in Mitsuba, see MITSUBA_RR_DEPTH); a medium boundary with index 1 is a Mitsuba "null" BSDF.

The spectral scenes (names sp_*) run with the spectral build and Mitsuba's scalar_spectral variant. Every colored input is a spectrum
on a 1 nm grid given to both renderers (lights as tables scaled to the asked luminance, reflectance and n, k as tables, the D65
environment as GoLD_rayt's own table); gray values are flat spectra. The sphere is also compared in XYZ, whose channels stay positive.

  python tools/test/mitsuba_compare.py --exe out/cli-dev/GoLD_rayt.exe --spp 256 [--scenes a,b] [--work <folder>] [--out <results.txt>]
  python tools/test/mitsuba_compare.py --exe out/cli-dev-spectral/GoLD_rayt.exe --variant scalar_spectral --spp 1024 [...]
  Without --scenes, all the scenes of the variant (RGB or spectral) run.
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


def scene(name, obj, env=0.0, light=True, floor=True, note="", spd=None):
    """spd: the light's spectrum for the spectral scenes, ("blackbody", K) or ("d65",); None = an RGB scene (white lights)."""
    return dict(name=name, camera=CAMERA, env=env, lights=[dict(LIGHT, spd=spd)] if light else [], spectral=spd is not None,
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
    # Spectral scenes (run with --variant scalar_spectral and the spectral build). Every colored input is a spectrum given to both
    # renderers on the same 1 nm grid: the two convert RGB to spectra by different methods, which would be a difference of the
    # inputs, not of the rendering. Gray values are flat spectra in both. The environment is D65 (the GoLD_rayt table) of luminance env.
    scene("sp_glassfurnace", dict(type="dielectric", ior=1.5, roughness=0.0), env=1.0, light=False, floor=False, spd=("d65",),
          note="spectral white furnace: smooth glass sphere in a D65 environment of luminance 1 (exactly the environment everywhere)"),
    scene("sp_diffuse", dict(type="diffuse", spectrum="orange"), spd=("blackbody", 2856),
          note="spectral: diffuse sphere with a smooth orange reflectance spectrum, tungsten light (blackbody 2856 K)"),
    scene("sp_gold", dict(type="conductor", metal="Au", roughness=0.0), env=0.2, spd=("blackbody", 1850),
          note="spectral: smooth gold (Johnson & Christy n, k at 1 nm), candle light (blackbody 1850 K) + D65 environment 0.2"),
    scene("sp_roughgold", dict(type="conductor", metal="Au", roughness=0.5), env=0.2, spd=("blackbody", 1850),
          note="spectral: rough gold (GGX alpha 0.25), candle light + D65 environment 0.2"),
    scene("sp_glass", dict(type="dielectric", ior=1.5, roughness=0.0), env=0.2, spd=("d65",), note="spectral: smooth glass, D65 light + environment"),
    scene("sp_roughglass", dict(type="dielectric", ior=1.5, roughness=0.3), env=0.2, spd=("blackbody", 2856),
          note="spectral: rough glass (alpha 0.09), tungsten light + D65 environment 0.2"),
    scene("sp_medium", dict(type="dielectric", ior=1.5, roughness=0.0, medium=dict(sigmaT=[20.0] * 3, albedo=[0.95] * 3, g=0.5)), env=0.2,
          spd=("blackbody", 2856), note="spectral: smooth glass boundary with a gray medium (sigma_t 20 /m, albedo 0.95, g 0.5), tungsten light"),
]

# ------------------------------------------------------------------------------------------------ spectra of the spectral scenes
SP_LAM = np.arange(360.0, 831.0)                             # 1 nm, the range of Mitsuba's CIE tables and of GoLD_rayt's sampling


def sp_light(spd):
    """Unscaled spectrum of a light on SP_LAM (each renderer then scales it to the luminance asked for)."""
    if spd[0] == "blackbody":
        h, c, k = 6.62607015e-34, 2.99792458e8, 1.380649e-23
        lm = SP_LAM * 1e-9
        return 1.0 / (lm ** 5 * (np.exp(h * c / (lm * k * spd[1])) - 1.0))
    return sp_d65()


def sp_d65():
    """GoLD_rayt's D65 table (380-780 nm, 5 nm, linear, 0 outside), read from its source by metal_accuracy.py."""
    sys.path.insert(0, os.path.join(ROOT, "tools", "studies"))
    import metal_accuracy as ma
    return np.asarray(ma.d65(SP_LAM), dtype=np.float64)


def sp_reflectance(name):
    if name == "orange":                                     # a smooth reflectance: low in the blue, high in the red
        return 0.12 + 0.72 / (1.0 + np.exp(-(SP_LAM - 585.0) / 18.0))
    raise ValueError(name)


def sp_metal_nk(metal):
    """n, k of the metal's Johnson & Christy table (the one GoLD_rayt uses) interpolated linearly to SP_LAM."""
    path = os.path.join(ROOT, "GoLD_rayt", "assets", "Metal_data", f"Johnson_{metal}.csv")
    blocks, cur = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        if line.startswith("wl,"):
            cur = line.split(",")[1]
            blocks[cur] = []
            continue
        a, b = (float(x) for x in line.split(","))
        blocks[cur].append((a * 1000.0, b))
    n, k = (np.array(sorted(blocks[key])) for key in ("n", "k"))
    return np.interp(SP_LAM, n[:, 0], n[:, 1]), np.interp(SP_LAM, k[:, 0], k[:, 1])


def write_csv(path, values):
    """A two-column table (wavelength in nm, value) for GoLD_rayt; returns the path with forward slashes."""
    with open(path, "w", newline="\n") as f:
        f.write("wavelength_nm,value\n" + "".join(f"{l:.0f},{v:.9g}\n" for l, v in zip(SP_LAM, values)))
    return path.replace(os.sep, "/")


# ------------------------------------------------------------------------------------------------ GoLD_rayt scene file
def write_gold(sc, folder):
    mats = {}
    for name, m in sc["materials"].items():
        if m["type"] == "diffuse" and "spectrum" in m:
            mats[name] = {"type": "diffuse", "reflectance": write_csv(os.path.join(folder, f"refl_{name}.csv"), sp_reflectance(m["spectrum"]))}
        elif m["type"] == "diffuse":
            mats[name] = {"type": "diffuse", "albedo": m["albedo"]}
        elif m["type"] == "conductor" and "metal" in m:
            n, k = sp_metal_nk(m["metal"])
            csv = os.path.join(folder, f"nk_{name}.csv")
            with open(csv, "w", newline="\n") as f:     # RefractiveIndex.info layout, micrometres
                f.write("wl,n\n" + "".join(f"{l / 1000:.4f},{v:.6g}\n" for l, v in zip(SP_LAM, n)))
                f.write("\nwl,k\n" + "".join(f"{l / 1000:.4f},{v:.6g}\n" for l, v in zip(SP_LAM, k)))
            mats[name] = {"type": "conductor", "table": csv.replace(os.sep, "/"), "roughness": m["roughness"], "multipleScattering": False}
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
        spectrum = {"rgb": [1, 1, 1]} if L.get("spd") is None else {"file": write_csv(os.path.join(folder, f"light{i}.csv"), sp_light(L["spd"]))}
        mats[f"light{i}"] = {"type": "emission", "spectrum": spectrum, "radiance": L["radiance"]}
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

    def spectrum(values):
        return {"type": "regular", "wavelength_min": float(SP_LAM[0]), "wavelength_max": float(SP_LAM[-1]),
                "values": ", ".join(f"{v:.9g}" for v in values)}

    def luminous(values, Y):
        """values scaled so that Mitsuba's film reads luminance Y (its CIE Y integral; a flat spectrum of 1 reads 1)."""
        ybar = np.array([float(mi.cie1931_y(float(l))) for l in SP_LAM])
        return spectrum(np.asarray(values) * Y / ((np.asarray(values) * ybar).sum() / 106.856895))

    spectral = sc.get("spectral", False)
    if spectral and "spectral" not in mi.variant():
        raise RuntimeError(f"{sc['name']} is a spectral scene: run it with --variant scalar_spectral and the spectral build")

    def gray(v):
        return {"type": "uniform", "value": float(v)} if spectral else rgb([v] * 3)

    has_medium = any("medium" in m for m in sc["materials"].values())
    cam = sc["camera"]
    d = {"type": "scene",
         "integrator": {"type": MITSUBA_VOLUME_INTEGRATOR if has_medium else "path", "max_depth": -1, "rr_depth": MITSUBA_RR_DEPTH},
         "sensor": {"type": "perspective", "fov": cam["fov"], "fov_axis": "y",
                    "to_world": mi.ScalarTransform4f().look_at(origin=cam["origin"], target=cam["target"], up=cam["up"]),
                    "sampler": {"type": "independent", "sample_count": spp},
                    "film": {"type": "hdrfilm", "width": w, "height": h, "pixel_format": "rgb", "rfilter": {"type": "box"}}}}
    if sc["env"] > 0:
        d["env"] = {"type": "constant", "radiance": luminous(sp_d65(), sc["env"]) if spectral else rgb([sc["env"]] * 3)}
    bsdfs = {}
    for name, m in sc["materials"].items():
        if m["type"] == "diffuse":
            if "spectrum" in m:
                refl = spectrum(sp_reflectance(m["spectrum"]))
            else:
                refl = gray(m["albedo"][0]) if spectral and len(set(m["albedo"])) == 1 else rgb(m["albedo"])
            bsdfs[name] = {"type": "diffuse", "reflectance": refl}
        elif m["type"] == "conductor":
            alpha = m["roughness"] ** 2
            if "metal" in m:
                n, k = sp_metal_nk(m["metal"])
                base = {"eta": spectrum(n), "k": spectrum(k)}
            else:
                base = {"eta": gray(m["eta"]), "k": gray(m["k"])}
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
                          "bsdf": {"type": "diffuse", "reflectance": gray(0.0)},
                          "emitter": {"type": "area", "radiance": luminous(sp_light(L["spd"]), L["radiance"]) if spectral else rgb([L["radiance"]] * 3)}}
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
            if spectral:
                assert len(set(md["sigmaT"])) == 1 and len(set(md["albedo"])) == 1, "spectral scenes use gray media"
                s["interior"] = {"type": "homogeneous", "sigma_t": gray(1.0), "scale": float(st), "albedo": gray(md["albedo"][0]),
                                 "phase": {"type": "hg", "g": float(md["g"])}}
            else:
                s["interior"] = {"type": "homogeneous", "sigma_t": rgb([v / st for v in md["sigmaT"]]), "scale": float(st),
                                 "albedo": rgb(md["albedo"]), "phase": {"type": "hg", "g": float(md["g"])}}
        d[f"obj{i}"] = s
    return d


def render_mitsuba(mi, sc, w, h, spp, seed):
    t = time.time()
    s = mi.load_dict(mitsuba_dict(mi, sc, w, h, spp))
    img = np.array(mi.render(s, spp=spp, seed=seed), dtype=np.float64)[..., :3]
    return img, time.time() - t


def mitsuba_spectral_scale(mi):
    """Factor that removes the offset of Mitsuba's spectral film: it divides by MI_CIE_Y_NORMALIZATION = 1 / 106.7503, while its own
    y-bar table integrates to 106.857 over [360, 830] nm, so every spectral render reads 1.001 times the luminance of its spectrum
    (a flat spectrum of 1 reads Y = 1.001). 1 for the RGB variants."""
    if "spectral" not in mi.variant():
        return 1.0
    lam = np.arange(float(mi.MI_CIE_MIN), float(mi.MI_CIE_MAX) + 0.5, 1.0)
    ybar = np.array([float(mi.cie1931_y(float(l))) for l in lam])
    integral = float(((ybar[1:] + ybar[:-1]) / 2).sum())
    return 1.0 / (integral * float(mi.MI_CIE_Y_NORMALIZATION))


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
    ap.add_argument("--reuse", action="store_true", help="read the renders already in the work folder instead of rendering again")
    a = ap.parse_args()
    import mitsuba as mi
    mi.set_variant(a.variant)
    m_scale = mitsuba_spectral_scale(mi)
    exe = os.path.abspath(a.exe)
    a.work = os.path.abspath(a.work)          # GoLD_rayt runs in its project folder, so its paths must be absolute
    os.makedirs(a.work, exist_ok=True)
    spectral_run = "spectral" in a.variant           # by default, the scenes that match the variant
    chosen = [s for s in SCENES if (s["name"] in a.scenes.split(",") if a.scenes else s["spectral"] == spectral_run)]
    lines = []

    def emit(s=""):
        print(s, flush=True)
        lines.append(s)

    emit(f"GoLD_rayt against Mitsuba {mi.__version__} ({a.variant}). {a.width}x{a.height}, {a.spp} spp, two renders per renderer (seeds 1, 2), blocks of {a.block}x{a.block} pixels.")
    emit("Written by tools/test/mitsuba_compare.py. ratio = Mitsuba / GoLD_rayt (image mean, R G B) with its z; |z|>3 = share of blocks (expected about 0.3 %);")
    emit("mean z = average block z (sign = which is brighter); |rel| = relative block difference (median, 95th percentile); noise = median relative noise of a block difference.")
    to_xyz = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])
    if spectral_run:
        emit(f"Spectral: Mitsuba's images are multiplied by {m_scale:.6f}, which removes the offset of its spectral film (it normalizes by 1/106.7503")
        emit("while its own y-bar integrates to 106.857, so every spectral render reads 1.001 times the luminance). Image means and blocks are compared in XYZ.")
    emit()
    for sc in chosen:
        folder = os.path.join(a.work, sc["name"])
        os.makedirs(folder, exist_ok=True)
        sf = write_gold(sc, folder)
        g, m, tg, tm = [], [], 0.0, 0.0
        for seed in (1, 2):
            gp, mp = os.path.join(folder, f"gold_{seed}.pfm"), os.path.join(folder, f"mitsuba_{seed}.npy")
            if a.reuse and os.path.exists(gp) and os.path.exists(mp):
                g.append(read_pfm(gp).astype(np.float64)); m.append(np.load(mp).astype(np.float64))
                continue
            img, t, _ = render_gold(exe, sf, gp, a.width, a.height, a.spp, seed)
            g.append(img); tg += t
            img, t = render_mitsuba(mi, sc, a.width, a.height, a.spp, seed)
            m.append(img); tm += t
            np.save(mp, img.astype(np.float32))                # as Mitsuba wrote it (the spectral scale is applied below)
        g, m = np.stack(g), np.stack(m) * m_scale
        r = compare(g @ to_xyz.T, m @ to_xyz.T, a.block) if spectral_run else compare(g, m, a.block)
        mask = sphere_mask(sc, a.width, a.height)
        sg, sm = g.mean(0)[mask].mean(0), m.mean(0)[mask].mean(0)
        sse = np.sqrt((((g[0] - g[1]) ** 2 / 2)[mask].mean(0) + ((m[0] - m[1]) ** 2 / 2)[mask].mean(0)) / (2 * mask.sum()))
        emit(f"--- {sc['name']}: {sc['note']} (GoLD_rayt {tg:.0f} s, Mitsuba {tm:.0f} s) ---")
        emit(f"  image mean  GoLD_rayt ({', '.join(f'{v:.4f}' for v in r['mean_g'])})  Mitsuba ({', '.join(f'{v:.4f}' for v in r['mean_m'])})")
        emit(f"  ratio ({', '.join(f'{v:.4f}' for v in r['ratio'])})  z ({', '.join(f'{v:+.1f}' for v in r['z_mean'])})")
        emit(f"  sphere ({int(mask.sum())} px): ratio ({', '.join(f'{v:.4f}' for v in sm / sg)}) +- ({', '.join(f'{v:.4f}' for v in sse / sg)})  z ({', '.join(f'{v:+.1f}' for v in (sm - sg) / sse)})")
        # the same in XYZ, whose channels stay positive for colors outside sRGB (a linear sRGB channel of such a color can be near 0 or negative)
        gx, mx = np.stack([x[mask] @ to_xyz.T for x in g]), np.stack([x[mask] @ to_xyz.T for x in m])
        xg, xm = gx.mean(axis=(0, 1)), mx.mean(axis=(0, 1))
        xse = np.sqrt((((gx[0] - gx[1]) ** 2 / 2).mean(0) + ((mx[0] - mx[1]) ** 2 / 2).mean(0)) / (2 * mask.sum()))
        emit(f"  sphere XYZ: ratio ({', '.join(f'{v:.4f}' for v in xm / xg)}) +- ({', '.join(f'{v:.4f}' for v in xse / xg)})  z ({', '.join(f'{v:+.1f}' for v in (xm - xg) / xse)})")
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
