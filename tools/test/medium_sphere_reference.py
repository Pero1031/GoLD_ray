"""
An independent reference for the furnace scenes of mitsuba_compare.py with an index-matched medium sphere: neither renderer's code.

In a uniform environment of radiance 1 an index-matched sphere of a homogeneous medium shows, in each pixel, the probability that light
entering along the pixel's ray leaves the sphere again without being absorbed (reciprocity; there is no boundary to reflect or bend).
This script follows that light with a plain random walk (exponential free flights with sigma_t, Henyey-Greenstein directions, the weight
multiplied by the albedo at every scattering, Russian roulette below a small weight), with the camera, the pixel footprint (a uniform
position in the pixel) and the sphere of the scene, and prints the mean over the sphere pixels with its standard error.

  python tools/test/medium_sphere_reference.py [--scene densegrayfurnace] [--paths 64] [--width 320 --height 240]
                                              [--renders <the scene's work folder of mitsuba_compare.py>]

With --renders, the sphere means of the two renders of each renderer in that folder (gold_<seed>.pfm, mitsuba_<seed>.npy) are compared
with the reference.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mitsuba_compare as mc                                  # noqa: E402  (scene descriptions only; Mitsuba is not imported)


def camera_rays(cam, w, h, px, py, rng):
    o, t, up = (np.array(cam[k], float) for k in ("origin", "target", "up"))
    fwd = (t - o) / np.linalg.norm(t - o)
    right = np.cross(fwd, up); right /= np.linalg.norm(right)
    upv = np.cross(right, fwd)
    th = np.tan(np.radians(cam["fov"]) / 2)
    x = px + rng.random(px.shape)
    y = py + rng.random(py.shape)
    u = (2 * x / w - 1) * th * w / h
    v = (1 - 2 * y / h) * th
    d = fwd + u[:, None] * right + v[:, None] * upv
    return o, d / np.linalg.norm(d, axis=1, keepdims=True)


def exit_distance(p, d, c, r):
    """Distance along d from p (inside the sphere) to the sphere."""
    oc = p - c
    b = np.sum(oc * d, axis=1)
    disc = b * b - (np.sum(oc * oc, axis=1) - r * r)
    return -b + np.sqrt(np.maximum(disc, 0.0))


def hg_turn(d, g, rng):
    n = len(d)
    u1, u2 = rng.random(n), rng.random(n)
    if abs(g) < 1e-6:
        cos = 1 - 2 * u1
    else:
        s = (1 - g * g) / (1 - g + 2 * g * u1)
        cos = (1 + g * g - s * s) / (2 * g)
    sin = np.sqrt(np.maximum(0, 1 - cos * cos))
    phi = 2 * np.pi * u2
    a = np.where(np.abs(d[:, 0:1]) > 0.9, np.array([[0.0, 1.0, 0.0]]), np.array([[1.0, 0.0, 0.0]]))
    t1 = np.cross(d, a); t1 /= np.linalg.norm(t1, axis=1, keepdims=True)
    t2 = np.cross(d, t1)
    return sin[:, None] * (np.cos(phi)[:, None] * t1 + np.sin(phi)[:, None] * t2) + cos[:, None] * d


def walk(p, d, c, r, sigma_t, albedo, g, rng):
    """Escape weight of each path that starts at p (on the sphere) going into it along d."""
    n = len(p)
    w = np.ones(n)
    out = np.zeros(n)
    alive = np.arange(n)
    p, d = p.copy(), d.copy()
    while len(alive):
        pa, da = p[alive], d[alive]
        t = -np.log1p(-rng.random(len(alive))) / sigma_t
        texit = exit_distance(pa, da, c, r)
        esc = t >= texit
        out[alive[esc]] = w[alive[esc]]
        keep = ~esc
        idx = alive[keep]
        p[idx] = pa[keep] + da[keep] * t[keep, None]
        w[idx] *= albedo
        d[idx] = hg_turn(da[keep], g, rng)
        # Russian roulette on small weights (unbiased)
        low = w[idx] < 0.05
        kill = low & (rng.random(len(idx)) < 0.5)
        w[idx[low & ~kill]] *= 2.0
        alive = idx[~kill]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scene", default="densegrayfurnace")
    ap.add_argument("--paths", type=int, default=64, help="paths per sphere pixel")
    ap.add_argument("--width", type=int, default=320)
    ap.add_argument("--height", type=int, default=240)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--renders", default="")
    a = ap.parse_args()
    sc = next(s for s in mc.SCENES if s["name"] == a.scene)
    m = sc["materials"]["obj"]
    assert m["type"] == "dielectric" and m["ior"] == 1.0 and "medium" in m and sc["env"] > 0 and not sc["lights"] and len(sc["objects"]) == 1
    md = m["medium"]
    sigma_t, albedo, g = md["sigmaT"][0], md["albedo"][0], md["g"]
    assert len(set(md["sigmaT"])) == 1 and len(set(md["albedo"])) == 1, "gray media only"
    sp = sc["objects"][0]
    c, r = np.array(sp["center"], float), float(sp["radius"])
    mask = mc.sphere_mask(sc, a.width, a.height)
    py, px = np.nonzero(mask)
    rng = np.random.default_rng(a.seed)
    sums = np.zeros(len(px)); sq = np.zeros(len(px))
    for k in range(a.paths):
        o, d = camera_rays(sc["camera"], a.width, a.height, px.astype(float), py.astype(float), rng)
        oc = o - c
        b = d @ oc
        disc = b * b - (oc @ oc - r * r)
        hit = disc > 0
        val = np.full(len(px), sc["env"])            # a jittered ray that misses the sphere sees the environment
        th = -b[hit] - np.sqrt(disc[hit])
        val[hit] = sc["env"] * walk(o + d[hit] * th[:, None], d[hit], c, r, sigma_t, albedo, g, rng)
        sums += val; sq += val * val
        if (k + 1) % 16 == 0:
            print(f"  {k + 1} paths per pixel: sphere mean {sums.sum() / ((k + 1) * len(px)):.5f}", flush=True)
    n = a.paths * len(px)
    mean = sums.sum() / n
    se = np.sqrt((sq.sum() / n - mean * mean) / n)
    print(f"{a.scene}: reference sphere mean {mean:.5f} +- {se:.5f} ({len(px)} pixels, {a.paths} paths each; sigma_t {sigma_t}, albedo {albedo}, g {g})")
    if a.renders:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "tools", "lib"))
        from hdr_io import read_pfm
        for label, load in (("GoLD_rayt", lambda k: read_pfm(os.path.join(a.renders, f"gold_{k}.pfm"))),
                            ("Mitsuba", lambda k: np.load(os.path.join(a.renders, f"mitsuba_{k}.npy")))):
            x, y = (np.asarray(load(k), dtype=np.float64)[mask][:, 1] for k in (1, 2))
            m = (x.mean() + y.mean()) / 2
            e = np.sqrt(((x - y) ** 2 / 2).mean() / (2 * len(x)))
            print(f"  {label:9s} sphere mean {m:.5f} +- {e:.5f}: ratio to the reference {m / mean:.4f}, z {(m - mean) / np.hypot(e, se):+.1f}")


if __name__ == "__main__":
    main()
