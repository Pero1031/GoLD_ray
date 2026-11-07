"""Reference port of the engine's RGB -> spectrum conversion (include/Color/SpectralConversion.hpp, detail::RgbBasis).

Three smooth basis spectra P_r, P_g, P_b on 5 nm nodes (380..780 nm) are found once as the smoothest solution (second differences) of
RGB(P_c) = e_c under D65 / sRGB with P_r + P_g + P_b = 1; the spectrum of an RGB color is  max(0, r P_r + g P_g + b P_b).

    import rgb_to_spectrum as r2s
    rho = r2s.reflectance_spectrum([0.8, 0.3, 0.1], wl_nm)       # an albedo (white -> 1 everywhere)

Run it to check the port against the engine: the spectral build renders a diffuse sphere inside a glowing sphere (a convex sphere in a
uniform light: radiance = albedo(lambda) x light(lambda)), and the film color must equal the prediction made here.

    python tools/rgb_to_spectrum.py check <spectral GoLD_rayt.exe>

It is also the tool for asking how much a real material loses when it is described by an RGB color (the test of README "Not done").
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import metal_accuracy as ma  # noqa: E402  (the film sensor, light spectra and D65 of the engine)

_NODES = np.arange(380.0, 781.0, 5.0)      # 81 nodes
_P = None


def basis():
    """P[channel][node] (3 x 81)."""
    global _P
    if _P is not None:
        return _P
    n = len(_NODES)
    lam = np.arange(380, 781)
    trap = np.where((lam == 380) | (lam == 780), 0.5, 1.0)
    sens = ma.SENSOR[(ma.LAM >= 380) & (ma.LAM <= 780)]                 # (401, 3): s_c(lambda), D65 -> (1,1,1)
    w = sens * (ma.d65(lam.astype(float)) * trap)[:, None]            # (401, 3)
    x = (lam - 380.0) / 5.0
    j = np.minimum(x.astype(int), n - 2)
    f = x - j
    A = np.zeros((3, n))
    for ch in range(3):
        np.add.at(A[ch], j, w[:, ch] * (1 - f))
        np.add.at(A[ch], j + 1, w[:, ch] * f)
    L = np.zeros((n - 2, n))
    for i in range(n - 2):
        L[i, i:i + 3] = (1, -2, 1)
    Q = L.T @ L + 1e-6 * np.eye(n)
    X = np.linalg.solve(Q, A.T)                                       # (n, 3)
    S = A @ X                                                         # (3, 3)
    Sinv = np.linalg.inv(S)
    P = np.zeros((3, n))
    for ch in range(3):
        b = np.full(3, -1.0 / 3.0)
        b[ch] += 1.0
        P[ch] = 1.0 / 3.0 + X @ (Sinv @ b)
    _P = P
    return P


def reflectance_spectrum(rgb, wl_nm):
    """The spectrum the engine uses for an RGB albedo, at the given wavelengths (nm)."""
    P = basis()
    s = sum(float(rgb[c]) * np.interp(np.clip(wl_nm, 380, 780), _NODES, P[c]) for c in range(3))
    return np.maximum(s, 0.0)


def _check(exe):
    """Engine vs this port: a diffuse sphere of albedo `rgb` in a tungsten glow."""
    cases = [(0.8, 0.8, 0.8), (0.8, 0.3, 0.1), (0.2, 0.5, 0.9), (0.1, 0.6, 0.2), (0.35, 0.1, 0.12), (0.05, 0.04, 0.04)]
    spd, rgb_light = ma.make_light(ma.blackbody(ma.LAM, 2856.0))
    lam = ma.LAM
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join([os.path.join(root, "librarys", "embree-4.4.0.x64.windows", "bin"),
                                   os.path.join(root, "librarys", "oidn-2.5.1.x64.windows", "bin"), env["PATH"]])
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from hdr_io import read_pfm
    print("albedo RGB            engine (spectral build)        this port                 max diff")
    worst = 0.0
    for c in cases:
        scene = {"options": {"width": 64, "height": 64, "spp": 256},
                 "camera": {"from": [0, 0, 3], "at": [0, 0, 0], "fov": 14},
                 "materials": {"m": {"type": "diffuse", "albedo": list(c)},
                               "glow": {"type": "emission", "spectrum": "tungsten", "radiance": 1, "twoSided": True}},
                 "objects": [{"shape": "sphere", "center": [0, 0, 0], "radius": 0.5, "material": "m"},
                             {"shape": "sphere", "center": [0, 0, 0], "radius": 30, "material": "glow", "light": False}]}
        import json
        with tempfile.TemporaryDirectory() as d:
            sj, out = os.path.join(d, "s.json"), os.path.join(d, "o.pfm")
            json.dump(scene, open(sj, "w"))
            subprocess.run([exe, sj, "--out", out], cwd=os.path.join(root, "GoLD_rayt"), env=env, capture_output=True, check=True)
            im = read_pfm(out).astype(float)
        engine = im[28:36, 28:36].reshape(-1, 3).mean(0)
        mine = ma.sensor_rgb(reflectance_spectrum(c, lam) * spd)
        worst = max(worst, np.abs(engine - mine).max())
        print("%-20s  %.4f %.4f %.4f        %.4f %.4f %.4f        %.4f" % (c, *engine, *mine, np.abs(engine - mine).max()))
    print("worst difference: %.4f" % worst)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "check":
        _check(sys.argv[2])
    else:
        print(__doc__)
