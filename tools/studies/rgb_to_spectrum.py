"""Reference port of the engine's RGB -> spectrum conversion (include/Color/SpectralConversion.hpp, detail::RgbBasis).

Three smooth basis spectra P_r, P_g, P_b on 5 nm nodes (380..780 nm) are found once as the smoothest solution (second differences) of
RGB(P_c) = e_c under D65 / sRGB with P_r + P_g + P_b = 1; the spectrum of an RGB color is  max(0, r P_r + g P_g + b P_b).

    import rgb_to_spectrum as r2s
    rho = r2s.reflectance_spectrum([0.8, 0.3, 0.1], wl_nm)       # an albedo (white -> 1 everywhere)

Run it to check the port against the engine: the spectral build renders a diffuse sphere inside a glowing sphere (a convex sphere in a
uniform light: radiance = albedo(lambda) x light(lambda)), and the film color must equal the prediction made here.

    python tools/studies/rgb_to_spectrum.py check <spectral GoLD_rayt.exe>

It is also the tool for asking how much a real material loses when it is described by an RGB color (the test of README "Not done").
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
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


# ---------------------------------------------------------------------------------------------------------------------------------
# Jakob & Hanika (2019), "A Low-Dimensional Function Space for Efficient Spectral Upsampling": a reflectance
#     s(lambda) = S(c0 x^2 + c1 x + c2),   S(t) = 1/2 + t / (2 sqrt(1 + t^2)),   x = (lambda - 360) / 470
# whose color under D65 (the engine's film sensor) is the given RGB. The coefficients are found here by a damped Gauss-Newton fit in
# CIELAB, walking from gray to the color; the paper tabulates them instead. Reflectances are in [0, 1]: an RGB channel above 1 is not
# representable (it is clamped to 0.999).
# ---------------------------------------------------------------------------------------------------------------------------------
_X = (ma.LAM - 360.0) / 470.0
_W3 = ma.M_RGB_TO_XYZ @ np.ones(3)
_D65 = None


def _sigmoid(t):
    return 0.5 + t / (2.0 * np.sqrt(1.0 + t * t))


def _jh_spectrum(c):
    return _sigmoid(c[0] * _X ** 2 + c[1] * _X + c[2])


def _lab(rgb):
    t = (ma.M_RGB_TO_XYZ @ np.asarray(rgb, dtype=float)) / _W3
    f = np.where(t > 216 / 24389, np.cbrt(t), (24389 / 27 * t + 16) / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def jakob_hanika_coefficients(rgb):
    """(c0, c1, c2) of the sigmoid spectrum whose color under D65 is `rgb`; returns (c, residual in CIELAB units)."""
    global _D65
    if _D65 is None:
        _D65 = ma.lights()["D65"][0]
    rgb = np.clip(np.asarray(rgb, dtype=float), 1e-3, 0.999)

    def model(c):
        return _lab(ma.sensor_rgb(_jh_spectrum(c) * _D65))

    m = float(np.clip(rgb.mean(), 1e-3, 0.999))
    c = np.array([0.0, 0.0, (2 * m - 1) / np.sqrt(1 - (2 * m - 1) ** 2)])      # flat spectrum of the mean value
    gray = np.full(3, m)
    for step in range(1, 13):                                                    # continuation: gray -> the color
        target = _lab(gray + (rgb - gray) * step / 12.0)
        lam = 1e-3
        for _ in range(60):
            r = model(c) - target
            J = np.zeros((3, 3))
            for k in range(3):
                d = np.zeros(3)
                d[k] = 1e-4
                J[:, k] = (model(c + d) - model(c - d)) / 2e-4
            try:
                dc = np.linalg.solve(J.T @ J + lam * np.eye(3), -J.T @ r)
            except np.linalg.LinAlgError:
                break
            new = model(c + dc) - target
            if np.linalg.norm(new) < np.linalg.norm(r):
                c, lam = c + dc, max(lam * 0.3, 1e-9)
                if np.linalg.norm(new) < 1e-6:
                    break
            else:
                lam *= 10.0
    return c, float(np.linalg.norm(model(c) - _lab(rgb)))


def jakob_hanika_spectrum(rgb, wl_nm):
    """The Jakob-Hanika reflectance of an RGB albedo at the given wavelengths (nm)."""
    c, _ = jakob_hanika_coefficients(rgb)
    return _sigmoid(c[0] * ((np.asarray(wl_nm, dtype=float) - 360.0) / 470.0) ** 2 + c[1] * ((np.asarray(wl_nm, dtype=float) - 360.0) / 470.0) + c[2])


def _check(exe):
    """Engine vs this port: a diffuse sphere of albedo `rgb` in a tungsten glow."""
    cases = [(0.8, 0.8, 0.8), (0.8, 0.3, 0.1), (0.2, 0.5, 0.9), (0.1, 0.6, 0.2), (0.35, 0.1, 0.12), (0.05, 0.04, 0.04)]
    spd, rgb_light = ma.make_light(ma.blackbody(ma.LAM, 2856.0))
    lam = ma.LAM
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join([os.path.join(root, "librarys", "embree-4.4.0.x64.windows", "bin"),
                                   os.path.join(root, "librarys", "oidn-2.5.1.x64.windows", "bin"), env["PATH"]])
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
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
