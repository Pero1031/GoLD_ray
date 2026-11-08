"""
Clear-sky environment map (latitude-longitude, Radiance .hdr) from single scattering in the atmosphere, plus the sun for the scene file.

Model (Nishita et al. 1993 / Bruneton & Neyret 2008 style, the simplest version that gives the colors of a clear sky):
  * Rayleigh scattering of air:  beta_R(lambda) = 13.5e-6 (550 / lambda)^4  per metre at sea level, density exp(-h / 8 km),
    phase 3/(16 pi) (1 + cos^2).
  * Mie scattering of aerosols:  beta_M = 4.4e-6 per metre (extinction; single-scattering albedo 0.9), density exp(-h / 1.2 km),
    Henyey-Greenstein phase with g = 0.76.
  * Ozone is left out (it only matters for the blues of twilight). Only light scattered once reaches the eye, so the sky near the
    horizon is a little too dark and too blue compared with the real one.
  * The observer stands on the ground (radius 6360 km, atmosphere to 6460 km). The sun is a point at infinity; the sky is computed
    for the 16 wavelengths 400 ... 700 nm, then turned into linear sRGB with the same CMF fit as the renderer.
  * The sun disc is NOT in the map (it would be a 46000 W/m^2/sr spike in two pixels): the scene file gets a small sphere light.
    This script prints its color and the numbers for the JSON.
  * Units: the unattenuated sun gives a white surface facing it a radiance of 1 (illuminance pi), the usual exposure of the scenes.

The map follows the engine's convention (y up, u = (atan2(z, x) + pi) / 2 pi, top row = zenith). Below the horizon: a gray ground
(albedo 0.25) lit by the sun and the sky.

  python tools/data/make_sky.py --sun-el 12 --sun-dir -0.15 -0.97 --out GoLD_rayt/assets/env/sky_el12.hdr
"""
import argparse
import os

import numpy as np

R_EARTH, R_ATMO = 6360e3, 6460e3
H_R, H_M = 8000.0, 1200.0
BETA_M_EXT, MIE_ALBEDO, MIE_G = 4.4e-6, 0.9, 0.76
LAMBDAS = np.linspace(400.0, 700.0, 16)
SUN_SOLID_ANGLE = 6.8e-5            # 0.53 degrees across

M_XYZ_TO_RGB = np.array([[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415], [0.0557, -0.2040, 1.0570]])


def cmf_wyman(l):
    def g(x, mu, s1, s2):
        s = np.where(x < mu, s1, s2)
        return np.exp(-0.5 * ((x - mu) / s) ** 2)
    x = 1.056 * g(l, 599.8, 37.9, 31.0) + 0.362 * g(l, 442.0, 16.0, 26.7) - 0.065 * g(l, 501.1, 20.4, 26.2)
    y = 0.821 * g(l, 568.8, 46.9, 40.5) + 0.286 * g(l, 530.9, 16.3, 31.1)
    z = 1.217 * g(l, 437.0, 11.8, 36.0) + 0.681 * g(l, 459.0, 26.0, 13.8)
    return np.stack([x, y, z], axis=-1)


def blackbody(l_nm, kelvin):
    c2 = 1.438777e7
    return 1e15 / (l_nm ** 5 * np.expm1(c2 / (l_nm * kelvin)))


def sun_spectrum():
    """Extraterrestrial solar irradiance (5778 K blackbody), scaled so that its luminance is pi (a white surface facing it: radiance 1)."""
    e = blackbody(LAMBDAS, 5778.0)
    y = (cmf_wyman(LAMBDAS)[:, 1] * e).sum() * (LAMBDAS[1] - LAMBDAS[0])
    return e * (np.pi / y)


def to_linear_rgb(spec):
    """spec: (..., 16) spectral radiance -> linear sRGB (..., 3)."""
    xyz = (spec[..., :, None] * cmf_wyman(LAMBDAS)).sum(axis=-2) * (LAMBDAS[1] - LAMBDAS[0])
    return xyz @ M_XYZ_TO_RGB.T


def ray_sphere_exit(o, d, radius):
    """Distance along d (unit) from o to the exit point of the sphere of `radius` around the origin (o inside)."""
    b = (o * d).sum(-1)
    c = (o * o).sum(-1) - radius ** 2
    return -b + np.sqrt(np.maximum(b * b - c, 0.0))


def hits_earth(o, d):
    b = (o * d).sum(-1)
    c = (o * o).sum(-1) - R_EARTH ** 2
    disc = b * b - c
    return (disc > 0) & (-b - np.sqrt(np.maximum(disc, 0.0)) > 0)


def beta(h, lam):
    """Scattering coefficients at height h (m): (Rayleigh (..., L), Mie scattering (..., 1), Mie extinction (..., 1))."""
    r = 13.5e-6 * (550.0 / lam) ** 4 * np.exp(-h / H_R)[..., None]
    m_ext = (BETA_M_EXT * np.exp(-h / H_M))[..., None]
    return r, m_ext * MIE_ALBEDO, m_ext


def optical_depth(o, d, length, n, lam):
    """Integral of the extinction along the segment: (..., L)."""
    t = np.zeros(o.shape[:-1] + (len(lam),))
    ds = length / n
    for i in range(n):
        p = o + d * ((i + 0.5) * ds)[..., None]
        h = np.linalg.norm(p, axis=-1) - R_EARTH
        r, _, me = beta(h, lam)
        t += (r + me) * ds[..., None]
    return t


def phase_hg(c, g):
    return (1 - g * g) / (4 * np.pi * (1 + g * g - 2 * g * c) ** 1.5)


def sky_radiance(dirs, sun_dir, e_sun, n_view=24, n_sun=8):
    """Single-scattering radiance (..., 16) of the sky in directions `dirs` (unit, y up), observer 1 m above the ground."""
    o = np.zeros(dirs.shape[:-1] + (3,))
    o[..., 1] = R_EARTH + 1.0
    length = ray_sphere_exit(o, dirs, R_ATMO)
    ds = length / n_view
    cos_t = (dirs * sun_dir).sum(-1)
    p_r = 3.0 / (16.0 * np.pi) * (1.0 + cos_t ** 2)
    p_m = phase_hg(cos_t, MIE_G)
    out = np.zeros(dirs.shape[:-1] + (len(LAMBDAS),))
    od_view = np.zeros_like(out)
    for i in range(n_view):
        p = o + dirs * ((i + 0.5) * ds)[..., None]
        h = np.linalg.norm(p, axis=-1) - R_EARTH
        r, ms, me = beta(h, LAMBDAS)
        od_view += (r + me) * ds[..., None]
        sd = np.broadcast_to(sun_dir, p.shape)
        l_sun = ray_sphere_exit(p, sd, R_ATMO)
        od_sun = optical_depth(p, sd, l_sun, n_sun, LAMBDAS)
        lit = ~hits_earth(p, sd)
        t = np.exp(-(od_view + od_sun))
        scat = (r * p_r[..., None] + ms * p_m[..., None]) * ds[..., None]
        out += np.where(lit[..., None], t * scat, 0.0)
    return out * e_sun


def sun_transmittance(sun_dir):
    o = np.array([0.0, R_EARTH + 1.0, 0.0])
    l = ray_sphere_exit(o, sun_dir, R_ATMO)
    return np.exp(-optical_depth(o, sun_dir, np.array(l), 64, LAMBDAS))


def write_hdr(path, img):
    """Radiance .hdr with new-style RLE scanlines. img: (H, W, 3) float, non-negative."""
    img = np.maximum(img, 0.0)
    m = img.max(axis=-1)
    mant, exp = np.frexp(m)
    scale = np.where(m > 1e-32, mant * 256.0 / np.maximum(m, 1e-32), 0.0)
    rgbe = np.zeros(img.shape[:2] + (4,), dtype=np.uint8)
    rgbe[..., :3] = np.clip(img * scale[..., None], 0, 255).astype(np.uint8)
    rgbe[..., 3] = np.where(m > 1e-32, exp + 128, 0).astype(np.uint8)
    h, w = img.shape[:2]
    with open(path, "wb") as f:
        f.write(b"#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n" + f"-Y {h} +X {w}\n".encode())
        for y in range(h):
            f.write(bytes([2, 2, w >> 8, w & 255]))
            for c in range(4):
                row = rgbe[y, :, c].tobytes()
                i = 0
                while i < w:
                    run = 1
                    while i + run < w and run < 127 and row[i + run] == row[i]:
                        run += 1
                    if run >= 4:
                        f.write(bytes([128 + run, row[i]]))
                        i += run
                    else:
                        n = 1
                        while i + n < w and n < 128 and not (i + n + 3 < w and row[i + n] == row[i + n + 1] == row[i + n + 2] == row[i + n + 3]):
                            n += 1
                        f.write(bytes([n]) + row[i:i + n])
                        i += n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sun-el", type=float, default=12.0, help="sun elevation, degrees")
    ap.add_argument("--sun-dir", type=float, nargs=2, default=(-0.15, -0.97), help="horizontal direction of the sun (x z)")
    ap.add_argument("--width", type=int, default=512)
    ap.add_argument("--ground", type=float, default=0.25, help="ground albedo")
    ap.add_argument("--out", default="sky.hdr")
    a = ap.parse_args()

    el = np.radians(a.sun_el)
    hz = np.array([a.sun_dir[0], a.sun_dir[1]]) / np.hypot(*a.sun_dir)
    sun_dir = np.array([hz[0] * np.cos(el), np.sin(el), hz[1] * np.cos(el)])
    e_sun = sun_spectrum()
    w, h = a.width, a.width // 2

    # directions of the pixel centers (engine convention)
    u = (np.arange(w) + 0.5) / w
    v = 1.0 - (np.arange(h) + 0.5) / h                       # v = 1 at the top row (zenith)
    phi = u * 2 * np.pi - np.pi
    theta = np.pi * (1.0 - v)
    d = np.stack([np.sin(theta)[:, None] * np.cos(phi)[None, :], np.cos(theta)[:, None] * np.ones(w)[None, :], np.sin(theta)[:, None] * np.sin(phi)[None, :]], axis=-1)

    img = np.zeros((h, w, 3))
    up = d[..., 1] >= 0.0
    rows = np.where(up.any(axis=1))[0]
    dirs = d[rows]
    chunk = 8
    for s in range(0, len(rows), chunk):
        sl = slice(s, s + chunk)
        rad = sky_radiance(np.maximum(dirs[sl], [-2, 1e-3, -2]) if False else np.where((dirs[sl][..., 1] < 1e-3)[..., None], np.stack([dirs[sl][..., 0], np.full_like(dirs[sl][..., 1], 1e-3), dirs[sl][..., 2]], -1), dirs[sl]), sun_dir, e_sun)
        img[rows[sl]] = to_linear_rgb(rad)
    # the ground: lit by the sun and by the sky (mean sky radiance of the upper hemisphere as a flat ambient)
    sky_mean = img[up].mean(axis=0)
    t_sun = to_linear_rgb(e_sun * sun_transmittance(sun_dir))
    ground = a.ground / np.pi * (t_sun * max(np.sin(el), 0.0) + np.pi * sky_mean)
    img[~up] = ground
    write_hdr(a.out, img)

    # the sun for the scene file
    spec = e_sun * sun_transmittance(sun_dir) / SUN_SOLID_ANGLE
    sun_rgb = to_linear_rgb(spec)
    lum = 0.2126 * sun_rgb[0] + 0.7152 * sun_rgb[1] + 0.0722 * sun_rgb[2]
    print(f"wrote {a.out} ({w}x{h}); sky luminance: zenith {0.2126 * img[0, 0, 0] + 0.7152 * img[0, 0, 1] + 0.0722 * img[0, 0, 2]:.3f}, "
          f"horizon row {0.2126 * img[h // 2 - 1, 0, 0] + 0.7152 * img[h // 2 - 1, 0, 1] + 0.0722 * img[h // 2 - 1, 0, 2]:.3f}")
    print(f"sun direction (x, y, z) = ({sun_dir[0]:.4f}, {sun_dir[1]:.4f}, {sun_dir[2]:.4f})")
    print(f"sun disc: color (linear sRGB, luminance 1) = [{sun_rgb[0] / lum:.3f}, {sun_rgb[1] / lum:.3f}, {sun_rgb[2] / lum:.3f}], radiance = {lum:.0f}")
    print(f'  scene file: {{ "type": "emission", "spectrum": {{ "rgb": [{sun_rgb[0] / lum:.3f}, {sun_rgb[1] / lum:.3f}, {sun_rgb[2] / lum:.3f}] }}, "radiance": {lum:.0f} }}')
    print(f"  sphere of radius {20 * np.tan(np.radians(0.265)):.4f} at distance 20: center = [{20 * sun_dir[0]:.3f}, {20 * sun_dir[1]:.3f}, {20 * sun_dir[2]:.3f}]")


if __name__ == "__main__":
    main()
