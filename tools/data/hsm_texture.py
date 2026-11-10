"""Turns a Topcon SR-5000/5100 ".hsm" measurement (a spectral image) into a spectral reflectance texture
(".sptex", read by Core/SpectralTexture.hpp), dividing by a white-plate measurement.

An .hsm file is a tar archive. Its "spectral.image" holds a small header (int32: id, width, height, 4,
bands) in a first block, then one block per wavelength (380..780 nm, 1 nm) holding the width x height
float32 radiance image row by row; blocks are 32768 floats or a multiple of it. This layout was checked
against the instrument's own spot CSV export (agreement within 2 %).

The texture stores a mean spectrum plus K principal-component spectra (default 6) and K coefficients
per texel, fitted over a wavelength range where the lamp has power (default 420..700 nm; outside it the
ratio to the white is noise). The engine holds the end values beyond that range.

Usage:
  python tools/data/hsm_texture.py sample.hsm --white white.hsm --out cloth.sptex [--k 6] [--range 420 700]
         [--white-reflectance 0.98] [--crop X Y W H] [--light-out lamp.csv] [--preview cloth.png]
  python tools/data/hsm_texture.py sample.hsm --info
"""
import argparse
import io
import os
import re
import struct
import sys
import tarfile

import numpy as np

BLOCK = 32768


def read_hsm(path):
    """Returns (cube (401, H, W) float32 radiance, info dict)."""
    with tarfile.open(path) as tar:
        raw = tar.extractfile("spectral.image").read()
        meta = tar.extractfile("measurement_data.xml").read().decode("utf-8", "ignore")
        recipe = tar.extractfile("recipe.xml").read().decode("utf-8", "ignore") if "recipe.xml" in tar.getnames() else ""
    hdr = np.frombuffer(raw[:20], dtype=np.int32)
    w, h, bands = int(hdr[1]), int(hdr[2]), int(hdr[4])
    stride = ((w * h + BLOCK - 1) // BLOCK) * BLOCK
    a = np.frombuffer(raw, dtype=np.float32)
    if a.size < stride * (bands + 1):
        raise ValueError(f"{path}: spectral.image too short for {w}x{h}x{bands}")
    cube = np.stack([a[stride * (k + 1):stride * (k + 1) + w * h].reshape(h, w) for k in range(bands)])
    info = {"width": w, "height": h, "bands": bands, "file": os.path.basename(path)}
    for key, pat in (("datetime", r'FileInfo datetime="([^"]+)"'), ("device", r'device_type="([^"]+)"'),
                     ("distance_mm", r'MeasureDistance mm="([^"]+)"'), ("integration", r'SpectralFilter[^>]*integration_time="([^"]+)"'),
                     ("saturated", r'is_saturated="([^"]+)"')):
        m = re.search(pat, meta) or re.search(pat, recipe)
        info[key] = m.group(1) if m else "?"
    return cube, info


def white_spectrum(white_cube, threshold=0.6):
    """Mean radiance spectrum over the bright plate (pixels above `threshold` of the maximum Y-like sum)."""
    y = white_cube[100:300].sum(axis=0)            # 480..680 nm, a luminance stand-in
    disc = y > threshold * y.max()
    spd = white_cube[:, disc].mean(axis=1)
    return spd, disc, (y[disc].std() / y[disc].mean())


def fit_basis(refl, k):
    """refl: (nL, N) spectra of the texels -> mean (nL,), basis (k, nL), coef (N, k)."""
    S = refl.T.astype(np.float64)
    mean = S.mean(axis=0)
    U, s, Vt = np.linalg.svd(S - mean, full_matrices=False)
    basis = Vt[:k]
    coef = U[:, :k] * s[:k]
    rec = mean + coef @ basis
    err = np.abs(rec - S)
    var = (s[:k] ** 2).sum() / (s ** 2).sum()
    return mean, basis, coef, {"variance": var, "max_err": err.max(), "p99_err": np.percentile(err, 99), "mean_err": err.mean()}


def write_sptex(path, w, h, lambda0, step, mean, basis, coef):
    k, nl = basis.shape
    with open(path, "wb") as f:
        f.write(b"SPTX")
        f.write(struct.pack("<5I", 1, w, h, k, nl))
        f.write(struct.pack("<2f", float(lambda0), float(step)))
        f.write(mean.astype("<f4").tobytes())
        f.write(basis.astype("<f4").tobytes())
        f.write(coef.reshape(h, w, k).astype("<f4").tobytes())


def write_curve(path, wavelengths, values, label, comments=()):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        for c in comments:
            f.write(f"# {c}\n")
        f.write(f"wavelength_nm,{label}\n")
        for wl, v in zip(wavelengths, values):
            f.write(f"{wl:g},{v:.6g}\n")


def preview_png(path, refl, lam):
    """sRGB under D65 of the reflectance cube (nL, H, W) on `lam` (nm), CIE 1931 2 deg via a compact fit."""
    try:
        from PIL import Image
    except ImportError:
        print("Pillow is not installed: no preview", file=sys.stderr)
        return
    # Wyman, Sloan, Shirley 2013 multi-lobe fit of the CIE 1931 2-degree observer
    def g(x, mu, s1, s2):
        s = np.where(x < mu, s1, s2)
        return np.exp(-0.5 * ((x - mu) / s) ** 2)
    xbar = 1.056 * g(lam, 599.8, 37.9, 31.0) + 0.362 * g(lam, 442.0, 16.0, 26.7) - 0.065 * g(lam, 501.1, 20.4, 26.2)
    ybar = 0.821 * g(lam, 568.8, 46.9, 40.5) + 0.286 * g(lam, 530.9, 16.3, 31.1)
    zbar = 1.217 * g(lam, 437.0, 11.8, 36.0) + 0.681 * g(lam, 459.0, 26.0, 13.8)
    # D65 approximated by a 6504 K black body normalised at 560 nm (a preview only)
    l_m = lam * 1e-9
    d65 = 1.0 / (l_m ** 5 * (np.exp(1.438776e-2 / (l_m * 6504.0)) - 1.0))
    d65 /= d65[np.argmin(np.abs(lam - 560))]
    cmf = np.stack([xbar, ybar, zbar]) * d65
    xyz = np.einsum("cl,lhw->chw", cmf, refl) / (ybar * d65).sum()
    M = np.array([[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415], [0.0557, -0.2040, 1.0570]])
    rgb = np.clip(np.einsum("ij,jhw->ihw", M, xyz), 0, 1)
    rgb = np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * rgb ** (1 / 2.4) - 0.055)
    img = (rgb.transpose(1, 2, 0) * 255 + 0.5).astype(np.uint8)
    Image.fromarray(img).resize((img.shape[1] * 3, img.shape[0] * 3), Image.NEAREST).save(path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sample", help="the sample's .hsm")
    ap.add_argument("--white", help="the white plate's .hsm (same lamp and distance)")
    ap.add_argument("--white-reflectance", type=float, default=1.0, help="reflectance of the white plate (default 1)")
    ap.add_argument("--out", help="output .sptex")
    ap.add_argument("--k", type=int, default=6, help="number of basis spectra (default 6)")
    ap.add_argument("--range", nargs=2, type=float, default=(420.0, 700.0), metavar=("LO_NM", "HI_NM"))
    ap.add_argument("--crop", nargs=4, type=int, metavar=("X", "Y", "W", "H"), help="texel rectangle to keep (default: all)")
    ap.add_argument("--light-out", help="also write the lamp's relative spectrum (white radiance / white reflectance) as an engine CSV")
    ap.add_argument("--preview", help="write an sRGB (D65) preview PNG of the reflectance")
    ap.add_argument("--info", action="store_true", help="print the file's header and statistics only")
    a = ap.parse_args(argv)

    cube, info = read_hsm(a.sample)
    lam = np.arange(380, 380 + info["bands"])
    print(f"{info['file']}: {info['width']} x {info['height']} texels, {info['bands']} bands, {info['device']}, {info['datetime']}, "
          f"distance {info['distance_mm']} mm, integration {info['integration']}, saturated {info['saturated']}")
    if a.info or not a.out:
        print(f"  radiance at 550 nm: median {np.median(cube[170]):.4g}, max {cube[170].max():.4g}")
        return
    if not a.white:
        ap.error("--white is required to make a reflectance texture")
    white_cube, winfo = read_hsm(a.white)
    spd, disc, nonuni = white_spectrum(white_cube)
    print(f"white {winfo['file']}: {winfo['datetime']}, distance {winfo['distance_mm']} mm; plate = {int(disc.sum())} px, "
          f"non-uniformity {nonuni * 100:.1f} % (std/mean over the plate); lamp peak {lam[spd.argmax()]} nm")
    refl = cube / np.maximum(spd[:, None, None], 1e-12) * a.white_reflectance
    if a.crop:
        x, y, w, h = a.crop
        refl = refl[:, y:y + h, x:x + w]
    nl, h, w = refl.shape
    sel = (lam >= a.range[0]) & (lam <= a.range[1])
    lam_sel = lam[sel]
    mean, basis, coef, stats = fit_basis(refl[sel].reshape(sel.sum(), -1), a.k)
    print(f"reflectance over {a.range[0]:g}..{a.range[1]:g} nm: median {np.median(refl[sel]):.3f}, p99 {np.percentile(refl[sel], 99):.3f}, max {refl[sel].max():.3f}; "
          f"{int((refl[sel] > 1).sum())} values above 1")
    print(f"basis k={a.k}: variance {stats['variance']:.4f}, error max {stats['max_err']:.4f}, p99 {stats['p99_err']:.4f}, mean {stats['mean_err']:.4f}")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    write_sptex(a.out, w, h, lam_sel[0], 1.0, mean, basis, coef)
    print(f"wrote {a.out} ({w} x {h} x {a.k}, {os.path.getsize(a.out) / 1024:.0f} KB)")
    if a.light_out:
        write_curve(a.light_out, lam, spd / spd.max(), "relative_power",
                    [f"lamp spectrum from the white plate {winfo['file']} ({winfo['datetime']}), SR-5100 ROI {winfo['width']} x {winfo['height']}",
                     f"relative spectral power (peak = 1); the plate's reflectance ({a.white_reflectance:g}) is assumed flat"])
        print(f"wrote {a.light_out}")
    if a.preview:
        preview_png(a.preview, refl[sel], lam_sel)
        print(f"wrote {a.preview}")


if __name__ == "__main__":
    main()
