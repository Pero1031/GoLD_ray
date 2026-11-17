"""Readers for the renderer's float images: Radiance .hdr (RGBE, the files stb_image_write writes: flat or new-style RLE scanlines),
.pfm, and the ENVI spectral image of --spectral-image.

    img = read_hdr(path)              # float32 array (H, W, 3), linear values
    cube, wl = read_envi("out/x")     # float32 (H, W, bands) and the band centres in nm (out/x_spectral.hdr + .img)
"""
import os
import re

import numpy as np


def read_pfm(path):
    """Portable float map, 3 channels (what `GoLD_rayt --out x.pfm` writes): exact 32-bit floats, negative values kept."""
    with open(path, "rb") as f:
        if f.readline().strip() != b"PF":
            raise ValueError("not a color PFM file: " + path)
        w, h = (int(v) for v in f.readline().split())
        scale = float(f.readline())
        a = np.fromfile(f, dtype="<f4" if scale < 0 else ">f4", count=w * h * 3)
    return a.reshape(h, w, 3)[::-1].astype(np.float32)   # rows are stored bottom to top


def read_hdr(path):
    if path.lower().endswith(".pfm"):
        return read_pfm(path)
    with open(path, "rb") as f:
        data = f.read()
    pos = 0

    def line():
        nonlocal pos
        end = data.index(b"\n", pos)
        s = data[pos:end]
        pos = end + 1
        return s

    if not line().startswith(b"#?"):
        raise ValueError("not a Radiance HDR file: " + path)
    while line().strip():           # header lines until the blank line
        pass
    res = line().split()            # b"-Y", H, b"+X", W
    h, w = int(res[1]), int(res[3])
    rgbe = np.zeros((h, w, 4), dtype=np.uint8)
    for y in range(h):
        if w >= 8 and w < 32768 and data[pos] == 2 and data[pos + 1] == 2 and (data[pos + 2] << 8 | data[pos + 3]) == w:
            pos += 4
            for c in range(4):
                x = 0
                while x < w:
                    n = data[pos]
                    pos += 1
                    if n > 128:
                        n -= 128
                        rgbe[y, x:x + n, c] = data[pos]
                        pos += 1
                    else:
                        rgbe[y, x:x + n, c] = np.frombuffer(data, dtype=np.uint8, count=n, offset=pos)
                        pos += n
                    x += n
        else:                       # flat scanline
            rgbe[y] = np.frombuffer(data, dtype=np.uint8, count=4 * w, offset=pos).reshape(w, 4)
            pos += 4 * w
    e = rgbe[..., 3].astype(np.int32)
    scale = np.where(e > 0, np.ldexp(1.0, e - (128 + 8)), 0.0)
    return (rgbe[..., :3].astype(np.float64) * scale[..., None]).astype(np.float32)


ENVI_DTYPES = {1: "u1", 2: "i2", 3: "i4", 4: "f4", 5: "f8", 12: "u2", 13: "u4", 14: "i8", 15: "u8"}


def envi_paths(path):
    """(header, data) paths of an ENVI image given its .hdr, its .img or the stem a render was written to (x for x_spectral.hdr)."""
    stem, ext = os.path.splitext(path)
    if ext.lower() in (".hdr", ".img"):
        return stem + ".hdr", stem + ".img"
    if ext.lower() in (".pfm", ".png", ".exr", ".bmp", ".jpg"):
        path = stem                   # the --out file of the render
    return path + "_spectral.hdr", path + "_spectral.img"


def read_envi_header(path):
    """The header of an ENVI image as a dict of strings; brace lists ({a, b, ...}, possibly over several lines) as lists of strings."""
    header, _ = envi_paths(path)
    with open(header, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()
    if not text.lstrip().startswith("ENVI"):
        raise ValueError("not an ENVI header: " + header)
    out = {}
    for m in re.finditer(r"^\s*([^=\n]+?)\s*=\s*(\{[^}]*\}|[^\n]*)", text, re.M):
        key, value = m.group(1).strip().lower(), m.group(2).strip()
        if value.startswith("{"):
            inner = value[1:-1].strip()
            out[key] = inner if key == "description" else [v.strip() for v in inner.split(",") if v.strip()]
        else:
            out[key] = value
    return out


def read_envi(path):
    """An ENVI image as (cube, wavelengths): float32 (H, W, bands) and the band centres (None if the header lists none).

    `path` is the .hdr, the .img or the stem of a render (`--out x.pfm --spectral-image ...` gives x_spectral.hdr / .img). Band sequential,
    band interleaved by line or by pixel, either byte order. The renderer's bands are the mean spectral radiance over [centre - fwhm/2,
    centre + fwhm/2] (docs/CAMERA.md)."""
    h = read_envi_header(path)
    _, data = envi_paths(path)
    w, ht, b = int(h["samples"]), int(h["lines"]), int(h["bands"])
    dtype = np.dtype(ENVI_DTYPES[int(h.get("data type", 4))]).newbyteorder(">" if int(h.get("byte order", 0)) else "<")
    a = np.fromfile(data, dtype=dtype, count=w * ht * b, offset=int(h.get("header offset", 0)))
    if a.size != w * ht * b:
        raise ValueError("%s holds %d values, the header says %d x %d x %d" % (data, a.size, w, ht, b))
    order = h.get("interleave", "bsq").lower()
    if order == "bsq":
        cube = a.reshape(b, ht, w).transpose(1, 2, 0)
    elif order == "bil":
        cube = a.reshape(ht, b, w).transpose(0, 2, 1)
    elif order == "bip":
        cube = a.reshape(ht, w, b)
    else:
        raise ValueError("unknown ENVI interleave: " + order)
    wl = np.array([float(v) for v in h["wavelength"]]) if "wavelength" in h else None
    if wl is not None and h.get("wavelength units", "nanometers").lower() in ("micrometers", "microns", "um"):
        wl = wl * 1000.0
    return np.ascontiguousarray(cube, dtype=np.float32), wl
