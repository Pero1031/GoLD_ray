"""Checks of the output gamuts (--gamut srgb|p3|rec2020).

  python tools/test/gamut_check.py matrices
      recomputes the sRGB -> Display P3 / Rec.2020 matrices from the primaries and the D65 white, and compares them
      with the numbers in GoLD_rayt/include/Color/Gamut.hpp (exact to the 7 digits written there).

  python tools/test/gamut_check.py png FILE.png [primaries]
      reads the PNG chunks: is there a cICP chunk, with which code points, and is the CRC right?
      (primaries: 1 = sRGB, 12 = Display P3, 9 = Rec.2020; with it the check fails when the file differs.)

  python tools/test/gamut_check.py png16 FILE16.png FILE8.png
      FILE16 was written with --bits 16, FILE8 with --bits 8 (an independent render of the same scene): IHDR says 16 bits and the CRCs are right,
      the pixel data decodes (own decoder; every row filter is "none"), a standard decoder (Pillow) opens it, and the mean color agrees with
      FILE8 (noise of two renders: a few percent at most). The 16-bit values must use more than 256 distinct levels per channel.

  python tools/test/gamut_check.py pfm A.pfm B.pfm p3|rec2020
      A was written with --gamut srgb and B with --gamut p3 (rec2020), two renders with the same settings. The renderer's random
      numbers are not seeded (and threads are scheduled dynamically), so the pixels differ by noise; the check compares the mean
      color of the image: mean(B) must equal the matrix applied to mean(A), within a few percent (noise of the two renders).
"""
import os, re, struct, sys, zlib

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools", "lib"))

from gamut import from_srgb   # noqa: E402  (the primaries and the matrices live in tools/lib/gamut.py)


def header_matrices():
    text = open(os.path.join(ROOT, "GoLD_rayt", "include", "Color", "Gamut.hpp"), encoding="utf-8").read()
    out = {}
    for key, var in (("p3", "p3"), ("rec2020", "rec2020")):
        start = re.search(r"static const Matrix3 " + var + r" = matrixFromRows\(", text).end()
        body = text[start:text.index(");", start)]
        nums = [float(v) for v in re.findall(r"Real\((-?[0-9.]+)\)", body)]
        assert len(nums) == 9, key
        out[key] = np.array(nums).reshape(3, 3)
    return out


def check_matrices():
    ok = True
    hm = header_matrices()
    for name in ("p3", "rec2020"):
        ref = from_srgb(name)
        err = np.abs(ref - hm[name]).max()
        white = hm[name] @ np.ones(3)
        print(f"{name:8s} max |reference - header| = {err:.2e}   header matrix on white = {np.round(white, 6)}")
        ok &= err < 2e-6 and np.allclose(white, 1.0, atol=2e-6)
    # the gold of the document: a candle-lit color outside sRGB is negative there and inside P3
    c = np.array([1.0, 0.45, -0.011])
    print("example (1, 0.45, -0.011) in sRGB ->", {n: np.round(hm[n] @ c, 4).tolist() for n in hm})
    print("PASS" if ok else "FAIL")
    return ok


def png_chunks(path):
    data = open(path, "rb").read()
    assert data[:8] == b"\x89PNG\r\n\x1a\n", "not a PNG"
    pos, chunks = 8, []
    while pos < len(data):
        n = struct.unpack(">I", data[pos:pos + 4])[0]
        typ = data[pos + 4:pos + 8]
        body = data[pos + 8:pos + 8 + n]
        crc = struct.unpack(">I", data[pos + 8 + n:pos + 12 + n])[0]
        chunks.append((typ.decode("latin1"), body, zlib.crc32(typ + body) & 0xFFFFFFFF == crc))
        pos += 12 + n
    return chunks


def check_png(path, want=None):
    chunks = png_chunks(path)
    print("chunks:", " ".join(t for t, _, _ in chunks[:6]), "...")
    bad = [t for t, _, good in chunks if not good]
    cicp = [(b, g) for t, b, g in chunks if t == "cICP"]
    order_ok = [t for t, _, _ in chunks].index("cICP") < [t for t, _, _ in chunks].index("IDAT") if cicp else True
    if cicp:
        prim, tc, mx, rng = cicp[0][0]
        print(f"cICP: primaries {prim}, transfer {tc}, matrix {mx}, full range {rng}, CRC {'ok' if cicp[0][1] else 'BAD'}, before IDAT: {order_ok}")
    else:
        print("no cICP chunk")
    ok = not bad and order_ok
    if want is not None:
        ok &= (cicp and cicp[0][0][0] == want) or (not cicp and want == 1)
    print("PASS" if ok else "FAIL")
    return ok


def decode_png16(path):
    """Own decoder for the files Film writes at 16 bits: RGB, no interlace, all rows filtered "none". Returns an (h, w, 3) uint16 array."""
    chunks = png_chunks(path)
    ihdr = [b for t, b, _ in chunks if t == "IHDR"][0]
    w, h, depth, ctype, _, _, interlace = struct.unpack(">IIBBBBB", ihdr)
    assert (depth, ctype, interlace) == (16, 2, 0), (depth, ctype, interlace)
    raw = zlib.decompress(b"".join(b for t, b, _ in chunks if t == "IDAT"))
    rows = np.frombuffer(raw, dtype=np.uint8).reshape(h, 1 + 6 * w)
    assert (rows[:, 0] == 0).all(), "a row uses a PNG filter"
    return rows[:, 1:].reshape(h, w, 3, 2).astype(np.uint16) @ np.array([256, 1], dtype=np.uint16) if False else         (rows[:, 1:].reshape(h, w, 3, 2)[..., 0].astype(np.uint16) << 8) | rows[:, 1:].reshape(h, w, 3, 2)[..., 1]


def check_png16(path16, path8):
    from PIL import Image
    ok = all(good for _, _, good in png_chunks(path16))
    a16 = decode_png16(path16)
    im = Image.open(path16)
    a8 = np.asarray(Image.open(path8).convert("RGB"))
    print(f"{path16}: {a16.shape[1]}x{a16.shape[0]}, Pillow opens it as {im.mode} {im.size}, CRCs {'ok' if ok else 'BAD'}")
    levels = [len(np.unique(a16[..., c])) for c in range(3)]
    m16, m8 = a16.reshape(-1, 3).mean(axis=0) / 65535.0, a8.reshape(-1, 3).mean(axis=0) / 255.0
    rel = np.abs(m16 - m8) / np.maximum(m8, 1e-3)
    print(f"distinct levels per channel: {levels} (8 bit would give at most 256)")
    print(f"mean color: 16 bit {np.round(m16, 4)}, 8 bit {np.round(m8, 4)}, relative difference {np.round(rel * 100, 2)} %")
    ok &= im.size == (a16.shape[1], a16.shape[0]) and max(levels) > 256 and bool((rel < 0.05).all())
    print("PASS" if ok else "FAIL")
    return ok


def read_pfm(path):
    from hdr_io import read_pfm as rp
    return rp(path)


def check_pfm(a_path, b_path, name, tol=0.03):
    a, b = read_pfm(a_path), read_pfm(b_path)
    m = from_srgb(name)
    ma, mb = a.reshape(-1, 3).mean(axis=0), b.reshape(-1, 3).mean(axis=0)
    pred = m @ ma
    rel = np.abs(pred - mb) / np.maximum(np.abs(pred), 1e-9)
    print(f"{a.shape[1]}x{a.shape[0]}  mean A {np.round(ma, 4)}  M*mean A {np.round(pred, 4)}  mean B {np.round(mb, 4)}  relative difference {np.round(rel * 100, 2)} %")
    print(f"negative values: A {(a < 0).mean() * 100:.2f} %, B {(b < 0).mean() * 100:.2f} %")
    ok = bool((rel < tol).all())
    print("PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "matrices":
        ok = check_matrices()
    elif len(sys.argv) >= 3 and sys.argv[1] == "png":
        ok = check_png(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else None)
    elif len(sys.argv) >= 4 and sys.argv[1] == "png16":
        ok = check_png16(sys.argv[2], sys.argv[3])
    elif len(sys.argv) >= 5 and sys.argv[1] == "pfm":
        ok = check_pfm(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print(__doc__)
        ok = False
    sys.exit(0 if ok else 1)
