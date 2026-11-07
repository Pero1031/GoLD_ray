"""Checks of the output gamuts (--gamut srgb|p3|rec2020).

  python tools/test/gamut_check.py matrices
      recomputes the sRGB -> Display P3 / Rec.2020 matrices from the primaries and the D65 white, and compares them
      with the numbers in GoLD_rayt/include/Color/Gamut.hpp (exact to the 7 digits written there).

  python tools/test/gamut_check.py png FILE.png [primaries]
      reads the PNG chunks: is there a cICP chunk, with which code points, and is the CRC right?
      (primaries: 1 = sRGB, 12 = Display P3, 9 = Rec.2020; with it the check fails when the file differs.)

  python tools/test/gamut_check.py pfm A.pfm B.pfm p3|rec2020
      A was written with --gamut srgb and B with --gamut p3 (rec2020), two renders with the same settings. The renderer's random
      numbers are not seeded (and threads are scheduled dynamically), so the pixels differ by noise; the check compares the mean
      color of the image: mean(B) must equal the matrix applied to mean(A), within a few percent (noise of the two renders).
"""
import os, re, struct, sys, zlib

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools", "lib"))

D65 = (0.3127, 0.3290)
PRIMARIES = {
    "srgb":    [(0.640, 0.330), (0.300, 0.600), (0.150, 0.060)],
    "p3":      [(0.680, 0.320), (0.265, 0.690), (0.150, 0.060)],
    "rec2020": [(0.708, 0.292), (0.170, 0.797), (0.131, 0.046)],
}


def rgb_to_xyz(prims, white=D65):
    """Linear RGB -> XYZ for primaries given as xy chromaticities and a white xy (Y = 1)."""
    p = np.array([[x / y, 1.0, (1 - x - y) / y] for x, y in prims]).T   # columns: XYZ of the primaries at Y = 1
    w = np.array([white[0] / white[1], 1.0, (1 - white[0] - white[1]) / white[1]])
    s = np.linalg.solve(p, w)
    return p * s


def from_srgb(name):
    """3x3 matrix: linear sRGB -> linear RGB of `name`."""
    return np.linalg.inv(rgb_to_xyz(PRIMARIES[name])) @ rgb_to_xyz(PRIMARIES["srgb"])


def header_matrices():
    text = open(os.path.join(ROOT, "GoLD_rayt", "include", "Color", "Gamut.hpp"), encoding="utf-8").read()
    out = {}
    for key, case in (("p3", "case Gamut::DisplayP3:"), ("rec2020", "case Gamut::Rec2020:")):
        start = re.search(re.escape(case) + r"\s*return Vector3\(", text).start()   # the matrix, not the name / code-point switches
        block = text[start:]
        block = block[:block.index("case Gamut::", 10)]
        m = np.zeros((3, 3))
        # rows are written one per line; each term is  Real(v) * c.x|y|z
        terms = re.findall(r"Real\((-?[0-9.]+)\)\s*\*\s*c\.([xyz])", block)
        nums = terms
        row = -1
        for line in block.splitlines():
            t = re.findall(r"Real\((-?[0-9.]+)\)\s*\*\s*c\.([xyz])", line)
            if not t:
                continue
            row += 1
            for v, ax in t:
                m[row, "xyz".index(ax)] = float(v)
        assert row == 2 and len(terms) == len(nums), key
        out[key] = m
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
    elif len(sys.argv) >= 5 and sys.argv[1] == "pfm":
        ok = check_pfm(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print(__doc__)
        ok = False
    sys.exit(0 if ok else 1)
