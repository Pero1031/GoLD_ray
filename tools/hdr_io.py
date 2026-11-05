"""Minimal Radiance .hdr (RGBE) reader, enough for the files written by stb_image_write (flat or new-style RLE scanlines).

    img = read_hdr(path)    # float32 array (H, W, 3), linear values
"""
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
