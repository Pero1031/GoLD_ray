"""Checks of the sensor's raw file (--sensor, <out>.dng) with programs that are not the renderer: LibRaw (rawpy) reads and develops the
DNG, colour-demosaicing demosaics its mosaic. Needs `pip install rawpy colour-demosaicing`.

    out\\cli-dev-spectral\\GoLD_rayt.exe --scene spd --width 480 --spp 64 --seed 1 --denoise --sensor default ^
        --camera-response GoLD_rayt\\assets\\CameraResponse\\Canon_5DMarkII.csv --out x.pfm
    python tools/test/sensor_dng_check.py x_denoised          (reads x_denoised.dng and x_denoised.pfm)

1. The file as LibRaw reads it: the pattern R G / G B, the black and white levels, the colour matrix and the white balance as this
   script's own TIFF reader finds them, and every value of the mosaic between black and white.
2. The picture again from the DNG alone: the mosaic minus black, white-balanced with 1 / AsShotNeutral, demosaiced by colour-demosaicing's
   Malvar 2004, turned into XYZ with the inverse of ColorMatrix1 and into linear sRGB. It must be the renderer's picture (the .pfm, which
   holds radiance: compared up to one scale) inside 4 pixels from the edges, within the rounding of the matrix to 4 decimals.
3. LibRaw's own development (its AHD demosaicing and its colour pipeline, linear sRGB) against the renderer's picture: the mean colour of
   the picture and the colours of smooth regions (2 %). Expose so that little clips (the log's share at the top of the range): raw
   developers treat clipped highlights each in their own way.
"""
import os
import struct
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "lib"))
from hdr_io import read_pfm                                             # noqa: E402

XYZ_TO_SRGB = np.array([[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415], [0.0557, -0.2040, 1.0570]])   # the renderer's


def dng_tags(path):
    """The tags of the first directory of a little-endian TIFF: {tag: tuple of values} (rationals as floats)."""
    b = open(path, "rb").read()
    assert b[:2] == b"II" and struct.unpack_from("<H", b, 2)[0] == 42, "not a little-endian TIFF"
    ifd = struct.unpack_from("<I", b, 4)[0]
    sizes = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 10: 8}
    tags = {}
    for k in range(struct.unpack_from("<H", b, ifd)[0]):
        tag, typ, count = struct.unpack_from("<HHI", b, ifd + 2 + 12 * k)
        size = sizes[typ] * count
        at = ifd + 2 + 12 * k + 8 if size <= 4 else struct.unpack_from("<I", b, ifd + 2 + 12 * k + 8)[0]
        if typ in (1, 2):
            v = tuple(b[at:at + count])
        elif typ == 3:
            v = struct.unpack_from("<%dH" % count, b, at)
        elif typ == 4:
            v = struct.unpack_from("<%dI" % count, b, at)
        else:
            raw = struct.unpack_from("<%d%s" % (2 * count, "i" if typ == 10 else "I"), b, at)
            v = tuple(raw[2 * i] / raw[2 * i + 1] for i in range(count))
        tags[tag] = v
    return tags


def main(stem):
    import rawpy
    from colour_demosaicing import demosaicing_CFA_Bayer_Malvar2004

    tags = dng_tags(stem + ".dng")
    picture = read_pfm(stem + ".pfm").astype(np.float64)
    black, white = tags[50714][0], tags[50717][0]
    cm = np.array(tags[50721]).reshape(3, 3)
    neutral = np.array(tags[50728])
    failed = []

    # 1. LibRaw reads the same file. Its camera-to-sRGB matrix (rgb_cam, applied after its white balance) times the balance is
    #    XYZ_TO_SRGB inv(ColorMatrix1), each row scaled to sum 1 (LibRaw normalizes its matrix so that the balanced white stays white).
    with rawpy.imread(stem + ".dng") as raw:
        mosaic = raw.raw_image_visible.astype(np.float64).copy()
        pattern = raw.raw_pattern.tolist()
        libraw_black = list(raw.black_level_per_channel)
        libraw_white = raw.white_level
        libraw_wb = np.array(raw.camera_whitebalance)[:3]
        libraw_matrix = np.array(raw.color_matrix)[:, :3] @ np.diag(libraw_wb)
    ours = XYZ_TO_SRGB @ np.linalg.inv(cm)
    matrix_diff = np.abs(libraw_matrix / libraw_matrix.sum(1, keepdims=True) - ours / ours.sum(1, keepdims=True)).max()
    ok1 = (pattern == [[0, 1], [3, 2]] and libraw_black == [black] * 4 and libraw_white == white and mosaic.shape == picture.shape[:2]
           and mosaic.min() >= 0 and mosaic.max() <= white and matrix_diff < 3e-3
           and np.allclose(libraw_wb / libraw_wb[1], neutral[1] / neutral, rtol=1e-4))
    print("1. LibRaw %s reads %d x %d, pattern %s, black %d, white %d, as-shot multipliers %s, colour matrix off by %.1e: %s"
          % (".".join(map(str, rawpy.libraw_version)), mosaic.shape[1], mosaic.shape[0], pattern, black, white,
             np.round(libraw_wb / libraw_wb[1], 4), matrix_diff, "OK" if ok1 else "FAIL"))
    if not ok1:
        failed.append(1)

    # 2. the renderer's picture again from the DNG's numbers
    lin = (mosaic - black) / (white - black)
    rows, cols = np.indices(lin.shape)
    site = np.where(rows % 2 == 0, np.where(cols % 2 == 0, 0, 1), np.where(cols % 2 == 0, 1, 2))
    balanced = lin / neutral[site]
    rgb_raw = demosaicing_CFA_Bayer_Malvar2004(balanced, "RGGB") * neutral
    srgb = rgb_raw @ (XYZ_TO_SRGB @ np.linalg.inv(cm)).T
    inner = (slice(4, -4), slice(4, -4))
    a, b = srgb[inner].reshape(-1, 3), picture[inner].reshape(-1, 3)
    bright = a[:, 1] > 0.05 * a[:, 1].max()
    scale = np.median(b[bright, 1] / a[bright, 1])
    rel = np.abs(a[bright] * scale - b[bright]) / np.maximum(b[bright].mean(1, keepdims=True), 0.05 * b[bright].mean())
    ok2 = np.percentile(rel, 99) < 2e-3
    print("2. picture from the DNG alone (colour-demosaicing Malvar 2004, ColorMatrix1, AsShotNeutral) against the .pfm: relative "
          "difference median %.2e, 99th percentile %.2e, largest %.2e: %s" % (np.median(rel), np.percentile(rel, 99), rel.max(),
                                                                                "OK" if ok2 else "FAIL"))
    if not ok2:
        failed.append(2)

    # 3. LibRaw's own development
    with rawpy.imread(stem + ".dng") as raw:
        dev = raw.postprocess(use_camera_wb=True, output_color=rawpy.ColorSpace.sRGB, gamma=(1, 1), no_auto_bright=True,
                              output_bps=16, user_flip=0).astype(np.float64)
    a, b = dev[inner], picture[inner]
    s = (b.reshape(-1, 3).mean(0) / a.reshape(-1, 3).mean(0))
    mean_colour = s / s[1] - 1                       # the mean colour's channel ratios against the renderer's (0 = same)
    gy, gx = np.gradient(b.mean(-1))
    smooth = (np.hypot(gx, gy) < 0.002 * b.mean()) & (b.mean(-1) > 0.05 * b.mean(-1).max())
    ratio = (a[smooth] * s[1]) / b[smooth]
    ok3 = np.all(np.abs(mean_colour) < 0.02) and np.all(np.abs(np.median(ratio, 0) - 1) < 0.02)
    print("3. LibRaw's development (AHD, its colour pipeline) against the .pfm: mean colour channel ratios %s, smooth regions (%d pixels) "
          "median ratios %s: %s" % (np.round(mean_colour + 1, 4), smooth.sum(), np.round(np.median(ratio, 0), 4), "OK" if ok3 else "FAIL"))
    if not ok3:
        failed.append(3)

    print("RESULT:", "PASS" if not failed else "FAIL %s" % failed)
    return 0 if not failed else 2


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1]))
