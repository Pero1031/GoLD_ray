"""
A dark forest behind a close-up, as an environment map (latitude-longitude, Radiance .hdr): the background of a photograph of
leaves taken with an open aperture, where the trees far behind are only dark green and the sunlit leaves among them become discs.

  python tools/data/make_forest_env.py out.hdr [--width 4096] [--seed 3] [--spots 260] [--gain 1]

  - Ground below the horizon: dark brown. Canopy above: dark green, made uneven by large soft blobs.
  - Middle-sized soft patches of dim, sunlit green.
  - Small bright spots (0.1 to 0.5 degrees): sunlit leaves and gaps to the sky. They are smaller than the blur of the lens,
    so each becomes a disc in the shape of the aperture. Most lie in the half of the sky behind the subject (-z), the camera's
    usual direction in the scene files; a few lie elsewhere.

The map follows the engine's convention (y up, u = (atan2(z, x) + pi) / 2 pi, top row = zenith); linear sRGB. The sun is not in
the map: give the scene a small sphere light.
"""
import argparse
import math

import numpy as np


def write_hdr(path, img):
    """Flat (uncompressed) RGBE."""
    m = img.max(axis=2)
    exp = np.where(m > 1e-32, np.ceil(np.log2(np.maximum(m, 1e-32))), 0)
    scale = np.where(m > 1e-32, 256.0 / np.exp2(exp), 0)
    rgbe = np.zeros(img.shape[:2] + (4,), dtype=np.uint8)
    rgbe[..., :3] = np.clip(img * scale[..., None], 0, 255).astype(np.uint8)
    rgbe[..., 3] = np.where(m > 1e-32, exp + 128, 0).astype(np.uint8)
    h, w = img.shape[:2]
    with open(path, "wb") as f:
        f.write(b"#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n" + f"-Y {h} +X {w}\n".encode())
        f.write(rgbe.tobytes())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--width", type=int, default=4096)
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--spots", type=int, default=260)
    ap.add_argument("--gain", type=float, default=1.0, help="multiplies the spots and patches")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    w = args.width
    h = w // 2

    theta = math.pi * (np.arange(h) + 0.5) / h            # from the zenith
    phi = 2 * math.pi * (np.arange(w) + 0.5) / w - math.pi
    st, ct = np.sin(theta), np.cos(theta)
    dirs = np.stack([st[:, None] * np.cos(phi)[None, :], np.repeat(ct[:, None], w, 1), st[:, None] * np.sin(phi)[None, :]], -1)

    img = np.zeros((h, w, 3))
    up = dirs[..., 1]
    ground = np.array([0.004, 0.003, 0.002])
    canopy = np.array([0.006, 0.010, 0.003])
    t = np.clip(up / 0.08 + 0.5, 0, 1)[..., None]
    img += ground * (1 - t) + canopy * t

    def random_dir(n, behind):
        """Directions, `behind` of them around -z (azimuth +-70 deg, elevation -8 .. 40 deg), the rest anywhere above -5 deg."""
        out = []
        for i in range(n):
            if i < behind:
                az = math.radians(rng.uniform(-70, 70))
                el = math.radians(rng.triangular(-8, 6, 40))
            else:
                az = rng.uniform(-math.pi, math.pi)
                el = math.radians(rng.uniform(-5, 60))
            out.append((math.cos(el) * math.sin(az), math.sin(el), -math.cos(el) * math.cos(az)))
        return np.array(out)

    def add_blobs(centres, radii_deg, colours, soft):
        for c, r, col in zip(centres, radii_deg, colours):
            r = math.radians(r)
            reach = r * (3.0 if soft else 1.0) + 2 * math.pi / w * 2
            el = math.asin(np.clip(c[1], -1, 1))
            th = math.pi / 2 - el
            y0 = max(0, int((th - reach) / math.pi * h) - 1)
            y1 = min(h, int((th + reach) / math.pi * h) + 2)
            sub = dirs[y0:y1]
            ang = np.arccos(np.clip(sub @ c, -1, 1))
            if soft:
                wgt = np.exp(-0.5 * (ang / r) ** 2)
            else:
                px = 2 * math.pi / w
                wgt = np.clip((r - ang) / px + 0.5, 0, 1)
            img[y0:y1] += wgt[..., None] * col

    # large dark unevenness of the canopy
    n = 40
    add_blobs(random_dir(n, 30), rng.uniform(6, 18, n), rng.uniform(0.3, 1.0, (n, 1)) * np.array([0.008, 0.014, 0.003]), True)
    # middle-sized soft patches of sunlit green
    n = 60
    hues = np.array([[0.035, 0.060, 0.010], [0.050, 0.055, 0.012], [0.025, 0.050, 0.015]])
    cols = hues[rng.integers(0, 3, n)] * rng.lognormal(0, 0.6, (n, 1)) * args.gain
    add_blobs(random_dir(n, 50), rng.uniform(1.0, 4.0, n), cols, True)
    # small bright spots
    n = args.spots
    hues = np.array([[0.55, 0.80, 0.15], [0.75, 0.80, 0.20], [0.95, 0.80, 0.35], [0.80, 0.88, 0.95]])
    pick = rng.choice(4, n, p=[0.45, 0.30, 0.15, 0.10])
    cols = hues[pick] * np.exp(rng.uniform(math.log(0.4), math.log(6.0), (n, 1))) * args.gain
    add_blobs(random_dir(n, int(n * 0.85)), rng.uniform(0.1, 0.5, n), cols, False)

    write_hdr(args.out, img.astype(np.float32))
    print(f"{w} x {h} -> {args.out}; mean {img.mean():.4f}, max {img.max():.2f}")


if __name__ == "__main__":
    main()
