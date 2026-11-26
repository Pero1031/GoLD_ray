"""Compares two PFM renders pixel by pixel (docs/GPU.md): the share of bit-identical pixels, the largest difference, and the means."""
import struct
import sys

import numpy as np


def read_pfm(path):
    with open(path, "rb") as f:
        kind = f.readline().strip()
        w, h = map(int, f.readline().split())
        scale = float(f.readline())
        data = np.frombuffer(f.read(), dtype="<f4" if scale < 0 else ">f4")
    return data.reshape(h, w, 3 if kind == b"PF" else 1).astype(np.float32)


def main():
    a, b = read_pfm(sys.argv[1]), read_pfm(sys.argv[2])
    same = np.all(a == b, axis=2)
    d = np.abs(a.astype(np.float64) - b.astype(np.float64))
    print(f"{a.shape[1]} x {a.shape[0]}: {100 * same.mean():.2f} % of the pixels bit-identical, "
          f"{100 * (d.max(axis=2) <= 1e-6 * np.maximum(1, np.abs(a).max(axis=2))).mean():.2f} % within 1e-6 relative")
    print(f"mean A {a.mean(axis=(0, 1))}, mean B {b.mean(axis=(0, 1))}, ratio B / A {b.mean() / a.mean():.6f}")
    print(f"largest |A - B| {d.max():.4g}; mean |A - B| {d.mean():.4g}")


if __name__ == "__main__":
    main()
