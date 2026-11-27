"""The kernel renderer (--device gpu or cpu-kernel, float) against the engine (--device cpu, double) on a scene file (docs/GPU.md).

The two renderers use different random numbers, so their pictures can only agree within the Monte Carlo noise. As in
mitsuba_compare.py, each renderer draws the scene twice with independent seeds; the noise of a pixel comes from the difference
of its two renders, and the means of the two renderers are compared:

  ratio = kernel / engine for the image mean (R G B), with its z (difference / standard error)
  per block of BxB pixels: z = (kernel - engine) / sqrt(var_kernel + var_engine); if both converge to the same picture,
  |z| > 3 happens in about 0.3 % of the blocks, and a systematic difference shows as many blocks of one sign.

It also renders the kernel scene once on the CPU and once on the GPU with the same seed and reports the share of bit-identical
pixels (the GPU module built without fast math gives 100 %).

With --aov the denoiser's features (<out>_albedo.pfm, <out>_normal.pfm) are compared the same way: most blocks of a feature image
are free of noise (a pixel inside one surface records the same albedo and normal for every sample), and there the two renderers must
agree up to float rounding; the blocks at edges, silhouettes and behind smooth glass are noisy, and are compared by z as the picture.
With --denoise as well, the means of the denoised pictures are printed (informative only: a denoised picture is a learned estimate).

    python tools/test/gpu_compare.py [--scene scenes/gpu_box.json] [--spp 1024] [--kernel gpu] [--aov] [-- --denoise]
"""
import argparse
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools", "lib"))
from hdr_io import read_pfm  # noqa: E402


def render(exe, scene, device, spp, seed, out, extra):
    cmd = [exe, scene, "--device", device, "--spp", str(spp), "--seed", str(seed), "--out", out] + extra
    t = time.time()
    r = subprocess.run(cmd, cwd=os.path.join(ROOT, "GoLD_rayt"), capture_output=True, text=True, errors="replace")
    if r.returncode != 0:
        sys.exit(f"render failed ({' '.join(cmd)}):\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    return read_pfm(out).astype(np.float64), time.time() - t


def block_means(img, b):
    h, w = img.shape[0] // b * b, img.shape[1] // b * b
    return img[:h, :w].reshape(h // b, b, w // b, b, -1).mean(axis=(1, 3))


def compare(e, k, b):
    """e, k: (2, H, W, 3), two renders of the engine and of the kernel renderer."""
    E, K = e.mean(0), k.mean(0)
    var_e = block_means((e[0] - e[1]) ** 2 / 2, b) / (2 * b * b)
    var_k = block_means((k[0] - k[1]) ** 2 / 2, b) / (2 * b * b)
    be, bk = block_means(E, b), block_means(K, b)
    lum = lambda x: x @ np.array([0.2126, 0.7152, 0.0722])
    var = var_e + var_k
    exact = var.max(axis=-1) <= (1e-7 * np.maximum(be, bk).max(axis=-1)) ** 2   # noiseless blocks (the environment alone)
    keep = (lum(be) > 0.02 * lum(be).max()) & ~exact
    z = (bk - be) / np.sqrt(var + 1e-30)
    rel = (bk - be) / np.maximum(be, 1e-12)
    se = np.sqrt((((e[0] - e[1]) ** 2 / 2).reshape(-1, 3).mean(0) + ((k[0] - k[1]) ** 2 / 2).reshape(-1, 3).mean(0)) / (2 * E.shape[0] * E.shape[1]))
    mean_e, mean_k = E.reshape(-1, 3).mean(0), K.reshape(-1, 3).mean(0)
    return dict(ratio=mean_k / mean_e, z_mean=(mean_k - mean_e) / se, blocks=int(keep.sum()), frac_z3=float((np.abs(z[keep]) > 3).mean()),
                mean_z=float(z[keep].mean()), med_rel=float(np.median(np.abs(rel[keep]))), noise_rel=float(np.median(np.sqrt(var)[keep] / np.maximum(be[keep], 1e-12))),
                exact=int(exact.sum()), exact_rel=float(np.abs(rel[exact]).max()) if exact.any() else 0.0)


def compare_feature(e, k, b):
    """e, k: (2, H, W, 3) feature images (albedo or normal) of the engine and of the kernel renderer, two renders each."""
    E, K = e.mean(0), k.mean(0)
    var_e = block_means((e[0] - e[1]) ** 2 / 2, b) / (2 * b * b)
    var_k = block_means((k[0] - k[1]) ** 2 / 2, b) / (2 * b * b)
    be, bk = block_means(E, b), block_means(K, b)
    var = (var_e + var_k).max(axis=-1)
    quiet = var <= 1e-14   # every sample of every pixel recorded the same features (up to float rounding)
    z = (bk - be) / np.sqrt(var_e + var_k + 1e-30)
    se = np.sqrt((((e[0] - e[1]) ** 2 / 2).reshape(-1, 3).mean(0) + ((k[0] - k[1]) ** 2 / 2).reshape(-1, 3).mean(0)) / (2 * E.shape[0] * E.shape[1]))
    d = K.reshape(-1, 3).mean(0) - E.reshape(-1, 3).mean(0)
    return dict(diff=d, z_mean=d / np.maximum(se, 1e-30), quiet=int(quiet.sum()), quiet_max=float(np.abs(bk - be)[quiet].max()) if quiet.any() else 0.0,
                noisy=int((~quiet).sum()), frac_z3=float((np.abs(z[~quiet]) > 3).mean()) if (~quiet).any() else 0.0,
                mean_z=float(z[~quiet].mean()) if (~quiet).any() else 0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scene", default="scenes/gpu_box.json")
    ap.add_argument("--spp", type=int, default=1024)
    ap.add_argument("--kernel", default="gpu", choices=["gpu", "cpu-kernel"])
    ap.add_argument("--exe", default=os.path.join(ROOT, "out", "cli-dev", "GoLD_rayt.exe"))
    ap.add_argument("--block", type=int, default=16)
    ap.add_argument("--seed-base", type=int, default=1, help="engine seeds base, base + 1; kernel seeds base + 2, base + 3")
    ap.add_argument("--folder", default=os.path.join(ROOT, "out", "gpu_compare"))
    ap.add_argument("--aov", action="store_true", help="also compare the albedo and normal features")
    ap.add_argument("extra", nargs="*", help="further options for every render (after --)")
    a = ap.parse_args()
    os.makedirs(a.folder, exist_ok=True)
    name = os.path.splitext(os.path.basename(a.scene))[0]
    extra = a.extra + (["--aov"] if a.aov else [])

    e, k, te, tk, outs_e, outs_k = [], [], [], [], [], []
    for seed in (a.seed_base, a.seed_base + 1):
        outs_e.append(os.path.join(a.folder, f"{name}_engine_{seed}.pfm"))
        img, t = render(a.exe, a.scene, "cpu", a.spp, seed, outs_e[-1], extra)
        e.append(img); te.append(t)
        outs_k.append(os.path.join(a.folder, f"{name}_{a.kernel}_{seed + 2}.pfm"))
        img, t = render(a.exe, a.scene, a.kernel, a.spp, seed + 2, outs_k[-1], extra)
        k.append(img); tk.append(t)
    s = compare(np.array(e), np.array(k), a.block)
    h, w = e[0].shape[:2]
    print(f"{a.scene}: {w}x{h}, {a.spp} spp, two renders each; blocks of {a.block}x{a.block} pixels")
    print(f"  time per render (whole program): engine {np.mean(te):.1f} s, kernel ({a.kernel}) {np.mean(tk):.1f} s")
    print(f"  image mean kernel / engine: R {s['ratio'][0]:.5f}  G {s['ratio'][1]:.5f}  B {s['ratio'][2]:.5f}   z {np.round(s['z_mean'], 2)}")
    print(f"  {s['blocks']} noisy blocks: |z| > 3 in {100 * s['frac_z3']:.2f} % (about 0.3 % expected), mean z {s['mean_z']:+.3f}, "
          f"median |relative difference| {100 * s['med_rel']:.3f} % against a noise of {100 * s['noise_rel']:.3f} %")
    if s["exact"]:
        print(f"  {s['exact']} noiseless blocks (environment only): largest relative difference {s['exact_rel']:.2e}")
    if a.aov:
        for feature in ("albedo", "normal"):
            fe = np.array([read_pfm(p[:-4] + f"_{feature}.pfm").astype(np.float64) for p in outs_e])
            fk = np.array([read_pfm(p[:-4] + f"_{feature}.pfm").astype(np.float64) for p in outs_k])
            f = compare_feature(fe, fk, a.block)
            print(f"  {feature}: image mean kernel - engine {np.array2string(f['diff'], precision=6)}, z {np.round(f['z_mean'], 2)}; "
                  f"{f['quiet']} noiseless blocks agree within {f['quiet_max']:.2e}; {f['noisy']} noisy blocks: |z| > 3 in {100 * f['frac_z3']:.2f} %, mean z {f['mean_z']:+.3f}")
    if "--denoise" in extra:
        de = np.array([read_pfm(p[:-4] + "_denoised.pfm").astype(np.float64) for p in outs_e]).mean(0)
        dk = np.array([read_pfm(p[:-4] + "_denoised.pfm").astype(np.float64) for p in outs_k]).mean(0)
        print(f"  denoised: image mean kernel / engine {np.array2string(dk.reshape(-1, 3).mean(0) / de.reshape(-1, 3).mean(0), precision=5)} (informative)")

    # the same seed on the CPU and the GPU: the share of bit-identical pixels
    p_cpu, _ = render(a.exe, a.scene, "cpu-kernel", 16, 7, os.path.join(a.folder, f"{name}_pair_cpu.pfm"), extra)
    p_gpu, _ = render(a.exe, a.scene, "gpu", 16, 7, os.path.join(a.folder, f"{name}_pair_gpu.pfm"), extra)
    same = np.all(p_cpu == p_gpu, axis=2).mean()
    print(f"  same seed, 16 spp, CPU kernel against GPU: {100 * same:.3f} % of the pixels bit-identical, largest difference {np.abs(p_cpu - p_gpu).max():.3g}")


if __name__ == "__main__":
    main()
