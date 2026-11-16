"""The figures of docs/LVC_BPT.md, section 7 (vertex merging), from the development cases of the self-check.

Input: the pictures that `--scene lvc` writes for one case with an --out of one's own (the reference, the render without merging and
the render with merging, as <stem>_ref.pfm, <stem>_m7.pfm, <stem>_m15.pfm), run from GoLD_rayt/:
    ..\\out\\cli-dev\\GoLD_rayt.exe --scene lvc --lvc-mask 1019 --spp 8 --width 384 --out <dir>/pic_defocus.png   (about 6.5 minutes)
    ..\\out\\cli-dev\\GoLD_rayt.exe --scene lvc --lvc-mask 1016 --spp 8 --width 384 --out <dir>/pic_focus.png     (about 1.5 minutes)

    python tools/figures/lvc_figures.py <dir> docs/img

Writes lvc_merging_through_glass.png (the caustic seen through the ball: reference, without merging, with merging) and
lvc_merging_focal_bias.png (the floor in focus: reference, with merging, and the signed difference).
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib"))
from hdr_io import read_pfm          # noqa: E402


def load(path):
    return np.asarray(read_pfm(path), dtype=np.float64)


def lum(img):
    return 0.2126 * img[..., 0] + 0.7152 * img[..., 1] + 0.0722 * img[..., 2]


def tone(img, white):
    return np.clip(img / white, 0, 1) ** (1 / 2.2)


def strip(n):
    fig, ax = plt.subplots(1, n, figsize=(12, 2.75))
    for a in ax:
        a.axis("off")
    return fig, ax


def main():
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)

    ref = load(os.path.join(src, "pic_defocus_ref.pfm"))
    m7 = load(os.path.join(src, "pic_defocus_m7.pfm"))
    m15 = load(os.path.join(src, "pic_defocus_m15.pfm"))
    white = np.percentile(lum(ref), 99.5)
    fig, ax = strip(3)
    for a, im, t in zip(ax, (ref, m7, m15), ("Reference (path tracing + photon map)", "Without merging", "With merging")):
        a.imshow(tone(im, white))
        a.set_title(t, fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "lvc_merging_through_glass.png"), dpi=100)
    plt.close(fig)
    print("through glass: mean with merging / reference = %.3f" % (lum(m15).mean() / lum(ref).mean()))

    ref = load(os.path.join(src, "pic_focus_ref.pfm"))
    m15 = load(os.path.join(src, "pic_focus_m15.pfm"))
    white = np.percentile(lum(ref), 99.0)
    d = lum(m15) - lum(ref)
    dl = np.percentile(np.abs(d), 99.7)
    fig, ax = strip(3)
    ax[0].imshow(tone(ref, white)); ax[0].set_title("Reference", fontsize=10)
    ax[1].imshow(tone(m15, white)); ax[1].set_title("With merging (r = 3 mm)", fontsize=10)
    ax[2].imshow(d, cmap="RdBu_r", vmin=-dl, vmax=dl); ax[2].set_title("Difference (red: too bright, blue: too dark)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "lvc_merging_focal_bias.png"), dpi=100)
    plt.close(fig)
    print("in focus: mean with merging / reference = %.3f" % (lum(m15).mean() / lum(ref).mean()))


if __name__ == "__main__":
    main()
