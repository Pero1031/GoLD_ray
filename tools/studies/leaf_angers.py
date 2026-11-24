"""
PROSPECT against the ANGERS leaves (docs/LEAF.md): how well the measured pigments give the measured colour.

For each of the 308 leaves the structure parameter N is fitted to the measured reflectance and transmittance at 750-900 nm (where
only water and dry matter absorb, both measured); the pigments are the measured ones (anthocyanins: the estimates of Feret et al. 2017).
Then PROSPECT (alpha 40, the measurement geometry) is compared with the measurement in the visible: RMSE and CIEDE2000 of the colour
of R and of T under D65.

Usage (after python tools/data/leaf_optics.py fetch):  python tools/studies/leaf_angers.py
"""
import os
import sys

import colour
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "tools", "lib"))
from prospect import load_angers, load_table, prospect  # noqa: E402


def main():
    t = load_table()
    chem, wl, Rm, Tm = load_angers()
    m = len(wl)
    lam = t["lambda"][:m]
    assert np.allclose(wl, lam)
    vis = (lam >= 400) & (lam <= 780)
    nir = (lam >= 750) & (lam <= 900)
    cmf = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"]
    d65 = colour.SDS_ILLUMINANTS["D65"]
    W = np.array([cmf[x] for x in lam[vis]]) * np.array([d65[x] for x in lam[vis]])[:, None]
    W /= W[:, 1].sum()
    white = colour.XYZ_to_xy(W.sum(0))

    def lab(s):
        return colour.XYZ_to_Lab(s[..., vis] @ W, white)

    def run(N, c):
        R, T = prospect(t, N=N, alpha=40, **c)
        return R[:m], T[:m]

    rows = []
    for i in range(len(chem["chl"])):
        c = {k: chem[k][i] for k in chem}
        err = lambda N: np.mean((run(N, c)[0][nir] - Rm[i, nir]) ** 2 + (run(N, c)[1][nir] - Tm[i, nir]) ** 2)
        N = min(np.arange(1.0, 4.0, 0.1), key=err)
        N = min(np.arange(max(1.0, N - 0.1), N + 0.1, 0.005), key=err)
        R, T = run(N, c)
        rows.append((N, np.sqrt(np.mean((R[vis] - Rm[i, vis]) ** 2)), np.sqrt(np.mean((T[vis] - Tm[i, vis]) ** 2)),
                     colour.delta_E(lab(R), lab(Rm[i])), colour.delta_E(lab(T), lab(Tm[i])), c["ant"]))
    r = np.array(rows)
    print("leaves                      n    N   RMSE R  RMSE T   dE00 R (median / 90 %)   dE00 T (median / 90 %)")
    for name, sel in (("all", np.ones(len(r), bool)), ("anthocyanins > 3 ug/cm2", r[:, 5] > 3), ("anthocyanins <= 3", r[:, 5] <= 3)):
        q = r[sel]
        print(f"{name:26s} {sel.sum():4d}  {np.median(q[:, 0]):.2f}  {np.median(q[:, 1]):.4f}  {np.median(q[:, 2]):.4f}"
              f"   {np.median(q[:, 3]):4.1f} / {np.percentile(q[:, 3], 90):4.1f}            {np.median(q[:, 4]):4.1f} / {np.percentile(q[:, 4], 90):4.1f}")


if __name__ == "__main__":
    main()
