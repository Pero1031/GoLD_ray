"""Fits a Kubelka-Munk fading law to every patch of a measured fading series, so that the renderer can ask for any time:
between the rounds, after the last one (extrapolated) and before the first one (negative hours: the patch as it was
before the exposure, a "renewal" of an already faded colour).

The law (per wavelength, one patch):
    Q = K/S = (1 - r)^2 / (2 r),  r = R / s
    Q(t) = Qinf + (Q0 - Qinf) exp(-k t)
    R(t) = s * (1 + Q - sqrt(Q^2 + 2 Q))
The colorant's absorption decays exponentially toward what is left when it is gone (Qinf: the paper, the binder, a
residue). Q0 comes from the first round (lightly smoothed), k(lambda) > 0 and Qinf(lambda) > 0 are smooth curves (one
value every 10 nm, linear between) fitted by least squares to all rounds. s >= 1 is a fixed scale per patch for
reflectance factors above 1 (bright dyes measured against a white reference that is darker than an ideal white): K/S
needs r < 1.

On most measured colours the law reproduces the rounds within about 1 CIEDE2000 (D65), while a model that predicts a
colorant's fading from its first colour alone (shared rates, the shared principal-component law of Kataoka and Okajima,
CSAJ 2024, or this law with shared k and Qinf) misses by several: how fast a colorant fades is a property of the
colorant, not of its colour. So the fit is per patch, and only for patches that were measured over time. A change that is complete before the first
measured time cannot tell how fast it was: k is capped at 4 / (the first interval), which keeps the law finite when it
is run backward. Before the first round and after the last one the law is an extrapolation.

Input : a series folder of tools/data/fading_refangle.py (hours.csv, r1/<color>.csv, r2/..., 1 nm curves).
Output: <series>/km/<color>.csv  "wavelength_nm,R0,Rinf,k_per_hour,scale" (R0, Rinf are reflectances, the scale is
        repeated on every row), and <series>/km/report.csv with the CIEDE2000 of every round.

Usage:  python tools/data/fading_km_fit.py GoLD_rayt/assets/Fading_measured/dye_paper [more series ...]
"""
import csv
import glob
import os
import sys

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import savgol_filter

try:
    import colour
except ImportError:  # the report needs colour-science; the fit does not
    colour = None


def read_curve(path):
    wl, r = [], []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            a, b = line.split(",")[:2]
            try:
                wl.append(float(a)); r.append(float(b))
            except ValueError:
                continue
    return np.array(wl), np.array(r)


def read_hours(path):
    rounds = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            a, b = line.split(",")[:2]
            try:
                rounds.append((int(a), float(b)))
            except ValueError:
                continue
    return rounds


def ks(r):
    r = np.clip(r, 1e-3, 0.999)
    return (1 - r) ** 2 / (2 * r)


def ks_inv(q):
    # 1 + q - sqrt(q^2 + 2q), written so that a huge q gives a small R instead of inf - inf
    return 1 / (1 + q + np.sqrt(q * q + 2 * q))


def predict(R0, s, k, Qinf, hours):
    Q0 = ks(R0 / s)
    return s * ks_inv(Qinf + (Q0 - Qinf) * np.exp(-np.outer(hours, k)))


class Colorimeter:
    """CIEDE2000 under D65 (CIE 1931 2 deg); curves held at their ends outside their range, as the renderer does."""

    def __init__(self):
        shape = colour.SpectralShape(380, 780, 1)
        self.wl = np.arange(380, 781)
        self.cmf = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"].copy().align(shape).values
        self.ill = colour.SDS_ILLUMINANTS["D65"].copy().align(shape).values
        self.k = 1 / (self.ill * self.cmf[:, 1]).sum()
        self.white = colour.XYZ_to_xy(self.k * (self.ill[:, None] * self.cmf).sum(0))

    def lab(self, wl, R):
        Rf = np.array([np.interp(self.wl, wl, r) for r in np.atleast_2d(R)])
        return colour.XYZ_to_Lab(self.k * (Rf * self.ill) @ self.cmf, self.white)

    def de(self, wl, A, B):
        return colour.delta_E(self.lab(wl, A), self.lab(wl, B), method="CIE 2000")


def fit_patch(wl, curves, hours, knot_step=10.0):
    R0 = savgol_filter(curves[0], 21, 3) if len(wl) > 21 else curves[0].copy()
    s = max(1.0, 1.02 * max(savgol_filter(c, 21, 3).max() if len(wl) > 21 else c.max() for c in curves))
    knots = np.arange(wl[0], wl[-1] + 1e-6, knot_step)
    if knots[-1] < wl[-1]:
        knots = np.append(knots, wl[-1])
    nk = len(knots)
    T = max(hours)
    last = savgol_filter(curves[-1], 21, 3) if len(wl) > 21 else curves[-1]

    def unpack(p):
        return np.exp(np.interp(wl, knots, p[:nk])), np.exp(np.interp(wl, knots, p[nk:]))

    def residual(p):
        k, Qinf = unpack(p)
        return (predict(R0, s, k, Qinf, hours) - curves).ravel()

    # A change that is over before the first measurement after 0 h only bounds k from below; any faster k fits as well
    # and would explode when the law is run backward. So the time constant is kept at least a quarter of that interval.
    first = min(h for h in hours if h > 0)
    k_hi = np.log(4.0 / first)
    p0 = np.r_[np.full(nk, np.log(1.0 / T)), np.log(np.interp(knots, wl, ks(last / s)) * 0.7 + 1e-4)]
    lo = np.r_[np.full(nk, np.log(1e-7)), np.full(nk, np.log(1e-6))]
    hi = np.r_[np.full(nk, k_hi), np.full(nk, np.log(1e3))]
    sol = least_squares(residual, np.clip(p0, lo + 1e-9, hi - 1e-9), bounds=(lo, hi))
    k, Qinf = unpack(sol.x)
    Rinf = s * ks_inv(Qinf)
    return R0, Rinf, k, s


def main(series_dirs):
    meter = Colorimeter() if colour else None
    for d in series_dirs:
        rounds = read_hours(os.path.join(d, "hours.csv"))
        hours = np.array([h for _, h in rounds])
        colors = sorted(os.path.splitext(os.path.basename(p))[0] for p in glob.glob(os.path.join(d, f"r{rounds[0][0]}", "*.csv")))
        out_dir = os.path.join(d, "km")
        os.makedirs(out_dir, exist_ok=True)
        report = []
        print(f"{d}: {len(colors)} patches, rounds at {', '.join(f'{h:g}' for h in hours)} h")
        for c in colors:
            wl, _ = read_curve(os.path.join(d, f"r{rounds[0][0]}", c + ".csv"))
            curves = np.array([read_curve(os.path.join(d, f"r{r}", c + ".csv"))[1] for r, _ in rounds])
            R0, Rinf, k, s = fit_patch(wl, curves, hours)
            with open(os.path.join(out_dir, c + ".csv"), "w", newline="", encoding="utf-8") as f:
                f.write(f"# Kubelka-Munk fading law of {os.path.basename(os.path.normpath(d))}/{c} (tools/data/fading_km_fit.py)\n")
                f.write("# R(l, h) = scale * ksinv(Qinf + (Q0 - Qinf) exp(-k h)), Q0 = ks(R0 / scale), Qinf = ks(Rinf / scale), h in hours of the series\n")
                w = csv.writer(f)
                w.writerow(["wavelength_nm", "R0", "Rinf", "k_per_hour", "scale"])
                for i in range(len(wl)):
                    w.writerow([f"{wl[i]:g}", f"{R0[i]:.6f}", f"{Rinf[i]:.6f}", f"{k[i]:.6e}", f"{s:.4f}"])
            if meter:
                P = predict(R0, s, k, ks(Rinf / s), hours)
                de = meter.de(wl, P, curves)
                report.append([c] + [f"{x:.2f}" for x in de])
                print(f"  {c:10s} dE00 per round " + " ".join(f"{x:4.2f}" for x in de) + f"   (mean after 0 h {de[1:].mean():.2f}, scale {s:.2f})")
        if report:
            with open(os.path.join(out_dir, "report.csv"), "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(["color"] + [f"dE00_{h:g}h" for h in hours])
                w.writerows(report)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    main(sys.argv[1:])
