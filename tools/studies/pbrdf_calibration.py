"""Spectralon as a transfer standard for the KAIST pBRDF gold table.

    python tools/studies/pbrdf_calibration.py     -> out/studies/measured_gold/calibration.txt

The Spectralon sample (dataset material 1, reflectance about 0.99, close to Lambertian up to about 60 deg) should give
fr = 0.99 / pi = 0.315 1/sr everywhere. What the table gives instead (scale, dependence on the incidence angle, excess
near theta_h = 0, per-wavelength differences) is the systematic error of the measurement, and the gold table shares it.
So the gold is corrected bin by bin: fr_gold_corrected = fr_gold * (0.99 / pi) / S(theta_h, theta_d), where S is the
Spectralon table median over phi_d (variant A: as measured; variant B: the excess at theta_h < 5 deg replaced by the
value at 5-8 deg, in case that excess is real Spectralon gloss rather than an instrument effect). Then the hemispherical
reflectance and the spectral shape of the gold are recomputed with both variants.
"""
import os
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("tools", 1)[0] + "tools/studies")
from measured_gold import directions, half_diff, fresnel_conductor, load_au   # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REF = os.path.join(ROOT, "references", "kaist_pbrdf")
OUT = os.path.join(ROOT, "out", "studies", "measured_gold")
LAMBERT = 0.99 / np.pi
lines = []


def say(s=""):
    print(s)
    lines.append(s)


def main():
    os.makedirs(OUT, exist_ok=True)
    g = np.load(os.path.join(REF, "gold_unpolarized.npz"))
    s = np.load(os.path.join(REF, "spectralon_unpolarized.npz"))
    gold = np.transpose(g["fr"], (2, 1, 0, 3)).astype(np.float64)        # [theta_h, theta_d, phi_d, wavelength]
    spec = np.transpose(s["fr"], (2, 1, 0, 3)).astype(np.float64)
    theta_h, theta_d, phi_d = (g[k].astype(np.float64) for k in ("theta_h", "theta_d", "phi_d"))
    wvls = g["wavelengths_nm"].astype(float)
    eta, kk = load_au(wvls)
    TH, TD, PD = np.meshgrid(theta_h, theta_d, phi_d, indexing="ij")
    wi, wo = directions(TH, TD, PD)
    ci, co = wi[..., 2], wo[..., 2]
    above = (ci > 0) & (co > 0)

    # 1. Spectralon: scale and spectral flatness near normal incidence
    near = above & (ci > 0.94) & (co > 0.94) & np.all(np.isfinite(spec), axis=-1)
    med = np.median(spec[near], axis=0)
    say(f"[Calib] Spectralon near normal (theta_i, theta_o < 20 deg): median {np.round(med, 4).tolist()} 1/sr, "
        f"expected {LAMBERT:.4f}; scale {np.round(med / LAMBERT, 3).tolist()}; relative to 650 nm {np.round(med / med[-1], 3).tolist()}")
    # by incidence angle, at theta_o 20-40 deg
    say("[Calib] Spectralon median (550 nm) by theta_i, theta_o in 20-40 deg; a Lambertian gives 0.315 at every angle:")
    ti = np.degrees(np.arccos(np.clip(ci, 0, 1)))
    to = np.degrees(np.arccos(np.clip(co, 0, 1)))
    for a in range(0, 90, 10):
        m = above & np.isfinite(spec[..., 2]) & (ti >= a) & (ti < a + 10) & (to >= 20) & (to < 40) & (TH > np.radians(5))
        say(f"     theta_i {a:2d}-{a + 10:2d}: {np.median(spec[m][:, 2]):.3f}  (x{np.median(spec[m][:, 2]) / LAMBERT:.2f}; cos theta_i {np.cos(np.radians(a + 5)):.2f})")

    # 2. the correction field S(theta_h, theta_d, wavelength): Spectralon median over phi_d, holes filled along theta_d
    S = np.nanmedian(spec, axis=2)                                        # [theta_h, theta_d, wavelength]
    for li in range(S.shape[-1]):
        for ih in range(S.shape[0]):
            col = S[ih, :, li]
            bad = ~np.isfinite(col)
            if bad.all():
                continue
            col[bad] = np.interp(np.flatnonzero(bad), np.flatnonzero(~bad), col[~bad])
    S_A = S
    S_B = S.copy()
    ridge = theta_h < np.radians(5)
    ref = np.nanmedian(S[(theta_h >= np.radians(5)) & (theta_h < np.radians(8))], axis=0)   # [theta_d, wavelength]
    S_B[ridge] = ref[None]
    excess = np.nanmedian(S_A[theta_h < np.radians(0.5)] / ref[None], axis=(0, 1))
    say(f"[Calib] Spectralon excess at theta_h < 0.5 deg relative to theta_h 5-8 deg (median over theta_d): {np.round(excess, 2).tolist()}")

    # 3. corrected gold: hemispherical reflectance by the same integration as measured_gold.py
    def lookup_idx(th_, td_, pd_):
        ih = np.clip(np.rint(90.0 * np.sqrt(np.clip(th_ / (np.pi / 2), 0, 1))).astype(int), 0, 90)
        id_ = np.clip(np.rint(np.degrees(td_)).astype(int), 0, 90)
        ip = np.clip(np.rint(np.degrees(pd_)).astype(int) + 180, 0, 360)
        return ih, id_, ip

    n_to, n_po = 360, 720
    t_o = (np.arange(n_to) + 0.5) / n_to * (np.pi / 2)
    p_o = (np.arange(n_po) + 0.5) / n_po * (2 * np.pi)
    TO, PO = np.meshgrid(t_o, p_o, indexing="ij")
    wo_g = np.stack([np.sin(TO) * np.cos(PO), np.sin(TO) * np.sin(PO), np.cos(TO)], -1)
    dw_g = np.sin(TO) * (np.pi / 2 / n_to) * (2 * np.pi / n_po)
    n_fine, cone = 1201, np.radians(3.0)
    u = np.linspace(-cone, cone, n_fine)
    U, V = np.meshgrid(u, u, indexing="ij")
    ang = np.sqrt(U ** 2 + V ** 2)
    inside = ang <= cone
    dw_f = (2 * cone / (n_fine - 1)) ** 2 * np.sinc(ang / np.pi)

    def integrate(table, Sfield, deg):
        t = np.radians(deg)
        wi_v = np.array([np.sin(t), 0.0, np.cos(t)])
        mirror = np.array([-np.sin(t), 0.0, np.cos(t)])
        t1 = np.array([np.cos(t), 0.0, np.sin(t)])
        t2 = np.array([0.0, 1.0, 0.0])
        wo_f = mirror * np.cos(ang)[..., None] + (t1 * U[..., None] + t2 * V[..., None]) * np.sinc(ang / np.pi)[..., None]
        wo_f /= np.linalg.norm(wo_f, axis=-1, keepdims=True)
        total = np.zeros(len(wvls))
        for wo_, dw_, mask in ((wo_f, dw_f, inside), (wo_g, dw_g, np.arccos(np.clip(wo_g @ mirror, -1, 1)) > cone)):
            ih, id_, ip = lookup_idx(*half_diff(np.broadcast_to(wi_v, wo_.shape), wo_))
            vals = table[ih, id_, ip, :]
            if Sfield is not None:
                vals = vals * LAMBERT / Sfield[ih, id_, :]
            w = (wo_[..., 2] * dw_ * mask)[..., None]
            total += np.nansum(np.where(np.isfinite(vals), vals * w, 0.0), axis=(0, 1))
        return total

    say("[Calib] gold hemispherical reflectance per wavelength (450 500 550 600 650 nm):")
    say("     theta_i   as measured                 corrected A (Spectralon as is)   corrected B (ridge removed)      Fresnel(theta_i), Johnson Au")
    rows = []
    for deg in (15, 30, 45, 60, 75):
        raw = integrate(gold, None, deg)
        ca = integrate(gold, S_A, deg)
        cb = integrate(gold, S_B, deg)
        F = fresnel_conductor(np.cos(np.radians(deg)), eta, kk)
        rows.append((deg, raw, ca, cb, F))
        say(f"     {deg:3d} deg   {np.round(raw, 3).tolist()}   {np.round(ca, 3).tolist()}   {np.round(cb, 3).tolist()}   {np.round(F, 3).tolist()}")
    say("[Calib] corrected / Fresnel(Johnson): A, then B")
    for deg, raw, ca, cb, F in rows:
        say(f"     {deg:3d} deg   {np.round(ca / F, 2).tolist()}   {np.round(cb / F, 2).tolist()}")

    # 4. spectral shape of the corrected spike and haze (650 nm = 1), against Fresnel
    fin = np.all(np.isfinite(gold), axis=-1) & above
    spike = fin & (TH < np.radians(0.5)) & (TD > np.radians(4)) & (TD < np.radians(40))
    haze = fin & (TH > np.radians(10)) & (TH < np.radians(35)) & (TD > np.radians(4)) & (TD < np.radians(40))
    F0 = fresnel_conductor(1.0, eta, kk)
    say(f"[Calib] Fresnel F0 Johnson, relative to 650 nm: {np.round(F0 / F0[-1], 3).tolist()}")
    for name, region in (("spike", spike), ("haze", haze)):
        Fr = fresnel_conductor(np.cos(TD[region])[:, None], eta[None, :], kk[None, :])
        for label, Sf in (("as measured", None), ("corrected A", S_A), ("corrected B", S_B)):
            vals = gold[region]
            if Sf is not None:
                ih = np.broadcast_to(np.arange(91)[:, None, None], TH.shape)[region]
                id_ = np.broadcast_to(np.arange(91)[None, :, None], TH.shape)[region]
                vals = vals * LAMBERT / Sf[ih, id_, :]
            r = np.exp(np.mean(np.log(vals / Fr), axis=0))
            say(f"[Calib] {name} {label:12s}: measured / Fresnel relative to 650 nm {np.round(r / r[-1], 3).tolist()}; absolute at 650 nm {r[-1]:.3f}")

    with open(os.path.join(OUT, "calibration.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote", os.path.join(OUT, "calibration.txt"))


if __name__ == "__main__":
    main()
