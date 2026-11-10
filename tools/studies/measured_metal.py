"""A measured metal BRDF table (KAIST pBRDF 2020 or hpBRDF 2025, unpolarized) against the engine's rough conductor.

    python tools/studies/measured_metal.py <table.npz> <metal.csv> <out_dir> <name>
    e.g.   ... references/kaist_pbrdf/gold_unpolarized.npz GoLD_rayt/assets/Metal_data/Johnson_Au.csv out/studies/measured_gold gold
           ... references/kaist_pbrdf/aluminum_unpolarized.npz GoLD_rayt/assets/Metal_data/Rakic_Al.csv out/studies/measured_aluminum aluminum

The engine's model (Materials/RoughConductor): GGX distribution with alpha = roughness^2, height-correlated Smith
masking, conductor Fresnel from the metal's (n, k) table (assets/Metal_data, the engine's two-block CSV), plus the
Kulla-Conty multiple-scattering lobe. Everything is re-implemented here in numpy; Beckmann is added for the lobe shape.

Table: fr[phi_d, theta_d, theta_h, wavelength] (Mueller element [0, 0]) in Rusinkiewicz half-difference coordinates,
made by references/kaist_pbrdf/extract_unpolarized.py (2020 data) or extract_hpbrdf.py (2025 data). Directions follow
the MERL BRDFRead.cpp convention (phi_h = 0 for an isotropic BRDF).

Measured:
  1. the roughness that fits the specular lobe (shape only), for GGX and Beckmann, and a two-lobe GGX fit;
  2. the spectral shape of the spike and of the haze, measured / Fresnel, relative to the reference wavelength;
  3. the directional albedo per wavelength (fine integration around the mirror direction) against the model.
"""
import csv
import os
import sys

import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.special import erf

lines = []


def say(s=""):
    print(s)
    lines.append(s)


# ---------------------------------------------------------------------------------------------
# The engine's model, in numpy
# ---------------------------------------------------------------------------------------------


def load_metal(path, wavelengths_nm):
    """eta, k at the wavelengths, as ComplexIorTable::fromCurves (union grid, linear interpolation, clamped)."""
    sections, cur = {}, None
    for row in csv.reader(open(path, encoding="utf-8-sig")):
        if not row:
            continue
        if row[0] == "wl":
            cur = row[1]
            sections[cur] = []
        else:
            sections[cur].append((float(row[0]) * 1000.0, float(row[1])))
    n, k = np.array(sorted(sections["n"])), np.array(sorted(sections["k"]))
    grid = np.union1d(n[:, 0], k[:, 0])
    eta_g = np.interp(grid, n[:, 0], n[:, 1])
    k_g = np.interp(grid, k[:, 0], k[:, 1])
    return np.interp(wavelengths_nm, grid, eta_g), np.interp(wavelengths_nm, grid, k_g)


def load_au(wavelengths_nm):   # kept for the older scripts
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    return load_metal(os.path.join(root, "GoLD_rayt", "assets", "Metal_data", "Johnson_Au.csv"), wavelengths_nm)


def fresnel_conductor(cos_i, eta, k):
    """Unpolarized conductor Fresnel (port of fresnel::fresnelConductor). cos_i broadcasts against eta, k."""
    cos_i = np.clip(np.abs(cos_i), 0.0, 1.0)
    c2 = cos_i ** 2
    s2 = 1.0 - c2
    eta2, k2 = eta ** 2, k ** 2
    t0 = eta2 - k2 - s2
    a2b2 = np.sqrt(t0 * t0 + 4.0 * eta2 * k2)
    t1 = a2b2 + c2
    a = np.sqrt(np.maximum(0.5 * (a2b2 + t0), 0.0))
    t2 = 2.0 * cos_i * a
    rs = (t1 - t2) / (t1 + t2)
    t3 = c2 * a2b2 + s2 * s2
    t4 = t2 * s2
    rp = rs * (t3 - t4) / (t3 + t4)
    return 0.5 * (rs + rp)


def D_ggx(cos_h, alpha):
    cos_h = np.maximum(cos_h, 1e-9)
    tan2 = (1.0 - cos_h ** 2) / cos_h ** 2
    return 1.0 / (np.pi * alpha ** 2 * cos_h ** 4 * (1.0 + tan2 / alpha ** 2) ** 2)


def lambda_ggx(cos_t, alpha):
    cos_t = np.maximum(np.abs(cos_t), 1e-9)
    tan2 = (1.0 - cos_t ** 2) / cos_t ** 2
    return 0.5 * (np.sqrt(1.0 + alpha ** 2 * tan2) - 1.0)


def D_beckmann(cos_h, alpha):
    cos_h = np.maximum(cos_h, 1e-9)
    tan2 = (1.0 - cos_h ** 2) / cos_h ** 2
    return np.exp(-tan2 / alpha ** 2) / (np.pi * alpha ** 2 * cos_h ** 4)


def lambda_beckmann(cos_t, alpha):
    cos_t = np.maximum(np.abs(cos_t), 1e-9)
    tan_t = np.sqrt(np.maximum(1.0 - cos_t ** 2, 0.0)) / cos_t
    a = 1.0 / np.maximum(alpha * tan_t, 1e-9)
    return 0.5 * (erf(a) - 1.0) + np.exp(-a * a) / (2.0 * a * np.sqrt(np.pi))


MODELS = {"GGX": (D_ggx, lambda_ggx), "Beckmann": (D_beckmann, lambda_beckmann)}


def single_scatter(kind, alpha, cos_h, cos_i, cos_o, F):
    """D G F / (4 cos_i cos_o); F carries a trailing wavelength axis."""
    D, lam = MODELS[kind]
    G = 1.0 / (1.0 + lam(cos_o, alpha) + lam(cos_i, alpha))
    geo = D(cos_h, alpha) * G / (4.0 * np.maximum(cos_i, 1e-9) * np.maximum(cos_o, 1e-9))
    return geo[..., None] * F


def directional_albedo_table(kind, alpha, n_mu=64, n_to=256, n_po=256):
    """E(mu) of the single-scattering lobe with F = 1 (GGXAlbedoTable), by quadrature over the outgoing hemisphere."""
    D, lam = MODELS[kind]
    mus = (np.arange(n_mu) + 0.5) / n_mu
    to = (np.arange(n_to) + 0.5) / n_to * (np.pi / 2)
    po = (np.arange(n_po) + 0.5) / n_po * (2 * np.pi)
    TO, PO = np.meshgrid(to, po, indexing="ij")
    wo = np.stack([np.sin(TO) * np.cos(PO), np.sin(TO) * np.sin(PO), np.cos(TO)], -1)
    dw = np.sin(TO) * (np.pi / 2 / n_to) * (2 * np.pi / n_po)
    E = np.zeros(n_mu)
    for i, mu in enumerate(mus):
        wi = np.array([np.sqrt(1 - mu * mu), 0.0, mu])
        h = wi + wo
        h /= np.linalg.norm(h, axis=-1, keepdims=True)
        G = 1.0 / (1.0 + lam(wo[..., 2], alpha) + lam(mu, alpha))
        f = D(h[..., 2], alpha) * G / (4.0 * mu * wo[..., 2])
        E[i] = np.sum(f * wo[..., 2] * dw)
    E_avg = 2.0 * np.sum(E * mus) / n_mu
    return mus, np.minimum(E, 1.0), E_avg


def average_fresnel(eta, k, n=64):
    mu = (np.arange(n) + 0.5) / n
    return 2.0 * np.sum(fresnel_conductor(mu[:, None], eta[None, :], k[None, :]) * mu[:, None], axis=0) / n


def ms_lobe(mus, E, E_avg, F_avg, cos_i, cos_o):
    """Kulla-Conty compensation lobe (RoughConductor::msBrdf), per wavelength."""
    F_ms = F_avg ** 2 * E_avg / (1.0 - F_avg * (1.0 - E_avg))
    Ei = np.interp(cos_i, mus, E)
    Eo = np.interp(cos_o, mus, E)
    scale = (1.0 - Eo) * (1.0 - Ei) / (np.pi * (1.0 - E_avg))
    return scale[..., None] * F_ms


# ---------------------------------------------------------------------------------------------
# Half-difference coordinates (MERL convention, phi_h = 0)
# ---------------------------------------------------------------------------------------------


def directions(th, td, pd):
    """(theta_h, theta_d, phi_d) -> wi, wo (arrays [..., 3]); the half vector lies in the x-z plane."""
    dx, dy, dz = np.sin(td) * np.cos(pd), np.sin(td) * np.sin(pd), np.cos(td)
    ix = dx * np.cos(th) + dz * np.sin(th)
    iy = dy
    iz = -dx * np.sin(th) + dz * np.cos(th)
    wi = np.stack([ix, iy, iz], -1)
    half = np.stack([np.sin(th), np.zeros_like(th), np.cos(th)], -1)
    wo = 2.0 * np.sum(wi * half, -1, keepdims=True) * half - wi
    return wi, wo


def half_diff(wi, wo):
    """wi, wo [..., 3] -> theta_h, theta_d, phi_d (MERL std_coords_to_half_diff_coords)."""
    h = wi + wo
    h /= np.maximum(np.linalg.norm(h, axis=-1, keepdims=True), 1e-12)
    th = np.arccos(np.clip(h[..., 2], -1, 1))
    ph = np.arctan2(h[..., 1], h[..., 0])
    c, s = np.cos(-ph), np.sin(-ph)
    x1 = wi[..., 0] * c - wi[..., 1] * s
    y1 = wi[..., 0] * s + wi[..., 1] * c
    z1 = wi[..., 2]
    c, s = np.cos(-th), np.sin(-th)
    x2 = x1 * c + z1 * s
    z2 = -x1 * s + z1 * c
    td = np.arccos(np.clip(z2, -1, 1))
    pd = np.arctan2(y1, x2)
    return th, td, pd


# ---------------------------------------------------------------------------------------------


def main():
    if len(sys.argv) < 5:
        print(__doc__)
        sys.exit(1)
    data_path, metal_path, out_dir, name = sys.argv[1:5]
    os.makedirs(out_dir, exist_ok=True)
    d = np.load(data_path)
    fr = np.transpose(d["fr"], (2, 1, 0, 3)).astype(np.float64)   # [theta_h, theta_d, phi_d, wavelength]
    theta_h, theta_d, phi_d = (d[k].astype(np.float64) for k in ("theta_h", "theta_d", "phi_d"))
    wvls = d["wavelengths_nm"].astype(float)
    n_lam = len(wvls)
    eta, kk = load_metal(metal_path, wvls)
    # the wavelengths that are printed and drawn (all are used in the computations)
    show = sorted(set(int(np.argmin(np.abs(wvls - w))) for w in (450, 500, 550, 600, 650, 700, 800, 900) if wvls.min() - 1 <= w <= wvls.max() + 1))
    i_ref = int(np.argmin(np.abs(wvls - 650)))
    i550 = int(np.argmin(np.abs(wvls - 550)))

    def fmt(arr, nd=3):
        return "[" + " ".join(f"{int(round(wvls[i]))}:{arr[i]:.{nd}f}" for i in show) + "]"

    say(f"[Measured:{name}] table {fr.shape} (theta_h, theta_d, phi_d, wavelength), {n_lam} wavelengths {wvls.min():.0f}-{wvls.max():.0f} nm; shown: {[int(round(wvls[i])) for i in show]}")
    say(f"[Measured:{name}] metal table {os.path.basename(metal_path)}: eta {fmt(eta)}, k {fmt(kk)}, F0 {fmt(fresnel_conductor(1.0, eta, kk))}")

    TH, TD, PD = np.meshgrid(theta_h, theta_d, phi_d, indexing="ij")
    wi, wo = directions(TH, TD, PD)
    cos_i, cos_o = wi[..., 2], wo[..., 2]
    cos_h = np.cos(TH)
    cos_d = np.cos(TD)
    finite = np.all(np.isfinite(fr), axis=-1)
    above = (cos_i > 0.0) & (cos_o > 0.0)
    say(f"[Measured:{name}] cells with both directions above the surface: {above.mean() * 100:.1f} %; of those measured (finite): {finite[above].mean() * 100:.1f} %")

    # 0. reciprocity: phi_d -> phi_d + pi
    j = np.arange(len(phi_d))
    a, b = fr[..., i550], fr[:, :, (j + 180) % 360, i550]
    for thr in (1e-3, 0.5):
        ok = np.isfinite(a) & np.isfinite(b) & (a > thr) & (b > thr) & above & (TH < np.radians(85))
        ratio = np.log(a[ok] / b[ok])
        say(f"[Measured:{name}] reciprocity (phi_d vs phi_d + pi), cells above {thr} 1/sr at 550 nm: median |log ratio| {np.median(np.abs(ratio)):.3f}, 90 % {np.percentile(np.abs(ratio), 90):.3f}  ({ok.sum()} pairs)")

    # the far tail
    tail = finite & (cos_i > 0.3) & (cos_o > 0.3) & (TH > np.radians(50))
    floor = np.nanmedian(fr[tail], axis=0)
    say(f"[Measured:{name}] far tail (theta_h > 50 deg) median: {fmt(floor, 4)} 1/sr = Lambertian albedo {fmt(np.pi * floor)}")

    # 1. single-lobe fits on the specular lobe (shape only)
    sel = finite & (cos_i > 0.1) & (cos_o > 0.1) & (TH < np.radians(35)) & (fr[..., i550] > 20 * floor[i550])
    say(f"[Measured:{name}] fit region: theta_h < 35 deg, value > 20 x floor: {sel.sum()} cells")
    F = fresnel_conductor(cos_d[sel], eta[i550], kk[i550])
    meas = np.log(fr[sel][:, i550])
    dth = np.gradient(theta_h)
    w = (np.sin(TH) * dth[:, None, None] * np.sin(TD))[sel]
    w /= w.sum()
    fits = {}
    for kind in MODELS:
        def cost(alpha):
            geo = single_scatter(kind, alpha, cos_h[sel], cos_i[sel], cos_o[sel], F[:, None])[:, 0]
            r = meas - np.log(np.maximum(geo, 1e-12))
            m = np.sum(w * r)
            return np.sum(w * (r - m) ** 2)
        res = minimize_scalar(cost, bounds=(0.001, 0.5), method="bounded", options={"xatol": 1e-5})
        alpha = res.x
        geo = single_scatter(kind, alpha, cos_h[sel], cos_i[sel], cos_o[sel], F[:, None])[:, 0]
        m = np.sum(w * (meas - np.log(np.maximum(geo, 1e-12))))
        fits[kind] = (alpha, np.sqrt(res.fun), np.exp(m))
        say(f"[Measured:{name}] {kind}: fitted alpha {alpha:.4f} (roughness {np.sqrt(alpha):.3f}), residual std of log(measured/model) {np.sqrt(res.fun):.3f} (x{np.exp(np.sqrt(res.fun)):.2f}), scale measured/model at 550 nm {np.exp(m):.3f}")

    alpha = fits["GGX"][0]
    mus, E, E_avg = directional_albedo_table("GGX", alpha)
    F_avg = average_fresnel(eta, kk)
    say(f"[Measured:{name}] GGX alpha {alpha:.4f}: E_avg {E_avg:.4f} (single scattering loses {100 * (1 - E_avg):.2f} %)")
    Fall = fresnel_conductor(cos_d[sel][:, None], eta[None, :], kk[None, :])
    model_all = single_scatter("GGX", alpha, cos_h[sel], cos_i[sel], cos_o[sel], Fall) + ms_lobe(mus, E, E_avg, F_avg, cos_i[sel], cos_o[sel])
    scale = np.exp(np.sum(w[:, None] * np.log(fr[sel] / np.maximum(model_all, 1e-12)), axis=0))
    say(f"[Measured:{name}] scale measured / model over the fitted cells: {fmt(scale)}; relative to {wvls[i_ref]:.0f} nm: {fmt(scale / scale[i_ref])}")

    # the theta_d = 5 deg cut and the two-lobe fit
    i_d = 5
    peak_meas = np.nanmean(fr[:, i_d, :, :], axis=1)      # [theta_h, wavelength]
    wi_c, wo_c = directions(theta_h, np.full_like(theta_h, theta_d[i_d]), np.zeros_like(theta_h))
    cuts = {}
    for kind, (a_fit, _, _) in fits.items():
        Fc = fresnel_conductor(np.cos(theta_d[i_d]), eta[None, :], kk[None, :])
        cuts[kind] = single_scatter(kind, a_fit, np.cos(theta_h), wi_c[..., 2], wo_c[..., 2], Fc)
        if kind == "GGX":
            cuts[kind] = cuts[kind] + ms_lobe(mus, E, E_avg, F_avg, wi_c[..., 2], wo_c[..., 2])
    Fc550 = fresnel_conductor(np.cos(theta_d[i_d]), eta[i550], kk[i550])
    ok_c = np.isfinite(peak_meas[:, i550]) & (theta_h < np.radians(60))

    def two_lobe(p):
        a1, w1, a2, w2 = np.exp(p)
        return (w1 * single_scatter("GGX", a1, np.cos(theta_h), wi_c[..., 2], wo_c[..., 2], np.array([Fc550]))[:, 0]
                + w2 * single_scatter("GGX", a2, np.cos(theta_h), wi_c[..., 2], wo_c[..., 2], np.array([Fc550]))[:, 0])

    def resid(p):
        return np.log(np.maximum(two_lobe(p), 1e-12))[ok_c] - np.log(peak_meas[ok_c, i550])
    best = least_squares(resid, x0=np.log([0.01, 0.3, 0.1, 0.1]), bounds=(np.log([1e-3] * 4), np.log([1.0, 10.0, 1.0, 10.0])))
    a1, w1, a2, w2 = np.exp(best.x)
    if a1 > a2:
        a1, w1, a2, w2 = a2, w2, a1, w1
    two_cut = two_lobe(np.log([a1, w1, a2, w2]))
    say(f"[Measured:{name}] two GGX lobes on the theta_d = 5 deg cut (550 nm): spike alpha {a1:.4f} (roughness {np.sqrt(a1):.3f}) weight {w1:.3f}, "
        f"haze alpha {a2:.4f} (roughness {np.sqrt(a2):.3f}) weight {w2:.3f}; residual std {np.std(resid(best.x)):.3f} (x{np.exp(np.std(resid(best.x))):.2f}); weights are multiples of the Fresnel term")
    for deg in (0, 5, 10, 20, 40, 60):
        i_h = int(np.argmin(np.abs(np.degrees(theta_h) - deg)))
        say(f"[Measured:{name}] lobe at theta_h {np.degrees(theta_h[i_h]):5.1f} deg, theta_d 5 deg, 550 nm: measured {peak_meas[i_h, i550]:9.4f}, "
            f"GGX {cuts['GGX'][i_h, i550]:9.4f} (x{peak_meas[i_h, i550] / cuts['GGX'][i_h, i550]:.2f}), two lobes {two_cut[i_h]:9.4f} (x{peak_meas[i_h, i550] / two_cut[i_h]:.2f})")

    # 2. the color of the spike and the haze
    spike = finite & (TH < np.radians(0.5)) & (TD > np.radians(4)) & (TD < np.radians(40))
    haze = finite & (TH > np.radians(10)) & (TH < np.radians(35)) & (TD > np.radians(4)) & (TD < np.radians(40)) & (fr[..., i550] > 5 * floor[i550])
    spectra = {}
    for label, region in (("spike (theta_h < 0.5 deg)", spike), ("haze (theta_h 10-35 deg)", haze)):
        Fr = fresnel_conductor(cos_d[region][:, None], eta[None, :], kk[None, :])
        r = np.exp(np.mean(np.log(fr[region] / Fr), axis=0))
        spectra[label] = r
        say(f"[Measured:{name}] {label}: measured / Fresnel per wavelength, relative to {wvls[i_ref]:.0f} nm: {fmt(r / r[i_ref])}  ({region.sum()} cells)")

    # 3. directional albedo: fine cone around the mirror direction + coarse hemisphere
    def lookup(th_, td_, pd_):
        ih = np.clip(np.rint(90.0 * np.sqrt(np.clip(th_ / (np.pi / 2), 0, 1))).astype(int), 0, 90)
        id_ = np.clip(np.rint(np.degrees(td_)).astype(int), 0, 90)
        ip = np.clip(np.rint(np.degrees(pd_)).astype(int) + 180, 0, 360)
        return fr[ih, id_, ip, :]

    n_to, n_po = 360, 720
    to = (np.arange(n_to) + 0.5) / n_to * (np.pi / 2)
    po = (np.arange(n_po) + 0.5) / n_po * (2 * np.pi)
    TO, PO = np.meshgrid(to, po, indexing="ij")
    wo_g = np.stack([np.sin(TO) * np.cos(PO), np.sin(TO) * np.sin(PO), np.cos(TO)], -1)
    dw_g = np.sin(TO) * (np.pi / 2 / n_to) * (2 * np.pi / n_po)
    n_fine, cone = 1201, np.radians(3.0)
    u = np.linspace(-cone, cone, n_fine)
    U, V = np.meshgrid(u, u, indexing="ij")
    ang = np.sqrt(U ** 2 + V ** 2)
    inside = ang <= cone
    dw_f = (2 * cone / (n_fine - 1)) ** 2 * np.sinc(ang / np.pi)

    def model_on(wi_v, wo_):
        h = wi_v + wo_
        h /= np.linalg.norm(h, axis=-1, keepdims=True)
        Fg = fresnel_conductor(np.sum(wi_v * h, -1)[..., None], eta[None, None, :], kk[None, None, :])
        ci = np.full(wo_.shape[:-1], wi_v[2])
        s = single_scatter("GGX", alpha, h[..., 2], ci, wo_[..., 2], Fg)
        return s, s + ms_lobe(mus, E, E_avg, F_avg, ci, wo_[..., 2])

    say(f"[Measured:{name}] directional albedo; model = GGX alpha {alpha:.4f} of the single-lobe fit:")
    albedo_rows = []
    for deg in (15, 30, 45, 60, 75):
        ti = np.radians(deg)
        wi_v = np.array([np.sin(ti), 0.0, np.cos(ti)])
        mirror = np.array([-np.sin(ti), 0.0, np.cos(ti)])
        t1 = np.array([np.cos(ti), 0.0, np.sin(ti)])
        t2 = np.array([0.0, 1.0, 0.0])
        wo_f = mirror * np.cos(ang)[..., None] + (t1 * U[..., None] + t2 * V[..., None]) * np.sinc(ang / np.pi)[..., None]
        wo_f /= np.linalg.norm(wo_f, axis=-1, keepdims=True)
        vals = lookup(*half_diff(np.broadcast_to(wi_v, wo_f.shape), wo_f))
        wgt_f = (wo_f[..., 2] * dw_f)[..., None] * inside[..., None]
        E_cone = np.nansum(np.where(np.isfinite(vals), vals * wgt_f, 0.0), axis=(0, 1))
        outside = np.arccos(np.clip(wo_g @ mirror, -1, 1)) > cone
        vals_g = lookup(*half_diff(np.broadcast_to(wi_v, wo_g.shape), wo_g))
        wgt_g = (wo_g[..., 2] * dw_g)[..., None] * outside[..., None]
        E_out = np.nansum(np.where(np.isfinite(vals_g), vals_g * wgt_g, 0.0), axis=(0, 1))
        E_meas = E_cone + E_out
        s_f, f_f = model_on(wi_v, wo_f)
        s_g, f_g = model_on(wi_v, wo_g)
        E_single = np.sum(s_f * wgt_f, axis=(0, 1)) + np.sum(s_g * wgt_g, axis=(0, 1))
        E_full = np.sum(f_f * wgt_f, axis=(0, 1)) + np.sum(f_g * wgt_g, axis=(0, 1))
        Fi = fresnel_conductor(np.cos(ti), eta, kk)
        albedo_rows.append((deg, E_meas, E_single, E_full, Fi))
        say(f"     {deg:3d} deg   measured {fmt(E_meas)}   single+ms {fmt(E_full)}   measured/(single+ms) {fmt(E_meas / E_full, 2)}")

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("turbo")
    colors = [cmap((wvls[i] - 400) / 550) for i in range(n_lam)]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    deg_h = np.degrees(theta_h)
    deg_h[0] = 0.005
    for li in show:
        ax[0].plot(deg_h, peak_meas[:, li], "o", ms=3, color=colors[li], label=f"measured {int(round(wvls[li]))} nm")
    ax[0].plot(deg_h, cuts["GGX"][:, i550], "-", color="k", label=f"GGX + ms, alpha {fits['GGX'][0]:.3f} (550 nm)")
    ax[0].plot(deg_h, cuts["Beckmann"][:, i550], "--", color="gray", label=f"Beckmann, alpha {fits['Beckmann'][0]:.3f} (550 nm)")
    ax[0].plot(deg_h, two_cut, "-", color="m", label=f"two GGX lobes, alpha {a1:.3f} + {a2:.3f} (550 nm)")
    ax[0].set_yscale("log")
    ax[0].set_xscale("log")
    ax[0].set_xlim(0.005, 80)
    ax[0].set_ylim(1e-4, 300)
    ax[0].set_xlabel("theta_h [deg]  (theta_d = 5 deg)")
    ax[0].set_ylabel("BRDF [1/sr]")
    ax[0].set_title(f"Lobe shape: measured {name} vs fitted models")
    ax[0].legend(fontsize=7)
    ax[0].grid(alpha=0.3)

    F0 = fresnel_conductor(1.0, eta, kk)
    ax[1].plot(wvls, F0, "k-", label=f"Fresnel F0, {os.path.basename(metal_path)}")
    for (label, r), mk in zip(spectra.items(), ("r-", "b-")):
        ax[1].plot(wvls, r / r[i_ref] * F0[i_ref], mk, label=f"measured {label}, scaled to F0 at {wvls[i_ref]:.0f} nm")
    ax[1].set_ylim(0, 1.2)
    ax[1].set_xlabel("wavelength [nm]")
    ax[1].set_ylabel("reflectance")
    ax[1].set_title("Spectral shape of the spike and the haze vs Fresnel")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    degs = [r[0] for r in albedo_rows]
    for li in show:
        ax[2].plot(degs, [r[1][li] for r in albedo_rows], "o-", color=colors[li], label=f"measured {int(round(wvls[li]))} nm")
        ax[2].plot(degs, [r[3][li] for r in albedo_rows], "--", color=colors[li], alpha=0.7)
    ax[2].set_xlabel("theta_i [deg]")
    ax[2].set_ylabel("directional albedo")
    ax[2].set_ylim(0, 1.05)
    ax[2].set_title("Directional albedo: measured (solid), GGX + ms (dashed)", fontsize=9)
    ax[2].legend(fontsize=7)
    ax[2].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"measured_{name}.png"), dpi=130)
    with open(os.path.join(out_dir, "results.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote", os.path.join(out_dir, "results.txt"), f"and measured_{name}.png")


if __name__ == "__main__":
    main()
