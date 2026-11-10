"""Directional albedo of a KAIST pBRDF sample from its unpolarized table, by the same integration as measured_gold.py.

    python tools/studies/pbrdf_albedo.py references/kaist_pbrdf/<sample>_unpolarized.npz [theta_i deg ...]

The table fr[phi_d, theta_d, phi_h ...] is read as written by references/kaist_pbrdf/extract_unpolarized.py. For each
incidence angle the outgoing hemisphere is integrated with a fine exponential-map grid (0.005 deg) within 3 deg of the
mirror direction and a 0.25 deg grid elsewhere; cells that were not measured (NaN) are skipped, and the cosine-weighted
solid angle they cover is reported. For a white diffuse standard (Spectralon, reflectance about 0.99) the albedo should
come out near 0.99 and the BRDF near 0.99 / pi = 0.315 1/sr if the table is an absolute BRDF.
"""
import sys
import numpy as np

sys.path.insert(0, __file__.rsplit("tools", 1)[0] + "tools/studies")
from measured_gold import half_diff   # noqa: E402  (the MERL half-difference convention)


def main():
    path = sys.argv[1]
    degs = [float(a) for a in sys.argv[2:]] or [15, 30, 45, 60, 75]
    d = np.load(path)
    fr = np.transpose(d["fr"], (2, 1, 0, 3)).astype(np.float64)   # [theta_h, theta_d, phi_d, wavelength]
    wvls = d["wavelengths_nm"].astype(int).tolist()
    print(f"[Albedo] {path}: wavelengths {wvls} nm, finite {np.isfinite(fr).mean() * 100:.1f} %, "
          f"median of finite values {np.nanmedian(fr):.4f} 1/sr (0.99 / pi = {0.99 / np.pi:.4f})")

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

    print("     theta_i   albedo per wavelength                      within 3 deg of mirror           NaN share of the hemisphere (cos-weighted)")
    for deg in degs:
        ti = np.radians(deg)
        wi_v = np.array([np.sin(ti), 0.0, np.cos(ti)])
        mirror = np.array([-np.sin(ti), 0.0, np.cos(ti)])
        t1 = np.array([np.cos(ti), 0.0, np.sin(ti)])
        t2 = np.array([0.0, 1.0, 0.0])
        wo_f = mirror * np.cos(ang)[..., None] + (t1 * U[..., None] + t2 * V[..., None]) * np.sinc(ang / np.pi)[..., None]
        wo_f /= np.linalg.norm(wo_f, axis=-1, keepdims=True)
        vals = lookup(*half_diff(np.broadcast_to(wi_v, wo_f.shape), wo_f))
        wgt_f = (wo_f[..., 2] * dw_f) * inside
        E_cone = np.nansum(np.where(np.isfinite(vals), vals * wgt_f[..., None], 0.0), axis=(0, 1))
        outside = np.arccos(np.clip(wo_g @ mirror, -1, 1)) > cone
        vals_g = lookup(*half_diff(np.broadcast_to(wi_v, wo_g.shape), wo_g))
        wgt_g = (wo_g[..., 2] * dw_g) * outside
        E_out = np.nansum(np.where(np.isfinite(vals_g), vals_g * wgt_g[..., None], 0.0), axis=(0, 1))
        nan_share = (np.sum(np.where(~np.isfinite(vals_g[..., 0]), wgt_g, 0.0)) + np.sum(np.where(~np.isfinite(vals[..., 0]), wgt_f, 0.0))) / np.pi
        print(f"     {deg:5.1f} deg  {np.round(E_cone + E_out, 3).tolist()}   {np.round(E_cone, 3).tolist()}   {nan_share * 100:5.1f} %")


if __name__ == "__main__":
    main()
