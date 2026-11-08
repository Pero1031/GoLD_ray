"""
Mie phase function of a cloud-droplet population (water spheres in air) -> table for the renderer.

What it does
  1. Mie scattering of a sphere (Bohren & Huffman, BHMIE): amplitudes S1, S2 for every angle, Q_sca, Q_ext.
  2. A droplet size distribution (modified gamma, effective radius / effective variance), averaged over
     the scattering cross section.
  3. An average over the visible wavelengths (the engine's phase function does not depend on the wavelength,
     so the fogbow / glory colors are averaged out; the shape of the lobes stays).
  4. Writes  theta_deg, phase  (phase = density per steradian, integrates to 1 over the sphere), a log-scale
     comparison with a Henyey-Greenstein lobe of the same mean cosine, and the quantities a cloud needs
     (extinction per liquid water content).

Checks (run with --check): Rayleigh limit, energy conservation for a non-absorbing sphere, Q_ext -> 2 for large x,
the normalization, the mean cosine of cloud droplets (about 0.85), the rainbow / glory angles.

Usage:
  python tools/studies/mie_phase.py --check
  python tools/studies/mie_phase.py --reff 10 --veff 0.1 --out GoLD_rayt/assets/Cloud/mie_reff10.csv --plot out/mie_reff10.png
"""
import argparse
import os
import sys

import numpy as np

# Refractive index of water (real part; the imaginary part is ~1e-9 .. 1e-7 in the visible: no absorption)
WATER_N = {400: 1.3435, 450: 1.3403, 500: 1.3374, 550: 1.3351, 600: 1.3332, 650: 1.3317, 700: 1.3304}


def mie_s1s2(x, m, mu):
    """Bohren & Huffman BHMIE: amplitude functions S1, S2 at cos(theta) = mu (array), Q_sca, Q_ext.
    x = 2 pi r / lambda (size parameter), m = relative refractive index (complex)."""
    nstop = int(x + 4.0 * x ** (1.0 / 3.0) + 2.0)
    y = m * x
    nmx = int(max(nstop, abs(y)) + 15)

    # Logarithmic derivative D_n(mx) by downward recursion (stable)
    d = np.zeros(nmx + 1, dtype=complex)
    for n in range(nmx, 0, -1):
        d[n - 1] = n / y - 1.0 / (d[n] + n / y)

    psi0, psi1 = np.cos(x), np.sin(x)
    chi0, chi1 = -np.sin(x), np.cos(x)
    xi1 = psi1 - 1j * chi1

    pi0 = np.zeros_like(mu)
    pi1 = np.ones_like(mu)
    s1 = np.zeros(mu.shape, dtype=complex)
    s2 = np.zeros(mu.shape, dtype=complex)
    qsca = 0.0
    qext = 0.0

    for n in range(1, nstop + 1):
        fn = (2.0 * n + 1.0) / (n * (n + 1.0))
        psi = (2.0 * n - 1.0) * psi1 / x - psi0
        chi = (2.0 * n - 1.0) * chi1 / x - chi0
        xi = psi - 1j * chi

        an = ((d[n] / m + n / x) * psi - psi1) / ((d[n] / m + n / x) * xi - xi1)
        bn = ((m * d[n] + n / x) * psi - psi1) / ((m * d[n] + n / x) * xi - xi1)

        # tau_n from pi_n and pi_{n-1}
        tau1 = n * mu * pi1 - (n + 1.0) * pi0
        s1 += fn * (an * pi1 + bn * tau1)
        s2 += fn * (an * tau1 + bn * pi1)

        qsca += (2.0 * n + 1.0) * (abs(an) ** 2 + abs(bn) ** 2)
        qext += (2.0 * n + 1.0) * (an.real + bn.real)

        # advance to n + 1
        psi0, psi1 = psi1, psi
        chi0, chi1 = chi1, chi
        xi1 = psi1 - 1j * chi1
        pi_next = ((2.0 * n + 1.0) * mu * pi1 - (n + 1.0) * pi0) / n
        pi0, pi1 = pi1, pi_next

    return s1, s2, 2.0 / x ** 2 * qsca, 2.0 / x ** 2 * qext


def phase_of_sphere(x, m, theta):
    """Unpolarized phase function (per steradian, integrates to 1) of one sphere, plus Q_sca and Q_ext."""
    mu = np.cos(theta)
    s1, s2, qsca, qext = mie_s1s2(x, m, mu)
    inten = 0.5 * (np.abs(s1) ** 2 + np.abs(s2) ** 2)       # integrates to pi x^2 Q_sca over the sphere (C_sca k^2)
    return inten / (np.pi * x * x * qsca), qsca, qext


def gamma_distribution(reff, veff, n=48):
    """Modified gamma size distribution n(r) ~ r^((1-3b)/b) exp(-r / (a b)) (Hansen & Travis 1974), a = reff, b = veff.
    Returns radii (um) and weights (number fraction)."""
    a, b = reff, veff
    r = np.linspace(a * 0.15, a * (1.0 + 7.0 * np.sqrt(veff) + 1.0), n)
    nr = r ** ((1.0 - 3.0 * b) / b) * np.exp(-r / (a * b))
    return r, nr / nr.sum()


def cloud_phase(reff, veff, wavelengths, theta, nsizes=48):
    """Phase function of a droplet population, averaged over the cross section and over the wavelengths.
    Returns (phase per steradian, per-wavelength phases, sigma_ext per liquid water content in m^2/g)."""
    r, w = gamma_distribution(reff, veff, nsizes)
    per_wl = []
    ext_per_lwc = []
    for lam in wavelengths:
        n_w = np.interp(lam, list(WATER_N.keys()), list(WATER_N.values()))
        acc = np.zeros_like(theta)
        csca_total = 0.0
        cext_total = 0.0
        vol_total = 0.0
        for ri, wi in zip(r, w):
            x = 2.0 * np.pi * ri / (lam * 1e-3)
            p, qsca, qext = phase_of_sphere(x, complex(n_w, 0.0), theta)
            csca = np.pi * ri ** 2 * qsca
            acc += wi * csca * p
            csca_total += wi * csca
            cext_total += wi * np.pi * ri ** 2 * qext
            vol_total += wi * 4.0 / 3.0 * np.pi * ri ** 3
        per_wl.append(acc / csca_total)
        # m^2 per gram of liquid water: (cext / V) / 1 g/cm^3 -> um^2/um^3 = 1/um ; 1 g/cm^3 = 1e6 g/m^3 ; 1/um = 1e6 /m
        ext_per_lwc.append((cext_total / vol_total) * 1e6 / 1e6)
    per_wl = np.array(per_wl)
    return per_wl.mean(axis=0), per_wl, float(np.mean(ext_per_lwc))


def integrate_sphere(theta, p):
    return 2.0 * np.pi * np.trapezoid(p * np.sin(theta), theta)


def mean_cosine(theta, p):
    return 2.0 * np.pi * np.trapezoid(p * np.sin(theta) * np.cos(theta), theta)


def hg(cos_t, g):
    d = 1 + g * g - 2 * g * cos_t
    return (1 - g * g) / (4 * np.pi * d ** 1.5)


def run_checks():
    ok = True

    def report(name, cond, detail):
        nonlocal ok
        ok &= bool(cond)
        print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")

    theta = np.linspace(0, np.pi, 3601)

    # 1. Rayleigh limit (x << 1, m = 1.33): the phase function is 3/(16 pi) (1 + cos^2)
    p, qsca, qext = phase_of_sphere(0.01, complex(1.33, 0), theta)
    ref = 3.0 / (16.0 * np.pi) * (1.0 + np.cos(theta) ** 2)
    err = np.max(np.abs(p - ref)) / ref.max()
    report("Rayleigh limit (x = 0.01)", err < 1e-3, f"max deviation {err:.2e} of the peak")

    # 2. Energy conservation: no absorption -> Q_ext = Q_sca
    worst = 0.0
    for x in (0.5, 3.0, 20.0, 114.0, 400.0):
        _, qs, qe = phase_of_sphere(x, complex(1.333, 0), theta)
        worst = max(worst, abs(qs - qe) / qe)
    report("Q_ext = Q_sca without absorption", worst < 1e-6, f"worst relative difference {worst:.2e}")

    # 3. Large sphere: Q_ext -> 2 (extinction paradox)
    _, _, qe = phase_of_sphere(400.0, complex(1.333, 0), theta)
    report("Q_ext -> 2 for x = 400", abs(qe - 2.0) < 0.1, f"Q_ext = {qe:.3f}")

    # 4. Normalization of the phase function (integral over the sphere = 1)
    p, _, _ = phase_of_sphere(114.0, complex(1.333, 0), theta)
    nrm = integrate_sphere(theta, p)
    report("Normalization (x = 114)", abs(nrm - 1.0) < 2e-3, f"integral = {nrm:.5f}")

    # 5. Cloud droplets: mean cosine about 0.85 (literature: 0.85 .. 0.88 in the visible)
    pc, _, ext = cloud_phase(10.0, 0.1, [450, 550, 650], theta, nsizes=32)
    g = mean_cosine(theta, pc)
    report("Mean cosine of cloud droplets (reff 10 um)", 0.82 < g < 0.90, f"g = {g:.3f}")
    report("Extinction per liquid water (reff 10 um)", abs(ext - 0.15) < 0.02,
           f"{ext:.3f} m^2/g   (3 / (2 rho reff) = {3.0 / (2.0 * 1.0 * 10.0) * 1.0:.3f} m^2/g)")

    # 6. Rainbow (a primary bow near 138 deg for water) and glory (a peak at 180 deg) in one 20 um droplet
    p, _, _ = phase_of_sphere(2 * np.pi * 20.0 / 0.55, complex(1.333, 0), theta)
    deg = np.degrees(theta)
    window = (deg > 125) & (deg < 150)
    bow = deg[window][np.argmax(p[window])]
    report("Primary rainbow of one droplet", 135.0 < bow < 142.0, f"maximum at {bow:.1f} deg")
    glory = p[-1] > 3.0 * np.median(p[(deg > 170) & (deg < 178)]) * 0.5 and p[-1] > p[(deg > 150) & (deg < 160)].mean()
    report("Backward peak (glory) of one droplet", glory, f"p(180) / p(150..160 mean) = {p[-1] / p[(deg > 150) & (deg < 160)].mean():.2f}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="run the self checks and exit")
    ap.add_argument("--reff", type=float, default=10.0, help="effective radius, um")
    ap.add_argument("--veff", type=float, default=0.1, help="effective variance")
    ap.add_argument("--step", type=float, default=0.05, help="angle step, degrees")
    ap.add_argument("--wavelengths", type=str, default="420,470,520,570,620,670", help="nm, comma separated")
    ap.add_argument("--out", type=str, default="")
    ap.add_argument("--plot", type=str, default="")
    args = ap.parse_args()

    if args.check:
        sys.exit(0 if run_checks() else 1)

    wls = [float(s) for s in args.wavelengths.split(",")]
    theta = np.radians(np.arange(0.0, 180.0 + 0.5 * args.step, args.step))
    p, per_wl, ext = cloud_phase(args.reff, args.veff, wls, theta)
    p = p / integrate_sphere(theta, p)
    g = mean_cosine(theta, p)
    print(f"reff = {args.reff} um, veff = {args.veff}, {len(theta)} angles")
    print(f"mean cosine g = {g:.4f}; extinction per liquid water = {ext:.4f} m^2/g")
    print(f"p(0) = {p[0]:.1f}, p(90 deg) = {p[len(p) // 2]:.4f}, p(180 deg) = {p[-1]:.4f}  (per steradian; isotropic = {1 / (4 * np.pi):.4f})")

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w", newline="\n") as f:
            f.write(f"# Mie phase function of cloud droplets: modified gamma distribution, reff = {args.reff} um, veff = {args.veff},\n")
            f.write(f"# water n(lambda), averaged over {args.wavelengths} nm. Generated by tools/studies/mie_phase.py.\n")
            f.write(f"# mean cosine = {g:.4f}; extinction per liquid water content = {ext:.4f} m^2/g\n")
            f.write("# theta_deg, phase_per_steradian (integrates to 1 over the sphere)\n")
            for t, v in zip(np.degrees(theta), p):
                f.write(f"{t:.4f}, {v:.8e}\n")
        print("wrote", args.out)

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        deg = np.degrees(theta)
        fig, ax = plt.subplots(figsize=(7, 4.2))
        ax.semilogy(deg, p * 4 * np.pi, label=f"Mie (r_eff {args.reff:g} um), g = {g:.3f}")
        ax.semilogy(deg, hg(np.cos(theta), g) * 4 * np.pi, "--", label=f"Henyey-Greenstein, g = {g:.3f}")
        ax.semilogy(deg, (0.7 * hg(np.cos(theta), 0.95) + 0.3 * hg(np.cos(theta), -0.3)) * 4 * np.pi, ":", label="two lobes (0.95 / -0.3)")
        ax.set_xlabel("scattering angle (deg)")
        ax.set_ylabel("4 pi p(theta)")
        ax.set_xlim(0, 180)
        ax.legend()
        fig.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(args.plot)), exist_ok=True)
        fig.savefig(args.plot, dpi=130)
        print("wrote", args.plot)


if __name__ == "__main__":
    main()
