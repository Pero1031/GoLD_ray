"""A random-walk reference for the radial profile of the subsurface scattering (docs/SSS.md, section 9.2).

The engine's BSSRDF uses Burley's normalized diffusion profile  P(r; d) = (exp(-r/d) + exp(-r/(3d))) / (8 pi d r)  scaled by the diffuse color A, with
d = l / s(A),  s(A) = 1.85 - A + 7 |A - 0.8|^3  (l: the scattering distance). This script measures what a real random walk gives, so that the error of that
approximation can be stated:

  * a semi-infinite homogeneous medium (extinction 1: lengths are in mean free paths), single-scattering albedo alpha, Henyey-Greenstein anisotropy g;
  * an index-matched boundary (the engine treats the Fresnel factors outside the profile), a pencil beam at the origin, normal incidence;
  * every photon is followed until it leaves the surface or is absorbed (analog: absorbed with probability 1 - alpha at each collision);
    where it leaves, its distance r from the origin goes into a histogram: R(r) 2 pi r dr (per incident photon).

For each case it prints the total diffuse reflectance A, then compares the measured radial distribution with the engine's profile in three ways:
  1. the engine's own mapping:  A from the walk,  l = 1 / (1 - g alpha)  (the reduced mean free path: the engine's `--sss` takes the scattering distance
     as an input; for an isotropic medium it is 1),  d = l / s(A);
  2. the best fit of d alone (a scale of the radius): how much the shape is wrong even when the scale is right;
  3. the mean radius (the engine's sampler has 2.5 d).
The errors are given as the Kolmogorov distance (the largest difference of the radial cumulative distributions, 0 = identical, 1 = disjoint) and as
ratios of the profile at a few radii (walk / profile).

Run:  python tools/studies/sss_random_walk.py [--walkers N] [--seed S] [--out FILE] [--fig FILE]
"""
import argparse
import time

import numpy as np


def s_of_A(A):
    return 1.85 - A + 7.0 * np.abs(A - 0.8) ** 3


def burley_cdf(r, d):
    return 1.0 - 0.25 * np.exp(-r / d) - 0.75 * np.exp(-r / (3.0 * d))


def burley_profile(r, d, A=1.0):
    return A * (np.exp(-r / d) + np.exp(-r / (3.0 * d))) / (8.0 * np.pi * d * np.maximum(r, 1e-12))


def sample_hg_cos(g, xi):
    if abs(g) < 1e-3:
        return 1.0 - 2.0 * xi
    s = (1.0 - g * g) / (1.0 - g + 2.0 * g * xi)
    return np.clip((1.0 + g * g - s * s) / (2.0 * g), -1.0, 1.0)


def random_walk(alpha, g, n, rng, edges, max_steps=400000, chunk=500000):
    """Histogram (weights per radial bin) of the photons that leave the surface, and the number of photons that never finished."""
    hist = np.zeros(len(edges) - 1)
    unfinished = 0
    done = 0
    while done < n:
        m = min(chunk, n - done)
        done += m
        x = np.zeros(m); y = np.zeros(m); z = np.zeros(m)
        u = np.zeros(m); v = np.zeros(m); w = np.ones(m)
        for step in range(max_steps):
            if x.size == 0:
                break
            s = -np.log1p(-rng.random(x.size))
            z = z + w * s
            esc = z < 0.0
            if esc.any():
                t0 = -(z[esc] - w[esc] * s[esc]) / w[esc]    # the distance from the last collision to the plane z = 0
                xe = x[esc] + u[esc] * t0
                ye = y[esc] + v[esc] * t0
                hist += np.histogram(np.hypot(xe, ye), bins=edges)[0]
            keep = ~esc
            x = x[keep] + u[keep] * s[keep]; y = y[keep] + v[keep] * s[keep]; z = z[keep]
            u = u[keep]; v = v[keep]; w = w[keep]
            alive = rng.random(x.size) < alpha        # absorbed with probability 1 - alpha
            x = x[alive]; y = y[alive]; z = z[alive]; u = u[alive]; v = v[alive]; w = w[alive]
            if x.size == 0:
                break
            cost = sample_hg_cos(g, rng.random(x.size))
            sint = np.sqrt(np.maximum(0.0, 1.0 - cost * cost))
            phi = 2.0 * np.pi * rng.random(x.size)
            cp, sp = np.cos(phi), np.sin(phi)
            near = np.abs(w) > 0.99999
            den = np.sqrt(np.maximum(1e-30, 1.0 - w * w))
            un = np.where(near, sint * cp, sint * (u * w * cp - v * sp) / den + u * cost)
            vn = np.where(near, sint * sp, sint * (v * w * cp + u * sp) / den + v * cost)
            wn = np.where(near, np.sign(w) * cost, -sint * cp * den + w * cost)
            nrm = np.sqrt(un * un + vn * vn + wn * wn)
            u, v, w = un / nrm, vn / nrm, wn / nrm
        else:
            pass
        unfinished += x.size
    return hist, unfinished


def analyse(hist, n, edges, g, alpha):
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    prob = hist / n                                   # the fraction of the incident photons per annulus
    A = prob.sum()
    cdf = np.cumsum(prob) / A                         # at the upper edge of each bin
    r_hi = edges[1:]
    centers = np.sqrt(edges[1:] * edges[:-1])
    ell = 1.0 / (1.0 - g * alpha)                     # the reduced mean free path
    d_map = ell / s_of_A(A)
    ks_map = float(np.max(np.abs(cdf - burley_cdf(r_hi, d_map))))
    scales = np.exp(np.linspace(np.log(0.2), np.log(5.0), 400))
    ks_fit = [float(np.max(np.abs(cdf - burley_cdf(r_hi, d_map * k)))) for k in scales]
    k_best = float(scales[int(np.argmin(ks_fit))])
    mean_r = float(np.sum(centers * prob) / A)
    ratios = {}
    for r0 in (0.1, 0.3, 1.0, 3.0, 10.0):                # in units of ell
        sel = (centers > 0.7 * r0 * ell) & (centers < 1.4 * r0 * ell)
        walk = prob[sel].sum() / area[sel].sum()
        model = (A * (burley_cdf(edges[1:][sel], d_map) - burley_cdf(edges[:-1][sel], d_map))).sum() / area[sel].sum()
        ratios[r0] = walk / model
    return dict(A=A, d_map=d_map, ell=ell, ks_map=ks_map, k_best=k_best, ks_fit=min(ks_fit), mean_r=mean_r,
                mean_model=2.5 * d_map, ratios=ratios, centers=centers, area=area, prob=prob, cdf=cdf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--walkers", type=int, default=1000000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=None, help="write the table to this file as well")
    ap.add_argument("--fig", default=None, help="save a figure (needs matplotlib)")
    args = ap.parse_args()

    edges = np.concatenate([[0.0], np.logspace(-3, 3.7, 90)])     # in mean free paths
    cases = [(0.5, 0.0), (0.8, 0.0), (0.9, 0.0), (0.95, 0.0), (0.99, 0.0), (0.999, 0.0), (0.9999, 0.0),
             (0.99, 0.5), (0.99, 0.8), (0.999, 0.8)]
    rng = np.random.default_rng(args.seed)
    lines = []
    head = ("alpha     g    |   A      s(A)  d_map  | KS(map)  d_fit/d_map  KS(fit) | mean r: walk  model(2.5d) | walk/model at r = 0.1 0.3 1 3 10 (x l)")
    lines.append(head)
    results = []
    t0 = time.time()
    for alpha, g in cases:
        hist, unfinished = random_walk(alpha, g, args.walkers, rng, edges)
        R = analyse(hist, args.walkers, edges, g, alpha)
        results.append((alpha, g, R))
        rr = " ".join("%5.2f" % R["ratios"][k] for k in (0.1, 0.3, 1.0, 3.0, 10.0))
        lines.append("%-8g %4.1f  | %6.3f %6.2f %6.2f | %6.3f   %6.2f     %6.3f  | %9.2f  %9.2f      | %s%s" % (
            alpha, g, R["A"], s_of_A(R["A"]), R["d_map"], R["ks_map"], R["k_best"], R["ks_fit"], R["mean_r"], R["mean_model"], rr,
            ("   (%d photons unfinished)" % unfinished) if unfinished else ""))
        print(lines[-1], "   [%.0f s]" % (time.time() - t0), flush=True)
    text = "\n".join(lines)
    print("\n" + text)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write("# random walk (%d walkers per case, seed %d) vs Burley normalized diffusion; lengths in mean free paths\n" % (args.walkers, args.seed))
            f.write(text + "\n")
    if args.fig:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sel = [(0.9, 0.0), (0.99, 0.0), (0.9999, 0.0)]
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
        for a, (alpha, g, R) in zip(ax, [r for r in results if (r[0], r[1]) in sel]):
            rc = R["centers"]
            walk = R["prob"] / R["area"]
            a.loglog(rc, walk, "k.", ms=4, label="random walk")
            a.loglog(rc, burley_profile(rc, R["d_map"], R["A"]), "r-", label="engine mapping  d = l / s(A)")
            a.loglog(rc, burley_profile(rc, R["d_map"] * R["k_best"], R["A"]), "b--", label="best-fit d (x%.2f)" % R["k_best"])
            a.set_title("albedo %g,  A = %.3f" % (alpha, R["A"]))
            a.set_xlabel("radius r (mean free paths)")
            a.set_ylim(max(walk.max() * 1e-9, 1e-14), walk.max() * 3)
        ax[0].set_ylabel("R(r)  (per area, per incident photon)")
        ax[0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(args.fig, dpi=130)
        print("figure:", args.fig)


if __name__ == "__main__":
    main()
