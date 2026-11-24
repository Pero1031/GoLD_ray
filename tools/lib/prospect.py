"""
PROSPECT in numpy (the leaf material's model, docs/LEAF.md), for studies and as an independent check of Materials/Prospect.cpp.

    table = load_table()                       # GoLD_rayt/assets/Leaf/prospect/dataSpec_PRO_v2.txt (tools/data/leaf_optics.py fetch)
    R, T = prospect(table, N=1.5, chl=40, car=8, ant=0, brown=0, ewt=0.01, lma=0.009, alpha=40)

E1 comes from scipy here (the engine has its own series / continued fraction).
"""
import os

import numpy as np
from scipy.special import exp1

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TABLE = os.path.join(ROOT, "GoLD_rayt", "assets", "Leaf", "prospect", "dataSpec_PRO_v2.txt")
ANGERS = os.path.join(ROOT, "GoLD_rayt", "assets", "Leaf", "angers")


def load_table(path=TABLE):
    with open(path) as f:
        head = f.readline().split()
    data = np.loadtxt(path, skiprows=1)
    return {name: data[:, i] for i, name in enumerate(head)}


def tav(alpha, n):
    """Transmissivity of a flat face of index n for isotropic light within alpha degrees of the normal (Stern 1964)."""
    n = np.asarray(n, dtype=float)
    n2 = n * n
    np_ = n2 + 1
    nm = n2 - 1
    a = (n + 1) ** 2 / 2
    k = -((n2 - 1) ** 2) / 4
    sa = np.sin(np.radians(alpha))
    b2 = sa * sa - np_ / 2
    b1 = 0 * b2 if alpha == 90 else np.sqrt(b2 * b2 + k)
    b = b1 - b2
    ts = (k * k / (6 * b ** 3) + k / b - b / 2) - (k * k / (6 * a ** 3) + k / a - a / 2)
    tp1 = -2 * n2 * (b - a) / np_ ** 2
    tp2 = -2 * n2 * np_ * np.log(b / a) / nm ** 2
    tp3 = n2 * (1 / b - 1 / a) / 2
    tp4 = 16 * n2 ** 2 * (n2 ** 2 + 1) * np.log((2 * np_ * b - nm ** 2) / (2 * np_ * a - nm ** 2)) / (np_ ** 3 * nm ** 2)
    tp5 = 16 * n2 ** 3 * (1 / (2 * np_ * b - nm ** 2) - 1 / (2 * np_ * a - nm ** 2)) / np_ ** 3
    return (ts + tp1 + tp2 + tp3 + tp4 + tp5) / (2 * sa * sa)


def prospect(t, N=1.5, chl=40.0, car=8.0, ant=0.0, brown=0.0, ewt=0.01, lma=0.009, alpha=40.0):
    """Reflectance and transmittance at the table's wavelengths (t["lambda"])."""
    n = t["nrefrac"]
    k = (chl * t["sac_chl"] + car * t["sac_car"] + ant * t["sac_ant"] + brown * t["sac_brown"]
         + ewt * t["sac_ewt"] + lma * t["sac_lma"]) / N
    with np.errstate(divide="ignore", invalid="ignore"):
        tau = np.where(k > 0, (1 - k) * np.exp(-k) + k * k * exp1(np.maximum(k, 1e-300)), 1.0)
    talf = tav(alpha, n)
    ralf = 1 - talf
    t12 = tav(90, n)
    r12 = 1 - t12
    t21 = t12 / n ** 2
    r21 = 1 - t21
    denom = 1 - r21 * r21 * tau * tau
    Ta = talf * tau * t21 / denom
    Ra = ralf + r21 * tau * Ta
    tt = t12 * tau * t21 / denom
    r = r12 + r21 * tau * tt
    D = np.sqrt((1 + r + tt) * (1 + r - tt) * (1 - r + tt) * (1 - r - tt))
    a = (1 + r * r - tt * tt + D) / (2 * r)
    b = (1 - r * r + tt * tt + D) / (2 * tt)
    bNm1 = b ** (N - 1)
    bN2 = bNm1 ** 2
    den = a * a * bN2 - 1
    Rsub = a * (bN2 - 1) / den
    Tsub = bNm1 * (a * a - 1) / den
    zero = r + tt >= 1
    Tsub = np.where(zero, tt / (tt + (1 - tt) * (N - 1)), Tsub)
    Rsub = np.where(zero, 1 - Tsub, Rsub)
    denom = 1 - Rsub * r
    return Ra + Ta * Rsub * tt / denom, Ta * Tsub / denom


def load_angers(path=ANGERS):
    """Measured pigments (dict of arrays: chl = a + b, car, ant, ewt, lma) and spectra (lambda, R[sample, lambda], T[...])."""
    with open(os.path.join(path, "DataBioch.txt")) as f:
        head = f.readline().split()
    b = np.loadtxt(os.path.join(path, "DataBioch.txt"), skiprows=1)
    col = {h: b[:, i] for i, h in enumerate(head)}
    chem = {"chl": col["CHLa"] + col["CHLb"], "car": col["CAR"], "ant": col["ANT_estimated"], "ewt": col["EWT"], "lma": col["LMA"]}
    R = np.loadtxt(os.path.join(path, "ReflectanceData.txt"), skiprows=1)
    T = np.loadtxt(os.path.join(path, "TransmittanceData.txt"), skiprows=1)
    return chem, R[:, 0], R[:, 1:].T, T[:, 1:].T
