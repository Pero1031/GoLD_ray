"""
The RIT bispectral (fluorescence) database for the fluorescent material (docs/FLUORESCENCE.md).

Gonzalez and Fairchild measured 144 materials (the files at the top of the folder: mostly prints, papers, inks and fluorescent
pigments; the bcraday*/ folders hold repeated measurements of BCRA ceramic reference tiles) with a Labsphere BFC-450
bispectral colorimeter ("Evaluation of bispectral spectrophotometry for accurate colorimetry of printing materials", Color Imaging
Conference 2000) and offered the matrix files on the Munsell Color Science Laboratory site. No license is stated, so the files are
not part of this repository: this script downloads them to GoLD_rayt/assets/Fluorescence/rit_bispectral/ (git-ignored).

  python tools/data/rit_bispectral.py fetch          download Labsphere.zip (6.8 MB) and unpack it
  python tools/data/rit_bispectral.py list [N]       the N most fluorescent samples (default 40), then the rest in the folder

Each .BFC file is text: rows = emission wavelengths 380..780 nm, columns = excitation wavelengths 300..780 nm, 10 nm steps; the
diagonal is the reflectance, the rest the fluorescence (radiance factor per 10 nm excitation band).
"""
import os
import sys
import urllib.request
import zipfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEST = os.path.join(ROOT, "GoLD_rayt", "assets", "Fluorescence", "rit_bispectral")
URL = "http://www.rit-mcsl.org/StudentResearch/sgonzalez/Labsphere.zip"


def fetch():
    os.makedirs(DEST, exist_ok=True)
    zpath = os.path.join(DEST, "Labsphere.zip")
    if not os.path.exists(zpath):
        print("downloading", URL)
        urllib.request.urlretrieve(URL, zpath)
    with zipfile.ZipFile(zpath) as z:
        bad = [n for n in z.namelist() if n.startswith(("/", "\\")) or ".." in n.replace("\\", "/").split("/")]
        if bad:
            sys.exit(f"refusing to unpack: unsafe names {bad[:3]}")
        z.extractall(DEST)
    print("unpacked to", DEST)


def load(path):
    lines = open(path, encoding="latin-1").read().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("r:c:"))
    ex = np.array([float(x) for x in lines[i].split()[1:]])
    rows = [l.split() for l in lines[i + 1:] if l.strip() and l.split()[0].replace(".", "").isdigit()]
    em = np.array([float(r[0]) for r in rows])
    m = np.array([[float(x) for x in r[1:]] for r in rows])
    title = lines[3].lstrip(";").strip() if len(lines) > 3 else ""
    return em, ex, np.clip(m, 0, None), title


def listing(n):
    base = os.path.join(DEST, "Labsphere")
    if not os.path.isdir(base):
        sys.exit("not found: run 'python tools/data/rit_bispectral.py fetch' first")
    seen, rows = set(), []
    for dirpath, _, files in os.walk(base):
        for f in files:
            if not f.lower().endswith(".bfc"):
                continue
            path = os.path.join(dirpath, f)
            key = os.path.relpath(path, base).lower()
            if key in seen:
                continue
            seen.add(key)
            em, ex, m, title = load(path)
            off = ex[None, :] - em[:, None]
            diag = [m[r, list(ex).index(e)] for r, e in enumerate(em) if e in ex]
            rows.append((float(m[off < 0].sum()), os.path.relpath(path, base).replace("\\", "/"), title, float(np.mean(diag))))
    rows.sort(reverse=True)
    print(f"{len(rows)} matrices in {base}; fluorescence = sum of the entries with excitation below emission")
    print("(IXCAXORA.BFC is wrong: Hua et al. 2023 found it a measurement error; it tops the list for that reason)")
    for F, f, t, R in rows[:n]:
        print(f"  {f:28s} fluorescence {F:6.2f}   mean reflectance {R:.2f}   {t}")


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "fetch":
        fetch()
    elif len(sys.argv) >= 2 and sys.argv[1] == "list":
        listing(int(sys.argv[2]) if len(sys.argv) > 2 else 40)
    else:
        print(__doc__)
