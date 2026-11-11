"""Converts the lab's fading-experiment reflectance tables ("<color>_refangle.csv") into the two-column
"wavelength_nm, reflectance" CSVs that GoLD_rayt reads ("type": "diffuse", "reflectance": "...csv").

Input layout (the robot-arm + SR-5000 measurements, 2023-2024):
  <series>/<N>回目.../<color>_refangle.csv          or
  <series>/<N>回目.../反射率/<color>_refangle.csv
One file holds ~26 rows: data_num, in_theta, in_phi, ref_theta, reh_phi, 380nm ... 780nm. The rows are the
measurement spots of one patch (the angles differ by about a degree), and the values are already the ratio
to the white reference measured in the same geometry. This script takes the median over the rows per
wavelength (robust against the few spots where the white was near zero) and writes one curve per color
and round, plus an index with the geometry and the spread across the spots.

Usage:
  python tools/data/fading_refangle.py <series dir> --out <output dir> [--name pigment]
  python tools/data/fading_refangle.py --all <sample root> --out <output dir>
"""
import argparse
import csv
import glob
import os
import re
import statistics
import sys

ROUND_RE = re.compile(r"(\d+)回目")
H_RE = re.compile(r"h\s*=\s*(\d+)")

SERIES = {  # --all: sub-folder of the lab's "sample" folder -> output name
    "顔料まとめ": "pigment",
    os.path.join("染料まとめ", "画用紙"): "dye_paper",
    os.path.join("染料まとめ", "折り紙"): "dye_origami",
}


HOURS_STEP = {"pigment": 150.0, "dye_paper": 30.0}  # hours of pseudo-sunlight between two rounds (Kataoka's thesis, chapter 2); round 1 = before exposure


def round_hours(name, rounds):
    """{round: exposure hours}: a fixed step per series, or the "(h=30)" in the folder names (the first round is then 0 h)."""
    if name in HOURS_STEP:
        return {rnd: (rnd - 1) * HOURS_STEP[name] for rnd, _, _, _ in rounds}
    named = {rnd: h for rnd, _, _, h in rounds if h is not None}
    if not named:
        return {}
    return {rnd: float(named.get(rnd, 0 if rnd == rounds[0][0] else float("nan"))) for rnd, _, _, _ in rounds}


def write_hours(out_dir, name, rounds):
    hours = round_hours(name, rounds)
    if not hours or any(h != h for h in hours.values()):
        print(f"{name}: exposure hours unknown, no hours.csv", file=sys.stderr)
        return
    with open(os.path.join(out_dir, name, "hours.csv"), "w", encoding="utf-8", newline="\n") as f:
        f.write("# exposure hours of each round (the time axis of \"fading\" in a scene file)\n")
        f.write("round,hours\n")
        for rnd in sorted(hours):
            f.write(f"{rnd},{hours[rnd]:g}\n")


def read_refangle(path):
    with open(path, encoding="utf-8", errors="ignore", newline="") as f:
        rows = [r for r in csv.reader(f) if r]
    header = rows[0]
    wl_cols = [(i, float(h.replace("nm", ""))) for i, h in enumerate(header) if h.strip().endswith("nm")]
    data = []
    for r in rows[1:]:
        if len(r) < len(header) or r[0].strip() == "":
            continue
        try:
            angles = [float(x) for x in r[1:5]]
            spectrum = [float(r[i]) if r[i].strip() != "" else float("nan") for i, _ in wl_cols]
        except ValueError:
            continue
        data.append((angles, spectrum))
    return [w for _, w in wl_cols], data


def summarize(wavelengths, data):
    """Median and standard deviation over the spots, per wavelength."""
    med, sd = [], []
    for k in range(len(wavelengths)):
        v = [s[k] for _, s in data if s[k] == s[k]]  # drop NaN
        if not v:
            med.append(float("nan")); sd.append(float("nan")); continue
        med.append(statistics.median(v))
        sd.append(statistics.pstdev(v) if len(v) > 1 else 0.0)
    return med, sd


def angle_summary(data):
    cols = list(zip(*[a for a, _ in data]))
    return [(min(c), max(c)) for c in cols]  # in_theta, in_phi, ref_theta, ref_phi


def find_rounds(series_dir):
    out = []
    for d in sorted(os.listdir(series_dir)):
        full = os.path.join(series_dir, d)
        if not os.path.isdir(full):
            continue
        m = ROUND_RE.search(d)
        if not m:
            continue
        files = glob.glob(os.path.join(full, "*_refangle.csv")) or glob.glob(os.path.join(full, "反射率", "*_refangle.csv"))
        if not files:
            continue
        h = H_RE.search(d)
        out.append((int(m.group(1)), d, sorted(files), int(h.group(1)) if h else None))
    return sorted(out)


def convert_series(series_dir, out_dir, name, index_rows, wl_range=(440.0, 680.0)):
    rounds = find_rounds(series_dir)
    if not rounds:
        print(f"{series_dir}: no rounds found", file=sys.stderr)
        return 0
    n = 0
    for rnd, folder, files, h in rounds:
        for path in files:
            color = os.path.basename(path)[: -len("_refangle.csv")]
            wavelengths, data = read_refangle(path)
            if not data:
                print(f"skip {path}: no rows", file=sys.stderr)
                continue
            med, sd = summarize(wavelengths, data)
            ang = angle_summary(data)
            lo, hi = wl_range
            keep = [i for i, w in enumerate(wavelengths) if lo <= w <= hi]
            wavelengths, med, sd = [wavelengths[i] for i in keep], [med[i] for i in keep], [sd[i] for i in keep]
            dest = os.path.join(out_dir, name, f"r{rnd}", f"{color}.csv")
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            above = sum(1 for v in med if v == v and v > 1.0)
            with open(dest, "w", encoding="utf-8", newline="\n") as out:
                out.write(f"# {name} / round {rnd}" + (f" (h={h})" if h is not None else "") + f" / {color}\n")
                out.write(f"# source: {path}\n")
                out.write(f"# median over {len(data)} spots; reflectance factor = sample / white in the same geometry\n")
                out.write(f"# kept {lo:g}..{hi:g} nm (outside, the lamp is too weak: the ratio is zeros and spikes); "
                          f"the engine extends the end values beyond this range\n")
                out.write(f"# incidence theta {ang[0][0]:.1f}..{ang[0][1]:.1f} deg, phi {ang[1][0]:.1f}..{ang[1][1]:.1f}; "
                          f"viewing theta {ang[2][0]:.1f}..{ang[2][1]:.1f}, phi {ang[3][0]:.1f}..{ang[3][1]:.1f}\n")
                if above:
                    out.write(f"# note: {above} wavelengths above 1.0 (left as measured)\n")
                out.write("wavelength_nm,reflectance\n")
                for w, v in zip(wavelengths, med):
                    if v == v:
                        out.write(f"{w:g},{v:.6g}\n")
            k550 = min(range(len(wavelengths)), key=lambda i: abs(wavelengths[i] - 550))
            index_rows.append({
                "series": name, "round": rnd, "h": "" if h is None else h, "color": color,
                "spots": len(data), "R550_median": f"{med[k550]:.4f}", "R550_sd": f"{sd[k550]:.4f}",
                "in_theta": f"{ang[0][0]:.1f}-{ang[0][1]:.1f}", "ref_theta": f"{ang[2][0]:.1f}-{ang[2][1]:.1f}",
                "above_1": above, "file": os.path.relpath(dest, out_dir).replace(os.sep, "/"), "source_folder": folder,
            })
            n += 1
    write_hours(out_dir, name, rounds)
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("series", nargs="?", help="folder holding the <N>回目 sub-folders")
    ap.add_argument("--all", metavar="SAMPLE_ROOT", help="convert the known series under the lab's sample folder")
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", help="output sub-folder name (default: the series folder's name)")
    ap.add_argument("--range", nargs=2, type=float, default=(440.0, 680.0), metavar=("LO_NM", "HI_NM"),
                    help="wavelengths to keep (default 440 680: outside, the lamp is weak and the ratio to the white turns into zeros and spikes; 380 780 keeps everything)")
    a = ap.parse_args(argv)
    index_rows, total = [], 0
    wl_range = tuple(a.range)
    if a.all:
        for sub, name in SERIES.items():
            d = os.path.join(a.all, sub)
            if os.path.isdir(d):
                total += convert_series(d, a.out, name, index_rows, wl_range)
            else:
                print(f"missing: {d}", file=sys.stderr)
    elif a.series:
        total += convert_series(a.series, a.out, a.name or os.path.basename(os.path.normpath(a.series)), index_rows, wl_range)
    else:
        ap.error("give a series folder or --all")
    if index_rows:
        os.makedirs(a.out, exist_ok=True)
        with open(os.path.join(a.out, "index.csv"), "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(index_rows[0].keys()))
            w.writeheader(); w.writerows(index_rows)
    print(f"{total} curves written to {a.out}")


if __name__ == "__main__":
    main()
