"""Converts the LED-cube's per-LED spectra (the lab's LEDCube_LEDoutput/1.csv .. 18.csv: "300nm,value" rows, 1 nm,
the maximum spectral radiance of each of the 18 LEDs) into "wavelength_nm, radiance" CSVs for the scene file's light mix:

  "spectrum": { "mix": [ { "file": "assets/Fading_measured/lights/led01.csv", "weight": 0.4 }, ... ] }

The scale is kept, so the weight of LED i is its output setting w_i (0..1) of the cube. The files are written for 360..830 nm
(the range of the engine's emission table; the LEDs emit nothing visible outside it).

Usage:
  python tools/data/ledcube_export.py <LEDCube_LEDoutput folder> --out <output folder>
"""
import argparse
import os
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", help="folder with 1.csv .. 18.csv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--range", nargs=2, type=float, default=(360.0, 830.0), metavar=("LO_NM", "HI_NM"))
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    n = 0
    for i in range(1, 19):
        src = os.path.join(a.folder, f"{i}.csv")
        if not os.path.exists(src):
            print(f"missing: {src}", file=sys.stderr)
            continue
        rows = []
        with open(src, encoding="utf-8", errors="ignore") as f:
            for line in f:
                parts = line.strip().split(",")
                if len(parts) < 2:
                    continue
                try:
                    rows.append((float(parts[0].lower().replace("nm", "")), float(parts[1])))
                except ValueError:
                    continue
        rows = [(w, v) for w, v in rows if a.range[0] <= w <= a.range[1]]
        peak = max(rows, key=lambda r: r[1])
        dest = os.path.join(a.out, f"led{i:02d}.csv")
        with open(dest, "w", encoding="utf-8", newline="\n") as out:
            out.write(f"# LED-cube LED {i} at full output (peak {peak[0]:g} nm), maximum spectral radiance as measured\n")
            out.write(f"# source: {src}\n")
            out.write("wavelength_nm,radiance\n")
            for w, v in rows:
                out.write(f"{w:g},{v:.6g}\n")
        n += 1
    print(f"{n} LED spectra written to {a.out}")


if __name__ == "__main__":
    main()
