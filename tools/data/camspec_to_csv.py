"""Camera spectral sensitivities of Jiang, Liu, Gu and Suesstrunk (WACV 2013) as one CSV per camera, for --camera-response.

The database (28 cameras, 400..720 nm in 10 nm steps) is not part of this repository: it is licensed CC BY-NC-SA 4.0 (Rochester
Institute of Technology). Download camspec_database.txt from http://www.gujinwei.org/research/camspec/db.html into
GoLD_rayt/assets/CameraResponse/ (git-ignored) and run

    python tools/data/camspec_to_csv.py GoLD_rayt/assets/CameraResponse/camspec_database.txt

Each camera becomes <folder>/<name>.csv with the columns wavelength_nm, R, G, B (the raw sensitivities as measured, each channel
scaled as in the database).
"""
import os
import re
import sys


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    src = sys.argv[1]
    out_dir = os.path.dirname(os.path.abspath(src))
    lines = [ln.strip() for ln in open(src, encoding="utf-8", errors="replace") if ln.strip()]
    if len(lines) % 4:
        sys.exit("%s: expected blocks of 4 lines (a name, then R, G and B)" % src)
    wavelengths = list(range(400, 721, 10))
    names = []
    for k in range(0, len(lines), 4):
        name = lines[k]
        curves = [[float(v) for v in lines[k + c].split()] for c in (1, 2, 3)]
        if any(len(c) != len(wavelengths) for c in curves):
            sys.exit("%s: %s: expected %d values per channel" % (src, name, len(wavelengths)))
        file = re.sub(r"[^A-Za-z0-9.+-]+", "_", name).strip("_") + ".csv"
        with open(os.path.join(out_dir, file), "w", encoding="utf-8", newline="\n") as f:
            f.write("# %s: camera spectral sensitivities (raw R, G, B) from Jiang, Liu, Gu and Suesstrunk, WACV 2013.\n" % name)
            f.write("# CC BY-NC-SA 4.0, Rochester Institute of Technology; http://www.gujinwei.org/research/camspec/db.html\n")
            f.write("wavelength_nm,R,G,B\n")
            for i, nm in enumerate(wavelengths):
                f.write("%d,%.7g,%.7g,%.7g\n" % (nm, curves[0][i], curves[1][i], curves[2][i]))
        names.append(file)
    print("%d cameras written to %s:" % (len(names), out_dir))
    print("  " + "  ".join(names))


if __name__ == "__main__":
    main()
