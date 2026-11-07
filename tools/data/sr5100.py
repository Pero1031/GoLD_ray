"""Reader and converter for CSV exports of the Topcon SR-5100 (2D spectroradiometer).

Two kinds of files are read:
  * spot files: a [Parameters] block, an optional [Measurement] block (position X, Y, luminance, chromaticity of every measured
    point) and a [Spectral] block with one spectrum per point (header ",X,Y,380,381,..." or ",X,Y,380nm,...");
  * image files: a [Parameters] block and one [Spectrum (N nm)] block per wavelength holding a H x W image.

Command line (every command works on the exports as they come from the instrument; nothing is written next to them):

    python tools/data/sr5100.py info  FILE...                          what is in the files
    python tools/data/sr5100.py spots FILE --out DIR                   every spot as "wavelength_nm,value" CSV (what the renderer reads)
    python tools/data/sr5100.py roi   IMAGE.csv --box X0 Y0 X1 Y1 --out FILE      mean spectrum of a box of an image
    python tools/data/sr5100.py roi   IMAGE.csv --peak 2 --out FILE               mean spectrum around the brightest pixel (radius 2)
    python tools/data/sr5100.py reflectance SAMPLE.csv WHITE.csv --out FILE       sample / white (times the white's own reflectance)

The values stay in the unit of the instrument (spectral radiance): ratios such as a reflectance need no unit.
Measurements are lab data and are not part of the repository (assets/Metal_measured/ is ignored by Git).
"""
import argparse
import os
import re
import sys
from dataclasses import dataclass, field

import numpy as np


@dataclass
class Sr5100:
    path: str
    params: dict = field(default_factory=dict)
    wavelengths: np.ndarray = None            # (n_wl,) in nm
    kind: str = "spots"                       # "spots" or "image"
    spectra: np.ndarray = None                # (n_spots, n_wl), kind == "spots"
    spots: list = field(default_factory=list)  # [{"index", "x", "y", ...measurement columns}], kind == "spots"
    cube: np.ndarray = None                   # (n_wl, H, W), kind == "image"


def _read_text(path):
    raw = open(path, "rb").read()
    for enc in ("utf-8-sig", "cp932"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", "replace")


def _wavelength(token):
    m = re.fullmatch(r"\s*([0-9]+(?:\.[0-9]+)?)\s*(?:nm)?\s*", token)
    return float(m.group(1)) if m else None


def _floats(parts):
    out = []
    for p in parts:
        p = p.strip()
        if p == "":
            continue
        out.append(float(p))
    return out


def read(path):
    """Reads one export. Raises ValueError when the file is not an SR-5100 export."""
    text = _read_text(path)
    f = Sr5100(path=path)
    section = None
    wl = []                 # wavelengths of a spot file (from its header) or of an image file (one per section)
    images = []
    meas_header, meas_rows = [], []
    spec_rows = []
    for line in text.splitlines():
        if not line.strip().strip(","):
            continue
        m = re.match(r"\s*\[([^\]]+)\]", line)
        if m:
            section = m.group(1).strip()
            im = re.fullmatch(r"Spectrum \((\d+(?:\.\d+)?) nm\)", section)
            if im:
                wl.append(float(im.group(1)))
                images.append([])
            continue
        parts = line.rstrip().split(",")
        if section == "Parameters":
            f.params[parts[0].strip()] = ",".join(p for p in parts[1:] if p.strip())
        elif section == "Measurement":
            if parts[0].strip() == "" and len(parts) > 2 and parts[1].strip() == "X":
                meas_header = [p.strip() for p in parts[3:] if p.strip()]
            else:
                meas_rows.append(parts)
        elif section == "Spectral":
            if parts[0].strip() == "" and len(parts) > 2 and parts[1].strip() == "X":
                wl = [w for w in (_wavelength(p) for p in parts[3:]) if w is not None]
            else:
                spec_rows.append(parts)
        elif section and section.startswith("Spectrum ("):
            images[-1].append(_floats(parts))
    if images:
        f.kind = "image"
        f.wavelengths = np.array(wl)
        f.cube = np.array(images, dtype=np.float32)
        return f
    if not spec_rows:
        raise ValueError(f"{path}: no [Spectral] or [Spectrum (N nm)] block: not an SR-5100 export?")
    f.wavelengths = np.array(wl)
    n = len(wl)
    f.spectra = np.array([_floats(r[3:3 + n]) for r in spec_rows], dtype=np.float64)
    meas = {int(r[0]): r for r in meas_rows if r[0].strip().isdigit()}
    for r in spec_rows:
        idx = int(r[0])
        spot = {"index": idx, "x": int(r[1]), "y": int(r[2])}
        if idx in meas:
            for name, v in zip(meas_header, _floats(meas[idx][3:3 + len(meas_header)])):
                spot[name] = v
        f.spots.append(spot)
    return f


def reflectance(sample, white, white_reflectance=1.0, min_white=0.02):
    """sample / white * white_reflectance, per wavelength. NaN where the white is below `min_white` of its peak: there the instrument sees
    (almost) no light of the source and the ratio is noise (an LED has nothing in the ultraviolet and the far red)."""
    sample, white = np.asarray(sample, dtype=np.float64), np.asarray(white, dtype=np.float64)
    ok = white >= min_white * white.max()
    out = np.full_like(sample, np.nan)
    out[ok] = sample[ok] / white[ok] * white_reflectance
    return out


def write_curve(path, wavelengths, values, label="value", comments=()):
    """A two-column CSV that IO/SpectralCurveLoader.hpp reads ('#' comment lines, a header line, then nm,value).
    Wavelengths whose value is NaN are left out."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as out:
        for c in comments:
            out.write(f"# {c}\n")
        out.write(f"wavelength_nm,{label}\n")
        for w, v in zip(wavelengths, values):
            if np.isfinite(v):
                out.write(f"{w:g},{v:.8g}\n")


def _describe(f):
    p = f.params
    keys = ("Measurement distance[mm]", "Integration time", "ND filter", "White plate correction", "Spectrum correction")
    info = "  ".join(f"{k.split('[')[0]}={p[k]}" for k in keys if k in p)
    if f.kind == "image":
        return f"image  {f.cube.shape[2]}x{f.cube.shape[1]} px, {len(f.wavelengths)} wavelengths {f.wavelengths[0]:g}-{f.wavelengths[-1]:g} nm, " \
               f"max {f.cube.max():.4g}\n      {info}"
    peaks = ", ".join(f"#{s['index']}({s['x']},{s['y']}) peak {v.max():.3g}@{f.wavelengths[v.argmax()]:g}nm"
                      for s, v in zip(f.spots, f.spectra))
    return f"spots  {len(f.spots)} point(s), {len(f.wavelengths)} wavelengths {f.wavelengths[0]:g}-{f.wavelengths[-1]:g} nm: {peaks}\n      {info}"


def _stem(path):
    return os.path.splitext(os.path.basename(path))[0]


def _comments(f, what):
    d = f.params.get("Measurement distance[mm]", "?")
    return [f"{what} from {os.path.basename(f.path)} (Topcon SR-5100; distance {d} mm; integration {f.params.get('Integration time', '?')})",
            "values are in the unit of the instrument (spectral radiance) unless stated otherwise"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("info"); a.add_argument("files", nargs="+")
    a = sub.add_parser("spots"); a.add_argument("file"); a.add_argument("--out", required=True, help="output folder")
    a = sub.add_parser("roi"); a.add_argument("file"); a.add_argument("--out", required=True, help="output CSV")
    a.add_argument("--box", nargs=4, type=int, metavar=("X0", "Y0", "X1", "Y1"), help="pixels X0..X1, Y0..Y1 (both included)")
    a.add_argument("--peak", type=int, metavar="R", help="the (2R+1)^2 pixels around the brightest pixel")
    a = sub.add_parser("reflectance"); a.add_argument("sample"); a.add_argument("white")
    a.add_argument("--out", required=True, help="output CSV")
    a.add_argument("--spot", type=int, default=None, help="spot index of the sample (default: the mean of all spots)")
    a.add_argument("--white-spot", type=int, default=None, help="spot index of the white (default: the mean of all spots)")
    a.add_argument("--white-reflectance", type=float, default=1.0, help="reflectance of the white standard (default 1)")
    a.add_argument("--min-white", type=float, default=0.02, help="leave out wavelengths where the white is below this fraction of its peak (default 0.02)")
    a.add_argument("--range", nargs=2, type=float, metavar=("NM_LO", "NM_HI"), help="keep only this wavelength range")
    args = ap.parse_args(argv)

    if args.cmd == "info":
        for p in args.files:
            try:
                print(f"{os.path.basename(p)}: {_describe(read(p))}")
            except Exception as e:
                print(f"{os.path.basename(p)}: {e}")
    elif args.cmd == "spots":
        f = read(args.file)
        if f.kind != "spots":
            sys.exit("this is an image file: use 'roi'")
        for s, v in zip(f.spots, f.spectra):
            out = os.path.join(args.out, f"{_stem(f.path)}_spot{s['index']}.csv")
            write_curve(out, f.wavelengths, v, "radiance", _comments(f, f"spot #{s['index']} at pixel ({s['x']}, {s['y']})"))
            print("wrote", out)
    elif args.cmd == "roi":
        f = read(args.file)
        if f.kind != "image":
            sys.exit("this is a spot file: use 'spots'")
        H, W = f.cube.shape[1:]
        if args.peak is not None:
            lum = f.cube.sum(axis=0)
            y, x = np.unravel_index(np.argmax(lum), lum.shape)
            x0, x1, y0, y1 = max(0, x - args.peak), min(W - 1, x + args.peak), max(0, y - args.peak), min(H - 1, y + args.peak)
        elif args.box:
            x0, y0, x1, y1 = args.box
        else:
            sys.exit("give --box or --peak")
        v = f.cube[:, y0:y1 + 1, x0:x1 + 1].reshape(len(f.wavelengths), -1).mean(axis=1)
        write_curve(args.out, f.wavelengths, v, "radiance", _comments(f, f"mean of pixels x {x0}-{x1}, y {y0}-{y1}"))
        print(f"wrote {args.out}  (pixels x {x0}-{x1}, y {y0}-{y1}; peak {v.max():.4g} at {f.wavelengths[v.argmax()]:g} nm)")
    elif args.cmd == "reflectance":
        s, w = read(args.sample), read(args.white)
        if s.kind != "spots" or w.kind != "spots":
            sys.exit("reflectance works on spot files")
        if len(s.wavelengths) != len(w.wavelengths) or not np.allclose(s.wavelengths, w.wavelengths):
            sys.exit("the sample and the white have different wavelengths")

        def pick(f, idx):
            if idx is None:
                return f.spectra.mean(axis=0)
            return next(v for sp, v in zip(f.spots, f.spectra) if sp["index"] == idx)
        r = reflectance(pick(s, args.spot), pick(w, args.white_spot), args.white_reflectance, args.min_white)
        if args.range:
            r[(s.wavelengths < args.range[0]) | (s.wavelengths > args.range[1])] = np.nan
        write_curve(args.out, s.wavelengths, r, "reflectance",
                    [f"reflectance of {os.path.basename(s.path)} relative to {os.path.basename(w.path)} (white reflectance {args.white_reflectance:g})"])
        ok = np.isfinite(r)
        print(f"wrote {args.out}  ({ok.sum()} wavelengths, {s.wavelengths[ok][0]:g}-{s.wavelengths[ok][-1]:g} nm; values {np.nanmin(r):.3g} - {np.nanmax(r):.3g})")


if __name__ == "__main__":
    main()
