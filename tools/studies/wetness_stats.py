"""Does a rendered wet surface change the image statistics the way real wet surfaces do?

Sawayama, Adelson & Nishida, "Visual wetness perception based on image color statistics", Journal of Vision 17(5):7, 2017,
photographed 32 materials dry and wet: wetting raised the u'v' saturation (1.43 times on average), left the hue alone, lowered the
mean luminance and raised the skew of the luminance histogram. An image operation with the same changes (their "wetness enhancing
transformation", WET: luminance Y -> Y^gamma, saturation raised) makes dry pictures look wetter, the more so the more hues a
picture has (hue entropy). This script measures those statistics on renders and applies the WET operator.

    python tools/studies/wetness_stats.py make OUTDIR            writes the pebble scenes (two palettes x three surfaces)
    (render every OUTDIR/*.json from GoLD_rayt with  --out FILE.pfm)
    python tools/studies/wetness_stats.py analyze OUTDIR         the statistics of each pair, a sheet (dry, wet, WET of dry) per pair
    python tools/studies/wetness_stats.py stats DRY WET          the statistics of any pair of pictures (.pfm linear, or .png sRGB)
    python tools/studies/wetness_stats.py wet IN OUT.png [gamma add factor]   the WET operator on one picture

The pebbles are spheres of random colour, half sunk into a grey ground, under a large dim light overhead and a lamp (a sphere of 30 cm at
about 1 m: its mirror image on a pebble is a few pixels wide).
Surfaces (all with the same body colour): "lambert" (a diffuse surface with no reflection at the top), "rough" (a rough clear coat of
index 1.5 over the diffuse body: the dry solid with its hazy surface reflection), "wet" (a smooth coat of index 1.33 over the body: a
film of water). The palettes: "many" (hues all around the circle) and "one" (one hue), for the hue entropy, and "natural"
(bright, pale, warm: closer to the dry stones of the paper's photographs; the other two are darker and more colourful).

What the numbers are: the saturation is s_uv = 13 sqrt((u' - u'n)^2 + (v' - v'n)^2) with the D65 white (the paper's equation 1), the
hue is the angle in the u'v' plane, the hue entropy is that of a 512-bin histogram (bits). The pictures of a pair share one exposure
(white = 1.25 times the 99.5th percentile of the dry picture's luminance) and are clipped there like a photograph, so a highlight
counts as a white pixel and not as a value of a few hundred.
"""
import colorsys
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from hdr_io import read_pfm  # noqa: E402

RGB2XYZ = np.array([[0.4124564, 0.3575761, 0.1804375], [0.2126729, 0.7151522, 0.0721750], [0.0193339, 0.1191920, 0.9503041]])
XYZ2RGB = np.linalg.inv(RGB2XYZ)
WHITE_UV = (0.19783, 0.46832)   # D65
BINS = 512


def load_linear(path):
    """Linear sRGB, float (H, W, 3)."""
    path = str(path)
    if path.lower().endswith(".pfm"):
        a = np.maximum(read_pfm(path).astype(np.float64), 0.0)
        h, w = a.shape[0] // 2 * 2, a.shape[1] // 2 * 2   # 2 x 2 pixels averaged: a quarter of the render noise, which would add to the SD and the saturation
        return a[:h, :w].reshape(h // 2, 2, w // 2, 2, 3).mean(axis=(1, 3))
    from PIL import Image
    v = np.asarray(Image.open(path).convert("RGB"), dtype=np.float64) / 255
    return np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4)


def to_uvY(lin):
    xyz = lin @ RGB2XYZ.T
    d = np.maximum(xyz[..., 0] + 15 * xyz[..., 1] + 3 * xyz[..., 2], 1e-12)
    return 4 * xyz[..., 0] / d, 9 * xyz[..., 1] / d, xyz[..., 1]


def from_uvY(u, v, Y):
    v = np.maximum(v, 1e-6)
    X = Y * 9 * u / (4 * v)
    Z = Y * (12 - 3 * u - 20 * v) / (4 * v)
    return np.stack([X, Y, Z], axis=-1) @ XYZ2RGB.T


def expose(lin, white):
    """The picture as a camera would record it: scaled so that `white` is 1, clipped per channel."""
    return np.clip(lin / white, 0.0, 1.0)


def white_of(lin):
    return 1.25 * np.percentile(lin @ RGB2XYZ[1], 99.5)


def statistics(img):
    """The paper's statistics of an exposed picture (values in [0, 1])."""
    u, v, Y = to_uvY(img)
    ok = Y > 1e-4
    du, dv = u[ok] - WHITE_UV[0], v[ok] - WHITE_UV[1]
    s = 13 * np.hypot(du, dv)
    hue = np.arctan2(dv, du)
    y = Y.ravel()
    m, sd = y.mean(), y.std()
    hist = np.histogram(hue, bins=BINS, range=(-np.pi, np.pi))[0].astype(float)
    p = hist[hist > 0] / hist.sum()
    L = 116 * np.cbrt(np.maximum(Y[ok], 1e-12)) - 16
    return {"sat": s.mean(), "hue": np.degrees(np.arctan2(np.sin(hue).mean(), np.cos(hue).mean())),
            "mean": m, "sd": sd, "skew": ((y - m) ** 3).mean() / sd ** 3, "kurt": ((y - m) ** 4).mean() / sd ** 4,
            "hueH": -(p * np.log2(p)).sum(), "chroma": (np.maximum(L, 0) * s).mean(), "clipped": (y >= 0.999).mean()}


def wet_operator(img, gamma=3.0, add=0.26, factor=1.0):
    """The WET operator on an exposed picture: Y -> Y^gamma, s_uv -> factor * s_uv + add (the hue is kept)."""
    u, v, Y = to_uvY(img)
    du, dv = u - WHITE_UV[0], v - WHITE_UV[1]
    r = np.hypot(du, dv)
    k = np.where(r > 1e-6, (factor * r + add / 13) / np.maximum(r, 1e-6), 1.0)
    return np.clip(from_uvY(WHITE_UV[0] + du * k, WHITE_UV[1] + dv * k, Y ** gamma), 0.0, 1.0)


def to_png(img, path):
    from PIL import Image
    v = np.where(img <= 0.0031308, 12.92 * img, 1.055 * np.maximum(img, 0) ** (1 / 2.4) - 0.055)
    Image.fromarray((np.clip(v, 0, 1) * 255 + 0.5).astype(np.uint8)).save(path)


HEAD = f"{'':28s} {'sat':>6s} {'hue':>6s} {'mean Y':>7s} {'SD':>6s} {'skew':>6s} {'kurt':>6s} {'hue H':>6s} {'chroma':>7s} {'white %':>8s}"


def line(name, st):
    return (f"{name:28s} {st['sat']:6.3f} {st['hue']:6.1f} {st['mean']:7.4f} {st['sd']:6.4f} {st['skew']:6.2f} {st['kurt']:6.2f}"
            f" {st['hueH']:6.2f} {st['chroma']:7.2f} {100 * st['clipped']:8.3f}")


def compare(dry_path, wet_path, label):
    dry, wet = load_linear(dry_path), load_linear(wet_path)
    white = 1.0 if str(dry_path).lower().endswith(".png") else white_of(dry)   # a PNG is already exposed
    a, b = statistics(expose(dry, white)), statistics(expose(wet, white))
    print(line(label + " dry", a))
    print(line(label + " wet", b))
    print(f"{'  wet / dry':28s} sat x{b['sat'] / a['sat']:.2f}, hue {b['hue'] - a['hue']:+.1f} deg, mean Y x{b['mean'] / a['mean']:.2f},"
          f" skew {b['skew'] - a['skew']:+.2f}, chroma x{b['chroma'] / a['chroma']:.2f}")
    return white, a, b


# --- the scenes

# The colours of the pebbles: (hue range, saturation range, value range) of HSV, and the albedo of the ground.
# "natural" is bright and pale like the dry pebbles, rock and wood of the paper's photographs; the first two are darker and more colourful.
PALETTES = {"many": ((0.0, 1.0), (0.3, 0.65), (0.25, 0.65), [0.3, 0.3, 0.3]),
            "one": ((0.055, 0.085), (0.3, 0.65), (0.25, 0.65), [0.3, 0.3, 0.3]),
            "natural": ((0.04, 0.13), (0.08, 0.35), (0.55, 0.9), [0.6, 0.55, 0.47])}
SURFACES = ("lambert", "rough", "wet")


def surface(kind, albedo):
    if kind == "lambert":
        return {"type": "diffuse", "albedo": albedo}
    coat = {"ior": 1.5, "roughness": 0.6} if kind == "rough" else {"ior": 1.33, "roughness": 0.03}
    return {"type": "layered", "coat": coat, "base": {"albedo": albedo}}


def make(out):
    out.mkdir(parents=True, exist_ok=True)
    for palette, (hues, sats, values, ground) in PALETTES.items():
        rng = np.random.default_rng(11)
        pebbles = []
        for _ in range(4000):
            if len(pebbles) == 240:
                break
            r = rng.uniform(0.012, 0.034)
            x, z = rng.uniform(-0.36, 0.36), rng.uniform(-0.42, 0.26)
            if all(np.hypot(x - px, z - pz) > 0.8 * (r + pr) for px, pz, pr, _ in pebbles):
                rgb = colorsys.hsv_to_rgb(rng.uniform(*hues), rng.uniform(*sats), rng.uniform(*values))
                pebbles.append((x, z, r, [round(c, 4) for c in rgb]))
        for kind in SURFACES:
            materials = {"ground": surface(kind, ground),
                         "lamp": {"type": "emission", "spectrum": "d65", "radiance": 24},
                         "sky": {"type": "emission", "spectrum": "d65", "radiance": 1.0}}
            objects = [{"shape": "quad", "material": "ground", "corner": [-1.5, 0, 1.5], "edge1": [3, 0, 0], "edge2": [0, 0, -3]},
                       {"shape": "sphere", "material": "lamp", "center": [0.35, 0.8, -0.25], "radius": 0.15},
                       {"shape": "quad", "material": "sky", "corner": [-1, 1.5, -1], "edge1": [2, 0, 0], "edge2": [0, 0, 2]}]
            for i, (x, z, r, rgb) in enumerate(pebbles):
                materials[f"p{i}"] = surface(kind, rgb)
                objects.append({"shape": "sphere", "material": f"p{i}", "center": [round(x, 4), round(0.35 * r, 4), round(z, 4)], "radius": round(r, 4)})
            scene = {"options": {"width": 1024, "spp": 128, "depth": 10}, "camera": {"from": [0, 0.55, 0.6], "at": [0, 0, -0.04], "fov": 32},
                     "materials": materials, "objects": objects}
            (out / f"pebbles_{palette}_{kind}.json").write_text(json.dumps(scene, indent=1) + "\n", encoding="utf-8")
    print("wrote", len(list(out.glob("pebbles_*.json"))), "scenes to", out)


def analyze(out):
    from PIL import Image, ImageDraw, ImageFont
    print("paper (32 photographed materials): saturation x1.43, hue unchanged, mean luminance lower, skew higher\n")
    print(HEAD)
    for palette in PALETTES:
        for dry_kind in ("lambert", "rough"):
            dry_path, wet_path = out / f"pebbles_{palette}_{dry_kind}.pfm", out / f"pebbles_{palette}_wet.pfm"
            if not (dry_path.exists() and wet_path.exists()):
                continue
            white, a, b = compare(dry_path, wet_path, f"{palette} / {dry_kind}:")
            dry, wet = expose(load_linear(dry_path), white), expose(load_linear(wet_path), white)
            op = wet_operator(dry)
            print(line("  WET operator on the dry", statistics(op)) + "\n")
            cells = [("dry (" + dry_kind + ")", dry), ("wet (rendered: a film of water)", wet), ("WET operator on the dry picture", op)]
            h, w = dry.shape[:2]
            pad, head = 8, 26
            try:
                font = ImageFont.truetype("arial.ttf", 17)
            except OSError:
                font = ImageFont.load_default()
            sheet = Image.new("RGB", (3 * (w + pad) + pad, h + head + 2 * pad), (255, 255, 255))
            draw = ImageDraw.Draw(sheet)
            for c, (label, img) in enumerate(cells):
                tmp = out / "_cell.png"
                to_png(img, tmp)
                draw.text((pad + c * (w + pad), pad + 2), label, fill=(0, 0, 0), font=font)
                sheet.paste(Image.open(tmp), (pad + c * (w + pad), pad + head))
                tmp.unlink()
            sheet.save(out / f"sheet_{palette}_{dry_kind}.png")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "make":
        make(Path(sys.argv[2]))
    elif cmd == "analyze":
        analyze(Path(sys.argv[2]))
    elif cmd == "stats":
        print(HEAD)
        compare(sys.argv[2], sys.argv[3], "")
    elif cmd == "wet":
        img = load_linear(sys.argv[2])
        args = [float(v) for v in sys.argv[4:7]]
        to_png(wet_operator(expose(img, white_of(img)), *args), sys.argv[3])
    else:
        print(__doc__)
