"""Scenes of skin under a layer of foundation, for stimuli where one control changes at a time (docs/COSMETICS.md).

    python tools/studies/cosmetic_stimuli.py make OUTDIR      writes the scenes
    (render every OUTDIR/*.json from GoLD_rayt with  --out FILE.png ; a fixed --seed makes a picture repeatable)

The picture is the close-up of scenes/skin_lamp_dry.json (a patch of skin 4.3 cm wide under a small lamp) without the dry flakes;
the skin has freckles and uneven pigment, which a layer has to hide. Names:
    bare                     no foundation
    PRESET_tT                one of the paper's four foundations (dewy1, dewy2, matte1, matte2) at optical thickness T
    cdC_aA_tT                a layer with diffuser fraction C and platelet roughness A (the other values of dewy1) at thickness T
"""
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def base_scene():
    src = REPO / "GoLD_rayt" / "scenes" / "skin_lamp_dry.json"
    scene = json.loads("".join(l for l in src.read_text(encoding="utf-8").splitlines(True) if not l.lstrip().startswith("//")))
    scene["materials"]["skin"].pop("scatter", None)
    scene["options"].update({"width": 640, "spp": 128, "denoise": False})
    return scene


def make(out):
    out.mkdir(parents=True, exist_ok=True)
    scenes = {"bare": None}
    for preset in ("dewy1", "matte2"):
        for t in (0.25, 0.5, 1.0, 2.0, 4.0):
            scenes[f"{preset}_t{t:g}"] = {"preset": preset, "thickness": t}
    for cd in (0.5, 0.8, 0.95, 1.0):
        for alpha in (0.05, 0.15, 0.4):
            scenes[f"cd{cd:g}_a{alpha:g}_t1"] = {"preset": "dewy1", "thickness": 1.0, "diffuserFraction": cd, "plateletRoughness": alpha}
    for name, foundation in scenes.items():
        scene = base_scene()
        if foundation:
            scene["materials"]["skin"]["foundation"] = foundation
        (out / f"{name}.json").write_text(json.dumps(scene, indent=1) + "\n", encoding="utf-8")
    print("wrote", len(scenes), "scenes to", out)


def sheet(out):
    """Two sheets of the rendered pictures (the thickness series, the grid of finishes) and a few numbers of each picture."""
    import numpy as np
    from PIL import Image, ImageDraw, ImageFont

    try:
        font = ImageFont.truetype("arial.ttf", 15)
    except OSError:
        font = ImageFont.load_default()

    def numbers(name):
        a = np.asarray(Image.open(out / f"{name}.png").convert("RGB"), dtype=float) / 255
        lin = np.where(a <= 0.04045, a / 12.92, ((a + 0.055) / 1.055) ** 2.4)
        y = lin @ np.array([0.2126, 0.7152, 0.0722])
        # the pattern of the skin (freckles, uneven pigment) without the pixel noise: the luminance averaged over 16 x 16 pixels, away from the glint
        h, w = y.shape[0] // 16 * 16, y.shape[1] // 16 * 16
        blocks = y[:h, :w].reshape(h // 16, 16, w // 16, 16).mean(axis=(1, 3))
        calm = blocks[blocks < np.percentile(blocks, 90)]
        return lin.reshape(-1, 3).mean(axis=0), calm.std() / calm.mean(), (y > 0.9).mean()

    def build(rows, row_labels, col_labels, path):
        first = Image.open(out / f"{rows[0][0]}.png")
        w, h = first.size[0] // 2, first.size[1] // 2
        pad, head, left = 6, 22, 110
        img = Image.new("RGB", (left + len(rows[0]) * (w + pad), head + len(rows) * (h + pad)), (255, 255, 255))
        draw = ImageDraw.Draw(img)
        for c, label in enumerate(col_labels):
            draw.text((left + c * (w + pad), 3), label, fill=(0, 0, 0), font=font)
        for r, names in enumerate(rows):
            draw.text((4, head + r * (h + pad) + h // 2 - 8), row_labels[r], fill=(0, 0, 0), font=font)
            for c, name in enumerate(names):
                img.paste(Image.open(out / f"{name}.png").convert("RGB").resize((w, h), Image.LANCZOS), (left + c * (w + pad), head + r * (h + pad)))
        img.save(path)
        print("wrote", path, img.size)

    thick = ("0.25", "0.5", "1", "2", "4")
    build([["bare"] + [f"{p}_t{t}" for t in thick] for p in ("dewy1", "matte2")], ["dewy 1", "matte 2"], ["bare skin"] + [f"thickness {t}" for t in thick],
          out / "sheet_thickness.png")
    cds, alphas = ("0.5", "0.8", "0.95", "1"), ("0.05", "0.15", "0.4")
    build([[f"cd{c}_a{a}_t1" for c in cds] for a in alphas], [f"roughness {a}" for a in alphas], [f"diffusers {float(c) * 100:.0f} %" for c in cds],
          out / "sheet_finish.png")
    print("picture                mean linear RGB           pattern contrast   bright pixels")
    for name in ["bare"] + [f"{p}_t{t}" for p in ("dewy1", "matte2") for t in thick] + [f"cd{c}_a{a}_t1" for a in alphas for c in cds]:
        rgb, contrast, bright = numbers(name)
        print(f"{name:20s}  ({rgb[0]:.3f}, {rgb[1]:.3f}, {rgb[2]:.3f})     {100 * contrast:6.2f} %        {100 * bright:6.3f} %")


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "make":
        make(Path(sys.argv[2]))
    elif len(sys.argv) == 3 and sys.argv[1] == "sheet":
        sheet(Path(sys.argv[2]))
    else:
        print(__doc__)
