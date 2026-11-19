"""
Apply real view transforms (OpenColorIO) to linear Rec.709 colors: a helper for view_transform_hiding.py.

Runs in any Python that has PyOpenColorIO, for example the one bundled with Blender (<blender>/<version>/python/bin/python.exe),
which also brings Blender's own OCIO config and LUTs. Input: an .npy array (..., 3) of scene-linear Rec.709 values (the renderer's
film; negative values allowed). Output: an .npz with one array per pipeline, the display-encoded values clipped to [0, 1] (what a PNG
would hold, before 8-bit rounding).

Pipelines (key = "<config>|<display>|<view>"):
  blender|sRGB|Standard, blender|sRGB|Filmic, blender|sRGB|AgX, blender|sRGB|Khronos PBR Neutral,
  blender|Display P3|Standard, blender|Display P3|AgX               Blender's config.ocio (input: "Linear Rec.709")
  aces|sRGB - Display|ACES 1.0 - SDR Video, aces|Display P3 - Display|ACES 1.0 - SDR Video
                                                                    OCIO's built-in ACES 1.3 CG config (input: "Linear Rec.709 (sRGB)")

  <ocio python> tools/studies/ocio_apply.py <in.npy> <out.npz> --blender-config <config.ocio> [--ev E] [--keys k1;k2;...]
"""
import argparse

import numpy as np
import PyOpenColorIO as ocio

ACES_CONFIG = "ocio://cg-config-v2.2.0_aces-v1.3_ocio-v2.4"
SCENE_LINEAR = {"blender": "Linear Rec.709", "aces": "Linear Rec.709 (sRGB)"}
DEFAULT_KEYS = ["blender|sRGB|Standard", "blender|sRGB|Filmic", "blender|sRGB|AgX", "blender|sRGB|Khronos PBR Neutral",
                "blender|Display P3|Standard", "blender|Display P3|AgX",
                "aces|sRGB - Display|ACES 1.0 - SDR Video", "aces|Display P3 - Display|ACES 1.0 - SDR Video"]


def processor(config, src, display, view):
    t = ocio.DisplayViewTransform(src=src, display=display, view=view)
    return config.getProcessor(t).getDefaultCPUProcessor()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--blender-config", required=True)
    ap.add_argument("--ev", type=float, default=0.0, help="exposure in stops, applied in scene-linear before the view (as Blender's Exposure)")
    ap.add_argument("--keys", default=";".join(DEFAULT_KEYS))
    a = ap.parse_args()

    configs = {"blender": ocio.Config.CreateFromFile(a.blender_config), "aces": ocio.Config.CreateFromFile(ACES_CONFIG)}
    x = np.load(a.inp).astype(np.float32) * np.float32(2.0 ** a.ev)
    shape = x.shape
    flat = np.ascontiguousarray(x.reshape(-1, 3))
    out = {}
    for key in a.keys.split(";"):
        cfg, display, view = key.split("|")
        p = processor(configs[cfg], SCENE_LINEAR[cfg], display, view)
        y = flat.copy()
        p.applyRGB(y)
        out[key] = np.clip(y, 0.0, 1.0).reshape(shape)
    np.savez(a.out, **out)
    print(f"OCIO {ocio.__version__}: {len(out)} pipelines, {flat.shape[0]} colors, EV {a.ev:+g} -> {a.out}")


if __name__ == "__main__":
    main()
