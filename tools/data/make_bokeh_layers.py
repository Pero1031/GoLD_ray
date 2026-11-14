"""Writes scenes/bokeh_layers.json: a portrait with four layers of depth, as in a photograph of a flower with blur in front of it and behind it.

    python tools/data/make_bokeh_layers.py GoLD_rayt/scenes/bokeh_layers.json

  front    big glossy balls 0.9 to 1.4 m from the camera at the edges of the picture: soft coloured masses (they overlap the subject's edge)
  subject  the gold Erato (a 3D scan, see assets/Meshes/erato/CREDITS.md), in focus, at 2.6 m
  behind   balls of gold, copper, glass and ceramic 1 to 2 m behind the subject, partly hidden by it: blurred, but still balls
  lights   small lights 6 to 10 m from the camera: discs in the shape of the aperture (7 blades)

The lens is 100 mm at f/1.6 (a portrait lens), the picture 3:2. The lights are placed by a script with a fixed seed.
"""
import json
import math
import random
import sys

out_path = sys.argv[1]
rng = random.Random(5)

cam_from = (0.19, 0.24, 2.58)
cam_at = (0.19, 0.25, 0.0)
focal_mm = 100.0
aspect = 3.0 / 2.0
sensor_h = 36.0 / aspect
vfov = 2 * math.atan(sensor_h / (2 * focal_mm))
hfov = 2 * math.atan(math.tan(vfov / 2) * aspect)


def at_screen(u, v, depth):
    """The point `depth` metres in front of the camera plane that appears at (u, v) of the picture (-0.5 .. 0.5)."""
    return (cam_from[0] + depth * math.tan(u * hfov), cam_from[1] + depth * math.tan(v * vfov), cam_from[2] - depth)


temps = [1800, 2200, 2700, 3200, 4000, 5000, 6500, 9000]

materials = {
    "gold":    {"type": "conductor", "metal": "Au", "roughness": 0.15},
    "copper":  {"type": "conductor", "metal": "Cu", "roughness": 0.05},
    "glass":   {"type": "dielectric", "ior": 1.5, "abbe": 30},
    "ceramic": {"type": "layered", "coat": {"ior": 1.5, "roughness": 0.05}, "base": {"albedo": [0.80, 0.78, 0.72]}},
    "rose":    {"type": "layered", "coat": {"ior": 1.5, "roughness": 0.05}, "base": {"albedo": [0.80, 0.18, 0.32]}},
    "amber":   {"type": "layered", "coat": {"ior": 1.5, "roughness": 0.05}, "base": {"albedo": [0.85, 0.50, 0.12]}},
    "floor":   {"type": "layered", "coat": {"ior": 1.5, "roughness": 0.10}, "base": {"albedo": 0.05}},
    "key":     {"type": "emission", "spectrum": "tungsten", "radiance": 20},
    "fill":    {"type": "emission", "spectrum": "tungsten", "radiance": 8},
    "rim":     {"type": "emission", "spectrum": "tungsten", "radiance": 14},
    "front":   {"type": "emission", "spectrum": "tungsten", "radiance": 5},
    "glow":    {"type": "emission", "spectrum": {"blackbody": 4200}, "radiance": 0.04, "twoSided": True},
}
for T in temps:
    materials["lamp%d" % T] = {"type": "emission", "spectrum": {"blackbody": T}, "radiance": 18}

objects = []


def sphere(c, r, m):
    objects.append({"shape": "sphere", "center": [round(x, 4) for x in c], "radius": r, "material": m})


# ---- subject: the gold Erato (the scene of scenes/erato_gold.json: 0.5 m tall, its base on the floor) ----
objects.append({"shape": "mesh", "file": "assets/Meshes/erato/erato_clean.obj", "scale": 0.0133, "rotate": {"axis": [0, 1, 0], "degrees": 90},
                "translate": [-0.031, 0.0169, 0.001], "material": "gold"})
objects.append({"shape": "sphere", "center": [0, -100, 0], "radius": 100, "material": "floor"})

# ---- front: big glossy balls at the edges, close to the camera; the first one's soft edge reaches the foot of the statue ----
sphere(at_screen(-0.44, -0.36, 1.0), 0.060, "rose")
sphere(at_screen(0.44, -0.14, 1.2), 0.040, "amber")
sphere(at_screen(-0.40, 0.40, 0.9), 0.018, "rose")

# ---- behind: balls on the floor behind the statue, the first one half hidden by it ----
for (u_x, depth, r, m) in [(0.054, 3.78, 0.080, "gold"), (0.42, 4.08, 0.060, "copper"), (0.62, 3.58, 0.070, "glass"), (-0.30, 3.98, 0.050, "ceramic")]:
    sphere((u_x, r, cam_from[2] - depth), r, m)

# ---- lights: the three lamps of the Erato scene, the rim one raised out of the picture ----
sphere((-0.6, 0.55, 0.7), 0.25, "key")
sphere((0.9, 0.30, 0.9), 0.20, "fill")
sphere((0.5, 0.95, -0.7), 0.20, "rim")
objects.append({"shape": "quad", "corner": [-0.10, 0.45, 2.20], "edge1": [0.40, 0, 0], "edge2": [0, 0, 0.30], "material": "front"})   # a soft fill for the front balls
objects.append({"shape": "sphere", "center": [0, 0, 0], "radius": 30, "material": "glow", "light": False})

# ---- light bokeh: small lights far behind, apart from each other, above the floor ----
screen = []
lights = []
tries = 0
while len(lights) < 22 and tries < 20000:
    tries += 1
    depth = rng.uniform(6.0, 10.0)
    u = rng.uniform(-0.5, 0.5) * 1.08
    v = rng.uniform(-0.5, 0.5) * 1.12
    p = at_screen(u, v, depth)
    if p[1] < 0.08:      # under the floor, or on it
        continue
    if any(math.hypot((u - a) * aspect, v - b) < 0.16 for a, b in screen):
        continue
    screen.append((u, v))
    lights.append((p, round(rng.uniform(0.02, 0.05), 4), rng.choice(temps)))
for p, r, T in lights:
    sphere(p, r, "lamp%d" % T)

lines = []
lines.append("{")
lines.append("  // Four layers of depth (a portrait, 100 mm, f/1.6, focus on the gold Erato): glossy balls close to the camera at the edges (blurred masses in front),")
lines.append("  // the statue (sharp), balls of gold, copper, glass and ceramic behind it (soft but still balls, partly hidden by the statue), and %d small lights far behind" % len(lights))
lines.append("  // (heptagonal discs). The mesh is not in the repository (assets/Meshes/erato/CREDITS.md). Written by tools/data/make_bokeh_layers.py.")
lines.append('  "options": { "width": 1200, "height": 800, "spp": 256, "depth": 24, "denoise": true },')
lines.append('  "camera": { "from": %s, "at": %s, "up": [0, 1, 0], "focalLength": %g, "fStop": 1.6, "focusAt": [0, 0.30, 0.06], "blades": 7 },'
             % (json.dumps(list(cam_from)), json.dumps(list(cam_at)), focal_mm))
lines.append('  "materials": {')
items = list(materials.items())
for k, (name, m) in enumerate(items):
    lines.append('    "%s": %s%s' % (name, json.dumps(m), "," if k + 1 < len(items) else ""))
lines.append("  },")
lines.append('  "objects": [')
for k, o in enumerate(objects):
    lines.append("    %s%s" % (json.dumps(o), "," if k + 1 < len(objects) else ""))
lines.append("  ]")
lines.append("}")
open(out_path, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
print("wrote", out_path, "with", len(lights), "lights")
