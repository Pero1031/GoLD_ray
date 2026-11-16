# Glass with water drops (2026-11-13)

**Status:** nested dielectrics (an outside refractive index per surface), a water drop generator, sample scenes.

## 1. What is needed

Whether a water drop looks "wet" on the glass is decided at the contact between the water and the glass. It is a boundary where light goes directly from water (1.333) into glass (1.5), and the reflection is only about 0.3%. With a trick that puts an air gap in between, two boundaries (4% each) appear, and the drop looks like a marble sitting on the glass.

The remedy is the same approach as Mitsuba 3: split the surface where the refractive index changes, and let each surface carry the indices of its inside and outside (`exteriorIor` of `dielectric`, `docs/SCENE_FILE.md`). A water drop is made of two surfaces:

| Surface | Material | Normal |
|---|---|---|
| Cap (the free surface in contact with air) | `ior` 1.333 | outward |
| Disc of the contact face (in contact with the glass) | `ior` 1.333, `exteriorIor` 1.5 | toward the glass (−z) |

And **the top face of the glass under the drop is removed**. A `cutout` (an RGBA mask; where the alpha is 0, every ray passes through) is attached to the mesh of the glass top face.

Apart from the priority approach of Arnold and RenderMan (Schmidt & Budge 2002: objects are made overlapping, and a ray carries an "interior list"), the contact face is written explicitly here.

## 2. Generator

```bat
python tools\data\make_drops.py GoLD_rayt\assets\Drops --plate 0.20 0.15 --coverage 0.25 --tilt 60 --seed 1
```

It writes `caps.ply` (caps), `contacts.ply` (contact faces), `plate_top.ply` (glass top face with uv), `plate_mask.png` (the holes) and `drops.csv`. The plate lies in the x-y plane at z = 0, the drops are on the +z side, and down is −y.

The model is phenomenological (an extension of the "part of a sphere" of Kaneda et al. 1993):
- The contact radius follows a log-normal distribution (median 1.2 mm, width 0.6, range 0.25-4 mm). The drops are placed without overlap, up to a fraction of the area (`--coverage`).
- The shape is a spherical cap with contact angle θ (≤ 90°; clean glass 20-40°, a used car window 50-70°). Large drops are flattened by gravity (the height is lowered using the capillary length of 2.7 mm as a guide).
- On a tilted plate (`--tilt`, the angle from horizontal), drops large enough to feel gravity become **pear-shaped**: a round head at the bottom and a tail tapering upward (the contact line is a star-shaped polygon from the origin; the head is a semicircle of radius a, the tail has length a(1 + 0.8 s b) and pointedness 0.6 s b, with b = min(1, (a/2.7 mm)²)). The apex shifts downward (the lower side is at the advancing angle and steep, the upper side at the receding angle and gentle). On a steep plate, some of the large drops become long, thin trails (a long pear shape with the apex near the lower end).
- **Mist of tiny droplets** (`--mist`, fraction of the area; default 2%): droplets of 0.05-0.25 mm laid in the gaps between the large drops. These are the fine grains around the large drops in photographs.

One drop has at most 40 segments × 10 rings of triangles (mist: 10 × 3). With 2370 drops + 12,000 mist droplets this is about 1.85 million triangles (50 MB). Generation takes 7 seconds. The mask paints the polygon of the contact line (any shape is fine).

## 3. Samples

`scenes/water_drops.json`: a street lamp at night (an emitting sphere) seen through the glass. The camera now has `aperture` (lens diameter) and `focus`. If the glass is in focus, the lamps become bokeh discs; if `focus` is set to the distance of the lamps, the drops blur and the lamps become sharp (this corresponds to the two conditions "focus on the drops / focus on the background" of Raindrop Clarity). Each drop forms a small inverted image of the lamps as a fisheye lens, and its rim glows. 900 × 600, 256 spp, about 4 minutes (development build).

`scenes/water_drops_soft.json`: the same plate seen under a soft light with the lamps placed outside the frame. Each drop mirrors upside down a colored sphere 1 m behind and the gradient of an area light diagonally above, and a small lamp on the camera side puts light into the rims. Use this one to look at the shape of the drops themselves.

`scenes/water_drops_overcast.json`: uses an overcast-sky HDR (CC0 from Poly Haven; not put in git because of its size) as the environment. It is the closest to photographs of real drops: because the sky is a large, soft light source, each drop shows an inverted image of the bright sky and the dark horizon, the upper rim mirrors the sky, and a dark ring appears near the contact line. In a clear-sky HDR the sun enters every drop as a hard point, which is not seen on a rainy day.

## 4. What was checked

- `--scene dispersion`, 4 (nested dielectrics. PASS in both RGB and spectral): the refraction angle for entry at 45° (air 28.12°, water 38.93°), total reflection (air) and transmission (water, 59.54°) on exit at 50°, the transmittance of a flat plate (1−R)/(1+R) (air 0.923, water 0.993), and the white furnace of a plate outside water, 1.000.
- Generator: that the hole in the mask coincides with the disc of the contact face, and that the rim of the cap lies at z = 0, are guaranteed by construction (the same ellipse is used). In rendering, it was confirmed by eye that the image of the lamp appears inside the drops and that the trails appear long and thin.
