# Camera and film output

What the camera simulates, what the written files contain, and what is verified. The keys of the scene file are listed in
[`SCENE_FILE.md`](SCENE_FILE.md) ("Camera keys"); the command-line options in `--help`.

## 1. The camera (`Render/Camera.hpp`)

| Part | What it does | Options |
|---|---|---|
| Projection | Perspective, from `from` / `at` / `up` and a vertical field of view, or a focal length on 35 mm film | `fov`, `--focal-length` |
| Thin lens | Defocus from an aperture (f-number), a plane of focus, a round or polygonal aperture (the shape of the bokeh) | `--fstop`, `--focus-dist`, `--focus-at`, `--blades` |
| Aberrations | Radial distortion r' = r (1 + k1 r² + k2 r⁴); lateral and axial chromatic aberration (spectral build only: they need the wavelength of the path) | `--distortion`, `--lateral-ca`, `--axial-ca` |
| Vignetting | The lens barrel cuts oblique bundles (dark corners, cat's-eye bokeh); the cos⁴ law | `--vignette`, `--natural-vignette` |
| Exposure time | Moving objects and a moving or panning camera blur over the shutter time | `--shutter`, `--cam-velocity`, `--pan` |

The lens of the scene file and of the command line are merged in `Cli/RenderPipeline.cpp` (`ResolveLens`, `BuildCamera`).

## 2. The written files (`Render/Film`, `IO/ImageWriter`)

| Output | Content |
|---|---|
| `.pfm`, `.hdr` | Linear radiance in the output gamut. `.pfm` keeps negative values and small channels: use it for measurements |
| `.png`, `.bmp`, `.jpg` | Display values: gamut (`--gamut`), exposure (`--exposure EV`), tone map (`--tonemap`), sRGB curve, 8 or 16 bit (`--bits`) |
| `--tonemap reinhard` | The default: x / (1 + x) per channel. Nothing clips, highlights are compressed |
| `--tonemap clip` | Linear up to white, clipped above it: brightness ratios below white are kept (stimuli). The log gives the share of pixels above white and how many stops the brightest is above it, i.e. the `--exposure` that avoids clipping |
| `--mask A,B` → `_mask.png` | The share of each pixel covered by objects of the named scene-file materials, as the camera first sees them. n × n camera rays per pixel (n = √spp, up to 16) through the lens and the shutter, so the edge is anti-aliased, defocused and blurred like the picture. Linear gray (no tone curve), one channel |
| `--depth-map` → `_depth.pfm` | Distance of the first surface along the view axis (not along the ray), through each pixel centre with the pinhole view; 0 = no surface |

The mask and the depth come from a separate pass after the render (`Render/MaskDepth.cpp`): it works with every render mode and leaves the
picture of a `--seed` run unchanged. Glass, water and mirrors are surfaces; an index-matched boundary (the surface of a fog) is not; a
translucent layer covers a pixel in proportion to its opacity.

## 3. Verified

`--scene film` (RGB and spectral builds):
- Display values: Reinhard and clip against their formulas, 8- and 16-bit PNG read back, `.pfm` unchanged by the display settings, the clip report.
- Mask and depth against closed forms: a sphere's silhouette area (0.01 %), its depth on the axis, a perpendicular wall at the same depth in every corner, a translucent layer, an index-matched boundary against glass, the files read back.
- On an Erato view, the mask agrees with an earlier mask made by rendering the statue black on white, once that render's tone curve is applied (area −0.07 %).

The lens, the aberrations, the vignetting and the motion blur have no self-check of their own yet.

## 4. Limits

- Light tracing (`--caustic-trace`, `--lvc`) connects light paths to a pinhole camera at time 0: its caustics are sharp, undistorted and not motion-blurred.
- The mask selects materials, not objects: two objects with the same material are masked together (give them copies of the material). A mix material is selected by the name of the mix.
- The depth has one ray per pixel (no anti-aliasing): an edge pixel has the depth of the surface at its centre. Gaussian splats are not surfaces for either output.
