# Scene files (JSON)

The CLI can take a scene as a JSON file. **The file goes in the first argument.**

```bat
GoLD_rayt.exe scenes\gold_bowl_tungsten.json
GoLD_rayt.exe scenes\gold_bowl_tungsten.json --spp 16 --width 400   :: options given after the file win over the values in the file
GoLD_rayt.exe --list-scenes                                          :: list of the built-in scenes and commands
```

In Visual Studio, put `scenes\gold_bowl_tungsten.json` in "Project Properties → Debugging → Command Arguments" and change only this file name to switch scenes. The working directory is `$(ProjectDir)`.

## Two forms

| Form | Content | Use |
|---|---|---|
| 1. Built-in scene + options | `"scene"` and `"options"` | Keep the contents of a command line in a file (`scenes/presets/`) |
| 2. Scene written as data | `"objects"` and the other keys | Write spheres, quads, meshes, lathes, lights and materials without recompiling (`scenes/*.json`) |

A form-2 file is rendered as the scene named `file`. `"scene"` and `"objects"` cannot be given together.

```jsonc
// Form 1
{ "scene": "goldbowl", "options": { "light": "tungsten", "width": 800, "spp": 256 } }
```

```jsonc
// Form 2 (excerpt)
{
  "options": { "width": 800, "spp": 256 },
  "camera": { "from": [0.02, 0.27, 0.33], "at": [0.035, 0.035, 0], "fov": 30 },
  "materials": {
    "gold": { "type": "conductor", "metal": "Au", "roughness": 0.12 },
    "bulb": { "type": "emission", "spectrum": "tungsten", "radiance": 90 }
  },
  "objects": [
    { "shape": "sphere", "center": [0, 0.03, 0], "radius": 0.028, "material": "gold" },
    { "shape": "sphere", "center": [-0.18, 0.2, 0.1], "radius": 0.03, "material": "bulb" }
  ]
}
```

## Writing conventions

- `//` comments are allowed. A misspelled key gives an **error that lists the allowed keys**. A top-level `"comment"` key is accepted and ignored.
- Relative paths inside the file are relative to the folder the program runs in (`GoLD_rayt/`), the same as on the command line.
- `"options"`: the command-line names without the `--`. A number or string is the value, `true` is a flag (`false` is left out). Example: `{ "width": 800, "denoise": true, "out": "x.png" }`. The output gamut and the PNG bit depth work the same way (`"gamut": "p3"`, `"bits": 16`), as do the exposure and tone map of the picture, the mask and depth outputs and the spectral image (`"exposure": -1, "tonemap": "clip", "mask": "statue", "depth-map": true, "spectral-image": "10,420,700"`; see [`CAMERA.md`](CAMERA.md)), and so does the random seed (`"seed": 7`: the same seed gives the same picture, noise included. It applies to path tracing only, not to the photon map or `--lvc`. If it is not given, a new seed is chosen at every run and printed on the `[Options]` line of the log).
- Units are metres. Vectors are `[x, y, z]`, colors are `[r, g, b]` (a single number is a gray).

## Key reference (form 2)

| Key | Content |
|---|---|
| `camera` | `from`, `at`, `up`, `fov` (degrees), and the lens, aberration, vignetting and motion-blur keys listed under "Camera keys" below. Default: looks from `[0, 0.5, 2.5]` toward `[0, 0, -1]` at 35° |
| `environment` | `"none"` (default), `"hdr"`, `"uniform"` (white), or `{ "type", "file", "scale" }`. The default file of `hdr` is `--env` |
| `materials` | name → material (below) |
| `lights` | array of point lights `{ "type": "point", "position", "color", "intensity" }` |
| `objects` | array of objects (below) |

**Camera keys**

- **Orthographic view.** `ortho` (the height of the view, m; 0 or omitted = perspective): parallel rays along the view direction, so the size of a thing does not depend on its distance; `fov` and `focalLength` then only set the aperture of an f-number. Not with light tracing (`--caustic-trace`, `--lvc`). See [`CAMERA.md`](CAMERA.md).
- **Lens** (either style is accepted). `aperture` (lens diameter, m; 0 = pinhole) and `focus` (focus distance, m; default: up to `at`). Or in photographic terms: `focalLength` (mm; replaces `fov`; 35 mm film equivalent, the long side of the picture is the 36 mm of the film), `fStop` (f-number; aperture diameter = focal length / f-number), `focusDistance` (distance to the plane of focus along the view direction, m) or `focusAt` (the point `[x, y, z]` to focus on), `blades` (number of aperture blades, 3 to 32; omitted or 0 = round), `bladeRotation` (degrees). If `fStop` and `focusDistance` / `focusAt` are written, they win over `aperture` and `focus`.
- **Aberrations.**
  - `distortion`: radial distortion, a number k1 or `[k1, k2]`. For the radius r from the center of the image (the half diagonal is 1), the pixel shows the direction that the ideal picture has at r(1 + k1 r² + k2 r⁴). k1 > 0 is barrel, k1 < 0 is pincushion.
  - `lateralCA`: lateral chromatic aberration. The size of the picture changes by this fraction per 100 nm. 0.004 makes red and blue differ in size by 0.8 %. Positive means that the red picture is smaller.
  - `axialCA`: axial chromatic aberration. The focal length changes by this fraction per 100 nm, so every wavelength has its own plane of focus. Together with `fStop` it colors the bokeh. 0.01 is a single lens, 0.0005 an achromat.
- **Vignetting.** `vignette` (f-number N of the opening of the lens barrel; the diameter of the opening = focal length / N; the smaller N is, the larger the opening and the weaker the darkening), `vignetteDist` (how many mm in front of the stop that opening is; default: half the focal length), `naturalVignette` (the cos⁴ law; 0 = none, 1 = as it is). The barrel opening clips the bundle of oblique rays, so the corners get darker and the bokeh balls there become cat's-eye (lemon) shaped. If the stop is closed down below the opening, the clipping disappears (`fStop` is required). `naturalVignette` works without a stop.
- **Motion blur.** `shutter` (time the shutter is open: a number of seconds, or a string such as `"1/60"`; 0 or omitted = no blur), `velocity` (camera velocity `[x, y, z]`, m/s), `pan` (speed at which the camera turns, degrees per second; a number turns to the right, `[yaw, pitch]` turns right and up. The pose at time 0 is `from` and `at`; the focus distance does not change). During the exposure the `velocity` of the objects (below) and the camera motion advance with time, and every ray sees the scene at its own time. Shadows and reflections move with the objects. The brightness does not depend on the length of the shutter (it is a time average). The pinhole cameras of `--lvc` and `--caustic-trace` stay at time 0 and do not blur. The chromatic aberrations need a different ray per wavelength, so they work in the spectral build only, and one path carries a single wavelength (more color noise).
- **Command line.** `--ortho`, `--focal-length`, `--fstop`, `--focus-dist`, `--focus-at`, `--blades`, `--blade-rotation`, `--distortion K1[,K2]`, `--lateral-ca`, `--axial-ca`, `--vignette`, `--vignette-dist`, `--natural-vignette`, `--shutter T` (seconds, or `1/60`), `--cam-velocity x,y,z` and `--pan YAW[,PITCH]` override the values in the file (they also work with built-in scenes).

**Materials (`type`)**

| type | Keys |
|---|---|
| `diffuse` | `albedo` (RGB), `reflectance` (CSV of a measured reflectance R(λ)), `fading` (measured reflectance that changes with the exposure time; see "Using spectroradiometer data" below), or `reflectanceTexture` (a `.sptex` with a reflectance spectrum per pixel; `tools/data/hsm_texture.py` makes it from a 2D spectroradiometer measurement; it is laid on a `quad`, with the `edge1` direction as u and the `edge2` direction as v). **Exactly one of these.** `texture` (PNG / JPG; laid on the uv of a mesh and multiplied with `albedo`; for the photographed texture of a scan or a 3DGS mesh: `docs/CAPTURE.md`; it combines with `albedo` only) |
| `conductor` | **Exactly one of** `metal` (`Au`, `Ag`, `Cu`: the Johnson & Christy tables), `table` (CSV of n, k from RefractiveIndex.info), `reflectance` (CSV of a measured normal-incidence reflectance R(λ); turned into n, k by Gulbrandsen's method). Also `roughness`, `anisotropy` (−1 to 1; > 0: rough along the tangent, < 0: hairline finish with grooves along the tangent), `rotation` (degrees; turns the tangent around the normal), `multipleScattering` |
| `dielectric` | `ior`, `exteriorIor` (**the index on the outside of the surface (the normal side)**, default 1 = air. It writes nested dielectrics in the same way as Mitsuba, by splitting surfaces where the index changes: for a water drop, the cap that touches air has `ior` 1.333, and the bottom disc that touches glass has `ior` 1.333 and `exteriorIor` 1.5, and the glass surface under it is not made. The BSDF uses only the ratio of inside to outside, and the dispersion is applied to the inside. Check 4 of `--scene dispersion` tests it), `roughness`, `abbe` (dispersion; 0 = none), `sellmeier` (gives the dispersion by the Sellmeier equation: a glass name `"diamond"`, `"sapphire"`, `"bk7"`, `"silica"`, or `{ "b": [b1, b2, b3], "c": [c1, c2, c3] }` (c in µm²). It cannot be combined with `abbe`. The RGB build uses the index at the d line), `distanceRoughness` (`{ "roughness", "start", "end" }`: the roughness grows to this value between the distances start and end from the y axis. It returns the slope of small distant waves that the mesh cannot resolve as a roughness. Used for the sea), `medium` (fills the inside with a scattering medium; below) |
| `thinFilm` | Interference of a thin transparent film on a smooth surface (soap bubble, oil film, coated metal). `ior` (film), `thickness` (a number in nm, or `{ "nm" }`, or a gradient `from` / `to` / `axis` / `range` / `exponent`, plus `warp` and `noise` (with `stretch`)), `coverage` (mask of where the film is), `substrate` (`"free"`, `{ "dielectric": n }`, `{ "metal": "Au" }`, `{ "table": "file.csv" }`). Theory and examples: `docs/THIN_FILM.md` |
| `layered` | A material with a transparent coat over a base (lacquered gold, car paint). `coat`: `{ "ior", "roughness", "tint", "samples" }` (index of the coat, roughness, transmittance RGB along the normal (1 = clear), number of traces per evaluation); instead of `tint`, `tintSpectrum` gives the transmittance along the normal per wavelength as a CSV (a varnish with a known absorption: the spectral build reads it at every wavelength, the RGB build uses its color under D65). `base`: either `{ "albedo" }` (diffuse), `{ "reflectanceTexture" }` (diffuse with a measured `.sptex` per texel, laid on the uv as for `diffuse`: a measured leather or wood under a lacquer) or `{ "metal" \| "table" \| "reflectance", "roughness", "anisotropy", "rotation", "multipleScattering" }` (metal). Theory and examples: `docs/GLOSS.md` |
| `skin` | Skin (`CoatedDiffuse`). Keys are listed under "Skin keys" below: `albedo`, `roughness`, `specular`, `scatterDistanceMm`, `sssScale`, `randomWalk`, `relief`, `gloss`, `pigment`, `film`, `scatter`, `lotion`, `foundation`, `primer`. Theory: `docs/SKIN.md`, `docs/SKIN_MOISTURE.md` |
| `hair` | Scattering of a single hair fiber (melanin coloring, the R, TT and TRT highlights; `docs/HAIR.md`). The absorption is given by **one of** `melanin` (`[eumelanin, pheomelanin]`; black `[8, 0]`, brown `[1.3, 0.1]`, red `[0.4, 1.4]`, blond `[0.25, 0.15]`, white (gray) hair `[0, 0]`), `color` (apparent color RGB, each in (0, 1]), `sigmaA` (absorption RGB per hair radius); the default is `melanin` [1.3, 0]. Also `eta` (default 1.55), `betaM` and `betaN` (longitudinal and azimuthal roughness, 0 to 1, default 0.3), `tilt` (tilt of the cuticle, degrees, default 2), `lobes` (draw only some of `["R", "TT", "TRT", "rest"]`; for explanation) |
| `measured` | A measured BRDF (an hpBRDF table compressed by PCA). `file` (`.hpbpca`), `alphaSample` (the GGX α used for importance sampling; it changes only the noise, not the picture; default 0.1). Polarization (`--stokes`) is supported |
| `mix` | Two materials of the file mixed over the surface: `a` and `b` (their names), and `mask` (an image laid on the uv; its red channel, as written in the file, is the weight of `a`: gray 128 = one half) or `weight` (a constant weight of `a`, default 0.5). Each hit behaves entirely as `a` or as `b`, chosen at random with that weight and fixed for the hit (as pbrt-v4's MixMaterial), so the surface reflects the mix on average; any reflecting material can be mixed (metal leaf on embossed leather, worn paint), but not another `mix`, a light, a `skin` with subsurface scattering or a `dielectric` with a `medium`. `Materials/MixMaterial.hpp` |
| `emission` | `spectrum`, `radiance`, `twoSided` |

**Skin keys**

- `albedo` (default 0.80, 0.52, 0.42), `roughness` (default 0.5), `specular` (weight of the sebum coat, default 0.4).
- `scatterDistanceMm` (RGB, mm; e.g. `[3.67, 1.37, 0.68]`): writing it turns on **subsurface scattering** (Burley's profile; `sssScale` is a multiplier).
- `randomWalk` (optional; `true`, or `{ "g": 0.8 }`): follows the subsurface light by **walking it inside the object** instead of using the profile. It needs a closed shape. Thick places keep the same color, and thin places let light through according to the thickness. `g` is the mean cosine of the phase function, default 0. `docs/SSS.md` section 12.
- `relief` (optional; bumps of pores and furrows; it bends the normal only): `cellMm`, `furrowDepthMm`, `furrowWidth`, `stretch`, `coarseScale`, `coarseWeight`, `poreDepthMm`, `poreRadiusMm`, `poreFraction`, `poreSizeVariation`, `poreDensityVariation`, `poreDensityFreq`, `poreDarkness` (darkens the pores, 0 to 1), `seed`.
- `wrinkles` inside `relief` (**furrows and pores made by graph synthesis**): `tileMm` (default 24), `resolution` (2048), `steps` (300), `depthScale`, `seed`, `levels` (optional; per level `spacingMm`, `angleDeg`, `base`, `distance`, `continuation`, `similarity`, `crossing`, `deposit`, `width`, `poreWidth`, `poreMin`, `poreBlend`, `poreSkew`, `wobble`, `wobbleFreq`, `depthUm`, `bulge`, `darkPores`). When it is written, the Worley furrows and pores are not used.
- `gloss` (optional; oily and dry places): `variation` (0 to 0.95), `frequency` (number of patches per metre), `seed`.
- `pigment` (optional; color variation): `freckleCellMm`, `freckleFraction`, `freckleRadiusMm`, `freckleStrength`, `clusterFreq`, `blotchFreq`, `blotchStrength`, `rednessFreq`, `rednessStrength`, `seed`.
- `film` (optional; **a liquid film that collects in the furrows**; needs `relief`): `thicknessUm` (mean thickness, µm; 0 = no film), `ior` (water 1.33, glycerin 1.47, oil about 1.45), `roughness` (roughness of the film surface, default 0.02), `softnessUm` (rounding of the film edge, default 1.5). The specular reflection follows the normal of the top surface of the film.
- `scatter` (optional; **white local scattering of dry skin (scales)**): `rho0` (reflectance where there are scales), `frequency`, `threshold` (raising it reduces the scales), `softness`, `hydration` (0 to 1), `kappaH` (the fraction by which hydration reduces the scattering), `gapUm` (depth to which a liquid fills the gaps), `seed`. If there is a `film`, the index matching reduces the scattering.
- `lotion` (optional; **the film and the hydration `timeS` seconds after a lotion was applied**; needs `relief`; not together with `film`):
  - `doseUm` (thickness of the applied liquid, default 10), `glycerin` (volume fraction of the humectant that remains, 0.05), `timeS`, `waterTauS` (evaporation of water, 120), `hydrationPeak` (0.8), `hydrationRiseS` (60), `hydrationDecayS` (1800), `waterIor`, `humectantIor`, `waterRoughness`, `residueRoughness`.
  - `bulk`: for the hydration H, multiplies the scatter distance of the subsurface scattering by 1 + bulk · H (default 0). The BSSRDF transmission is not reliable in thin places: `docs/SKIN_MOISTURE.md` section 5.4. If `scatter` is present, the hydration of the lotion also overrides its `hydration`.
  - `oil` (volume fraction of oil in the liquid; 0 = a watery lotion. Above 0 it is an **emulsion**: the milky layer of oil droplets disappears with time and an oil film remains), `oilIor` (1.45), `dropletUm` (droplet diameter, default 1), `coalesceTauS` (coalescence of the droplets, 600), `dropletG` (0.9), `occlusion` (fraction by which the merged oil film stops the loss of water from the stratum corneum, 0 to 1, default 0.3; 0 = the oil does not keep the water), `occlusionOilUm` (oil film thickness at which 63 % of it is reached, default 1).
  - **Layering** uses `layers` (a list of products; e.g. `[{ "atS": 0, "doseUm": 10 }, { "atS": 180, "doseUm": 10, "oil": 0.08 }]` = a lotion, then an emulsion 180 s later). Per layer: `atS` (time of application, on the same clock as `timeS`; a layer later than `timeS` is not applied yet), `doseUm`, `glycerin`, `hydrationPeak`, `oil`, `oilIor`, `dropletUm`, `coalesceTauS`, `dropletG`. When `layers` is written, these keys are not written directly under `lotion`.
- `foundation` (optional; **a layer of foundation** on top): `preset` (`dewy1`, `dewy2`, `matte1`, `matte2`: the values of the four products of the paper; the other keys override the preset), `thickness` (optical thickness; 0 = no layer; 0.25 to 0.5 is a thin coat, 4 nearly hides what is underneath), `diffuserFraction` (fraction of spherical particles, 0 to 1; the closer to 1, the more matte), `plateletRoughness` (spread of the orientation of the platelet particles, 0.001 to 1; the smaller, the sharper the highlight), `plateletTiltDeg`, `diffuserAlbedo` and `plateletAlbedo` (particle colors RGB), `g1`, `g2` and `lobeWeight` (how the spheres scatter), `leveling` (how much the layer levels the skin relief, 0 to 1, default 1. At 1 the furrows and pores do not show in the shine of the layer; at 0 the layer follows the relief. When there is a `primer`, the default is 0 and the surface of the layer is set by the primer). `docs/COSMETICS.md`.
- `primer` (optional; **the base under the foundation**: a film that fills the furrows and pores, with the foundation on top; needs `relief`; not together with `film` or `lotion`): `thicknessUm` (mean thickness, µm; 0 = none; the deepest places fill first), `ior` (default 1.40, silicone), `roughness` (roughness of the film surface, default 0.15), `softnessUm` (rounding of the edge, default 3). Without a foundation, the gloss of the skin follows the top surface of the film. `color` (color of the primer: the color RGB of an infinitely thick film; writing it makes a **colored primer**), `opacity` (opacity of the colored layer = scattering coefficient × thickness, default 0.5; 0 = transparent, above 1 it nearly hides the skin; needs `color`). `docs/COSMETICS.md` sections 4.1 and 4.2.

`medium` of `dielectric`: fills the inside of a closed shape (a sphere, an OBJ with correct orientation) with a scattering medium. Light that enters through the glass boundary scatters inside and then comes out (wax, milk, jade).

```jsonc
"jade": { "type": "dielectric", "roughness": 0.2,
          "medium": { "sigmaT": [250, 150, 250], "albedo": [0.85, 0.998, 0.95], "g": 0.3 } }
```

- `sigmaT`: extinction coefficient (1/m). The larger, the denser. `albedo`: the fraction of light that survives scattering (1 = no absorption). Both are a single number (gray) or `[r, g, b]` (colored). `g`: anisotropy of the phase function (−1 to 1; positive = forward scattering). Adding `g2` and `weight` gives a two-lobe (forward + backward) phase function `weight·HG(g) + (1−weight)·HG(g2)`. The mean cosine is `weight·g + (1−weight)·g2`. With the same σs′ = σs(1−mean cosine), a thick medium looks almost the same, and the differences appear in thin parts, in backlight and for light near the surface.
- `phaseTable`: gives the phase function as a table (not together with `g`, `g2`, `weight`). The CSV has rows "angle (degrees), value", equally spaced from 0° to 180°; `#` is a comment. The values are a density per steradian (the scale is normalized automatically). Mie scattering of cloud droplets is `assets/Cloud/mie_reff10.csv` (made by `tools/studies/mie_phase.py`; `docs/CLOUD.md`).
- `--media-nee` (`"media-nee": true` in `options`): samples a light at every collision in the medium and passes straight through boundaries of index 1. It is needed when a small light such as the sun lights a cloud or fog (`scenes/cloud_mie.json`).
- For a colored medium, a single path chooses the distance with one representative wavelength (one channel in the RGB build) and weights the other wavelengths by the ratio of the probability densities (`--scene volume` checks it against a white furnace and against gray). A colored medium has more noise than a gray one at the same spp.
- Adding `noise` (optional) makes a medium whose density varies in space (cloud, fog). `sigmaT` is the value where the density is maximal.

```jsonc
"cloud": { "type": "dielectric", "ior": 1.0,   // ior 1: the boundary is invisible
           "medium": { "sigmaT": 420, "albedo": 0.9995, "g": 0.6,
                       "noise": { "frequency": 22, "octaves": 5, "lo": 0.42, "hi": 0.62, "seed": 7,
                                  "maskCenter": [0, 0, 0], "maskRadius": 0.1, "maskStart": 0.45 } } }
```

  `frequency`: number of patterns per metre in the first layer (the size of a cloud is about 1/frequency). `octaves`: number of layers. `lo` / `hi`: the density is 0 where the noise is below lo and 1 where it is above hi (a narrow width gives a sharp outline, a wide one a haze). `seed`: another cloud. `maskCenter` / `maskRadius` / `maskStart`: sets the density to 0 toward the surface of a sphere (keeps the outline of a cloud inside a sphere away from the sphere). The density is determined in world coordinates.
  The mechanism is null collisions (the distance is chosen with the maximum density, and with probability 1−q nothing happens).
- The upper limit of the path length is `--volume-cap`.

Spectroradiometer measurement data can be converted with `tools/data/sr5100.py` into two-column CSV files (`wavelength_nm, value`) (below).

`spectrum`: a preset name `"d65"`, `"tungsten"`, `"candle"`, `"led"` (also written `{ "preset": "d65" }`), or `{ "blackbody": 2856 }`, `{ "rgb": [1, 1, 1] }`, `{ "file": "x.csv" }` (a table of wavelength and value), `{ "mix": [ { "file": "led1.csv", "weight": 0.4 }, { "file": "led2.csv", "weight": 1 } ] }` (the sum of several sources lit together, e.g. the color mixing of an LED cube. The weights are 0 or more and only their ratio matters; the brightness is set by `radiance`).

`--wb` (as an option, or `"options": { "wb": true }`) adapts the film to the **light of the first `emission` object** (CAT16; it removes the color cast of the light). Point lights and the environment map are not used for the adaptation.

**Objects (`shape`)**

| shape | Keys |
|---|---|
| `sphere` | `center`, `radius` |
| `quad` | `corner`, `edge1`, `edge2` (faces the direction of `cross(edge1, edge2)`) |
| `mesh` | `file` (OBJ, PLY or STL, decided by the extension. PLY is read in ascii and binary; vertex normals, uv and the per-face `texcoord` of MeshLab are used; vertex colors are not read. STL is read in ascii and binary; vertices at the same position are joined and smooth normals are made from the orientation (counter-clockwise is outside); normals, colors and uv in the file are not used. STL units are often mm and Z is up, so match them with `scale` and `rotate`), `cutout` (an RGBA image laid on the uv; where the alpha is below 1/2 there is a hole that every ray passes through. Used to cut out the glass surface under a water drop: `docs/WATERDROPS.md`. It uses its own BVH), `scale`, `rotate` (`{ "axis", "degrees" }`), `translate`. The order is scale → rotate → translate. How to prepare a scan or a 3DGS mesh: `docs/CAPTURE.md` |
| `lathe` | `profile` (an array of `[radius, height]`, or `{ "shell": { "radius", "wall", "endDegrees", "steps" } }`: the wall of a bowl), `segments`, `translate` |
| `ocean` | A sea surface (made from a wave spectrum; `docs/WATER.md`). A circular mesh around the origin. `windSpeed` (m/s, 1 to 40, default 8), `windDirection` (degrees, the direction the waves travel, `atan2(z, x)`), `spread` (spread s of the direction), `peakEnhance` (γ of JONSWAP, 1 = Pierson–Moskowitz), `chop` (0 to 1, sharpness of the crests), `time` (seconds), `minWavelength` (m), `seed`, `rInner` and `rOuter` (m), `rings` and `sectors` (fineness of the mesh) |
| `splats` | Draws a **3D Gaussian Splatting point cloud** (`point_cloud.ply`) as emitting and occluding particles (`docs/CAPTURE.md` stage 2). `file`, `scale` / `rotate` / `translate` (as for a mesh; the y-down of COLMAP is `"rotate": { "axis": [1, 0, 0], "degrees": 180 }`), `alphaMin` (responses below this are ignored, default 1/255), `intensity` (multiplier of the radiance), `linearize` (turn the sRGB values back to linear, default true), `degree` (order of the spherical harmonics to use, 0 to 3), `keep` and `remove` (a box `{ "min", "max" }`; keeps only / drops the particles whose center is inside it; used to remove the particles where an object is to be placed). **Do not write `material`** (the color is the capture's). One per scene. Path tracing only (the photon map, manifold NEE and LVC-BPT do not see it). The particles are not a target of next-event estimation (a diffuse surface receives the light of the particles through BSDF samples only) |
| `hair` | Procedurally made hair (`docs/HAIR.md`). `style`: `"lock"` (a lock hanging from above) or `"head"` (hair flowing over a spherical head). Common: `strands`, `segments`, `length` (m), `rootRadius` and `tipRadius` (radius of the fiber, m; a real one is about 0.00005), `clumps` (number of clumps), `seed`. `lock`: `root` (position of the root), `width`, `tipWidth` (width at the tip / width at the root), `waveAmp` and `waveLength`, `twist`, `frizz`, `flyaway` (fraction of fibers that leave the lock). `head`: `center`, `radius` (the head sphere), `droop`, `outward`, `sway` and `swayFreq`, `headFrizz`, `clumpStrength`, `maxPolar` (degrees), `hairline`. The material is `hair` |

Spheres, quads, meshes and lathes (`lathe`) can have a `velocity` (`[x, y, z]`, m/s). During the `shutter` of the camera, the object moves in a straight line at this velocity (time 0 is the written position). A light (a sphere or quad with an `emission` material whose `light` is true) cannot be moved (the light sampler knows only the position at time 0). If the light is only for the picture, write `"light": false`. The sea, hair and splats cannot be moved.

Every object has a `material` (required), except `splats`. A sphere or quad with an `emission` material is by default a **light** (a target of next-event estimation). With `"light": false` it is not a light (for uses such as a faint light around everything). An `emission` material can be used on a sphere or a quad only.

## What was checked

The three files in `scenes/` (`gold_spheres`, `gold_bowl_tungsten`, `studio`) are transcriptions of the code versions `gold`, `goldbowl --light tungsten` and `studio`. In both the RGB and the spectral build, the mean color of four renders agrees with the code version within three standard errors, and the pictures cannot be told apart.

## Adding a material or a shape

Change three places.

1. `include/Cli/SceneDescription.hpp`: add the item to the type.
2. Read and validate the JSON (add the keys to the list of allowed keys): the shapes and the scene in `src/Cli/SceneFile.cpp`, the materials in `src/Cli/SceneFileMaterials.cpp`.
3. `src/Cli/Scenes/FileScene.cpp`: make the engine object from that type.

## Using spectroradiometer (SR-5100) data

`tools/data/sr5100.py` reads the CSV exports of the SR-5100 (two formats, spots and images) and converts them to the two-column CSV that this renderer reads.

```bat
python tools\data\sr5100.py info measurement.csv                                :: list of the contents (number of points, peak, measurement distance, integration time)
python tools\data\sr5100.py spots measurement.csv --out output_folder           :: one "wavelength_nm, radiance" CSV per spot
python tools\data\sr5100.py roi image.csv --peak 2 --out light.csv              :: mean spectrum around the brightest pixel of the image (5x5)
python tools\data\sr5100.py reflectance sample.csv white.csv --out reflectance.csv   :: sample / white (wavelengths where the white is weak are not output)
python tools\data\hsm_texture.py sample.hsm --white white.hsm --out cloth.sptex --preview cloth.png   :: .hsm (spectral image) to a texture of per-pixel reflectance spectra (compressed to 6 principal components, 420 to 700 nm)
```

- The spectrum of an LED is given with `"spectrum": { "file": "light.csv" }` (the unit of the values does not matter; the brightness is normalized to 1 when it is used).
- A measured reflectance can be used as a metal with `"type": "conductor", "reflectance": "reflectance.csv"`, and as a diffuse surface (cloth etc.) with `"type": "diffuse", "reflectance": "reflectance.csv"`. For a diffuse surface, the spectral build uses the curve as it is, and the RGB build turns it into its color under D65 (outside the measured range the end values are extended; a value above 1 gives a warning).
- **Reflectance that changes with the exposure time (fading)** is used with `"type": "diffuse", "fading": { "series": "series folder", "color": "color", "hours": 45 }`. The series folder contains `hours.csv` (`round,hours`: the exposure time of each round, in increasing order) and `r1/color.csv`, `r2/color.csv`, … (two-column reflectance CSV files). If the time is between two rounds, the two measured curves are mixed linearly per wavelength (no extrapolation: before the first round and after the last round the end round is kept, with a warning). `--fade-hours H` overrides the `hours` of the file for all `fading` materials (used to make a sequence of pictures that differ only in time).
  - **Exposure that depends on the place:** between `"unexposedHours"` (time at places that received nothing, default 0) and `"hours"` (time at places that received everything), the time is set by an **exposure fraction f (0 to 1)**. f is the product of the following factors (only those that are written; the others are 1): `"map"` (a gray image; white = received everything; the top left of the image is uv = (0, 0)), `"mat"` (`{ "inset": 0.07 | [left, right, bottom, top], "soft": 0.003 }`: a band at the edge of the uv that was covered, like the mat of a picture frame; uv takes the length of the edge as 1), `"gradient"` (`{ "bright": [x,y,z], "dark": [x,y,z], "exponent": 1 }`: 1 at bright, 0 at dark, linear along the axis in between; from the window side to the back, etc. It works on any shape), `"noise"` (`{ "amount": 0 to 1, "frequency": number of patterns per metre, "octaves", "seed" }`: unevenness of the fading. It works on any shape). `map` and `mat` use the uv of the object (`quad`, a mesh with uv, `lathe`). `tools/data/fading_refangle.py` converts the table of the lab's fading experiment into this form. The check is 5 of `--scene metals`.
- **Fading from a color alone (dyed paper)**: a color without a measurement can also be faded, from an RGB color (or an image). `"fading": { "rgb": [0.8, 0.2, 0.15], "model": "dye_table.csv", "hours": 90, "fastness": 1 }`, or `"image": "picture.png"` instead of `"rgb"` (an sRGB image; laid on the uv of a `quad` etc.).
  - Mechanism: the color is turned into a smooth spectrum by the engine's RGB-to-spectrum conversion, and the dye part of the Kubelka–Munk K/S (`max(0, K/S(color) − K/S(paper))`) is multiplied by the concentration c(t): `K/S(t) = K/S(paper) + c(t) × dye part`. A color brighter than the paper (white) does not change. The RGB build is an approximation that treats each channel as a reflectance and does the same calculation (the difference from the spectral build is ΔE00 1.6 on average and 4.6 at most over the 9 colors of the dyed paper).
  - `model` is a table of `hours,c` (c(0) = 1; `tools/data/fading_dye_model.py` makes it from the measured series of dyed paper, as the median of the c fitted per color). Outside the measured range the end value is kept (no extrapolation). The reflectance of the paper used for the table is `paper` (default 0.9; use the same value as the table). `fastness` is a multiplier of the time (2 = a dye that fades twice as fast). `hours`, `unexposedHours`, `map`, `mat`, `gradient` and `noise` are the same as in the `fading` of a series.
  - **Accuracy (dyed paper, predicting one color with that color left out; the lab fading experiment):** the color predicted from the fresh RGB alone is ΔE00 4.9 on average from the measurement (6 to 7 in the last round). Without fading it is 9.5 (13.8 in the last round). Black, blue and green fit badly (mixed dyes, dyes that hardly fade). **It does not fit the dyes of pigments and origami paper** (for pigments, not fading is closer). For saturated colors, the RGB-to-spectrum conversion cuts negative components at 0, which adds another ΔE00 of 2 to 4. The speed of fading cannot be guessed from the color (weighting by hue did not improve it), so the median over all colors is used.
- The individual spectra of an LED cube can be made into this renderer's CSV files with `tools/data/ledcube_export.py` (spectral radiance at maximum output; weight = output setting). Listing them in a `mix` with weights gives the color mixing of the cube.
- A spectral image (`.hsm`) becomes a texture of per-pixel reflectance spectra with `"type": "diffuse", "reflectanceTexture": "cloth.sptex"` (laid on a `quad`. The spectral build uses the spectrum of the pixel, and the RGB build the D65 color of each pixel. A measured cloth rendered under the same light as the measurement agreed with the measured image to chromaticity 0.001 and a median relative pixel difference of 3 %).
- Measurement data are experimental data, so they are not put in the repository (`assets/Metal_measured/` is outside Git).
