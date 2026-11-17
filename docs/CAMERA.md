# Camera and film output

What the camera simulates, what the written files contain, and what is verified. The keys of the scene file are listed in
[`SCENE_FILE.md`](SCENE_FILE.md) ("Camera keys"); the command-line options in `--help`.

## 1. The camera (`Render/Camera.hpp`)

| Part | What it does | Options |
|---|---|---|
| Projection | Perspective, from `from` / `at` / `up` and a vertical field of view, or a focal length on 35 mm film | `fov`, `--focal-length` |
| Orthographic view | Parallel rays along the view axis from a rectangle H high, centred on `from`: the size of a thing does not depend on its distance (patch stimuli, measurement geometry). An aperture still defocuses (a telecentric lens) | `ortho`, `--ortho H` |
| Thin lens | Defocus from an aperture (f-number), a plane of focus, a round or polygonal aperture (the shape of the bokeh) | `--fstop`, `--focus-dist`, `--focus-at`, `--blades` |
| Aberrations | Radial distortion r' = r (1 + k1 r² + k2 r⁴); lateral and axial chromatic aberration (spectral build only: they need the wavelength of the path) | `--distortion`, `--lateral-ca`, `--axial-ca` |
| Vignetting | The lens barrel cuts oblique bundles (dark corners, cat's-eye bokeh); the cos⁴ law | `--vignette`, `--natural-vignette` |
| Diffraction | The Fraunhofer pattern of the aperture (`Render/Diffraction`): an Airy disc for a round one, a star for blades (n spikes for an even n, 2n for an odd one), its size λN, so red reaches farther than blue. Integrated over the pixels (36 mm over the long side of the picture) and convolved with every picture, `.pfm` included. At the real size (`1`) only very bright points show a faint star; `S` > 1 enlarges the pattern S times | `--diffraction S`, `diffraction` |
| Exposure time | Moving objects and a moving or panning camera blur over the shutter time | `--shutter`, `--cam-velocity`, `--pan` |

The lens of the scene file and of the command line are merged in `Cli/RenderPipeline.cpp` (`ResolveLens`, `BuildCamera`).

### A real lens (`--lens FILE`, camera `"lens"`; `Render/LensSystem`)

Instead of the thin lens, the prescription of a real lens: its spherical surfaces and its aperture stop, one line per surface from the
front, `radius thickness n_d diameter [Abbe number]` in millimetres (a radius of 0 is the stop; the layout of pbrt's lens files). A camera
ray starts on the film, crosses the surfaces from the back to the front by Snell's law and is lost at a rim, the stop or total internal
reflection (Kolb, Mitchell and Hanrahan 1995). What the thin lens takes as parameters follows from the prescription: distortion, field
curvature, natural and optical vignetting, cat's-eye bokeh at the edges, focus breathing, and, in the spectral build with Abbe numbers,
axial and lateral chromatic aberration.

- **Set-up:** the focal length, the principal planes and the entrance pupil come from rays traced a micrometre off the axis. The front
  principal plane sits at the camera position; the film (36 mm on the long side) is placed by Gauss's lens equation for the focus distance;
  `--fstop` closes the stop (no wider than the lens opens), `--blades` makes it a polygon. The lens and the film decide the field of view
  (`fov` and `focalLength` are not used), and the thin lens's aberration and vignetting settings are not used.
- **Sampling:** for 64 bands of film radius, the box on the rear vertex's plane through which rays get through (from a 128 × 128 grid of
  trial rays) is stored; a camera ray aims at a point of that box, turned to the film point's angle. Its weight is the film irradiance
  factor cos⁴θ × box area / distance², divided by its value at the film centre, so the centre of the picture keeps the scene's radiance and
  the corners darken as the lens makes them.
- **Lenses in `assets/Lenses/`:** `dgauss_50mm` (double-Gauss f/2, Tronnier), `wide_22mm` (Nakamura, f/2.8), `fisheye_10mm` (Muller,
  f/4), all from W. J. Smith, *Modern Lens Design*, as scaled in pbrt's lens files (they give n_d only, so the glass is not dispersive), and
  `singlet_bk7_100mm`, an uncorrected N-BK7 lens with Abbe number for the chromatic aberration. The set-up prints the focal length, the
  f-number and the field: 50.4 mm f/2.02, 22.0 mm f/2.79 and 9.99 mm f/3.99 for the three.

## 2. The written files (`Render/Film`, `IO/ImageWriter`)

| Output | Content |
|---|---|
| `.pfm`, `.hdr` | Linear radiance in the output gamut. `.pfm` keeps negative values and small channels: use it for measurements |
| `.png`, `.bmp`, `.jpg` | Display values: gamut (`--gamut`), exposure (`--exposure EV`), tone map (`--tonemap`), sRGB curve, 8 or 16 bit (`--bits`) |
| `--tonemap reinhard` | The default: x / (1 + x) per channel. Nothing clips, highlights are compressed |
| `--tonemap clip` | Linear up to white, clipped above it: brightness ratios below white are kept (stimuli). The log gives the share of pixels above white and how many stops the brightest is above it, i.e. the `--exposure` that avoids clipping |
| `--mask A,B` → `_mask.png` | The share of each pixel covered by objects of the named scene-file materials, as the camera first sees them. n × n camera rays per pixel (n = √spp, up to 16) through the lens and the shutter, so the edge is anti-aliased, defocused and blurred like the picture. Linear gray (no tone curve), one channel |
| `--depth-map` → `_depth.pfm` | Distance of the first surface along the view axis (not along the ray), through each pixel centre with the pinhole view; 0 = no surface |
| `--spectral-image STEP[,FIRST,LAST]` → `_spectral.img` + `.hdr` | The mean spectral radiance of every pixel in bands of STEP nm (default 380–780 nm; spectral build), as an ENVI image (float32, band sequential; the header lists the band centres). Units of the renderer: a surface the RGB image shows as white under D65 has the D65 spectrum (100 at 560 nm). Made from the same paths and weights as the RGB picture, which it leaves unchanged |

A band's value is the mean radiance over the band (a box of STEP nm); to compare with a measurement at 1 nm, average the measurement over the
same bands. Reading it in Python: `np.fromfile(stem + "_spectral.img", "<f4").reshape(bands, height, width)`, or `read_envi` of
`tools/lib/hdr_io.py`; `tools/lib/colorimetry.py` turns it into XYZ (any CIE observer, the bands integrated as boxes), CIELAB, CIEDE2000
and display RGB with colour-science, so colour science on a render can be done in Python after it. The renderer's RGB uses an analytic
fit of the CIE 1931 functions: against the tabulated ones its colours differ by CIEDE2000 0.4 (median) to 1.2 on the `spd` scene.

The mask and the depth come from a separate pass after the render (`Render/MaskDepth.cpp`): it works with every render mode and leaves the
picture of a `--seed` run unchanged. Glass, water and mirrors are surfaces; an index-matched boundary (the surface of a fog) is not; a
translucent layer covers a pixel in proportion to its opacity.

### The camera's spectral sensitivity (`--camera-response FILE`; `Color/CameraResponse`)

The film normally sees as the CIE 1931 observer: an ideal camera that meets the Luther condition. A real camera's R, G, B sensitivities are
not a linear mix of the colour matching functions, so spectra that match for the eye can differ for the camera (and the other way round),
and its colours drift from the eye's most under narrow-band light. With `--camera-response` (spectral build) the film does what a camera
does: raw = the light weighted with its sensitivities, white balance (D65 becomes (1, 1, 1)), and a 3 × 3 colour matrix to sRGB fitted by
least squares on 216 smooth reflectances under D65, white kept white. The three linear steps make one weight per wavelength in place of the
CIE one. `<out>_camraw.pfm` holds the camera's raw RGB before white balance and matrix (G of D65 = 1), for comparing with a real camera's
raw file (both linear; scale them by a white or grey patch).

- **Data:** a CSV of `wavelength_nm,R,G,B`. The 28 cameras of Jiang et al., "What is the Space of Spectral Sensitivity Functions for
  Digital Color Cameras?" (WACV 2013; 400–720 nm in 10 nm steps) are licensed CC BY-NC-SA 4.0 and so are not in this repository: download
  `camspec_database.txt` from http://www.gujinwei.org/research/camspec/db.html into `GoLD_rayt/assets/CameraResponse/` (git-ignored) and
  run `python tools/data/camspec_to_csv.py GoLD_rayt/assets/CameraResponse/camspec_database.txt` (one CSV per camera). A measured camera
  goes into a CSV of the same form.
- **Example:** the three SR-5100 cloth swatches under D65 and under the lab's LED, the ideal film against three cameras (CIEDE2000 of
  the mean colours): under D65 0.1–1.0 (the matrices are fitted there), under the LED 0.5–3.9 on the cloth and 2.4–6.7 on the grey table
  (the phone camera worst). The grey's error under the LED comes from the white balance, which stays at D65.

## 3. Verified

`--scene camera` (RGB and spectral builds): `project()` undoes the ray through a pixel centre; with an aperture, the rays of a film position
meet on the plane of focus (perspective and orthographic); the orthographic rays are parallel to the axis and start on the H × W rectangle;
a sphere covers the same pixel area, π R² over the pixel area, at 3 m and at 9 m, and its depth follows the sphere.

A real lens, also in `--scene camera`: the singlet's focal length matches the thick-lens formula and its principal planes lie symmetric;
focused at 1 m, rays near the axis from the axis point meet on the film (0.06 µm) and farther out they cross it nearer the lens by an amount
that grows as the square of the height in the stop (third-order spherical aberration: 2.76 and 11.0 µm at 0.4 and 0.8 mm); the f-number wide
open is the focal length over the stop's diameter and the stop closes to f/8; in the spectral build the F and C lines' focal lengths differ by
f / V_d (0.988 of it); in a camera the double-Gauss's weight averages 1 at the film centre and 0.21 in the corner, the picture is upright and
a wall 10 m away has the depth 10 in the middle and the corners.

Diffraction, also in `--scene camera`: the transform of the aperture (disc, 6 and 7 blades) against a direct sum over it (7 × 10⁻⁵ of the
area); the Airy disc holds 83.8 % of the light inside its first dark ring and disc and heptagon nearly all of it within 60 λN; an enlarged
pattern has 6 spikes for 6 blades and 14 for 7, across the edges (within 1°); a uniform picture stays uniform and one bright pixel spreads
exactly as the kernel, not mirrored. The kernel of a real-size pattern holds 99.6–99.9 % of the light within its 511 px (the rest is the
estimate's error and the far tail; it is renormalized to 1).

`--scene film` (RGB and spectral builds):
- Display values: Reinhard and clip against their formulas, 8- and 16-bit PNG read back, `.pfm` unchanged by the display settings, the clip report.
- Mask and depth against closed forms: a sphere's silhouette area (0.01 %), its depth on the axis, a perpendicular wall at the same depth in every corner, a translucent layer, an index-matched boundary against glass, the files read back.
- On an Erato view, the mask agrees with an earlier mask made by rendering the statue black on white, once that render's tone curve is applied (area −0.07 %).
- Camera response (spectral build): the CIE film's own curves given as a camera give the identity white balance and matrix; Gaussian sensitivities map D65 to (1, 1, 1) and a grey to the same grey, and the raw file's D65 equals the integral of the curves; a light plus a metameric black (no component along the CIE curves) keeps its CIE colour exactly but changes through the camera (by 1 %).
- Spectral image (spectral build): a uniform sky gives 0.5 × the band means of D65 (largest deviation 0.74 % at a fixed seed, the sampling noise); weighted with the sensor curves, the bands of a red ball give its RGB picture (0.02 %); the ENVI files read back.
- On measured data: the three cloth swatches of the SR-5100 (spectral reflectance textures) under D65, seen from above with the orthographic camera. Their reflectance recovered from the spectral image (0.3 × swatch / gray table) matches the textures within 0.001 RMS over 420–700 nm in 10 nm bands.

The aberrations, the polygonal aperture, the vignetting and the motion blur have no self-check of their own yet.

## 4. Limits

- Light tracing (`--caustic-trace`, `--lvc`) connects light paths to a pinhole camera at time 0: its caustics are sharp, undistorted and not motion-blurred. It refuses the orthographic view (its weights are those of a pinhole).
- Diffraction is the pattern of a point in focus, laid over the rendered picture (whose defocus comes from the lens samples); it needs an aperture. Each of X, Y and Z is spread as for a D65-like spectrum within it, since the picture keeps no spectrum per pixel. It is not applied to the AOVs, the mask, the depth or the spectral image, and the denoiser works on the picture before it. The kernel takes 7–10 s to compute.
- The camera response balances to D65 only (it refuses `--wb`; a camera would balance to the scene's light), its matrix is fitted on smooth reflectances under D65, and the sensor stays ideal: no noise, no colour filter mosaic, no saturation.
- A real lens has spherical surfaces only (no aspheres), no reflections between them (no ghosts or flare) and no absorption in the glass; a radius of 0 is always the stop (a flat glass surface cannot be written). It refuses light tracing and the orthographic view. The RGB build traces n_d only. The focus uses the paraxial principal planes, so at close range the plane of sharpest focus can sit slightly off the asked distance (spherical aberration).
- In the orthographic view the vignetting has no effect (every chief ray is on the axis), and a slanted view down at a floor shows no horizon.
- The mask selects materials, not objects: two objects with the same material are masked together (give them copies of the material). A mix material is selected by the name of the mix.
- The depth has one ray per pixel (no anti-aliasing): an edge pixel has the depth of the surface at its centre. Gaussian splats are not surfaces for either output.
- A path carries 4 wavelengths, so a narrow band gets few samples: the log gives the samples per pixel of the poorest band, and its noise goes as one over their square root. Average over a region, widen the bands or raise spp; 1 nm bands as the SR-5100 records them need thousands of spp. The spectral image is not denoised and not white-balanced (`--wb` acts on the RGB picture only). It is refused with `--lvc`, `--caustic-trace` and `--caustic-passes`, which add RGB pictures of their own. Its units are relative: comparing with a measurement in W·sr⁻¹·m⁻²·nm⁻¹ needs a common reference (e.g. a white or gray patch in both).
