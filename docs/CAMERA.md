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
| Physical exposure | The f-number, the shutter time and the ISO speed set the brightness as in a camera (below); off by default, when the picture keeps the scene's radiance whatever the f-number and shutter | `--iso`, `iso`; `--nits` |

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
- **Lenses in `assets/Lenses/`:** three patented designs with the patents' glasses (n_d and Abbe number), so all three are dispersive
  in the spectral build: `dgauss_50mm` (double-Gauss f/2, Tronnier, US 2,673,491), `wide_22mm` (f/2.8, 76°, for a 35 mm rangefinder
  camera, A. Nakamura, US 4,211,472, embodiment 4) and `fisheye_10mm` (f/4, 156°, R. Muller, US 4,647,161, embodiment 1); and
  `singlet_bk7_100mm`, an uncorrected N-BK7 lens. They came from W. J. Smith, *Modern Lens Design*, through pbrt's lens files, which
  round n_d to three decimals and leave out the Abbe numbers. The double-Gauss and the fisheye keep those radii (the patents', scaled);
  the wide-angle takes the patent's own, as Smith's is a variant with three other glasses. Their F and C foci differ by 0.19 %
  (double-Gauss, 93 µm), 0.28 % (wide-angle) and 0.48 % (fisheye) of f: the residual axial colour, seen as purple and green fringes on
  the bokeh of the fast double-Gauss wide open. The lateral colour stays under 0.1 % of the image height up to the edge of the field
  (11 µm for the fisheye). The set-up prints the focal length, the f-number and the field: 50.0 mm f/2.00, 22.0 mm f/2.79 and
  10.0 mm f/4.00 for the three.
- **Noise of a dispersive lens:** each wavelength takes its own path through the glass, so a camera path carries one wavelength instead
  of four. The colour noise is about three times that of the same lens without dispersion, and the denoiser does not keep the colour of
  such noise (on bright bokeh it shifted red by up to 20 %): give it more samples (four times or more) or compare without denoising.

### Physical exposure (`--iso S --nits K`, camera `"iso"`; `Render/Exposure.hpp`)

Without it, the picture holds the scene's radiance: f/16 is as bright as f/2. With an ISO speed the camera records what a real one does:

- **Film irradiance** at the centre of the picture: E = L G, with G = π / (4 N_w²) for the thin lens, N_w = N s / (s − f) the working
  f-number focused at the distance s (it grows at close range), and for a real lens (`--lens`) G from its traced rays. The corners darken
  as the lens makes them (a real lens) or with `--natural-vignette 1` (the thin lens).
- **Exposure and ISO:** H = E t (lux seconds). The ISO speed is the standard output sensitivity of ISO 12232: H = 10 / S is recorded as mid
  grey, 0.18 linear (118 in 8-bit sRGB). The scene's radiance has no unit of its own: `--nits K` says that a radiance whose Y is 1 has the
  luminance K cd/m². The written value is 0.018 S t K G times the radiance.
- **Needs** an f-number (no pinhole), a shutter time (`--shutter`, which also blurs what moves) and `--nits`; not with the orthographic view.
  The log gives the working f-number, EV100, the gain in stops and the luminance that comes out as mid grey.
- **Output:** the gain is added to `--exposure` (which stays as exposure compensation) on the `.png` / `.bmp` / `.jpg`; the `.pfm` keeps the
  radiance. With `--tonemap clip` over-exposure clips at white as on a sensor, and the clip report says by how many stops.

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
and display RGB with colour-science, so colour science on a render can be done in Python after it. The film uses the same tabulated
CIE 1931 functions (1 nm), so Python's defaults give the RGB picture back (`tools/test/spectral_image_check.py`).

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

### The sensor (`--sensor`; `Render/Sensor`)

Without it the film is ideal: each pixel holds its exact mean. With `--sensor` the written pictures are what a camera's sensor records of
that mean, in these steps (on the camera's raw values with `--camera-response`, else on the film's linear sRGB):

0. **Colour filter mosaic:** a single sensor sees colour through a filter on each pixel, the Bayer pattern R G / G B repeated over 2 × 2
   pixels (green twice, as it carries most of the luminance). Each pixel records one channel (`cfa=none`: all three, as three sensors
   would).
1. **Electrons:** the light becomes on average raw × (the white's electrons), and the count is drawn from a Poisson distribution, as
   photons arrive at random. Its standard deviation is the square root of the mean, so N electrons have a signal-to-noise ratio of √N
   (shot noise). This is the noise of the light itself; no sensor does better.
2. **Read noise:** reading the charge out adds Gaussian noise, `read` electrons before the amplifier and `adc` electrons (counted at the
   base ISO) after it. The ISO setting is the amplifier's gain, so at S times the base ISO the second part counts base / S of itself:
   σ = √(read² + (adc · base / S)²).
3. **Analog-to-digital converter:** the amplified signal is rounded to an integer (DN) of `bits` bits above a black level of 1/32 of the
   range, which keeps the noise below black (no light reads 0 on average). The range ends at the full well at the base ISO and at the
   white above it, as a camera sets its white level: a full pixel reads the top value, flat.
4. **Back** to the film's units, (DN − black) / K / (the white's electrons), then the camera's white balance.
5. **Demosaicing:** each pixel's two missing channels are interpolated from its neighbours, on white-balanced values (as raw developers
   do: the corrections below assume that R, G and B of a grey are equal), then the camera's colour matrix. `demosaic=bilinear` averages
   the nearest pixels of the missing colour. `demosaic=malvar` (the default; Malvar, He and Cutler, ICASSP 2004) adds a correction from
   the pixel's own colour that follows edges where R, G and B change together, as they mostly do in natural pictures: at a grey edge it
   halves the false colour (RMS 0.0068 against 0.0136) and keeps the edge sharper; where the colours change independently it errs more
   (30 against 34 dB PSNR on a blurred random colour field). Beyond the picture's edge the mosaic is mirrored about the edge pixel, which
   keeps the pattern.

- **Raw file:** `<out>.dng` (and `_denoised.dng`): the mosaic as recorded, 16-bit, in Adobe's Digital Negative format, with the black and
  white levels, the colour matrix (CIE XYZ to the camera's raw under D65) and the white balance, so that a raw developer can develop it
  like a camera's. `tools/test/sensor_dng_check.py` checks a render's DNG with LibRaw (rawpy) and colour-demosaicing: LibRaw reads the
  pattern, the levels, the balance and a matrix equal to the renderer's (0.1 %); the picture rebuilt from the DNG's numbers alone is the
  renderer's (4 × 10⁻⁵ median); LibRaw's own development (AHD) matches the renderer's colours within 1 % (Canon 5D Mark II response, and
  the ideal camera). Highlights that clip are developed each developer's own way.
- **Clipped highlights:** the white balance multiplies R and B (by 2.2 and 1.5 for the Canon 5D Mark II under D65) after the channels
  clipped at the top of the range, so a highlight whose G clipped turns pink, as on a camera, unless a developer recovers it.

- **ISO:** the same light gives the same electrons at any ISO. At S times the base ISO the white is reached with 1/S of the full well:
  the picture is brighter, the shot noise relative to white √S times larger, the read noise after the amplifier S times smaller in
  electrons. So raising the ISO beats brightening the picture afterwards until the read noise before the amplifier dominates. Below the
  base ISO the well is full before the white.
- **Exposure:** with `--iso` the sensor gets the physical exposure's light and `--exposure` pushes the recording afterwards, noise and
  all; without `--iso` the picture's white is the full well at the base ISO and `--exposure` sets the light.
- **Default:** `well=50000,base=100,read=1.5,adc=3,bits=14,cfa=rggb,demosaic=malvar`, round numbers of a 24-megapixel full-frame sensor (6 µm pixels): mid
  grey at ISO 100 holds 9000 electrons, a signal-to-noise ratio of 95. For a given camera use its measured full well and read noise.
  The scale is plausible: D65 brings 1.39 × 10¹⁶ photons per lux second and square metre, so a 6 µm pixel at the white of ISO 100
  (0.56 lux s) receives 2.8 × 10⁵ photons, of which a fifth become electrons (its green filter passes about a third of them).
- **Order:** after diffraction and after the denoiser. The main picture also carries the render's own (Monte Carlo) noise; the
  `_denoised` copy has the sensor's alone. `_camraw.pfm` holds the recorded raw. The spectral image, the AOVs, the mask, the depth and
  the polarizer views stay noiseless. The log gives the electrons of white and mid grey, the latter's signal-to-noise ratio and the share
  of pixels at the top of the range.

## 3. Verified

`--scene camera` (RGB and spectral builds): `project()` undoes the ray through a pixel centre; with an aperture, the rays of a film position
meet on the plane of focus (perspective and orthographic); the orthographic rays are parallel to the axis and start on the H × W rectangle;
a sphere covers the same pixel area, π R² over the pixel area, at 3 m and at 9 m, and its depth follows the sphere.

A real lens, also in `--scene camera`: the singlet's focal length matches the thick-lens formula and its principal planes lie symmetric;
focused at 1 m, rays near the axis from the axis point meet on the film (0.06 µm) and farther out they cross it nearer the lens by an amount
that grows as the square of the height in the stop (third-order spherical aberration: 2.76 and 11.0 µm at 0.4 and 0.8 mm); the f-number wide
open is the focal length over the stop's diameter and the stop closes to f/8; in the spectral build the F and C lines' focal lengths differ by
f / V_d (0.988 of it), and through the whole double-Gauss the real rays near the axis cross it where a paraxial y-u trace of the
prescription puts the F, d and C foci (within 10⁻⁵ µm; f = 50.0 mm, the patent's); in a camera the double-Gauss's weight averages 1 at the
film centre and 0.21 in the corner, the picture is upright and a wall 10 m away has the depth 10 in the middle and the corners.

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

Physical exposure, also in `--scene camera`: the thin lens's irradiance factor is π / (4 N²) at infinity and a quarter of it focused at twice the focal length; the double-Gauss's, from its traced rays, is π / (4 N²) within 0.04 % (it meets the sine condition) and drops by (8 / 2)² from f/2 to f/8 within 0.03 %; mid grey comes out at the luminance a reflected-light meter (K = 12.5) exposes for, within 2 % (12.7). A uniform sky (Y 1, `--nits 100`) at f/4, 1/100 s, ISO 100 is written as 0.0887 linear against 0.0884 (8-bit steps), and a quarter of it at f/8.

The sensor, also in `--scene film`: at the base ISO the electrons are the light's and their variance is their mean plus the read noise
squared plus the ADC's rounding, within 3 %, and its slope against the mean in DN is the conversion gain (the photon transfer curve); the
same light at ISO 100 and 1600 gives the same 20 electrons and the signal-to-noise ratios of the formula (3.53 and 4.24); above white all
pixels read the top value, and at ISO 50 the well is full at half the white; no light reads 0 on average; a film keeps its mean, gets the
shot noise's spread and draws the same noise for the same seed. Both demosaicings give a uniform colour back everywhere and a linear ramp
inside, and at a slanted grey edge the false colour that colour-demosaicing 0.2.6 (a separate implementation in Python) gives (RMS
0.013633 and 0.006779, and one pixel's three values to 10⁻⁶); the recorded mosaic holds R G / G B; the DNG file reads back.

The aberrations, the polygonal aperture, the vignetting and the motion blur have no self-check of their own yet.

## 4. Limits

- Light tracing (`--caustic-trace`, `--lvc`) connects light paths to a pinhole camera at time 0: its caustics are sharp, undistorted and not motion-blurred. It refuses the orthographic view (its weights are those of a pinhole).
- Diffraction is the pattern of a point in focus, laid over the rendered picture (whose defocus comes from the lens samples); it needs an aperture. Each of X, Y and Z is spread as for a D65-like spectrum within it, since the picture keeps no spectrum per pixel. It is not applied to the AOVs, the mask, the depth or the spectral image, and the denoiser works on the picture before it. The kernel takes 7–10 s to compute.
- The camera response balances to D65 only (it refuses `--wb`; a camera would balance to the scene's light), its matrix is fitted on smooth reflectances under D65.
- The sensor (`--sensor`) has no optical low-pass filter (fine detail aliases into false colour, as on cameras without one), no dark
  current and no fixed-pattern noise (the small differences of gain and offset between pixels). Its DNG carries no EXIF data (the ISO and
  the sensor are in its description).
- Physical exposure ignores the glass's transmittance (a real lens passes about 90 %: its T-stop is a little slower than its f-number); the sensor's noise and full well come with `--sensor`. It sets the exposure of the centre from the focus distance and the f-number at the d line.
- A real lens has spherical surfaces only (no aspheres), no reflections between them (no ghosts or flare) and no absorption in the glass; a radius of 0 is always the stop (a flat glass surface cannot be written). It refuses light tracing and the orthographic view. The RGB build traces n_d only. The focus uses the paraxial principal planes, so at close range the plane of sharpest focus can sit slightly off the asked distance (spherical aberration).
- In the orthographic view the vignetting has no effect (every chief ray is on the axis), and a slanted view down at a floor shows no horizon.
- The mask selects materials, not objects: two objects with the same material are masked together (give them copies of the material). A mix material is selected by the name of the mix.
- The depth has one ray per pixel (no anti-aliasing): an edge pixel has the depth of the surface at its centre. Gaussian splats are not surfaces for either output.
- A path carries 4 wavelengths, so a narrow band gets few samples: the log gives the samples per pixel of the poorest band, and its noise goes as one over their square root. Average over a region, widen the bands or raise spp; 1 nm bands as the SR-5100 records them need thousands of spp. The spectral image is not denoised and not white-balanced (`--wb` acts on the RGB picture only). It is refused with `--lvc`, `--caustic-trace` and `--caustic-passes`, which add RGB pictures of their own. Its units are relative: comparing with a measurement in W·sr⁻¹·m⁻²·nm⁻¹ needs a common reference (e.g. a white or gray patch in both).
