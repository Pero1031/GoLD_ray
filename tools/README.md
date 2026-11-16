# tools/

Scripts around the renderer. Python scripts need `numpy`, `Pillow` and (for some figures) `matplotlib`. Run them from the repository folder.

| Folder / file | What it is for |
|---|---|
| `build_cli.ps1` | builds the command-line renderer into `out\cli` (`-Spectral`: the spectral build, `-Dev`: no link-time code generation, a rebuild after one `.cpp` takes seconds) |
| `fetch_deps.ps1` | downloads the third-party libraries that the repository does not contain (see README, section 9) |
| **`test/`** | |
| `regress.py` | small renders and self-checks of a build, compared with an earlier build (use it when changing the renderer's structure) |
| `ab_test.py` | interleaved A/B runs of two builds: mean color (z score), time, how often a self-check fails |
| `gamut_check.py` | checks of `--gamut srgb\|p3\|rec2020`: the matrices against the primaries, the PNG's cICP chunk, the 16-bit PNG, `.pfm` of two gamuts against each other |
| **`figures/`** | figures for any scene |
| `compare_figure.py` | RGB build against spectral build under several lights, with a CIEDE2000 difference map (the README hero) |
| `thinfilm_figure.py` | the RGB build against the spectral build of a thin-film scene, side by side (`docs/THIN_FILM.md`) |
| `label_render.py` | tone maps a `.pfm` render and writes labels on it |
| `skin_figures.py` | the figures of `docs/SKIN.md` |
| `hair_figures.py` | the figures of `docs/HAIR.md`: a strip of the hair lobes, captions on the color comparison |
| `lvc_figures.py` | the vertex merging of the bidirectional path tracer: the caustic seen through a glass ball with and without merging, and the bias of the radius at a focal spot (`docs/LVC_BPT.md`) |
| `cosmetic_maps.py` | thick foundation layers as colored patches over the angles of incidence and observation, in the style of Lanza et al. (2024) (`docs/COSMETICS.md`) |
| **`studies/`** | the analyses that README and the `docs/` pages rest on |
| `metal_accuracy.py` | reference solution for the color of metals: spectral (1 nm integral) against the RGB approximation |
| `rgb_to_spectrum.py` | the engine's RGB-to-spectrum conversion as a reference port, and the Jakob-Hanika alternative |
| `gold_compare_linear.py`, `gold_compare_linear_fig.py` | RGB build against spectral build on the gold bowl, measured on linear renders |
| `metamer_mismatch.py` | how far can the color under another light be from the truth when only the D65 color (or m spectral moments) is known: linear programs over the bounded reflectances (`docs/SPECTRUM_FROM_RGB.md`) |
| `sss_random_walk.py` | a random-walk reference for the radial profile of subsurface scattering (a semi-infinite medium, index matched): the error of the engine's Burley profile, with the best-fit scale (`docs/SSS.md`, 9.1) |
| `erato_scielab.py` | the statue comparison with S-CIELAB: does the difference survive the eye's spatial filtering of color (several viewing distances)? |
| `gold_compare_gamut.py`, `erato_compare_gamut.py` | the same comparison on the bowl and on the statue, measured on what an sRGB, Display P3 or Rec.2020 display shows (with the noise floor) |
| `measured_gold.py`, `measured_metal.py` | a measured metal BRDF table (KAIST pBRDF 2020, hpBRDF 2025; unpolarized) against the engine's rough conductor (`docs/MEASURED_GOLD.md`) |
| `pbrdf_albedo.py`, `pbrdf_calibration.py` | the directional albedo of a KAIST pBRDF sample, and Spectralon as a transfer standard for its gold table (`docs/MEASURED_GOLD.md`) |
| `mie_phase.py` | the Mie phase function of a cloud-droplet population: the table the cloud scenes read (`docs/CLOUD.md`) |
| `cosmetic_stimuli.py` | scenes of skin under a layer of foundation where one control changes at a time (`docs/COSMETICS.md`) |
| `skin_moisture_cone.py`, `skin_moisture_layers.py` | which mechanism of a lotion changes how much light thin skin lets through; a lotion, then an emulsion (`docs/SKIN_MOISTURE.md`) |
| `wetness_stats.py` | whether a rendered wet surface changes the image statistics the way real wet surfaces do (`docs/WETNESS.md`) |
| **`data/`** | |
| `sr5100.py` | reads and converts CSV exports of the Topcon SR-5100 (spectral radiance) into the two-column CSV the renderer reads |
| `hsm_texture.py` | turns a Topcon SR-5000/5100 `.hsm` spectral image into a spectral reflectance texture (`.sptex`), dividing by a white-plate measurement |
| `ledcube_export.py` | converts the per-LED spectra of the lab's LED cube into the `wavelength_nm, radiance` CSVs of the scene file's light mix |
| `fading_refangle.py`, `fading_dye_model.py` | the fading experiment's reflectance tables as renderer CSVs; the dye concentration table c(t) for "fading from RGB" |
| `make_sky.py` | a clear-sky environment map (`.hdr`) from single scattering in the atmosphere, plus the sun for the scene file |
| `make_ring.py` | writes the meshes of the diamond ring scene: a round brilliant-cut diamond (flat facets, planar by construction) and a ring band with four prongs |
| `make_drops.py` | water drops on a glass plate: the caps, the contact discs, the plate's top face and the cut-out mask under the drops (`docs/WATERDROPS.md`) |
| `make_bokeh_layers.py` | writes `scenes/bokeh_layers.json`: a portrait with four layers of depth (glossy balls in front, the gold Erato in focus, balls behind it, small lights far behind: 100 mm, f/1.6, 7 blades) |
| `mesh_clean.py` | prepares an OBJ or PLY for the renderer: welds vertices, recomputes smooth normals, keeps the uv on request, turns the up axis, scales to a height in metres, drops floating fragments, closes the open boundaries (a scan or a mesh from 3D Gaussian Splatting: `docs/CAPTURE.md`) |
| **`lib/`** | shared code (not run directly): `colordiff.py` (sRGB to Lab, CIEDE2000, label font), `hdr_io.py` (Radiance `.hdr` and `.pfm` readers), `gamut.py` (sRGB / Display P3 / Rec.2020 matrices and display), `scielab.py` (S-CIELAB spatial filtering) |

Scripts that make a figure from renders expect the renders from `GoLD_rayt.exe ... --out x.pfm --denoise` (see the docstring at the top of each script). The denoiser clips negative values, so for colors outside sRGB (a candle-lit gold) use the raw `.pfm`, not `_denoised`.
