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
| **`studies/`** | the analyses that README and `docs/METAL_ACCURACY.md` rest on |
| `metal_accuracy.py` | reference solution for the color of metals: spectral (1 nm integral) against the RGB approximation |
| `rgb_to_spectrum.py` | the engine's RGB-to-spectrum conversion as a reference port, and the Jakob-Hanika alternative |
| `gold_compare_linear.py`, `gold_compare_linear_fig.py` | RGB build against spectral build on the gold bowl, measured on linear renders |
| `metamer_mismatch.py` | how far can the color under another light be from the truth when only the D65 color (or m spectral moments) is known: linear programs over the bounded reflectances (`docs/SPECTRUM_FROM_RGB.md`) |
| `erato_scielab.py` | the statue comparison with S-CIELAB: does the difference survive the eye's spatial filtering of color (several viewing distances)? |
| `gold_compare_gamut.py`, `erato_compare_gamut.py` | the same comparison on the bowl and on the statue, measured on what an sRGB, Display P3 or Rec.2020 display shows (with the noise floor) |
| **`data/`** | |
| `sr5100.py` | reads and converts CSV exports of the Topcon SR-5100 (spectral radiance) into the two-column CSV the renderer reads |
| `make_ring.py` | writes the meshes of the diamond ring scene: a round brilliant-cut diamond (flat facets, planar by construction) and a ring band with four prongs |
| `mesh_clean.py` | cleans an OBJ for the renderer: welds vertices, recomputes smooth normals, optionally closes the open boundaries |
| **`lib/`** | shared code (not run directly): `colordiff.py` (sRGB to Lab, CIEDE2000, label font), `hdr_io.py` (Radiance `.hdr` and `.pfm` readers), `gamut.py` (sRGB / Display P3 / Rec.2020 matrices and display), `scielab.py` (S-CIELAB spatial filtering) |

Scripts that make a figure from renders expect the renders from `GoLD_rayt.exe ... --out x.pfm --denoise` (see the docstring at the top of each script). The denoiser clips negative values, so for colors outside sRGB (a candle-lit gold) use the raw `.pfm`, not `_denoised`.
