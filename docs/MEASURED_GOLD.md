# Measured gold and analytic gold: the KAIST pBRDF gold sample against `RoughConductor` (2026-11-10)

**Status:** the comparisons of shape, spectral shape and albedo, and the scale check with Spectralon (sample 1, the white standard), have been done.
Conclusion: near normal incidence the absolute values of the dataset are 0.64 times the theory, and fall further with incidence angle (a systematic error of the measurement system). With Spectralon as a transfer standard for the correction, the hemispherical reflectance of the gold is **0.63 to 0.79 times** that of Johnson & Christy's gold (θi ≤ 45°).
The aluminum of the successor hpBRDF (2025, 68 wavelengths) was compared by the same procedure (section 6): the shape is the same "spike + haze", and the absolute value is again about 0.4 times the theory, falling with angle. **For both datasets the absolute values cannot be used as they are. The conclusion about shape (two GGX lobes are needed) agrees.**
Section 7: if the table is taken to be "BRDF × cos θi × an instrument constant of 0.60", the angle dependence of Spectralon, gold and aluminum is explained at once, and gold comes out at about 0.75 times the theory.
Reproduction: `python tools/studies/measured_gold.py` (about 30 seconds), `python tools/studies/pbrdf_calibration.py` (about 1 minute). Output goes to `out/studies/measured_gold/` (`results.txt`, `measured_gold.png`, `calibration.txt`). Copies of the raw logs: `docs/data/measured_gold_results.txt`, `measured_gold_calibration.txt`, `measured_aluminum_results.txt`.

This study works on unpolarized extractions of the tables. The renderer itself reads the 2025 hpBRDF, with its measured Mueller matrices (used by `--stokes`, see [POLARIZATION.md](POLARIZATION.md)), through a PCA-compressed table (`.hpbpca`) as the scene-file material `"type": "measured"` (`include/Materials/MeasuredHpBRDF.hpp`, example `scenes/hpbrdf_aluminum.json`). The tables are not in the repository.

## 1. Data

Sample 6, "Gold", of the KAIST polarimetric BRDF dataset (Baek, Zeltner, Ku, Hwang, Tong, Jakob, Kim, "Image-Based Acquisition and Modeling of Polarimetric
Reflectance", ACM TOG 39(4), SIGGRAPH 2020). What was confirmed from the materials table and README of the distribution page, and from supplementary material 1:

- The table `pbrdf_T[φd 361, θd 91, θh 91, wavelength 5, 4, 4]` (Mueller matrices, float32). θh is sampled quadratically (θh = (j/90)²·90°); θd and φd are in 1° steps. The wavelengths are 450 / 500 / 550 / 600 / 650 nm (narrow-band filters).
- Unmeasured bins are NaN (`pbrdf_T`). `pbrdf_T_f` is the filled-in version. θd < 5° is entirely unmeasured (the minimum angle between light and camera is 6.75°).
- **The unit of the values and the method of absolute calibration are written neither on the distribution page nor in the README** (the calibration file holds the irradiance of the light).
- The unpolarized component (element (0, 0) of the Mueller matrix) was used. The extracted table is in `references/kaist_pbrdf/` (not tracked by git).

## 2. How it was compared

The engine's `RoughConductor` (GGX, α = roughness², height-correlated Smith, the conductor Fresnel of Johnson & Christy's Au, Kulla-Conty multiple-scattering compensation) was
reimplemented in numpy and evaluated at the same (θh, θd, φd). The conversion of directions follows the convention of MERL's `BRDFRead.cpp`.

1. **Shape:** find the α that minimizes the variance of log(measured / model) over bins with θh < 35° and values at least 20 times the floor (stray light); the scale is obtained separately. GGX and Beckmann. A two-lobe GGX was also fitted to the θd = 5° slice.
2. **Spectral shape:** in the spike (θh < 0.5°) and the haze (θh 10 to 35°), divide the per-wavelength ratio measured / Fresnel (Johnson) by its value at 650 nm.
3. **Directional albedo:** integrate the measured table over the outgoing hemisphere (0.005° steps within 3° of the specular direction, 0.25° steps outside). Compare with the model's single scattering and with single + multiple scattering.

## 3. Results

**Shape.** The measurement has two components: a sharp spike about 1° wide (half maximum) and a broad haze beyond θh of 10°, and at θh > 50° it reaches a floor of 0.001 to 0.003 sr⁻¹ (0.4 to 1% in Lambert terms).

| Model | Fitted α (roughness) | Residual of log(measured / model) |
|---|---|---|
| GGX, 1 lobe | 0.049 (0.22) | ×1.9 |
| Beckmann, 1 lobe | 0.146 (0.38) | ×3.2 |
| GGX, 2 lobes (θd = 5° slice) | spike 0.022 (0.15) weight 0.39, haze 0.24 (0.49) weight 0.08 | ×1.23 |

GGX has long tails and so fits better than Beckmann, but one lobe cannot represent the spike and the haze together (the peak is off by half, and by a factor of 4 at 5°). Two lobes fit to ×1.2. The weights are multipliers on the Fresnel term, and their sum of 0.47 agrees with the ratio in the albedo below.

**Spectral shape (relative to 650 nm).**

| | 450 | 500 | 550 | 600 | 650 nm |
|---|---|---|---|---|---|
| Fresnel F0 (Johnson & Christy) | 0.43 | 0.50 | 0.83 | 0.95 | 1 |
| Measured spike / Fresnel | 0.66 | 0.73 | 0.86 | 0.90 | 1 |
| Measured haze / Fresnel | 1.13 | 1.18 | 1.13 | 1.05 | 1 |

The spike is **more saturated** than Johnson's gold (34% lower at 450 nm), and the haze is **close to neutral** (13 to 18% higher in blue). A natural reading is that dirt or a film on the gold surface dilutes the color of the haze, while the spike carries the color of the bulk metal.

**Directional albedo (measured / model single + multiple scattering).**

| θi | 450 | 500 | 550 | 600 | 650 nm |
|---|---|---|---|---|---|
| 15° | 0.38 | 0.41 | 0.44 | 0.45 | 0.46 |
| 30° | 0.34 | 0.37 | 0.39 | 0.40 | 0.42 |
| 45° | 0.36 | 0.38 | 0.38 | 0.39 | 0.40 |
| 60° | 0.41 | 0.40 | 0.38 | 0.37 | 0.38 |
| 75° | 0.39 | 0.38 | 0.34 | 0.31 | 0.32 |

The measured gold returns only **0.3 to 0.46 times the light that Johnson & Christy's gold should return**, independently of wavelength and almost independently of incidence angle. The unmeasured bins (NaN) hold no energy (the same values with the filled-in table).

## 4. Checking the scale: Spectralon (sample 1) as a transfer standard (`tools/studies/pbrdf_calibration.py`)

Spectralon is a diffuse standard with reflectance of about 0.99, so its BRDF should be 0.99 / π = 0.315 sr⁻¹ in all directions. The table of the same dataset gives:

| Configuration | Table value (550 nm) | Ratio to theory |
|---|---|---|
| θi, θo < 20° | 0.203 sr⁻¹ (±2% at all five wavelengths) | **0.64** |
| θi 30 to 40° (θo 20 to 40°) | 0.159 | 0.50 |
| θi 60 to 70° | 0.105 | 0.33 |
| θi 80 to 90° | 0.087 | 0.28 |

That is: (1) the absolute value is 0.64 times near normal incidence, (2) it falls further with incidence angle (more gently than cos θi), (3) it is flat across wavelength (the spectral calibration is good), (4) there is an excess of 1.2 times at θh ≈ 0 (θi ≈ θo) (several times at 60° and above).
The gold table was measured with the same instrument, so taking **"gold ÷ Spectralon × 0.315"** in the same geometric configuration cancels (1) and (2). The correction was made in two versions: one that keeps the excess at θh < 5° (A), and one that replaces it with the values at 5 to 8° (B).

**Hemispherical reflectance of gold after correction / Fresnel of Johnson & Christy** (450 / 500 / 550 / 600 / 650 nm):

| θi | Correction A | Correction B |
|---|---|---|
| 15° | 0.63 / 0.67 / 0.72 / 0.73 / 0.75 | 0.64 / 0.68 / 0.73 / 0.75 / 0.77 |
| 30° | 0.62 / 0.66 / 0.69 / 0.72 / 0.76 | 0.64 / 0.69 / 0.72 / 0.74 / 0.79 |
| 45° | 0.64 / 0.69 / 0.67 / 0.69 / 0.72 | 0.70 / 0.76 / 0.74 / 0.76 / 0.79 |
| 60° | 0.52 / 0.54 / 0.49 / 0.48 / 0.50 | 0.65 / 0.67 / 0.62 / 0.62 / 0.64 |
| 75° | 0.12 (the correction breaks down) | 0.17 (same) |

At 60° and above Spectralon departs from Lambertian (toward specular at grazing angles, and darker elsewhere), so it cannot be used as a transfer standard. For θi ≤ 45° **the corrected gold is 0.63 to 0.79 times the theory**, with blue lower than red (more saturated). The spectral ratios of the spike and the haze hardly change with the correction (because Spectralon is flat across wavelength).

## 5. How to read this

- The conclusion about shape does not depend on the scale: **to render this sample, two GGX lobes (a sharp specular lobe + haze) are needed.** The engine's `RoughConductor` has one lobe.
- **The absolute values of the dataset cannot be used as they are** (0.64 times + angle dependence). When using them, take the ratio to Spectralon. This agrees with the absence of a definition of the values and of the calibration method on the distribution page and in the README.
- Reciprocity (φd against φd + π) is off by a median of 25% inside the specular lobe of the gold, and 16% for Spectralon. The table has measurement noise and unevenness from interpolation.

## 6. Aluminum of the hpBRDF (Moon et al., SIGGRAPH Asia 2025) (`tools/studies/measured_metal.py`)

"Aluminum" of the successor dataset from the same laboratory (68 wavelengths, 414 to 950 nm, 8 nm steps, the same table format as the 2020 version, 13 GB per material, published on HuggingFace; there is no pure gold) was compared by the same procedure (`references/kaist_pbrdf/aluminum_unpolarized.npz`, extracted from the original file `Aluminum.hpbrdf`).
The optical constants are McPeak 2015 and Rakić 1995 (RefractiveIndex.info, CC0; added as `assets/Metal_data/McPeak_Al.csv` and `Rakic_Al.csv`). The two give almost the same results.
Aluminum has almost flat reflectance in the visible (F0 ≈ 0.92) and a shallow dip near 800 nm (0.86), so it is a check point for both scale and spectral shape.
Output goes to `out/studies/measured_aluminum/` (McPeak) and `measured_aluminum_rakic/`.

| Item | Result |
|---|---|
| Shape | The same "spike + haze" as gold. GGX 1 lobe ×2.3, Beckmann ×2.7, **GGX 2 lobes ×1.14** (spike α 0.048, weight 0.31; haze α 0.48, weight 0.15) |
| Reciprocity | Median 10% inside the specular lobe (> 0.5 sr⁻¹) (25% for the 2020 gold; the table is smooth because it is filled in) |
| Hemispherical reflectance / Fresnel | θi 15°: 0.38, 30°: 0.45, 45°: 0.31, 60°: 0.18, 75°: 0.12 (almost the same at all wavelengths) |
| Spectral shape (ratio to 650 nm) | Spike: fluctuates by ±20% from band to band, 1.0 to 1.25 over 446 to 702 nm, 0.96 at 798 nm, 0.86 at 902 nm. Haze: decreases monotonically from 1.27 at 446 nm to 0.84 at 902 nm, toward blue. The dip at 800 nm is buried in the fluctuation and cannot be identified |

**How to read this.** The 2025 version also has **absolute values of about 0.4 times the theory, falling with incidence angle**. It is natural to see this as a systematic error of the same type as the Spectralon of the 2020 version (scale + angle dependence), and
that two datasets show "the metal returns only 40% of the theory" is thought to be a problem of the measurement processing (the conversion from sphere images to BRDF), not of the sample. Across wavelength there are ±20% fluctuations from band to band and a monotonic slope of the haze.
The conclusion about shape (two lobes are needed) agrees in the two datasets.

## 7. What the values of the table are: taking them as "BRDF × cos θi × a constant" (2026-11-10)

The body of the paper (2025, 11 pages, read in full with `pypdf`) has no definition of the values; it says only that M is determined by the least squares of equation (5) (the difference between the measured intensity f and the prediction obtained by passing the light's Stokes vector [L^λ, 0, 0, 0] through the polarization elements and M).
The cosine of the incidence angle, the irradiance and the distance to the light do not appear in the equations. The 2020 calibration file has `M_air_inv` (the inverse matrix of the measurement looking directly at the light).
So it was assumed that "table = Stokes vector of the reflected light ÷ Stokes vector looking directly at the light = BRDF × cos θi × an instrument constant", and it was checked whether hemispherical integral ÷ expected value ÷ cos θi is constant independently of θi (550 nm):

| θi | 5° | 15° | 30° | 45° | 60° | 75° |
|---|---|---|---|---|---|---|
| Spectralon 2020 (÷ 0.99) | 0.60 | 0.60 | 0.60 | 0.64 | 0.76 | 1.33 |
| Gold 2020 (÷ Johnson's Fresnel) | 0.42 | 0.45 | 0.45 | 0.54 | 0.76 | 1.29 |
| Aluminum 2025 (÷ McPeak's Fresnel) | 0.41 | 0.39 | 0.55 | 0.44 | 0.36 | 0.47 |

Spectralon is exactly constant at 0.60 up to 45° (beyond that, the known behavior of Spectralon itself of going specular at grazing angles, and the excess common to the instrument). Gold is also constant up to 30°, and the rise at 60° and above is the same as for Spectralon. Aluminum fluctuates but has no trend.
**Reading: the table is not a BRDF but a "radiance factor" (with cos θi included), multiplied by an instrument constant (0.60 for the 2020 version).** If the constant does not depend on the material, gold is 0.45 / 0.60 = **0.75 times** Johnson's, and aluminum about 0.45 / 0.60 times McPeak's; both are natural values for "a real object 20% or so below the theory because of oxidation, dirt and plating on the surface".
The procedure for use: **table × (1 / cos θi) × (1 / 0.60)**.

**What was learned from the body of the 2020 paper (`polarbrdf_main.pdf`, read in full with pypdf).** The supplementary material of 2025 is not in the ACM DL either, and is not public. The body of the 2020 paper says:
- Radiometric calibration: "to correct for the polarization bias of the LED spectrum and for the differences in transmittance of the polarization elements, **a Spectralon sphere (Labsphere SRM-990, reflectance 99%) was photographed in five bands, and a constant multiplier was determined for each band, to guarantee color reproduction and energy conservation**". The authors intended absolute values, but since the measured table is fr × cos θi × 0.60, it is presumed that cos θi remained in the calculation that determined the multiplier (for example, a mean over the pixels of the sphere was used).
- **The gold sample is "a gold-plated metal sphere on a carbon steel base"**. For validation, a flat gold sample of the same kind (Thorlabs' unprotected gold mirror PF20-03-M03) was measured with a laser ellipsometer and compared with the ideal Mueller matrix calculated from the refractive index of gold at 633 nm. What is compared is **the shape of the Mueller matrix normalized to −1 to 1**, not absolute values. It states explicitly that "a discrepancy remains because of the difference in composition between the sphere and the flat sample".
- So that gold is 0.75 times the theory and is low in blue is natural for "plated gold", and the spike + haze can be read as "the specular reflection of the plated surface + the roughness of the base".
- In practice, the correction above (÷ cos θi ÷ 0.60) is sufficient.
