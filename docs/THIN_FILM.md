# Thin-film interference: theory and implementation

The iridescence of soap bubbles, oil films on puddles, heated titanium and coated metals appears because **light reflected from the top and bottom surfaces of a thin transparent film interferes**. The color is set by the film thickness, the angle and the wavelength. The spectral build of this project computes the reflectance of the film at the wavelength carried by each path and so produces this color. The RGB build uses the color under white light (the value integrated over wavelength). **Under white light the two look almost the same; when the light is not white, the colors differ** (section 6).

- Code: `include/Materials/ThinFilm.hpp`, `src/Materials/ThinFilm.cpp` (material `thinFilm`)
- Self-check: `GoLD_rayt.exe --scene thinfilm` (comparison with analytic formulas, sampling weights, white furnace; passes in both RGB and spectral builds)
- Scenes: `scenes/soap_bubble.json` (black background), `scenes/soap_bubble_hdr.json` (the same bubble in front of an HDR environment; `--env` and `environment` can add an HDR background to the first as well), `coated_metal.json`

## 1. Intuition: two reflections add up

Light hits a film of thickness d and index n₁ in air (index n₀). Wave ①, reflected at the top surface of the film, and wave ②, which enters the film, is reflected at the bottom surface and returns, overlap.

- ② is delayed relative to ① by its round trip inside the film. The optical path difference is **2 n₁ d cos θ₁** (θ₁ is the angle inside the film), and the phase delay is 2β.

      β = 2π n₁ d cos θ₁ / λ

- In a reflection that goes from the lower-index side to the higher-index side, the phase turns by π. Going from air into the film (soap bubble), ① turns by π, and from the film into air, ② does not. Therefore **for a soap bubble, 2 n₁ d cos θ₁ = (m + ½) λ gives constructive interference** and m λ gives destructive interference.
- Example (soap bubble, n₁ = 1.33, normal incidence): for d = 300 nm, the constructive wavelengths are 798/(m + ½) nm, and m = 1 gives 532 nm (green), which is strongly reflected. At d = 100 nm or less, the two waves are almost in antiphase and cancel, and the film looks black.
- The period of the oscillation against wavelength is Δλ ≈ λ² / (2 n₁ d cos θ₁): about 110 nm for d = 1000 nm and about 40 nm for d = 3000 nm. **The thicker the film, the more times it oscillates inside one color channel and is averaged out, so the color becomes paler.** This is why RGB cannot reproduce it.

## 2. One interface: the Fresnel amplitude reflection coefficient

For the boundary between two media of index ñ = n + iκ (complex; for metals κ is the extinction coefficient), crossed at angle θᵢ, the amplitude reflection coefficients for s and p polarization are as follows (ñᵢ and ñₜ are the incident and transmitting sides).

      r_s = (ñᵢ cos θᵢ − ñₜ cos θₜ) / (ñᵢ cos θᵢ + ñₜ cos θₜ)
      r_p = (ñₜ cos θᵢ − ñᵢ cos θₜ) / (ñₜ cos θᵢ + ñᵢ cos θₜ)

The refraction angle follows from Snell's law `ñ sin θ = constant` (the same in every layer), and `cos θ = √(1 − sin²θ)` is complex. With this, **a metal substrate enters without special treatment** (the values equal those of the existing `fresnelConductor`; the agreement at d = 0 is checked in the self-check).

## 3. Reflection of the whole film: sum of multiple reflections (Airy's formula)

There is a reflection r₀₁ at the top surface and a reflection r₁₂ at the bottom surface of the film. The light inside the film travels back and forth between the two surfaces many times, and a little of it leaves each time. Adding the outgoing waves gives a geometric series.

      r = r₀₁ + t₀₁ t₁₀ r₁₂ e^{2iβ} · Σₖ (r₁₀ r₁₂ e^{2iβ})ᵏ

Using the Stokes relations `r₁₀ = −r₀₁` and `t₀₁ t₁₀ = 1 − r₀₁²` and simplifying,

      r = (r₀₁ + r₁₂ e^{2iβ}) / (1 + r₀₁ r₁₂ e^{2iβ})

The intensity reflectance is `|r|²`, and for unpolarized light it is the mean of s and p.

      R = (|r_s|² + |r_p|²) / 2

If the film does not absorb and the substrate is a dielectric, **the transmittance is T = 1 − R**. For a metal substrate, T = 0.

**Branch of the complex square root**: inside the film and the substrate, choose the root for which the wave does **not grow** as it propagates (Im(ñ cos θ) ≥ 0). Beyond total internal reflection (the evanescent wave) cos θ is purely imaginary, and the result depends on how the sign is chosen. The implementation does not leave this to the principal branch of `std::sqrt`; it chooses by this condition (because the principal branch depends on the sign of zero).

## 4. What the self-check verifies

| Check | Reason |
|---|---|
| At d = 0, agrees with the ordinary Fresnel reflectance (dielectric and metal) | connection of the formulas |
| R = 0 for a quarter-wave antireflection film (n₁ = √n₂, d = λ/(4n₁)) | sign of the phase |
| Free film in air, normal incidence: R = F sin²β / (1 + F sin²β), F = 4R₀/(1 − R₀)² | independent analytic formula |
| For a thick film (40 µm), the mean over 400 to 700 nm agrees with the incoherent sum 2R₀/(1 + R₀) | interference washes out |
| R stays within 0 to 1. Total internal reflection from the denser side gives R = 1 | branch for the evanescent wave |
| The mean of the `sample()` weights agrees with R and T for each wavelength (Monte Carlo) | consistency of the implementation |
| White furnace (a sphere of film in a uniform white environment): the visible radiance is 1 | energy conservation |

## 5. Implementation

**Form of the material**: a smooth surface (delta-function reflection and transmission). `ThinFilm` (derived from `Material`) chooses between reflection and transmission.

- **Three kinds of substrate**: a free film (air on both sides: a soap bubble; transmission passes straight through), a dielectric substrate (oil on water; transmission bends with the index of the substrate) and a conductor substrate (coated metal; reflection only).
- **Film thickness field** (`FilmThickness`): a constant thickness, or a "dripping" model in which the liquid flows downward under gravity.
  - **Gradient**: increases along an axis (up-down) from `from` to `to`. The shape is `t^exponent` (with `exponent` > 1 the thin part at the top lasts longer). The top of a soap bubble becomes thin and turns into a black film (a few tens of nm or less; the reflections cancel and it looks black), and it thickens toward the bottom.
  - **Warp** (`warp`): bends the position with a smooth displacement field (± a few cm). The horizontal stripes undulate and swirl (the liquid circulating).
  - **Noise** (`noise`): adds a variation of ± several tens of nm with value noise (the same as the density field of media). With `stretch` set to (1, 2.5, 1), the stripes become fine vertically and long horizontally.
  - **Coverage** (`coverage`, optional): decides where the film exists with a separate noise (where there is no film, light passes straight through; irregular patches such as an oil film).
- **Reflectance per wavelength**: the spectral build computes `R(λᵢ)` at the 4 wavelengths λᵢ carried by the path. The direction does not depend on wavelength (the index is constant), so, unlike dispersive glass, the secondary wavelengths need not be dropped.
- **Choosing reflection or transmission**: one random number decides for all wavelengths. The probability of reflection is the mean pR of R over the 4 wavelengths, and the weights are Rᵢ/pR for reflection and Tᵢ/(1 − pR) for transmission. The expected values are Rᵢ and Tᵢ for each wavelength.
- **Radiance correction**: since the index changes on transmission, 1/η² is multiplied in radiance mode (as in `Dielectric`). For a free film it is 1.
- **RGB build**: the reflectance is integrated on a 5 nm grid with the film sensor (weights that make white (1,1,1) under D65). That is, "the color of the film under white light". The color when the light spectrum is changed, or the color after stacking reflections, cannot be reproduced. This is the meaning of the comparison with the spectral version.

## 6. Difference from RGB: when the light is not white

The thin film of the RGB build is the color under white light R_c(θ, d), the reflectance R(λ, θ, d) integrated with the D65 sensor weights (a table of angle x thickness; close to the best that RGB can do).

1. **Under white light, for a single reflection, the RGB and spectral versions are almost the same.** Even if the film is thick and R(λ) oscillates many times within one channel, the integrated color is the same (D65 area light in `scenes/soap_bubble.json`: mean ΔE00 of the render 2.6, median 1.5, including noise and the denoiser).
2. **When the light is not white, they differ.** The spectral version computes ∫ R(λ) L(λ) s_c(λ) dλ, but the RGB version computes R_c × L_c, and the two do not agree, because R oscillates with wavelength and the light spectrum is tilted. The same soap bubble rendered with a tungsten (2856 K) area light has a mean ΔE00 of **5.2** (95th percentile 15), and with a candle (1850 K) **7.2** (95th percentile 27). The RGB version is stronger in blue-green, and the spectral version is stronger in yellow and pink.
3. **Stacking reflections widens the difference.** For two reflections at the film, R_c² differs from ∫ R² (the same structure as F(λ)ⁿ for gold).
4. **The change with angle can be represented in the RGB version too** (the table is interpolated over angle). The shift of the stripes toward blue with angle appears in both.

Figure: `tools/figures/thinfilm_figure.py` (RGB version on the left, spectral version on the right).

## 7. Usage

```json
{ "type": "thinFilm", "ior": 1.33,
  "thickness": { "from": 15, "to": 900, "axis": [0, 1, 0], "range": [0.5, -0.5], "exponent": 1.7,
                 "warp": { "amplitude": 0.07, "frequency": 3.2, "seed": 5 },
                 "noise": { "amplitude": 70, "frequency": 5, "octaves": 4, "seed": 11, "stretch": [1, 2.5, 1] } },
  "substrate": "free" }
```

- `thickness`: a number (nm), or `{ "nm": 450 }`, or a gradient (`from` and `to` in nm, `axis` the axis, `range` the position along the axis [start, end] (m), `exponent`) + `warp` (`amplitude` (m), `frequency`, `seed`) + `noise` (`amplitude` (± nm), `frequency` (number of patterns per 1 m), `octaves`, `seed`, `stretch`). `coverage` (`frequency`, `octaves`, `lo`, `hi`, `seed`) is a mask for where the film exists.
- `substrate`: `"free"` (default), `{ "dielectric": 1.33 }`, `{ "metal": "Au" }` (`Ag`, `Cu`), `{ "table": "file.csv" }`.
- Guide values: soap solution n ≈ 1.33 (thickness a few tens to 1000 nm), oil ≈ 1.45 (on water, substrate 1.33), lacquer ≈ 1.5, titanium oxide ≈ 2.3 (on gold or copper the color change is large).

## 8. Oil film

- **An oil film (oil on water, `substrate: dielectric`) was rendered, but it did not look good (even with a black background the colors are pale and muddy), so it was removed from the scenes.** The feature remains, and `coverage` can make irregular patches. The color of the film turns white at grazing viewing angles (70° or more oblique) and is most vivid at 30 to 50° from the vertical.

## References

See [REFERENCES.md](REFERENCES.md), group "Dispersion and thin films".
