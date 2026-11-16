# Make-up (foundation) layers (2026-11-09)

**Status:** stage A (the layer alone: calculation and checks, section 3) and stage B (putting it on skin, `foundation` in JSON, section 4) are done. A primer (it fills furrows and pores, `primer` in JSON, section 4.1) was added on 2026-11-12. The aim is to **make stimuli with knobs that can be moved one at a time** (a picture to show is a bonus if one comes out). Only foundation is treated for now.

## 1. The paper

D. Lanza, J. Padrón-Griffe, A. Pranovich, A. Muñoz, J. R. Frisvad, A. Jarabo, "Practical Appearance Model for Foundation Cosmetics", *Computer Graphics Forum* 43(4), 2024 (EGSR). The full text and figures 4, 6 and 9 were read. **It was implemented from the description in the paper. The public code (non-commercial use only) was not used** (only one configuration file of the public repository was looked at, to check the units of the parameters).

**Model:** the foundation layer is treated as a thin medium with two kinds of suspended particles.

| Particle | Real counterpart | Scattering | Look |
|---|---|---|---|
| Spherical diffuser | titanium dioxide, silica | two Henyey-Greenstein lobes | matte |
| Platelet | talc, kaolin, mica | specular reflection; the platelet orientation follows an SGGX distribution | shiny |

- The layer is index-matched to air (no interface, no refraction). In the paper's fit, the refractive index came out at about 1.
- The layer is thin, so light leaves from the point where it entered (the lateral spread inside the layer is ignored).
- Colour comes from absorption only (at each collision, the albedo of the particle that was hit is multiplied in).
- Light is followed by a random walk inside the layer (the position-free Monte Carlo of Guo et al. 2018).

**Formulas:** the unit of length is the mean free path in the densest direction.

    σt(ω) = c_d + (1 − c_d) Ĉp(ω)                       extinction. Ĉp is the projected area of the platelets (1 along the mean platelet normal)
    Λ(ω)·f(ω → ω′) = [ Λd c_d f_d + Λp (1 − c_d) Ĉp(ω) f_p ] / σt(ω)     what a collision sends to ω′
    f_d = w_g HG(g1) + (1 − w_g) HG(g2)                 diffusers
    f_p = D_SGGX(ω_h) / (4 Ĉp(ω))                       platelets (no Fresnel; ω_h is the half vector)

## 2. Knobs (`CosmeticLayer::Params`)

| Knob | Symbol | Appearance it changes |
|---|---|---|
| `thickness` | t | optical thickness of the layer. Hiding power |
| `diffuserFraction` | c_d | fraction of spheres. The closer to 1, the more matte; the smaller, the more the platelets (shine) count |
| `plateletRoughness` | α_p | spread of the platelet orientation. The smaller, the sharper the highlight |
| `plateletTiltDeg` | θ_p | tilt of the mean platelet orientation (degrees; tilted towards +y) |
| `diffuserAlbedo`, `plateletAlbedo` | Λd, Λp | colour of the particles (RGB) |
| `g1`, `g2`, `lobeWeight` | g1, g2, w_g | scattering of the spheres (forward or backward) |

According to the paper, c_d changes the **strength** of the highlight and α_p its **shape** (two separate dimensions of gloss).

## 3. Results of stage A

**Implementation:** `Materials/CosmeticLayer.hpp` (phase functions and sampling of platelets and spheres, the random walk inside the layer `walk()`, an estimate with a prescribed outgoing direction `reflectionSample()`). The check is `--scene cosmetic` (PASS in both RGB and spectral, about 30 seconds).

**Checks:**

| Check | Result |
|---|---|
| Platelets: the distribution of normals gives the projected area, the phase function integrates to 1, samples follow the phase function (roughness 0.15 / 0.5 / 0.3 + tilt 35°) | difference 0.3% or less |
| Spheres: integral 1, mean cosine of the samples equals w_g g1 + (1 − w_g) g2 | 0.3001 against 0.3000 |
| Energy: with non-absorbing particles, reflection + transmission = 1 (thickness 0.25 / 1 / 4, three mixtures) | difference 0 (no terminated paths) |
| Closed form of single scattering (isotropic spheres, thickness 1, incidence 30°): light reflected after one collision, light that passes without collision | 0.1255 against 0.1256, 0.3151 against 0.3152 |
| Reciprocity: f(a → b) = f(b → a) (tilted platelets, coloured particles, 3 pairs) | difference 0.2% |
| Reflectance of the two estimates (the plain walk and the estimate with a prescribed outgoing direction) | agree (within 0.4%) |

**Properties of the layer alone for the 4 products fitted in the paper** (values of the paper's table 2. Green channel, incidence 30°):

| Product | Thickness | Reflected | Straight through | Scattered and transmitted | Absorbed |
|---|---|---|---|---|---|
| Dewy 1 | 0.25 | 0.06 | 0.75 | 0.18 | 0.01 |
| | 0.5 | 0.12 | 0.57 | 0.30 | 0.02 |
| | 1 | 0.22 | 0.32 | 0.42 | 0.04 |
| | 4 | 0.47 | 0.01 | 0.33 | 0.18 |
| | 16 | 0.56 | 0 | 0.03 | 0.41 |
| Matte 2 | 0.25 | 0.11 | 0.76 | 0.09 | 0.05 |
| | 0.5 | 0.18 | 0.57 | 0.15 | 0.10 |
| | 1 | 0.27 | 0.33 | 0.21 | 0.19 |
| | 4 | 0.39 | 0.01 | 0.10 | 0.49 |
| | 16 | 0.40 | 0 | 0 | 0.60 |

- **At the thicknesses the paper used for its pictures on skin (0.25 to 0.5), 60 to 80% of the light passes straight through without being scattered.** The skin below is still clearly visible in this setting.
- **At thickness 4, the straight-through fraction falls to 1%** (the skin below is hidden). The reflection levels off with thickness: 0.56 for dewy 1 and 0.40 for matte 2.
- At the same thickness, the transmission differs by product: at thickness 1, the light that is scattered and reaches the skin is 0.42 for dewy 1 and 0.21 for matte 2 (matte 2 absorbs more).

**Highlight of a thick layer (16)** (incidence 40°. Light entering a cone of 5° half angle in the plane of incidence. Per steradian):

| Product | Back −40° | Normal 0° | Mirror +40° | Brightest direction | Brightest / normal |
|---|---|---|---|---|---|
| Dewy 1 | 0.118 | 0.160 | 0.319 | +40° | 1.99 |
| Dewy 2 | 0.096 | 0.106 | 0.092 | +50° | 1.78 |
| Matte 1 | 0.124 | 0.147 | 0.262 | +40° | 1.78 |
| Matte 2 | 0.080 | 0.100 | 0.075 | 0° | 1.00 |

(Dewy 2 has an extremely small platelet roughness of 0.002. Its highlight is narrower than the cone of 5° half angle, so it comes out weak in the cone average.)

### Comparison with the paper's figure 4 (measurement) and figure 6 (fit)

The check also outputs the same quantity as the paper's figures (the reflectance π f relative to a white diffuse plate, for incidence angles 20 to 70° and viewing angles −40 to +70°). `tools/figures/cosmetic_maps.py` turns it into coloured patches. They were compared with the paper's figures (the figures are the paper's, so they are not reproduced here).

| Product | Paper's measurement and fit | This implementation (values of table 2) | Verdict |
|---|---|---|---|
| Dewy 1 | a broad white band in the specular direction | same | matches |
| Dewy 2 | a thin sharp white line near the specular direction | same (the line is about +10° from the specular direction) | roughly matches |
| Matte 1 | **no highlight** (bright only at grazing incidence of 70°) | **a white band appears in the specular direction** | **does not match** |
| Matte 2 | no highlight, almost uniform | same | matches |

- With the values of table 2, matte 1 does not match the paper's figure (its values are almost the same as dewy 1's, and a similar highlight appears). It is not used for stimuli; the knobs are chosen by hand.
- **The direction of the platelet tilt was chosen to match the measurement of matte 2.** The text of the paper does not fix the direction (sign, axis). If tilted towards the light source, matte 2 (tilt 52°) becomes brighter on the side of the light source, which the measurement does not show. With the opposite direction (tilted in the direction the light travels), it becomes uniform and matches the measurement.

## 4. Stage B: putting it on skin (2026-11-09)

**Usage:** write `foundation` in the JSON skin material (`docs/SCENE_FILE.md`).

```jsonc
"skin": { "type": "skin", ...,
          "foundation": { "preset": "dewy1", "thickness": 1.0 } }                               // from the values of the paper's 4 products
"skin": { "type": "skin", ...,
          "foundation": { "thickness": 1.0, "diffuserFraction": 0.95, "plateletRoughness": 0.1,  // the knobs directly
                          "diffuserAlbedo": [0.97, 0.93, 0.87] } }
```

`preset` is `dewy1`, `dewy2`, `matte1` or `matte2` (the paper's table 2). Keys written after it override single values. When run, a `[Foundation] ...` line prints how much the layer reflects and lets through along the normal.

**How it is put together (difference from the paper):** the paper solves the layer and the skin in one random walk. Here, **only the reflection of the layer is solved by a random walk, and the skin below is seen through a "table of what the layer lets through"** (the same way as the flakes and emulsion layers that the skin material already has).

| Component | Treatment |
|---|---|
| Reflection of the layer itself (matte scattering, platelet highlight) | one sample per evaluation, by a random walk inside the layer (`reflectionSampleT`). Noisy, but the mean is the BSDF of the layer |
| Specular lobe of the skin (shine of sebum) | seen only through the light that **passes through the layer without a collision**: multiplied by `straight(μi) · straight(μo)` |
| Diffuse lobe of the skin, white scatter of flakes | seen through **all the light that passes through the layer**: `through(μi) · through(μo) / (1 − Rback · A)`. The denominator is the sum of the light going back and forth between the skin and the layer (A is the albedo of the skin) |
| Subsurface scattering | `through(μo)` on the viewing side; on the side where light enters, the diffuse mean `throughDiffuse / (1 − Rback · A)` |

The table (`CosmeticCoat`) is built when the scene is read, from 24 bins of direction cosine × 100 000 random walks (the azimuth is averaged). The error of the table does not go away when more render trials are added (the same error applies to every trial), so it must be kept small: with 20 000 walks, the white furnace value moved by ±0.6% depending on the seed. Sampling chooses a lobe of the layer with a probability according to the size of the layer's reflection, and draws the direction from a mixture of "the specular reflection when the first collision hits a platelet" and a cosine distribution (the probability density can be written in closed form, so multiple importance sampling can be used as it is).

**Smoothing (`leveling`, 0 to 1, default 1):** the reflection of the layer is calculated not with the normal bent by the relief of the skin but with the **normal of the layer** (`SurfaceInteraction::nl`). The normal of the layer is the bent normal and the smoothed (relief-free) normal, mixed by `leveling`. At 1, the furrows and pores do not show in the shine of the layer. At 0, it is the same as the first version. **The lobes of the skin itself (shine, diffuse, subsurface scattering) keep the relief**, so in a thin layer the relief shine of the skin is seen in the light that passes through, and as the layer gets thicker it is replaced by the smoothed shine of the layer. The amount the layer lets through is also taken at the angle to the normal of the layer. The integrator weights by cos(normal of the skin, wi), so the lobe of the layer is multiplied by cos(normal of the layer, wi) / cos(normal of the skin, wi) (the same way as the specular lobe of a film). The thickness of the layer is the same everywhere, and `leveling` is a knob independent of thickness.

**Intermediate values of the smoothing do not give a natural in-between.** Tried with sharp platelets (roughness 0.05, 50% platelets), the network of furrows stays clearly visible in the shine even at `leveling` 0.7. Mixing the normals only weakens the lines of the furrows, and sharp platelets react to the slightest tilt. Use 1 (full smoothing). With full smoothing, a layer rich in platelets shines in an orderly way, like a resin plate.

**Direction of the platelet tilt:** when the shape has a tangent (sphere, plate, mesh with UV), the coordinate system of the layer is aligned to that tangent. With an arbitrary tangent, the direction of the tilt jumped across a line on the surface and the shine showed a seam (every product of the paper has a tilt of 1.5° or more).

**Checks** (7 of `--scene cosmetic`. PASS in both RGB and spectral):

| Check | Result |
|---|---|
| Table of a non-absorbing layer: reflection + transmission = 1 | difference 0 |
| Skin with a layer of thickness 10⁻⁶: evaluation and probability density equal those without a layer | difference 0 |
| Albedo of `sample()` and `eval()` agree (dewy 1 and matte 2, with and without subsurface scattering, thickness 1) | 0.565 against 0.567, 0.235 against 0.237, 0.457 against 0.455, 0.229 against 0.230 |
| Skin with relief + layer: the normal of the layer (tilt 0° at smoothing 1, equal to the skin's normal at 0), agreement of `sample()` and `eval()` (smoothing 1 / 0 / 0.5) | 0.555 against 0.555, 0.557 against 0.557, 0.561 against 0.561 |
| White furnace: a white sphere with a non-absorbing layer (thickness 1) | without layer 0.998 to 0.999, with layer 0.988 to 0.992 (3 seeds; about 0.8% darker) |
| White furnace: a skin sphere with subsurface scattering with dewy 1 (thickness 1) | (0.802, 0.526, 0.428) → (0.710, 0.494, 0.411). Does not exceed 1 |

The existing checks related to skin (`relief`, `sss`, `furnace`, `metals`, `aov`, `layered`) PASS in RGB. `relief` was fixed on 2026-11-11 (PASS in RGB and spectral × 4 seeds). 8 of `--scene cosmetic` (information) compares the white furnace of a non-absorbing layer with and without subsurface scattering (for the foundation layer the difference is 0.003 or less).

**Stimuli:** `tools/studies/cosmetic_stimuli.py make OUTDIR` writes 23 scenes on a close-up of skin (`skin_lamp_dry.json` without the flakes; with freckles and colour unevenness): no layer / the paper's 2 products × 5 thicknesses / 4 sphere fractions × 3 platelet roughnesses. With `--seed` fixed, the same pictures can be reproduced.

**Result of rendering them** (23 pictures. 640 px, 128 spp, seed 1. About 70 seconds per picture with a development build. `tools/studies/cosmetic_stimuli.py sheet OUTDIR` outputs two contact sheets and the numbers. The pictures are kept in `out/`):

(The following are pictures after smoothing was introduced, with smoothing 1.)

- **Thickness (the paper's 2 products):** matte 2 gets a darker brown with thickness; at thickness 2 the highlights and freckles vanish, and at thickness 4 the surface is uniform. Dewy 1 gets brighter (it was greyish before the platelet colour was fixed; see below), and a soft light of the layer spreads around the sharp highlight of the skin (shine broken up by the furrows). At thickness 2 the shine of the skin is no longer visible and only the light of the layer remains. The freckles almost vanish at thickness 4.
- **Finish (thickness 1, sphere fraction × platelet roughness):** lowering the sphere fraction (more platelets) makes the light stronger (at 50% it burns out to white). Raising the platelet roughness makes the light broader and weaker (0.05: small and sharp, 0.15: large and soft, 0.4: the whole surface hazes to whitish). At 100% spheres, the platelet knobs do not act and the result is matte. The two knobs act separately, as the paper says: "the fraction is the strength and the roughness is the shape".
- A thick layer of matte 2 shows oblique bands of light and dark. The platelets are tilted by 52° and their direction was aligned to the tangent of the shape (directional reflection of tilted platelets).
- **A problem seen, and its fix:** in the first version, the network of furrows appeared sharply in the shine of the shiny layer. The layer followed the bent normal of the skin relief as it was (real foundation fills the furrows and smooths the surface; the same kind of mistake as in the first prototype of moisturising). **It was fixed by giving the layer its own normal** (see "Smoothing" above).
- **Platelet colour of dewy 1 (fixed).** The first version used the platelet albedo (0.75, 0.85, 0.99) of the paper's table 2 as it is, and with a thick layer the skin shifted to grey. Compared with the measurement in the paper's figure 4 (dewy 1 is pinkish), this value is too bluish: against the background of the figure (about 232, 195, 166), the table's value gives a greyish yellow-green (191, 187, 170), while (0.99, 0.85, 0.75) with the channels reversed gives pink (209, 186, 159), and the ratios G/R and B/R match the measurement (reversing the diffuser colour Λd instead makes it blue, which does not match). I take the table's values to be a typo in the order blue, green, red, and the preset uses (0.99, 0.85, 0.75). **The basis is the comparison with the colours of the figure.** This is a point where I differ from the paper's values; the values of the other products are left as they are.

### 4.1 Primer (fills pores and furrows; 2026-11-12)

**Usage:** write `primer` in the skin material (`docs/SCENE_FILE.md`). It needs `relief`. Besides `thicknessUm`, the keys `ior`, `roughness` and `softnessUm` are accepted (defaults 1.40, 0.15 and 3 µm).

```jsonc
"skin": { "type": "skin", ..., "relief": {...},
          "primer": { "thicknessUm": 10 },
          "foundation": { "preset": "dewy1", "thickness": 1.0 } }
```

**Mechanism:** it uses the same form as the film that collects in furrows and pores (`SkinFilm`; `docs/SKIN_MOISTURE.md`). For the skin height h(p), the top surface of the film is z = smoothmax(h, L) (filled from the deep places up to the level L), and L is determined from the "mean thickness". The surface of the foundation layer follows the normal of this top surface (ridges stay as the skin, filled furrows are almost flat). `leveling` only flattens it further (the default is 0 when there is a `primer`). Without foundation, the shine of the skin follows the top surface of the film (index 1.40, roughness 0.15). **The thickness is the knob that decides "how much of the furrows is filled".**

**Check** (c3 of `--scene cosmetic`. PASS in both RGB and spectral): the tilt of the surface of the layer (default relief, mean) becomes flatter with thickness: 9.7 to 10.0° without primer → 3.7 to 4.6° at 5 µm → 0.1° at 20 µm → 0.01° at 60 µm. The albedos of `sample()` and `eval()` agree in all 4 cases (difference 0.002 or less).

**Result of rendering** (close-up of `skin_lamp_dry.json`, foundation dewy 1, thickness 1, 480 px, 64 spp. With this relief the furrows are filled completely at 10 µm): without primer (`leveling` 0), the network of furrows is visible in the highlight. At 10 µm and 30 µm the network disappears and the light source is reflected round and soft. With primer only and no foundation (30 µm), the freckles remain visible, the light source is reflected small and sharp, and the shine looks as if covered with glass.

**Basis:**
- The action of a primer (a silicone elastomer fills dents, a spherical powder scatters light, pigments correct the colour) comes from summaries of search results of manufacturers' descriptions, patents and reviews (no peer-reviewed paper was found).
- Index 1.40, roughness 0.15 and edge rounding 3 µm are assumed values.

### 4.2 Coloured primer (2026-11-12)

**Usage:** add `color` and `opacity` to `primer` (`opacity` without `color` is an error; the default of `opacity` is 0.5).

```jsonc
"primer": { "thicknessUm": 10, "color": [0.55, 0.72, 0.48], "opacity": 0.6 }   // a green primer
```

**Model:** the primer is made a Kubelka-Munk layer (two-flux) (`Materials/KubelkaMunkLayer.hpp`). `color` = the colour R∞ of an infinitely thick film (the "intrinsic colour" of Ikkai 2023), `opacity` = X = S·D (scattering coefficient × thickness). Per channel, with a = (R∞ + 1/R∞)/2 and b = √(a² − 1), the reflection of the layer alone is R = sinh(bX)/(a sinh bX + b cosh bX) and the transmission T = b/(a sinh bX + b cosh bX). The diffuse colour A of the skin below the layer becomes **R + T²A/(1 − RA)** (Kubelka 1954; the same as the two-layer formulas of Doi et al. 2005 and Yoshida and Okiyama 2021).

**How it enters the skin:** it uses the same path as the pigment (`pigment`): it becomes the factor `albedoTint` multiplied onto the diffuse colour of the skin (with subsurface scattering, a square root each on the way in and on the way out). The foundation above sees the new colour as the skin. **Only the diffuse colour changes**: the shine of the skin (sebum) still follows the normal of the top surface of the primer (section 4.1), and the spread of subsurface scattering stays that of the skin.

**Check** (c4 of `--scene cosmetic`. PASS in RGB): with no layer (X = 0), R = 0 and T = 1; at X = 60, R = R∞ and T = 0; with no absorption, R + T = 1; **a layer of thickness 0.7 on a layer of thickness 0.3 = a layer of thickness 1.0** (Kubelka additivity; difference 10⁻⁷ or less); the colour factor of the skin follows the formula.

**Result of rendering** (close-up of `skin_lamp_dry.json`, redness strong (`rednessStrength` 1.0), primer 10 µm, colour (0.55, 0.72, 0.48), 480 px, 64 spp. The numbers are log(R/G) over 16 × 16 blocks excluding the highlight):

| | Mean | Spread (standard deviation) |
|---|---|---|
| No primer | 0.386 | 0.0200 |
| Colourless primer | 0.402 | 0.0203 |
| Green, opacity 0.3 | 0.234 | 0.0110 |
| Green, opacity 0.6 | 0.121 | 0.0084 |
| Green, opacity 1.2 | −0.015 | 0.0066 |

- The thicker it is, the more the redness decreases and the more the unevenness of the redness (and the darkness of the freckles) is hidden. At opacity 0.6, the redness is about one third and the unevenness about 40%.
- **At opacity 1.2 the skin becomes grey-green** (R ≈ G). This is in the same direction as the warning in beauty articles that green on skin without redness looks unhealthy.
- A colourless primer changes neither the colour nor the unevenness (the furrows are filled and the shine only becomes round).

**Basis:**
- Two-flux is a theory of diffuse light. The colour of the layer is handled by the conversion from RGB to spectrum (with a raised floor).

## 5. Differences from the paper

- Colour uses 3 RGB channels (the paper's measurements are spectral, but the fit is linear RGB).
- **The assembly on skin is an approximation that differs from the paper** (section 4).
- When a white layer (flakes, lotion) is put on skin with subsurface scattering, red looked about 1 to 2% darker. **Resolved (2026-11-13):** the 1 to 2% was fluctuation of the check (26 000 paths); measured with 100 times as many paths, the difference in red is −0.05% ± 0.2%. The real components are two: (1) For red, the gain of the white layer is small, +0.7%, and with subsurface scattering it cancels to near 0 by the part where the entry and exit fall on different patches (about −0.7%). (2) The lobes around the bent normal discard samples below the horizon and lose about 1.2%, but the exit lobe of subsurface scattering uses the smooth normal and does not lose. A white layer moves light from the exit lobe to the flakes lobe, so it becomes darker by the moved fraction × 1.2% (this is why a residue remains even when the weight of the coat is set to 0). `docs/SKIN.md` section 6, `out/studies/flakes_sss`. For the foundation layer, the difference between with and without subsurface scattering is 0.003 or less and invisible (8 of `--scene cosmetic`). 11c and 13c of `--scene relief` were rebuilt from fixed ratios into three tests: bookkeeping (no change on a white base), local model (an exact predicted value), direction (lower bound = −ρ × horizon loss), each judged with a width of 4 standard errors.
- The diffuse lobe without subsurface scattering was aligned to a formula that counts the flakes and lotion layers up to the exit transmission (1 − ρ) and the back-and-forth 1/(1 − ρR) (2026-11-11. Before, only the entry (1 − ρ) was counted, and the brightness added by the white layer was estimated 2 to 4 times too large). The added brightness agrees with that of the path with subsurface scattering (output of 11c and 13c of `--scene relief`).

## 6. Next stage

| Stage | Content |
|---|---|
| C | Measure the stimuli: with a sweep of thickness × finish (c_d, α_p), on a patch lit from the front, turn the hiding power (how visible the colour unevenness is), the reduction of light returning from inside the skin, and the highlight into numbers. What to use as indicators will be chosen to match the quantities used in research on skin material appearance |

**Note on wording (2026-11-10):** the "toumeikan" of skin (a clear, fresh look of skin; a kansei quality, not physical transparency) is a kansei quality of skin seen from the front, not that light really passes through (the see-through look, "sukekan"). An earlier version put "light that passes through (a cone of backlight)" in stage C and wrote "compatibility of hiding power and toumeikan" as a candidate research question. A backlit condition that makes skin see-through cannot be measured on human skin and is not used in experiments on skin material appearance, so it was removed. What the layer calculation can say is limited to this: the thicker the layer, the more the colour unevenness is hidden, and the less light returns from inside the skin (the paper also writes "both shiny and matte reduce the effect of subsurface scattering of the skin a little", but gives no numbers).

## 7. References

See [REFERENCES.md](REFERENCES.md), group "Foundation make-up".
