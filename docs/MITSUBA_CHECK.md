# Cross-check against Mitsuba 3

The self-checks of the engine compare it with formulas and with small Monte Carlo codes of its own. This document compares it with another renderer: the same scenes rendered by GoLD_rayt and by Mitsuba 3 (version 3.7.1, the `scalar_rgb` and `scalar_spectral` variants), compared within the Monte Carlo noise. The comparison found one bug in GoLD_rayt (now fixed), one modelling choice that differs, and two places where Mitsuba is the one that differs (from an independent reference, and from its own color matching table).

## 1. Conclusions

- **Diffuse surfaces, a uniform environment, smooth and rough glass, and scattering media inside smooth, rough and index-matched boundaries agree** within the noise (0.03 to 0.4 % on the sphere, depending on the scene).
- **Rough glass was 0.5 % too bright in GoLD_rayt: a bug, fixed.** `Dielectric::eval` and `pdf` returned values for pairs of directions that no microfacet can connect (the generalized half vector back-facing to one of them). Only light sampling reaches such pairs, so they added light. After the fix the rough glass sphere agrees within 0.03 % (section 3.1).
- **Rough metal differs at the rim of the sphere by the form of the Smith masking term:** GoLD_rayt uses the height-correlated form, Mitsuba the separable one. A build of GoLD_rayt with the separable form agrees within 0.05 % (section 3.2).
- **Dense media agree with an independent random walk.** For a dense gray medium in an index-matched sphere, GoLD_rayt is within 0.06 % of a plain random walk written for this check; Mitsuba is 0.5 % low (section 3.3).
- **The spectral builds agree too** (section 4): measured gold (Johnson & Christy n, k) under a candle, a reflectance spectrum under tungsten, glass, rough glass and a medium, within 0.1 % once a 0.1 % offset of Mitsuba's spectral film is removed. Rough gold differs at the rim by the Smith form, as in RGB.

## 2. What was compared

`tools/test/mitsuba_compare.py` describes each scene once and writes it as a GoLD_rayt scene file and as a Mitsuba scene dictionary. A scene is a sphere (radius 0.3 m) on a diffuse floor (albedo 0.5) under a square area light (radiance 12) and, in most scenes, a uniform environment (radiance 0.2); the furnace scenes have the sphere alone in an environment of radiance 1. 320 × 240 pixels, the RGB builds of both renderers.

| Scene | Sphere |
|---|---|
| diffuse, env | diffuse (albedo 0.8, 0.5, 0.3) under the light, or under the environment alone |
| conductor | rough metal, GGX α 0.25, n 0.2, k 3.5 (the same in every channel) |
| glass, roughglass | glass n 1.5, smooth or GGX α 0.09 |
| glassfurnace, roughglassfurnace | the same glass alone in a uniform environment |
| medium, roughmedium, matchedmedium | a scattering medium (σt 20 /m, albedo 0.99, 0.9, 0.7, Henyey-Greenstein g 0.5) inside a smooth, rough (α 0.09) or index-matched (n 1) boundary |
| densemedium, roughdensemedium, densegray, densegraymatched | a dense medium like wax (σt 200 /m, albedo 0.995, 0.98, 0.9 or gray 0.98, g 0.3) |
| densewhitefurnace, densegrayfurnace | an index-matched sphere of the dense medium (albedo 1 or 0.98) alone in a uniform environment |

The conventions of the two renderers were matched: a vertical field of view; a box pixel filter; GGX with α = roughness² (the GoLD_rayt `roughness` is √α); the multiple-scattering compensation of the GoLD_rayt conductor off (Mitsuba has none); white lights; no limit on the length of a path in either renderer; a null BSDF in Mitsuba for an index-matched boundary. Mitsuba uses `volpathmis` for media and runs **without Russian roulette** (section 3.3).

Each renderer draws a scene twice with independent seeds. The noise of each pixel comes from the difference of the two renders, and the comparison uses the mean over the pixels of the sphere (with its standard error) and 16 × 16 pixel blocks (the share of blocks whose difference is more than 3 standard errors; about 0.3 % if the two images converge to the same one).

## 3. Results

| Scene | Mitsuba / GoLD_rayt on the sphere | |
|---|---|---|
| diffuse | 0.9999 ± 0.0002 | agrees |
| env | 1.0001 ± 0.0001 | agrees |
| conductor | 0.9987 ± 0.0005 | Smith form (3.2) |
| glass | 0.9993 ± 0.0008 | agrees |
| roughglass | 0.9997 ± 0.0008 | agrees (0.9951 before the fix) |
| glassfurnace | 1.0000 (both exactly 1) | agrees |
| roughglassfurnace | 0.9998 ± 0.0001 | agrees (0.9935 before the fix) |
| medium | 1.0006, 1.0024, 1.0054 ± 0.002 to 0.003 | agrees |
| roughmedium | 1.0011, 1.0003, 0.9991 ± 0.0018 | agrees |
| matchedmedium | 1.0013, 1.0008, 0.9994 ± 0.0015 | agrees |
| densemedium (256 spp) | 0.9987, 0.9977, 0.9979 ± 0.004 | agrees |
| roughdensemedium (256 spp) | 0.9970, 0.9936, 0.9937 ± 0.0034 | agrees within 2 standard errors |
| densegray (256 spp) | 0.9977 ± 0.0037 | agrees |
| densegraymatched (256 spp) | 0.9951 ± 0.0023 | Mitsuba low (3.3) |
| densewhitefurnace (256 spp) | 0.9997 ± 0.0003 (both 1) | agrees |
| densegrayfurnace (256 spp) | 0.9942 ± 0.0004 | Mitsuba low (3.3) |

The scenes without a dense medium use 1024 samples per pixel. The three numbers are R, G and B where they differ. Full output: [mitsuba_compare_results.txt](data/mitsuba_compare_results.txt).

### 3.1 Rough glass: back-facing microfacets (a bug in GoLD_rayt, fixed)

Before the fix the rough glass sphere was 0.5 % brighter in GoLD_rayt, over the whole sphere, and 0.65 % in the furnace. For a refraction, the microfacet that connects two directions has the generalized half vector `wo + η wi`. A microfacet is a facet of the surface, so this normal must face the outside of the surface, and both directions must meet its front side. Some pairs of directions give a half vector for which this is impossible. pbrt-v4 and Mitsuba return 0 for them; GoLD_rayt turned the half vector toward `wo` and evaluated the formula anyway (its GGX is even in the half vector). Integrated over the directions, these impossible pairs held 0.1 to 1.2 % of the light near normal incidence and 5 to 7 % near grazing incidence (glass, α 0.09). Sampling the BSDF never produces such a pair, but sampling a light does, and then the light was added. `Dielectric::eval` and `pdf` now orient the half vector to the side of the normal and return 0 when it is back-facing to either direction, as the layered coat (`LayeredCoat`) already did. The self-checks (`volume`, `furnace`, `dispersion`, `ocean`, `lights`) pass in both builds.

### 3.2 Rough metal: the form of the Smith masking term (a choice, not a bug)

The microfacet model needs the share of microfacets that both directions see. GoLD_rayt uses the height-correlated form `G = 1 / (1 + Λ(wo) + Λ(wi))`, as pbrt-v4 does; Mitsuba uses the separable form `G1(wo) G1(wi)`, which counts the two directions as independent and therefore shadows more. Both are common; the height-correlated form is the more accurate model of a random surface (Heitz 2014). The difference is under 0.1 % up to 50° and grows toward grazing (0.7 % at 70°, 2.5 % at 80° for α 0.25), so the rough metal sphere differs only at its rim (0.8 % in the outer tenth of the radius). A build of GoLD_rayt with the separable form agrees with Mitsuba within 0.05 % (1.0005 ± 0.0005). For rough glass at α 0.09 the two forms differ by under 0.05 %.

### 3.3 Dense media: Mitsuba's Russian roulette, and an independent reference

With its default settings Mitsuba rendered the dense media darker than GoLD_rayt. Two causes were separated:

- **Russian roulette.** Mitsuba applies the roulette at every interaction, scatterings in a medium included, and keeps a path with probability min(throughput, 0.95). In a dense medium of high albedo a path scatters thousands of times with a throughput near 1, so it survives with about 0.95ⁿ and the survivors carry a weight of 1/0.95ⁿ. The estimate stays unbiased in expectation, but the survivors are too rare to be sampled, and a finite render comes out dark: the sphere of `densewhitefurnace`, which must be exactly 1, came out at 0.906. GoLD_rayt applies its roulette at surface bounces only, not to collisions in a medium, and gives 1.0000. The comparison therefore runs Mitsuba without roulette (then 0.9997).
- **An independent reference.** In a uniform environment of radiance 1, an index-matched sphere of a medium shows in each pixel the probability that light entering along the pixel's ray leaves the sphere again. `tools/test/medium_sphere_reference.py` computes this with a plain random walk (exponential free flights, Henyey-Greenstein directions, the albedo as a weight), with the camera and the pixel footprint of the scene and neither renderer's code. For the dense gray medium (albedo 0.98, g 0.3) the reference is 0.70146 ± 0.00016; GoLD_rayt gives 0.70186 ± 0.00022 (ratio 1.0006) and Mitsuba without roulette 0.69778 ± 0.00022 (ratio 0.9947). A run by hand with g = 0 and g = −0.3 gave the same picture (GoLD_rayt 1.0006 and 1.0006, Mitsuba 0.9933 and 0.9927), and the deficit lies over the whole sphere, not at its rim. Where Mitsuba loses this half percent was not investigated; it showed with the index-matched (null) boundary, and with a glass boundary the comparison is within its noise of about 0.4 %.

## 4. Spectral builds

The spectral build of GoLD_rayt against Mitsuba's `scalar_spectral` variant, on seven scenes of the same layout (1024 samples per pixel). The two renderers turn RGB values into spectra by different methods, which would be a difference of the inputs rather than of the rendering, so every colored input is given as a spectrum on a 1 nm grid to both: the lights as tables (blackbodies of 1850 K and 2856 K, D65) scaled to the same luminance, a smooth orange reflectance spectrum, and the n, k of gold from the Johnson & Christy table that GoLD_rayt uses, interpolated linearly. Gray values are flat spectra in both, and the environment is GoLD_rayt's D65 table. The sphere is compared in XYZ, whose channels stay positive (the orange sphere under tungsten lies outside sRGB, where a linear sRGB channel is near 0).

| Scene | Mitsuba / GoLD_rayt on the sphere (X, Y, Z) | |
|---|---|---|
| white furnace, glass sphere in a D65 environment | 0.9999, 1.0000, 0.9999 ± 0.0002 | agrees |
| diffuse sphere, orange reflectance spectrum, tungsten | 0.9997, 1.0001, 0.9993 ± 0.0003 | agrees |
| smooth gold, candle + D65 environment | 1.0003, 1.0002, 0.9999 ± 0.0005 | agrees |
| rough gold (α 0.25), candle + D65 environment | 0.9979, 0.9979, 0.9971 ± 0.0005 | Smith form (as in 3.2) |
| smooth glass, D65 light + environment | 1.0000, 1.0001, 1.0001 ± 0.0009 | agrees |
| rough glass (α 0.09), tungsten + D65 environment | 0.9989, 0.9991, 0.9991 ± 0.0009 | agrees |
| gray medium (σt 20 /m, albedo 0.95) in glass, tungsten | 0.9995, 0.9995, 1.0001 ± 0.002 | agrees |

- **Mitsuba's spectral film reads 0.1 % high.** Before this correction every scene came out 0.1 % brighter in Mitsuba, including the environment seen directly (Y 1.00100 against 1.00000 for an environment of luminance 1 by Mitsuba's own color matching functions). Mitsuba divides the luminance by `MI_CIE_Y_NORMALIZATION` = 1 / 106.7503, while its own ȳ table integrates to 106.857 over 360 to 830 nm, so every spectral render reads 1.001000 times the luminance of its spectrum. The comparison multiplies Mitsuba's spectral images by 0.999001. GoLD_rayt normalizes with its own table, so a light of luminance 1 reads exactly 1.
- **Rough gold differs by the Smith form**, as the rough metal in RGB (3.2). The difference is larger than for the gray metal there, most in Z: gold's reflectance of blue light rises toward grazing angles, where the two forms differ most. A spectral build with the separable form agrees: 0.9995, 0.9997, 1.0002 ± 0.0005.
- The spectral build was not compared on chromatic media (both renderers would need spectral coefficients) or dispersive glass (Mitsuba's dielectric takes one index).

Full output: [mitsuba_compare_spectral_results.txt](data/mitsuba_compare_spectral_results.txt).

## 5. Limits

- One camera and one kind of scene. Textures, meshes, point lights, thin films, layered materials, hair, skin and the BSSRDF are not covered: Mitsuba has no direct counterpart for several of them.
- The dense media were compared at 256 samples per pixel (Mitsuba without roulette takes about 2.5 minutes per render at that count on this PC); their noise on the sphere is about 0.4 %.
- Mitsuba's GPU variant (`cuda_ad_rgb`) renders simple scenes on this PC's Quadro M4000 but failed to compile the medium scenes, so the CPU variant was used throughout.

## 6. Reproduction

Mitsuba 3 for Python (`pip install mitsuba`; version 3.7.1 here) and the RGB and spectral builds of GoLD_rayt.

```bash
# all scenes without a dense medium (about 20 minutes)
python tools/test/mitsuba_compare.py --exe out/cli-dev/GoLD_rayt.exe --spp 1024 --scenes diffuse,env,conductor,glass,roughglass,glassfurnace,roughglassfurnace,medium,roughmedium,matchedmedium
# the dense media (about 25 minutes)
python tools/test/mitsuba_compare.py --exe out/cli-dev/GoLD_rayt.exe --spp 256 --scenes densemedium,roughdensemedium,densegray,densegraymatched,densewhitefurnace,densegrayfurnace
# the independent reference, compared with the renders of the run above (about 10 minutes)
python tools/test/medium_sphere_reference.py --scene densegrayfurnace --paths 512 --renders out/studies/mitsuba_compare/densegrayfurnace
# the spectral scenes (about 20 minutes); --reuse reads renders already in the work folder instead of rendering again
python tools/test/mitsuba_compare.py --exe out/cli-dev-spectral/GoLD_rayt.exe --variant scalar_spectral --spp 1024 --work out/studies/mitsuba_compare_spectral
```

The renders and the per-block z scores are kept in `out/studies/mitsuba_compare/<scene>/` (spectral: `out/studies/mitsuba_compare_spectral/<scene>/`).
