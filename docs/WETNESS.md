# How wet things look: image statistics, and renders (2026-11-09)

**Status:** the existing renders were measured with the statistics of the paper. The model of wet skin is in `docs/SKIN_MOISTURE.md`. Reproduce with `tools/studies/wetness_stats.py`.

## 1. The paper

Sawayama, Adelson, Nishida, "Visual wetness perception based on image color statistics", *Journal of Vision* 17(5):7, 2017 (full text and figures read).

**Why wetting changes the appearance** (a summary of the paper's figure 1):

1. On a dry rough surface, the Fresnel reflection at the surface is a white haze. A water film is smooth, so the haze is replaced by sharp highlights.
2. Water has a refractive index closer to the material's than air has. Light enters the material more easily and leaves it more easily.
3. Light that is totally reflected inside the water surface enters the material again. The path through the material gets longer, and the surface becomes darker and more saturated.
4. In porous materials (cloth, sand, paper), water fills the gaps, and light passes through several particles.

**Dry and wet photographs of 32 materials (clay, soil, cloth, wood, vegetables)** (as far as read from figure 3):

| Statistic | Dry photograph | Wet photograph |
|---|---|---|
| Mean saturation (u′v′) | 0.1 to 0.4 (mostly 0.1 to 0.2) | 0.15 to 0.6 (1.43 times on average; about 1 to 3 times depending on the material) |
| Hue | - | almost unchanged |
| Mean luminance | 0.4 to 0.7 | 0.2 to 0.5 |
| Standard deviation of luminance | 0.04 to 0.1 | 0.04 to 0.1 (no systematic change) |
| Skew of luminance | -1 to +0.5 | -1 to +4 |

If a wet surface is photographed with the specular reflection removed, it only becomes darker: **it is physically wet but does not look wet.**

**Perception experiment:** an image operation (the WET operator) that sets the luminance to Y^γ (γ = 3) and raises the saturation makes dry photographs look wet. Saturation and skew act separately, and lowering the mean luminance alone does not make a picture look wet. The operator works better on images with more hues.

**What the paper does not cover:** water drops, puddles, changes of shape. The stimuli of the experiments are a photographic texture and an overlap of square plates seen from the front (called a "3D Mondrian", but the faces are flat and have no shading from curvature); there is no solid with a smooth curved surface.

## 2. What was measured

- **Statistics:** the u′v′ saturation s_uv = 13 √((u′ − u′ₙ)² + (v′ − v′ₙ)²) (white point D65), the hue (the angle in the u′v′ plane), the mean, standard deviation, skew and kurtosis of the luminance, and the hue entropy (512 bins, in bits).
- **Exposure:** white is 1.25 times the 99.5th percentile of the luminance of the dry picture. The dry and wet pictures of a pair share one exposure, and values above white are clipped (as in a photograph). So that render noise does not add to the standard deviation and the saturation, 2 × 2 pixels are averaged before measuring.
- **Scene:** 240 spheres (pebbles) half sunk into a grey ground. A large dim light overhead and a spherical light source of 30 cm diameter (at about 1 m). The surface is diffuse only (`lambert`), a transparent layer of index 1.5 and roughness 0.6 over a diffuse body (`rough`), or a layer of index 1.33 and roughness 0.03 over a diffuse body (`wet`, a water film). The palettes are hues all around the circle (`many`), a single hue (`one`), and bright, pale, warm colours (`natural`). 1024 px, 128 spp.

## 3. Results (change from dry to wet)

| Pair | Saturation | Hue | Mean luminance | Skew | Chroma C\*uv |
|---|---|---|---|---|---|
| Paper (32 materials) | ×1.43 on average | unchanged | lower | higher | - |
| Pebbles `many`: lambert → wet | ×1.08 | −0.3° | ×0.65 | +1.2 | ×0.89 |
| Pebbles `one`: lambert → wet | ×1.07 | −1.6° | ×0.64 | +1.2 | ×0.88 |
| Pebbles `natural`: lambert → wet | ×1.26 | −1.7° | ×0.73 | +0.6 | ×1.10 |
| Pebbles `many`: rough → wet | ×1.14 | −1.0° | **×1.29** | +0.05 | ×1.28 |
| Skin, small light source, water film 5 µm (`SKIN_MOISTURE.md` 5.7) | ×1.11 | +0.6° | ×0.98 | +2.1 | ×1.10 |
| Skin, area light, just after lotion (same, 5.5) | ×1.07 | +0.4° | ×1.00 | −0.2 | ×1.07 |

- **A water film on a purely diffuse surface moves all four statistics in the same direction as the paper's photographs.** The size depends on the brightness of the base: saturation rises 1.08 times for dark, strongly coloured pebbles and 1.26 times for bright, pale ones. In the model of the water film, light reflected inside the film enters the base again (mechanism 3), so the brighter the base, the more the saturation rises.
- **My first conclusion, "the saturation rise is too small; a mechanism is missing", was wrong.** The first pebbles (`many`) had saturation 0.32, mean luminance 0.31 and standard deviation 0.16: darker, more colourful and more contrasty than the paper's dry photographs. `natural`, with colours closer to the photographs, falls inside the paper's range (about 1 to 3 times depending on the material).
- **If the dry surface is represented by a "rough transparent layer", wetting makes it brighter** (×1.29, the opposite of reality). A layer of index 1.5 reflects more light internally than water (1.33) does. It is not suitable as a model of a dry surface.
- **Reducing the hues to one leaves the ratios almost unchanged** (hue entropy 8.7 against 4.8 bits), as the paper says: the physical change does not depend on the number of colours.
- **For skin under an area light, only the saturation moves and the skew does not rise.** This corresponds to the paper's condition that without highlights a surface does not look wet even if it is wet.

## 4. Image processing alone does not make a solid look wet

A picture of dry pebbles processed with the WET operator (γ = 3, saturation +0.26, the paper's values) lies, in the statistics, further on the "wet" side than the physically wetted picture (`natural`: saturation 0.46 against 0.26, skew 1.95 against 0.77). Yet it looks like **the lighting was only made darker**, not wet.

- A transformation applied to the whole picture does not know "where the light source should be reflected". On a smooth sphere, the brightest part is the broad lit side; raising the tone curve only strengthens the shading.
- In the paper the transformation worked on flat textures with fine bright points all over, where the emphasised points read as small highlights wherever they are. In the paper's regression analysis too, the more a picture varied in luminance at fine spatial scales, the better the transformation worked.
- Therefore, **for stimuli where the shape determines the positions of shading and highlights, wetness has to be made physically**. The statistics are used as a yardstick of whether the picture has "the cues people read as wetness" (not as a sufficient condition).

## 5. References

- M. Sawayama, E. H. Adelson, S. Nishida, "Visual wetness perception based on image color statistics", *Journal of Vision* 17(5):7, 2017. (Full text read.)
- J. Lekner, M. C. Dorf, "Why some things are darker when wet", *Applied Optics* 27(7), 1988. (Cited second-hand from the paper above. Not read.)
- S. A. Twomey, C. F. Bohren, J. L. Mergenthaler, "Reflectance and albedo differences between wet and dry surfaces", *Applied Optics* 25(3), 1986. (Same. Not read.)
- H. W. Jensen, J. Legakis, J. Dorsey, "Rendering of wet materials", *Rendering Techniques '99*, 1999. (Same. Not read.)
