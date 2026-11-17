# Fading: colours that change with the light they received

A diffuse colour can be shown as it looks after some hours of light, from the lab's fading measurements. Three ways exist, from the most to the least measured: the measured rounds, a fading law fitted to them, and a colour alone.

## 1. The measurements

The lab's fading experiment (Kataoka, master's thesis, Yokohama National University, 2025, chapter 2) exposed coloured samples to pseudo-sunlight (SXL-500V2 at full power, 371.76 W/m²) and measured their reflectance with a spectroradiometer every few rounds:

| Series | Samples | Rounds |
|---|---|---|
| Dyed drawing paper | 9 colours | 0 to 150 h, every 30 h |
| Dyed origami paper | 9 colours | 0 to 150 h, every 30 h |
| Watercolour pigments on canvas | 9 colours | 0 to 900 h, every 150 h |

The data are lab data and are not in the repository. `tools/data/fading_refangle.py` converts the lab's tables into a series folder (`hours.csv`, `r1/<colour>.csv`, ...).

## 2. Three ways to fade a colour

| Way | Scene file | Time | When |
|---|---|---|---|
| **Measured rounds** | `"fading": { "series", "color", "hours" }` | Linear blend of the two rounds around the time; held at the ends | The measured colour, inside the measured time |
| **Fitted law** | the same with `"fit": "km"` | Any time, also before the first round (negative hours) and after the last | The measured colour, smooth in time, extrapolated |
| **Colour alone** | `"fading": { "rgb" or "image", "model", "hours", "fastness" }` | The table of the dye concentration | A colour that was never measured (an approximation) |

All three take the same exposure over the surface (`map`, `mat`, `gradient`, `noise`), and `"outdoorDays"` / `"outdoorYears"` in place of `"hours"`. The keys are in [SCENE_FILE.md](SCENE_FILE.md).

## 3. The fitted law

Per wavelength, the Kubelka–Munk K/S of the colorant decays exponentially toward what is left when it is gone:

  K/S(t) = K/S∞ + (K/S₀ − K/S∞) e^(−k t),  R = 1 + K/S − √((K/S)² + 2 K/S)

K/S₀ comes from the first round. The rate k(λ) and the remainder K/S∞(λ) are smooth curves, fitted to all rounds by `tools/data/fading_km_fit.py <series>`, which writes `<series>/km/<colour>.csv`.

- **Fit to the rounds (CIEDE2000, D65):** 0.5 to 0.9 for the dyed paper, 0.5 to 1.1 for most pigments (red 2.5, vermilion 1.9), 1.0 to 2.7 for the origami paper (whose rounds scatter).
- **Per colour, not shared.** Predicting a colorant's fading from its first colour misses by 5 to 6 with any shared law: the thesis's principal-component law (5.4 as published), this law with shared rates, or coupled coefficients. How fast a colorant fades and how its spectrum changes depend on its chemistry, not on its colour, and 9 colours are too few to learn them.
- **Extrapolation.** After the last round the colour approaches its remainder. Before the first round the colorant gets stronger; going back as far as the whole measured span already oversaturates some colours (a yellow turns orange). A change finished before the first measured time cannot tell its rate, so k is capped at 4 / (first interval).

## 4. Time outdoors

One day in the device does the damage of **2.911 days of direct sunlight outdoors** (the thesis, section 6.1: the annual direct irradiation at Tsukuba in 2022, 2692.77 MJ/m², against the device, with ultraviolet weighted 0.95 and visible light 0.05). Heat, humidity, weather and coatings are left out: an order of magnitude, not a prediction.

## 5. Checks

`--scene metals` sections 5 (rounds and exposure), 5c (the fitted law: the first round at 0 h, the exponential in K/S, the remainder, a stronger colorant before 0 h, finite far back and far ahead, the sampled law for surfaces with varying exposure, the outdoor conversion) and 6 (colour alone). They pass in both builds.
