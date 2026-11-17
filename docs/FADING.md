# Fading: colours that change with the light they received

A diffuse colour can be shown as it looks after some hours of light, from measurements of samples that faded. Three ways exist, from the most to the least measured: the measured rounds, a fading law fitted to them, and a colour alone.

## 1. The measurements

The method follows the lab's fading experiment (T. Kataoka and K. Okajima, "Aging Simulation due to Sunlight Based on Spectral Reflectance", Color Science Association of Japan, 55th annual meeting, 2024): coloured papers were faded under pseudo-sunlight (SXL-500V2) and their spectral reflectance was measured every 30 h with a 2D spectroradiometer (SR-5100), as the ratio of the sample's radiance to that of a standard white. The data belong to the lab and are not in the repository.

A series folder holds `hours.csv` (the exposure time of each round) and `r1/<colour>.csv`, `r2/<colour>.csv`, ... (reflectance curves). `tools/data/fading_refangle.py` converts the lab's tables into this form.

## 2. Three ways to fade a colour

| Way | Scene file | Time | When |
|---|---|---|---|
| **Measured rounds** | `"fading": { "series", "color", "hours" }` | Linear blend of the two rounds around the time; held at the ends | The measured colour, inside the measured time |
| **Fitted law** | the same with `"fit": "km"` | Any time, also before the first round (negative hours) and after the last | The measured colour, smooth in time, extrapolated |
| **Colour alone** | `"fading": { "rgb" or "image", "model", "hours", "fastness" }` | The table of the dye concentration | A colour that was never measured (an approximation) |

All three take the same exposure over the surface (`map`, `mat`, `gradient`, `noise`). The keys are in [SCENE_FILE.md](SCENE_FILE.md).

## 3. The fitted law

Per wavelength, the Kubelka–Munk K/S of the colorant decays exponentially toward what is left when it is gone:

  K/S(t) = K/S∞ + (K/S₀ − K/S∞) e^(−k t),  R = 1 + K/S − √((K/S)² + 2 K/S)

K/S₀ comes from the first round. The rate k(λ) and the remainder K/S∞(λ) are smooth curves, fitted to all rounds by `tools/data/fading_km_fit.py <series>`, which writes `<series>/km/<colour>.csv`.

- **Fit to the rounds.** On most colours the law reproduces the measured rounds within about 1 CIEDE2000.
- **Per colour, not shared.** The 2024 study modelled all colours with shared principal components and shared time laws. Predicting a colour that was left out of the fit, that approach, this law with shared rates, and coupled coefficients all miss by several CIEDE2000: how fast a colorant fades, and how its spectrum changes, depend on its chemistry rather than its colour, and a handful of colours cannot teach that. So the law is fitted per colour, and only for colours measured over time.
- **Extrapolation.** After the last round the colour approaches its remainder. Before the first round the colorant gets stronger; going back as far as the whole measured span already oversaturates some colours (a yellow turns orange). A change finished before the first measured time cannot tell its rate, so k is capped at 4 / (first interval).

## 4. Checks

`--scene metals` sections 5 (rounds and exposure), 5c (the fitted law: the first round at 0 h, the exponential in K/S, the remainder, a stronger colorant before 0 h, finite far back and far ahead, the sampled law for surfaces with varying exposure) and 6 (colour alone). They pass in both builds.
