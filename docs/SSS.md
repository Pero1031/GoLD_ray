# Subsurface scattering (SSS): implementation notes

Notes on the implementation of materials such as skin, where light enters through the surface, scatters inside, and leaves at a **different point**.
They describe the current implementation of GoLD_rayt (as of 2026-11).

- [1. Overview](#1-overview)
- [2. Radial profile](#2-radial-profile)
- [3. Material (CoatedDiffuse)](#3-material-coateddiffuse)
- [4. Sampling the exit point (probe rays)](#4-sampling-the-exit-point-probe-rays)
- [5. Lobe at the exit point](#5-lobe-at-the-exit-point)
- [6. Flow of the integrator](#6-flow-of-the-integrator)
- [7. Spectral mode](#7-spectral-mode)
- [8. Usage and parameters](#8-usage-and-parameters)
- [9. Verification](#9-verification)
- [10. Computational cost](#10-computational-cost)
- [11. File list](#11-file-list)
- [12. Walking inside the object](#12-walking-inside-the-object-randomwalk-2026-11-10)

---

## 1. Overview

### 1.1 What is approximated

The process in which light enters at a point `po` and leaves at a point `pi` is expressed with the BSSRDF (bidirectional scattering-surface reflectance distribution function) `S`.

```
S(po, wo, pi, wi) = (1 − F(wo)) · Sp(|po − pi|) · Sw(wi)
```

| Term | Meaning | Implementation |
|---|---|---|
| `1 − F(wo)` | fraction that passes the surface coat and enters | lobe selection of `CoatedDiffuse` |
| `Sp(r)` | fraction of the light that comes out at distance `r` from the entry point (radial profile) | `Subsurface.hpp`, `SubsurfaceProbe.hpp` |
| `Sw(wi)` | directional distribution at the exit point (a cosine lobe including the Fresnel transmission) | `SubsurfaceExit` |

`Sp` is assumed to depend only on distance (not on the orientation or shape at the points). This is an approximation that assumes a "flat surface".

### 1.2 How it is handled inside path tracing

SSS is one **event** in the path tracing. When the event "entered the inside" is chosen at the surface,
the exit point is searched for stochastically, and the path continues from that point.

```
 camera → entry point po ──(SSS event)──→ exit point pi → light sample (NEE) / next direction → …
            │
            └─(specular lobe) reflection at the coat → next direction as usual
```

### 1.3 Why this method was chosen

| Method | Adopted? | Reason |
|---|---|---|
| Random walk (actually walking inside) | not at first | at the time there was no medium mechanism (`Medium`) yet. The skin of a PMX model is not a closed mesh. **There is now a mode for closed shapes (section 12)** |
| **BSSRDF + probe rays** | **adopted** | works even if the mesh is not closed. BVH and Embree can be used as they are. The standard method for skin |
| Image-space blur | rejected | inconsistent with path tracing |

---

## 2. Radial profile

### 2.1 Burley's normalized diffusion

For each channel `c` (RGB; 4 wavelengths in spectral mode), the area density of the exit point is expressed by

```
Sp_c(r) = A_c · P(r; d_c)

P(r; d) = ( e^(−r/d) + e^(−r/(3d)) ) / ( 8π · d · r )
```

`P` integrates to 1 over the whole plane (`∫ P · 2πr dr = 1`). Therefore the integral of `Sp_c` over the whole surface
is `A_c`, the diffuse color (albedo). Energy is conserved.

### 2.2 How the parameters are determined

```
d_c = ℓ_c / s(A_c)          s(A) = 1.85 − A + 7·|A − 0.8|³
```

| Symbol | Meaning | Origin of the value |
|---|---|---|
| `A_c` | diffuse color (texture color at the entry point x material color) | material. The value at the **entry point** is used |
| `ℓ_c` | scattering distance (mean free path) | default: red 3.67, green 1.37, blue 0.68 mm x `--sss-scale` |
| `s(A)` | Burley's approximation (sets the shape from the albedo) | fixed formula |

Because the texture color at the entry point is used, the picture (eyebrows, skin patterns) stays sharp, and **only the light is blurred**.

### 2.3 Sampling

The cumulative distribution function (CDF) of the radius has a closed form.

```
CDF(r) = 1 − ¼·e^(−r/d) − ¾·e^(−r/(3d))
```

The inverse is solved with Newton's method kept inside a bracket (at most 24 iterations in the code), giving the radius `r` from a uniform random number `u`.
The mean radius is `2.5·d` (used in the verification).

### 2.4 Upper limit of the radius

The radius `rMax` of the probe sphere is the largest of the radii at "CDF = 0.999" over all channels.
Samples with a larger radius are discarded (about 0.1% of the energy is lost).

---

## 3. Material (CoatedDiffuse)

When `setSubsurface(profile)` is called on `CoatedDiffuse`, the SSS event is used instead of the diffuse lobe.
The surface coat (GGX specular, weight `w`) remains a local reflection.

| Lobe | Probability of being chosen | Content |
|---|---|---|
| Specular (coat) | `pS = clamp( w·F(cosθo), 0.05, 0.95 )` | as usual. Sample a direction and go on |
| SSS event | `1 − pS` | "entered the inside". No direction is returned. The weight is `(1 − w·F) / (1 − pS)` |

- If `pS` is not clamped, the weight of the SSS event is exactly 1.
- `eval` and `pdf` return only the specular lobe (there is no diffuse lobe).
- A new flag `BxDFFlags::Subsurface` is attached to the BSDF sample result and is tested with `isSubsurface()`.
  Here `f` is "just a weight", and neither the cosine nor the pdf is multiplied in.

### Sharing the material

The skin materials (face, neck, body) share **the same profile pointer**. The candidates for the exit point are
only the surfaces whose pointer matches. Clothes and hair are not candidates.

---

## 4. Sampling the exit point (probe rays)

`SampleSubsurfaceExit` (`Render/SubsurfaceProbe.hpp`) decides the exit point `pi` from the entry point `po`.
It is the same method as the BSSRDF sampling of PBRT v3.

### 4.1 Procedure

1. **Channel selection**: choose one of all the channels uniformly.
2. **Radius and angle**: sample the radius `r` from the profile of the chosen channel, and the angle `φ` uniformly.
   If `r ≥ rMax`, it fails (the path ends).
3. **Axis selection**: choose the direction of the probe among the local axes of the entry point (tangent s, bitangent t, normal n).

   | Axis | Probability |
   |---|---|
   | normal n | 1/2 |
   | tangent s | 1/4 |
   | bitangent t | 1/4 |

4. **Probe ray**: from the position offset by `(r cosφ, r sinφ)` in the plane perpendicular to the chosen axis, shoot a segment
   that pierces a sphere of radius `rMax`. Its length is `2·sqrt(rMax² − r²)`.
   Even on curved or creased surfaces, it can find "the surface about `r` away from the entry point".
5. **Collecting hits**: follow the hits on the segment from the front (at most 32).
   Only **surfaces with the same profile** are added as candidates (at most 8). Other materials are passed through.
6. **Choosing a candidate**: choose one of the candidates uniformly. If none is found, it fails.

### 4.2 Probability density (pdf)

Take into account that the chosen exit point `pi` could have been found from "any axis and any channel".

```
pdf = (1 / number of candidates) · Σ_axis [ axis probability · |axis component of n_pi| · (1/number of channels) Σ_c P(rProj_axis; d_c) ]
```

- `rProj_axis`: the length of `po − pi` projected onto the plane perpendicular to that axis.
- `|axis component of n_pi|`: a correction for the exit surface being oblique to the probe axis (the area Jacobian).

### 4.3 Weight

```
weight_c = A_c · P(|po − pi|; d_c) / pdf
```

The actual 3D distance `|po − pi|` is used. This weight is multiplied into the path throughput `beta`.

### 4.4 Exit vertex

- The normal is aligned to the **outward** direction of the surface (`scene.hit()` flips the normal toward the direction of the probe ray, so it is flipped back).
- The material is replaced by `SubsurfaceExit`.

---

## 5. Lobe at the exit point

The directional distribution `Sw` at the exit point is a cosine lobe including the Fresnel transmission.

```
Sw(wi) = ( 1 − w·F(cosθi) ) / ( c·π )
c      = 1 − 2·w·∫₀¹ F(μ) μ dμ       (numerical integration with 128 points)
```

`c` is chosen so that `∫ Sw · cosθ dω = 1` over the hemisphere. Since it does not depend on the outgoing direction `wo`,
the path can continue as if it had "arrived" at the exit point.

| Operation | Content |
|---|---|
| `eval` | `Sw(wi)` (only when `wi` is in the outward hemisphere) |
| `pdf` | `cosθi / π` |
| `sample` | sample a direction from the cosine distribution. Weight `f·cos/pdf = (1 − w·F)/c` |

---

## 6. Flow of the integrator

The processing at each vertex of `Integrator::Li` (**at most 2 stages**).

```
Stage 0 (entry point)
 ├ light sample (NEE) …… only the specular lobe contributes
 └ lobe selection
    ├ specular lobe ……………… ordinary BSDF sample → next direction
    └ SSS event ………… SampleSubsurfaceExit
                             beta ×= (event weight) x (profile weight)
                             rec = exit vertex, lastRec = exit vertex (for MIS)
Stage 1 (exit point)
 ├ light sample (NEE) …… the Sw lobe of SubsurfaceExit
 └ next direction from the cosine distribution → continue the path
```

- If no exit point is found, the path ends.
- Since `lastRec` is the exit vertex, the MIS weight when a light is hit later is computed relative to the exit point.
- The light sample at the entry point is only for the specular lobe. All diffuse reflection is handled by the light sample at the exit point.

---

## 7. Spectral mode

- The number of channels changes from 3 (RGB) to 4 (wavelengths) (`kSssChannels`).
- The albedo is `RgbColor::reflectance(lambda)`, and the scattering distance uses the 3 RGB values smoothly converted to spectra.
- For the radius sample, one wavelength chosen uniformly from the 4 wavelengths is used, and the pdf is computed as the mean over the 4 wavelengths.

---

## 8. Usage and parameters

```bat
GoLD_rayt.exe --scene pmx --pmx-index 4 --pose-index 0 --expr shy --pmx-yaw -25 --pmx-studio ^
              --pmx-zoom 4.0 --pmx-focus 0.79 --width 420 --spp 128 --sss --denoise --out face.png
```

| Option | Content |
|---|---|
| `--sss` | enables SSS on the skin materials of a PMX model (names containing 体, 肌, 颜 or 首) |
| `--sss-scale S` | multiplies the scattering distance by `S` (specifying it also enables SSS) |

| `--sss-scale` | Look |
|---|---|
| 1 | close to real size. Reddish light wraps into the shadow side and the cheeks become soft (subtle) |
| 3 | for an anime look. Cheeks and neck glow clearly orange |

- The effect is visible in close-ups of the face. In a full-body composition, 1 px is about 3 mm, so it is hardly visible.
- With `--pmx-matte`, the skin preset is disabled and SSS has no effect either.
- To use it from code:

```cpp
auto profile = std::make_shared<SubsurfaceProfile>(Vector3(3.67, 1.37, 0.68) * 0.001, /*eta=*/1.4);
auto skin = std::make_shared<CoatedDiffuse>(baseColor, /*roughness=*/0.5, /*coat weight=*/0.4, /*metallic=*/0, texture);
skin->setSubsurface(profile);   // share the same profile among all materials of face, neck and body
```

---

## 9. Verification

Run with `GoLD_rayt.exe --scene sss` (passes in both RGB and spectral modes). On failure the exit code is 2.

| Check | Content | Result |
|---|---|---|
| Normalization of the profile | agreement of `∫ P · 2πr dr` and the CDF | agrees |
| Inverse CDF | `CDF(sampleRadius(u)) = u` | agrees (`--scene sss` within 1e-9 at 5 points. A Python replica of the same algorithm was run on 100,000 points: maximum error 3e-16, at most 7 iterations) |
| Mean radius | comparison with the theoretical value `2.5·d` (within 4σ) | agrees |
| White furnace | a flat skin lit by a uniform sky. The expected value is the reflectance 0.8 | 0.797 |
| Knife edge | the radial profile across a shadow edge, compared with the numerically integrated theory (6 points) | all points within 4σ |

The knife-edge check uses a point light (height 100 m) and a black occluding plate to light only one side of the edge uniformly.
The theoretical value at distance `x` from the edge is computed from the fraction of the profile mass on the far side of the straight line,
`D(t) = ∫ (e^(−s/d) + 3e^(−s/(3d))) / (8π) dφ` (`s = t/cosφ`).

During the verification, a bug in the inverse-CDF solver (it fell back to bisection after convergence) was found and fixed.

### 9.1 Comparison with a random walk (error of the profile shape)

How well the engine's Burley profile (`A · P(r; d)`, `d = ℓ / s(A)`) reproduces the spread of light in real scattering was measured against a random-walk reference (2026-11-08). The script is `tools/studies/sss_random_walk.py` (numpy; independent of the engine), and the raw table is [data/sss_random_walk.txt](data/sss_random_walk.txt).

**Conditions**: a semi-infinite homogeneous medium (extinction coefficient 1; lengths in mean free paths), single-scattering albedo α, Henyey-Greenstein anisotropy g, an **index-matched boundary** (the same assumption as the engine, which handles the Fresnel factors outside the profile), and a thin beam hitting the origin perpendicularly. Each photon is followed until it exits or is absorbed, and the radius of the exit point is counted (1 million photons per case).

**How they are compared**: take the total reflectance obtained by the walk as `A`, and build the profile with the engine's mapping (`ℓ` is 1 for isotropic, and the reduced mean free path `1/(1−gα)` for anisotropic). The measures are the maximum difference of the radial cumulative distribution (Kolmogorov distance, KS; 0 = identical), the KS when only `d` is optimized, the mean radius, and the ratio "walk / model" at several radii.

| α | g | A (walk) | KS (engine mapping) | best `d` / mapping | KS (best `d`) | mean radius walk / model (2.5d) |
|---|---|---|---|---|---|---|
| 0.5 | 0 | 0.115 | 0.021 | 1.05 | 0.008 | 0.68 / 0.63 |
| 0.8 | 0 | 0.285 | 0.013 | 1.02 | 0.007 | 1.03 / 0.99 |
| 0.9 | 0 | 0.415 | **0.010** | 0.97 | 0.006 | 1.33 / 1.36 |
| 0.95 | 0 | 0.535 | 0.014 | 0.95 | 0.007 | 1.66 / 1.73 |
| 0.99 | 0 | 0.752 | 0.025 | 1.06 | 0.014 | 2.57 / 2.28 |
| 0.999 | 0 | 0.913 | 0.081 | 1.30 | 0.036 | 4.11 / 2.64 |
| 0.9999 | 0 | 0.972 | 0.111 | 1.46 | 0.054 | 5.79 / 2.74 |
| 0.99 | 0.5 | 0.665 | 0.072 | 1.20 | 0.026 | 4.71 / 4.12 |
| 0.99 | 0.8 | 0.522 | 0.122 | 1.32 | 0.050 | 9.61 / 8.13 |
| 0.999 | 0.8 | 0.813 | 0.112 | 1.41 | 0.033 | 16.42 / 12.01 |

![Radial profile. Dots: random walk, red: the engine's mapping, blue dashed: best value of `d`](img/sss_profile_random_walk.png)

**What was found**
- **In the range of skin diffuse colors (isotropic, A = 0.3 to 0.75), they agree well.** KS is 0.010 to 0.025, and 1 to 2.5% of the exiting energy has its radius off. The mapping of `d` is also close to the best `d` (within ±6%).
- **When A is large (0.9 or more, α ≥ 0.999), the tail does not match.** The mean radius of the walk is 1.6 to 2.1 times that of the model. Even with `d` multiplied by 1.3 to 1.5, KS only falls to 0.04 to 0.05, so a difference in shape remains (Burley's `s(A)` is an approximation fitted around A of about 0.8).
- **With anisotropy (g = 0.5, 0.8), the region near the entry point does not match.** The similarity rule that takes the reduced mean free path as `ℓ` does not work near the entry point (radius up to 0.3 times ℓ). The walk's value is 0.5 to 0.7 times the model's, so the model overestimates the neighborhood. Far away (3 times ℓ or more) it is the other way round, 1.1 to 1.4 times, and the walk has the longer tail. The g of skin is around 0.8.
- The ratio in the tail (radius 10ℓ) is noisy, because few photons are counted (the 4.4 times at α = 0.5 is mostly the fluctuation of a few tens of photons).

**What this comparison can say**
- What it can say: the error of the profile shape for a **flat semi-infinite medium**. It is 1 to 2.5% for ordinary skin colors, and 7 to 12% when the albedo is high or the anisotropy is strong.

**Reproduction**: `python -I tools/studies/sss_random_walk.py --walkers 1000000 --out docs/data/sss_random_walk.txt --fig out/sss/profile_rw.png` (about 2 minutes).

### 9.2 Sphere and thin disc (BSSRDF and volumetric random walk)

Section 9.1 was a flat semi-infinite medium. To see how much the result shifts once curvature and thickness enter, the same medium was rendered in two ways inside the engine and compared (2026-11-08. `--scene sssrw`, `src/Cli/Tools/SssRandomWalk.cpp`, raw output in [data/sss_bssrdf_vs_volume.txt](data/sss_bssrdf_vs_volume.txt)).

**Conditions**
- **Volume**: the inside of a sphere or disc is filled with a scattering medium (gray, albedo α, anisotropy g). The boundary is a **null interface with index 1** (no reflection or refraction), lit by a point light with in-medium NEE.
- **BSSRDF**: the same shape with the SSS of `CoatedDiffuse` (no coat, weight 0). The scattering distance ℓ is the reduced mean free path (2 mm), and the diffuse color `A` is the total reflectance of the semi-infinite medium of 9.1 (0.752 for α = 0.99, 0.415 for 0.9, 0.522 for α = 0.99, g = 0.8).
- The only differences are what the profile model cannot represent (curvature, thickness, the shape of the profile).
- Lit by a point light with a narrow cone, 256 spp. The values are the film luminance (green), where 1 is the intensity of a white Lambertian surface facing the light.

**Sphere** (radius `R`, light 60° to the side. The terminator is at about `x/R ≈ −0.5` on the diameter)

| Condition | Whole-image total (BSSRDF / volume) | Shadow side (x/R = −0.6) | Center (0) | Light side (+0.6) | Edge of the light side (+0.96) |
|---|---|---|---|---|---|
| α = 0.99, R = 2ℓ | 0.88 | 0.64 | 0.84 | 0.96 | 1.40 |
| α = 0.99, R = 5ℓ | 0.91 | 0.54 | 0.83 | 0.98 | 1.17 |
| α = 0.99, R = 10ℓ | 0.91 | 0.63 | 0.87 | 0.99 | 1.09 |
| α = 0.9, R = 5ℓ | 0.78 | 0.34 | 0.76 | 0.93 | 1.00 |
| α = 0.99, g = 0.8, R = 5ℓ | 0.88 | 0.50 | 0.76 | 0.95 | 1.52 |

- The light side (where the light falls) agrees within 7% (except the edge).
- **The shadow side and the center are darker with the BSSRDF** (0.34 to 0.87 times). The light that travels through the inside of the object and wraps around beyond the terminator is missing. Even with a radius of 10ℓ it hardly changes (0.63), so it is not a problem of small spheres only.
- Over the whole image, it is 9 to 12% darker at α = 0.99 and 22% darker at α = 0.9.

**Disc** (radius 15ℓ, thickness `t`, mean over the central region. Reflection = viewed from the same side as the light, transmission = lit from the opposite side. The light is on the axis, in line with the camera)

| t / ℓ | Reflection: volume | Reflection: BSSRDF | Ratio | Transmission: volume | Transmission: BSSRDF | Ratio |
|---|---|---|---|---|---|---|
| 0.25 | 0.068 | 0.725 | **10.7** | 0.067 | 0.623 | **9.3** |
| 0.5 | 0.136 | 0.732 | 5.4 | 0.133 | 0.553 | 4.2 |
| 1 | 0.259 | 0.725 | 2.8 | 0.237 | 0.426 | 1.8 |
| 2 | 0.439 | 0.721 | 1.6 | 0.339 | 0.265 | 0.78 |
| 4 | 0.618 | 0.732 | 1.18 | 0.302 | 0.114 | **0.38** |
| 8 | 0.709 | 0.745 | 1.05 | 0.122 | 0.022 | **0.18** |

(α = 0.99, g = 0)

- **The reflection of the BSSRDF is constant regardless of thickness** (the semi-infinite value, about 0.73). In reality, when the disc is thin, light escapes from the back and the reflection becomes smaller. Therefore **for t ≤ 1ℓ the BSSRDF is 3 to 11 times brighter**, and **for t ≥ 4 to 8ℓ it is within 5 to 18%**.
- **Transmission** is **4 to 9 times too high** when thin (t ≤ 0.5ℓ), agrees around 1.5ℓ, and is **2.6 to 5.5 times too low when thick (t ≥ 4ℓ)**.
- With a low albedo (α = 0.9, A = 0.415), the tendency is the same: reflection 7.0, 3.6, 2.0, 1.3, 1.1, 1.1 times, transmission 5.6, 2.4, 1.0, 0.41, 0.19, 0.09 times.
- Anisotropy (g = 0.8) was measured with the light shifted 45° from the axis (on the axis, the camera looks straight into the forward-scattering peak, which makes the comparison meaningless). Reflection is 9.0, 3.8, 2.0, 1.3, 1.04, 0.97 times, the same tendency as isotropic. Transmission is 1.3, 0.72, 0.46, 0.36, 0.28, 0.17 times.
- With the light at 45° (α = 0.99, g = 0), reflection is 7.9, 4.1, 2.3, 1.5, 1.14, 1.02 times and transmission is 6.9, 3.3, 1.6, 0.81, 0.44, 0.21 times, the same tendency as normal incidence (the thinner, the larger the overestimate; the values are somewhat smaller).

**Meaning for skin**
- The default scattering distances are red 3.7 mm, green 1.4 mm and blue 0.7 mm. For an ear 3 mm thick, `t/ℓ` becomes 0.8 for red, 2.2 for green and 4.4 for blue.
- **Reflection is overestimated for red** (about 3 times), **about 1.6 times for green, and almost right for blue**. **For transmission, red is overestimated about 2 times, green is 0.8 times, and blue is underestimated at 0.4 times.** The color of thin parts may shift from channel to channel.

**What this can say**
- What it can say: the error of the BSSRDF for a **gray medium** with an index-matched boundary. On a sphere the shadow side is dark, and on a thin disc the thinner it is the larger the overestimate, and the thicker it is the more the transmission is underestimated.

**Reproduction**: `GoLD_rayt.exe --scene sssrw --spp 256 --sssrw-albedo 0.99 --sssrw-g 0 --sssrw-A 0.752 --sssrw-sphere 5` (about 90 seconds). `--sssrw-lamp-deg 45` shifts the light of the disc 45° from the axis.

---

## 10. Computational cost

Ratios to no SSS for a face close-up (width 420, 128 spp, same conditions).

| Setting | Time | Ratio to no SSS |
|---|---|---|
| no SSS | 77 s | 1.0x |
| `--sss` (scale 1) | 130 s | about 1.7x |
| `--sss-scale 3` | 185 s | about 2.4x |

- The larger the scattering distance, the wider the area the probe examines, so it is slower.
- These are numbers for a composition in which skin fills most of the frame. For a full body the increase is smaller.
- The numbers are a guide, since on this PC the speed varies by up to about 1.4 times from run to run.

---

## 11. File list

| File | Content |
|---|---|
| [`include/Materials/Subsurface.hpp`](../GoLD_rayt/include/Materials/Subsurface.hpp) | formulas of the profile, `SubsurfaceProfile`, `SubsurfaceExit` |
| [`include/Render/SubsurfaceProbe.hpp`](../GoLD_rayt/include/Render/SubsurfaceProbe.hpp) | sampling of the exit point (`SampleSubsurfaceExit`) |
| [`include/Materials/CoatedDiffuse.hpp`](../GoLD_rayt/include/Materials/CoatedDiffuse.hpp) | `setSubsurface`, lobe selection, computation of the parameters |
| [`include/Materials/Material.hpp`](../GoLD_rayt/include/Materials/Material.hpp) | `BxDFFlags::Subsurface`, hooks such as `subsurface()` |
| [`include/Render/Integrator.hpp`](../GoLD_rayt/include/Render/Integrator.hpp) | the 2-stage loop (NEE and BSDF sampling in `Li`) |
| [`include/Materials/SubsurfaceWalk.hpp`](../GoLD_rayt/include/Materials/SubsurfaceWalk.hpp) | conversion of the coefficients and the medium of the walk mode (section 12) |
| [`src/Cli/Options.cpp`](../GoLD_rayt/src/Cli/Options.cpp) | `--sss`, `--sss-scale` |
| [`src/Cli/PmxSupport.cpp`](../GoLD_rayt/src/Cli/PmxSupport.cpp) | assignment of the profile to the skin of a PMX model |
| [`src/Cli/Checks/SssChecks.cpp`](../GoLD_rayt/src/Cli/Checks/SssChecks.cpp) | `--scene sss` |
| [`src/Cli/Tools/SssRandomWalk.cpp`](../GoLD_rayt/src/Cli/Tools/SssRandomWalk.cpp) | `--scene sssrw` (section 9.2) |

---

## 12. Walking inside the object (`randomWalk`, 2026-11-10)

The BSSRDF "jumps" from the point of entry to another point on the surface, so it does not know the thickness of the object (9.2). This mode keeps the same skin material but **lets the light actually walk inside the object**. The controls are the same as for the BSSRDF (color and scattering distance), and the layers on the skin (pigment, flaking skin (scales), lotion, foundation) can be used as they are.

```jsonc
"skin": { "type": "skin", "albedo": [0.80, 0.52, 0.42], "scatterDistanceMm": [3.67, 1.37, 0.68],
          "randomWalk": true }            // or { "g": 0.8 }
```

It needs a closed shape (with outward normals) and works only with the path-tracing integrator (LVC-BPT and photon maps do not support it). When run, a `[Skin] ...` line prints the internal coefficients.

### 12.1 Treatment of the boundary

It is the same "diffusing interface" as in film-production renderers (Chiang et al. 2016, Wrenninge et al. 2017; the exit of Cycles is the same).

| Case | Treatment |
|---|---|
| Entering from outside | the light that passed the coat (1 − F) enters with an inward cosine distribution. Refraction is not computed |
| Inside | the walk of the existing medium (`PathIntegrator::advanceInMedium`. For colored coefficients, one channel chooses the distance and the others are reweighted) |
| Reaching the surface from inside | it always exits (no internal reflection). The exit direction and weight are those of the BSSRDF exit lobe (`SubsurfaceExit`). A light is aimed at directly there |
| Layers (pigment, flaking skin, lotion, foundation) | the same coefficients as for the BSSRDF are applied at the entrance and the exit |

Implementation: `Materials/SubsurfaceWalk.hpp` (conversion of the coefficients and the medium), `CoatedDiffuse::setSubsurfaceWalk()`. For an intersection reached from inside, `perturbShading` flips the normal outward (so that light selection and shadow rays face outward).

### 12.2 From color and scattering distance to the internal coefficients

    A = (1 − s)(1 − 0.139 s) / (1 + 1.17 s)        reflectance of a semi-infinite medium under diffuse light (van de Hulst's formula)
    s² = σa / (σa + σs (1 − g))                     solve this for s from A (a quadratic equation)
    1 / ℓ = σa + σs (1 − g)                         the scattering distance ℓ is the mean free path converted to isotropic

- It is set so that thick regions return the same color (A) as the BSSRDF. It is the same formula as Cycles.
- Changing `g` does not change s or ℓ. The color and spread of thick regions stay the same; only the regions decided by a small number of scatterings (thin regions, near the entrance) change.
- In the spectral build the conversion is done per wavelength (the conversion is nonlinear, so converting in RGB and then going to spectra shifts the color).

### 12.3 Checks (4 to 8 of `--scene sss`. PASS in both RGB and spectral)

| Check | Result |
|---|---|
| Error of the conversion formula (exact value: Chandrasekhar's H function. A = 0.1 to 0.95) | at most 0.0018 (the limit 0.004 is checked) |
| Boundary lobe: the fraction entering is 1 − w F, and the mean exit weight is 1 (`sample()` and `eval()`) | 0.9829 vs 0.9829, 1.000 |
| White-furnace sphere against an independent walk written inside the check (radius 2 to 10ℓ, color 0.4 and 0.8, g = 0 and 0.8, a colored skin) | within 0.4% in all 5 cases (the spectral colored skin is for reference only) |
| A thick sphere (radius 1000ℓ) returns the diffuse color (g = 0 and 0.8) | 0.8006, 0.8002 vs 0.8 |

- **A sphere of radius 40ℓ is still 3% brighter** (0.826). This is because light that spread sideways leaves early through the convex surface; it is not a defect (it agrees with the independent walk).
- White furnace of a thin sphere (radius 0.5ℓ): walk 0.99, BSSRDF 0.31. The BSSRDF loses the portion for which no exit is found.

### 12.4 Backlit finger (the composition of `scenes/skin_fingers.json`. Center of the finger, linear value of red)

| BSSRDF | This mode (diffusing interface) | Refracting interface + medium (`medium` of `dielectric`) |
|---|---|---|
| 0.35 | 0.46 | 0.82 |

**The treatment of the interface changes the transmitted light by about 1.8 times.** A refracting interface gathers the light toward the normal and keeps the light inside by internal reflection. These are single-image values at 768 spp.

### 12.5 References

See [REFERENCES.md](REFERENCES.md), group "Participating media and subsurface scattering".
