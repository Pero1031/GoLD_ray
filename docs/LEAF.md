# Leaves

A leaf is thin: light that enters it is scattered inside and leaves from both faces. Its colour comes from what absorbs inside
(chlorophyll, carotenoids, anthocyanins, brown pigments), and in autumn these change. Backlit, a leaf glows with the light it
transmits; this is what makes autumn foliage look the way it does.

The `leaf` material follows Kruk, Noizet and Lucas (EGSR 2026): the diffuse reflection and transmission come from **PROSPECT**, a
model from remote sensing that gives a leaf's reflectance R(λ) and transmittance T(λ) from its contents, and a GGX lobe stands for
the waxy cuticle on each face.

```json
"maple": { "type": "leaf", "chlorophyll": 5, "carotenoids": 8, "anthocyanins": 15, "N": 1.6, "roughness": 0.55 }
```

The coefficient table of PROSPECT and the ANGERS leaf database are not part of the repository: `python tools/data/leaf_optics.py fetch`
downloads them to `GoLD_rayt/assets/Leaf/` (git-ignored).

## How it works

Code: `Materials/Prospect` (the model), `Materials/Leaf` (the BSDF).

- **PROSPECT** (Jacquemoud and Baret 1990; PROSPECT-D, Féret et al. 2017) treats a leaf as a pile of N plates with flat faces. One
  plate lets through the fraction τ(k) = (1 − k) e^−k + k² E1(k) of the isotropic light that enters it, where the absorption
  k(λ) = Σ K_i(λ) C_i / N sums each constituent's specific absorption times its content. The faces reflect by Fresnel's law
  averaged over the incoming directions; Stokes' formulas stack the plates for any real N. It is evaluated per wavelength of the
  path (the spectral build) or once, as linear sRGB under D65 (the RGB build).
- **The BSDF.** On the side the light is on, a GGX lobe reflects F D G / (4 μo μi). What it does not reflect enters the leaf:
  - reflection Rd (1 − E(μo)) (1 − E(μi)) / (π (1 − E_avg)²),
  - transmission T (1 − E(μo)) (1 − E'(μi)) / (π (1 − E_avg)(1 − E'_avg)), where E' belongs to the other face.

  E(μ) is the directional albedo of the face's GGX lobe and E_avg its cosine average. Rd = R − (1 − tav(90°, n)) is PROSPECT's
  reflectance without the reflection of the first face, which the GGX lobe now provides.
- **One change from the paper.** The paper adds its GGX lobe to the whole of PROSPECT's reflectance. That reflectance already holds
  the first face's reflection, so the paper counts that light twice. That is about 0.09 of reflectance, more than the whole diffuse
  reflection of a green leaf in blue and red.
- Under light that is uniform over all directions, the lobes reflect exactly Rd and transmit exactly T. The pair is reciprocal, and
  no direction returns more light than it receives.
- The upper and the lower face have their own roughness (default 0.55 and 0.75). The diffuse lobes are the same on both faces, as in
  PROSPECT.

## Checks (`--scene leaf`, both builds, ~30 s)

| | check | result |
|---|---|---|
| a | E1, tav and R, T of a green, a red and a brown leaf (450–700 nm, α 90° and 40°) against an independent implementation in numpy (`tools/lib/prospect.py`, scipy's E1); no absorber gives R + T = 1 | ≤ 3e-16 |
| b | under uniform light: what the diffuse lobes reflect = Rd, what they transmit = T, per wavelength | ≤ 0.004 % |
| c | reciprocity of the reflection; a transmission from above = the reverse from below | 4e-15 |
| d | the directional albedo through `sample()` = through `eval()`, both faces, θo 0 / 50 / 80°; the pdf integrates to ≤ 1 | ≤ 2.2 σ |
| e | energy: a leaf without absorbers returns 0.95–0.97 (the GGX lobe and the cuticle factors lose the rest), never above 1 | holds |

## How good is PROSPECT at the colour of a leaf?

`tools/studies/leaf_angers.py` takes the 308 leaves of the ANGERS database (with 32 anthocyanin-rich sycamore maple leaves that
Féret et al. added). For each leaf it fits N in the near infrared and puts the measured pigments into PROSPECT. It then compares the
result with the measured R and T in the visible.

| leaves | n | RMSE of R | RMSE of T | ΔE00 of R (median / 90 %) | ΔE00 of T (median / 90 %) |
|---|---|---|---|---|---|
| all | 308 | 0.014 | 0.018 | 3.7 / 6.8 | 4.2 / 7.1 |
| anthocyanins > 3 µg/cm² | 35 | 0.008 | 0.017 | 2.6 / 4.1 | 6.5 / 20.2 |

The model gives the colour of a real leaf from its pigments to a few CIEDE2000 units. The worst cases are the transmitted colour of
red leaves, where the anthocyanin contents are estimates, not measurements.

## Limits

- PROSPECT is symmetric: real leaves are lighter underneath (the spongy mesophyll), and the cuticle may differ.
- One set of contents per material for now; contents that vary over the leaf (veins, margins, age) are the next step.
- The diffuse lobes are Lambertian. A leaf's transmission is slightly forward-peaked, and real veins and hairs scatter
  anisotropically.
