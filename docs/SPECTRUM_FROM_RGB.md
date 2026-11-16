# From RGB to spectrum: the unavoidable error, and the number of moments

An RGB renderer, and a conversion from RGB to a spectrum, know a surface only as "RGB under D65". This study examines how far the color under another light (tungsten, candle, LED) can be predicted from this information alone (2026-11-08). Using the reference code of the moment representation of Peters et al. (2019), it also measures "how many moments must be kept for the color to be guaranteed when the light changes".

The full output of the calculation: [metamer_mismatch_results.txt](data/metamer_mismatch_results.txt). Script: `tools/studies/metamer_mismatch.py`.

## 1. What was measured

The set of bounded reflectances (0 ≤ R(λ) ≤ 1) that have the same D65 RGB is a convex set, and the true reflectance is one member of it (the metamer set). Under a different light, each member of the set has a different color (metamer mismatch; prior work includes Wyszecki and Stiles, Finlayson and Morovic, and Logvinenko et al.). The spread of those colors is the **error that no conversion can avoid**, because whichever spectrum is chosen, the truth may be another member.

- The colors at the ends of the spread are found by a linear program per direction (scipy's HiGHS; 150 directions).
- **Radius**: the smallest, over single estimates, of the worst-case ΔE00 of that estimate (an approximation using the vertices and the centroid as candidates; a value that bounds it from above).
- Kinds of constraint: (a) the D65 color (three numbers), (b) in addition, smoothness (curvature at most f times the largest curvature of that spectrum itself), (c) m cosine moments of the phase (the version of Peters et al. with the "mirrored signal"; with and without the wavelength warp).
- Spectra: 11 measured cloth spectra (lab data; 430 to 700 nm used, clipped to [0, 1], 5 nm moving average) and 4 Fresnel reflectances of gold (normal incidence and 60°, one reflection and three reflections). Lights: A, candle, LED.

## 2. Results (radius, mean ΔE00)

| Constraint | Cloth: A | Cloth: candle | Cloth: LED | Gold: A | Gold: candle | Gold: LED |
|---|---|---|---|---|---|---|
| D65 color only | 16.8 | 22.0 | 20.7 | 8.8 | 11.7 | 8.3 |
| + smooth (f = 1) | 8.9 | 10.8 | 9.1 | 6.1 | 8.3 | 4.1 |
| 3 moments (with warp) | 17.2 | 20.4 | 23.0 | 8.3 | 10.0 | 8.8 |
| 4 | 9.4 | 9.1 | 21.9 | 5.0 | 3.8 | 7.9 |
| 6 | 5.7 | 6.3 | 13.3 | 1.8 | 1.7 | 3.9 |
| 10 | 4.6 | 4.1 | 8.6 | 1.2 | 0.9 | 2.5 |
| 6 / 10 moments without warp | 23.4 / 4.6 | 18.9 / 3.6 | 32.1 / 13.6 | 5.4 / 2.4 | 5.1 / 1.6 | 12.4 / 4.7 |
| (Reference) error of the engine's current conversion | 0.6 | 0.7 | 0.6 | 1.5 | 1.8 | 0.9 |

## 3. How to read this

1. **The lower bound with no prior knowledge is too large to serve as a yardstick for comparing conversions.** With the D65 color alone it is 17 to 22 for cloth and 9 to 12 for gold. The actual error of the current conversion (0.6 to 1.8) is small because real reflectances are smooth (the conversion chooses the smoothest spectrum). Even with the smoothness constraint added, the radius stays at 6 to 11. **To decide how good a conversion is, a set of real reflectances (statistical prior knowledge) is needed.**
2. **Three moments with warp carry almost the same information as the D65 color** (17.2 and 16.8 for cloth). This agrees with the claim of the paper and serves as a check of this implementation (a measure over the phase, the warp). With three to five moments without warp there is almost no guarantee (32 to 48 for cloth).
3. **Adding moments narrows the guarantee when the light changes.** For gold (a gentle one-step rise), six moments bring the worst-case difference under A and candle to 2 or less. For cloth, even ten moments leave 4 to 5. A sharp light such as an LED needs more (8.6 for cloth, 2.5 for gold). The sharper the light, the more moments are needed.
4. **This is a worst-case guarantee, not the error of a maximum-entropy reconstruction (the actual MESE).**

## Reproduction

```bash
python -I tools/studies/metamer_mismatch.py --smooth 1,2,4 --moments 3,4,5,6,8,10 --dirs 150
python -I tools/studies/metamer_mismatch.py --moments 3,4,5,6,8,10 --dirs 150 --warp references/peters2019/ReferenceImplementation.py
```

`--warp` loads, by path, the reference code of Peters et al. (BSD 3-clause, `references/peters2019/`; not included in the repository).
