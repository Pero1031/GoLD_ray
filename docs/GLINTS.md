# Glints: a GGX distribution made of a finite number of facets

**Status:** stage 1 of 5. The distribution itself and its checks are in (`Microfacet/GlintNDF`, `--scene glints`). It is not yet a
material: nothing renders with it. Next: a material seen directly by the camera (the pixel's footprint from the camera ray), importance
sampling, later bounces, example scenes.

A smooth microfacet model says how much of the surface faces each direction, as if every pixel held infinitely many facets. Metallic
paint, lamé and pearl pigments hold few enough for single facets to flash. The glinty distribution keeps the smooth model's average
and puts a finite number of facets under it.

## 1. The model

After P. Kemppinen, L. Paulin, T. Thonat, J.-M. Thiery, J. Lehtinen and T. Boubekeur, "Evaluating and Sampling Glinty NDFs in Constant
Time", *ACM Transactions on Graphics* 44(6), 2025 (SIGGRAPH Asia). I read the authors' preprint in full and wrote the code from its
text; their shader was not used.

- **A facet is a point in four dimensions**: two for its place on the surface (texture space), two for its normal. The normals live in
  a disc U inscribed in the unit square, where the map T(ω) = L(Mω / |Mω|) puts a micro-normal: M = diag(1/α_x, 1/α_y, 1) turns the
  GGX distribution into the normals of a hemisphere, and L(x, y, z) = (x, y) / √(1 + z) is Lambert's equal-area projection.
- **One facet per cell of a 4D grid**, at a random place in the cell (and no facet where that place falls outside the disc). The
  cell's index seeds the random numbers, so nothing is stored or precomputed.
- **A facet answers** the half vector with a Gaussian lobe N_a in U and the pixel with the pixel's filter f (a Gaussian on the
  surface):

  D\*(x, ω) = D(ω) × Σ_i N_a(T(ω) − μ_i) f(x_i − x)

  in units where the grid holds one facet per unit of (area on the surface × area in U). Averaged over the surface the sum is 1, so
  D\* averages to the smooth D.
- **The grid follows the pixel.** Level l has spatial cells of side 2^l and 2^l × 2^l angular cells: the same number of facets at
  every level, traded between place and direction. The level is chosen so that a spatial cell is 2 to 4 standard deviations of the
  filter wide; a few spatial cells then cover the pixel.
- **Near facets one by one, far facets by their expectation.** Within those spatial cells, the facets of the 2 × 2 angular cells
  nearest to T(ω) are summed. The others are replaced by what they give on average: the lobe's mass over the unit square minus that
  block (error functions). The cost per evaluation is therefore bounded, whatever the size of the pixel.
- **Between two levels** each facet is kept or dropped by a roulette with the level's weight instead of being blended, so a glint keeps
  its brightness.

Parameters (`GlintParams`): the roughness α_x, α_y of the whole distribution, the micro-roughness (the width of one facet's lobe), the
density (facets per unit area of texture space), a seed for another arrangement.

## 2. Where the code differs from the paper

- **Normalization.** The printed formulas do not fix the constants unambiguously (the Gaussian is written with a 1D normalization,
  the result is scaled by 1/π). Here the constants follow from the requirement that D\* average to D (the paper's Eq. 2).
- **The lobe's width.** The paper widens a lobe by √D as a stand-in for the map's Jacobian. Here it is the true Jacobian of T,
  det M / (8 |Mω|³), so the micro-roughness is the lobe's standard deviation in radians of the half vector (their geometric mean
  along the two axes; the map stretches a lobe, as in the paper).
- **The spatial neighbourhood** covers 4 standard deviations of the filter (up to 5 × 5 cells) instead of a fixed small block: a
  renderer that is compared with references should not lose the filter's tails. This choice, and the level being 2 to 4 standard
  deviations per cell, are mine; the paper does not state its own.
- **The roulette is exact by default** (keep with probability w). The paper softens it (ε = 0.1) for smooth animation;
  `rouletteSoftness` restores that.
- **The error function** is the standard library's, not a rational approximation.
- Not yet here: anisotropic facets, coloured facets, surfaces without texture coordinates (triplanar), importance sampling.

## 3. What the model does not do (the paper's own statements)

- **It does not conserve energy exactly.** Lobes of facets near the disc's rim reach outside it, and the expectation is taken over the
  square that holds the disc. The loss grows when the micro-roughness is not small against the roughness (section 4).
- **Roughness and micro-roughness are not independent.** A micro-roughness near the roughness itself blurs the highlight: physically
  sensible, but the two must be chosen together.
- **Path tracing is left as future work** in the paper (it would need the pixel's footprint carried along the path). That part will be
  this renderer's own.

## 4. Checks (`--scene glints`, both builds, about 3 minutes)

| Check | Result |
|---|---|
| a. The map T there and back; its Jacobian against finite differences; D against the engine's GGX (anisotropic α) | 4 × 10⁻¹⁵; 5 × 10⁻⁸; 9 × 10⁻¹⁶ |
| b. Levels 0 and 1, where every facet is summed, against a plain sum over all facets (which visits 6 standard deviations of the filter) | sum of differences / sum below 7 × 10⁻⁶ |
| c. The average over the surface is the smooth GGX: 8 cases, levels −2 to 9.2, footprints stretched up to 8 : 1, anisotropic roughness | within the standard error (0.002 to 0.3 %) in every case |
| d. Level 3: the far facets replaced by their expectation against the plain sum of all 64 angular cells | averages 1.006 and 1.000, each ± 0.002 |
| e. Energy against the smooth distribution (half vectors drawn from D cos θ, averaged over the surface): α 0.1 / 0.3 / 0.6, micro-roughness 0.01 / 0.05, levels 1.5 / 6.5 | within the error in 10 cases; −2.1 % and −0.8 % for α 0.1 with micro-roughness 0.05 (a lobe half as wide as the distribution) |
| f. Many facets per pixel: the spread of D\* / D over pixels | 0.04 % at level 12, 10⁻⁶ at level 16; exactly 1 beyond level 20 |
| g. (information) the paper's softened roulette, level 2.5 | average off by −0.2 ± 0.3 % |

The averages are unbiased estimates, so they are judged by their standard error (4 standard errors plus 0.3 %); the energy loss is an
intended approximation and is judged by a bound (2 %).
