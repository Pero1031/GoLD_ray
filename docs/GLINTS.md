# Glints: a GGX distribution made of a finite number of facets

**Status:** stage 3 of 5. The distribution (`Microfacet/GlintNDF`) is an option of the rough metal (`"glints"` on a `conductor`,
section 3) for surfaces the camera sees directly, with its own sampling where a pixel holds few facets; `--scene glints` checks
them. Next: later bounces, more scenes.

![Three gold spheres of the same roughness: smooth, 300 000 facets per square metre, 10 million](img/glints_spheres.png)

(`scenes/glints_spheres.json`, spectral build, 256 spp, 88 s. The same gold and the same roughness 0.3 three times, under two small
lamps: the smooth distribution, 3 × 10⁵ facets per m² and 10⁷. The grain on the lower halves is the floor's light in single facets, and in good part still noise: at 64 spp the noise there
is 78 % of the mean on the middle sphere and 40 % on the right one, so about half of that at 256 spp.)

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
- **Importance sampling is a different procedure** (section 3): the paper's has no density that can be written down exactly, and
  an unbiased path tracer needs one.
- Not yet here: anisotropic facets, coloured facets.

## 3. In the renderer

- **Scene file.** A `conductor` takes `"glints": { "density", "microRoughness", "filter", "seed" }`: facets per square unit of the
  scene (metres), the width of one facet's lobe in radians (default 0.01), the pixel filter's standard deviation in pixels (default
  0.5), another arrangement. The roughness is the conductor's own. Only D changes: Fresnel, masking and the multiple-scattering
  compensation are those of the smooth metal. `"sampling": false` switches the facets' own sampling off (below).
- **The pixel's size** comes from the camera (`PathIntegrator::setPixelFootprints`, switched on by a scene with glints): once per
  pixel, the rays of the pixel's centre and of its two neighbours through the same point of the lens give the spacing of neighbouring
  rays at the lens and its growth with distance (any camera: pinhole, real lens, orthographic). A camera ray carries the two numbers,
  the hit turns them into a length (`SurfaceInteraction::pixelSpacing`), and the material makes the filter from it: a disc across the
  ray, projected along the ray onto the surface (an ellipse, long in the plane of incidence; the cosine is held above 0.02).
- **No texture coordinates.** The facets are laid out in world space on the three coordinate planes: three arrangements, each seen
  through the footprint projected onto its plane, mixed by a roulette with the weights |n_k| / Σ|n| (the paper's triplanar extension;
  a plane the surface is nearly perpendicular to is left out). Spheres, quads and meshes without UVs work alike. The density is per
  area of the coordinate plane, so a surface tilted against all three carries somewhat fewer facets per area than one that faces an
  axis (by |n_k| for each arrangement).
- **Only where the camera looks directly.** A ray that has bounced does not know a pixel's size, and the material answers with the
  smooth distribution: the facets' average, so the picture's energy is the same and only the glints are missing in reflections
  (stage 4). A camera ray keeps its pixel through the boundary of a fog and through a see-through layer.
- **The pixel is filtered twice**: by the facets' filter and by the samples spread over the pixel. A glint is therefore a little
  wider than the filter says. A smaller `filter` gives sharper glints and needs more samples.
- **Sampling the facets.** Under a large light the smooth distribution's samples rarely land inside the narrow lobes, and the
  picture is noisy. Where a pixel holds few facets, half of the samples are therefore drawn from the lobes: (1) a level, with the
  blend's probabilities; (2) an angular cell, uniformly; (3) a facet of that cell under the pixel, in proportion to its weight (the
  pixel filter at the facet); (4) a point from the facet's lobe, cut to a square window of half width min(4 σ, half a cell). The
  window makes the density exact and as cheap as an evaluation: a point can only come from the 2 × 2 angular cells nearest to it.
  An empty cell or a point outside the disc gives no sample, which the density accounts for. The far facets' expectation is left
  to the smooth sampling, with which this one is mixed half and half.
- **Only below level 1.5** (about three facets per pixel). Measured on `glints_spheres.json` at 64 spp, the lower half of the middle
  sphere, lit by the floor (noise from two renders with different seeds, the brightest 0.5 % of the differences clipped):

  | Facets per m² (level) | Noise without / with | Time without / with |
  |---|---|---|
  | 3 × 10³ (about −1.2) | 357 % / 114 % of the mean | 12 s / 18 s |
  | 3 × 10⁴ (about 0.5; the right sphere) | 187 % / 92 % | the same renders |
  | 3 × 10⁵ (about 2.1) | 78 % / 62 % | 19 s / 44 s |

  With few facets the variance falls about tenfold for half as much time again. With many, the gain does not pay for the density
  (per time, the plain sampling wins), so the facets are sampled only below level 1.5. The weakness is step 2: the cell is chosen
  without knowing how much its facets weigh under the pixel. Under a small lamp the mixed sampling is a little noisier than the
  plain one (8 % against 4 % on the upper half at 3 × 10³), because the lamp is found by light sampling either way.
- **No denoiser.** It takes single bright pixels for noise.
- Pictures without glints are unchanged: three scene files rendered before and after with the same seed are identical byte for byte.

## 4. What the model does not do (the paper's own statements)

- **It does not conserve energy exactly.** Lobes of facets near the disc's rim reach outside it, and the expectation is taken over the
  square that holds the disc. The loss grows when the micro-roughness is not small against the roughness (section 5).
- **Roughness and micro-roughness are not independent.** A micro-roughness near the roughness itself blurs the highlight: physically
  sensible, but the two must be chosen together.
- **Path tracing is left as future work** in the paper (it would need the pixel's footprint carried along the path). That part will be
  this renderer's own.

## 5. Checks (`--scene glints`, both builds, about 6 minutes)

| Check | Result |
|---|---|
| a. The map T there and back; its Jacobian against finite differences; D against the engine's GGX (anisotropic α) | 4 × 10⁻¹⁵; 5 × 10⁻⁸; 9 × 10⁻¹⁶ |
| b. Levels 0 and 1, where every facet is summed, against a plain sum over all facets (which visits 6 standard deviations of the filter) | sum of differences / sum below 7 × 10⁻⁶ |
| c. The average over the surface is the smooth GGX: 8 cases, levels −2 to 9.2, footprints stretched up to 8 : 1, anisotropic roughness | within the standard error (0.002 to 0.3 %) in every case |
| d. Level 3: the far facets replaced by their expectation against the plain sum of all 64 angular cells | averages 1.006 and 1.000, each ± 0.002 |
| e. Energy against the smooth distribution (half vectors drawn from D cos θ, averaged over the surface): α 0.1 / 0.3 / 0.6, micro-roughness 0.01 / 0.05, levels 1.5 / 6.5 | within the error in 10 cases; −2.1 % and −0.8 % for α 0.1 with micro-roughness 0.05 (a lobe half as wide as the distribution) |
| f. Many facets per pixel: the spread of D\* / D over pixels | 0.04 % at level 12, 10⁻⁶ at level 16; exactly 1 beyond level 20 |
| g. (information) the paper's softened roulette, level 2.5 | average off by −0.2 ± 0.3 % |
| h. The material: rough gold with glints on a plane tilted against all three coordinate planes, its BRDF averaged over the plane against the smooth gold's, for pixels of 0.004, 0.04 and 0.6 units; and with no pixel size | −0.9 ± 0.3 %, 0.0 ± 0.3 %, −0.06 ± 0.01 %; exactly the smooth BRDF |
| i. The pixel's size at the hit: a plane 5 m in front of a pinhole camera, and an orthographic view, against closed forms | 0.2 % and 0.0 % (RGB build); 0.4 % and 0.2 % (spectral build: its noise) |
| j. Sampling the facets: six moments (1, x, y, x², y², xy) of the drawn half vectors against the same moments of the density, at levels −1.5, 0.6, 2.4 (stretched footprint, weight 0.6) and 4.3 | within 4 standard errors + 0.002; the share of draws that give a sample (100, 100, 55.6, 75.6 %) equals the density's integral |
| k. The material's reflectance under uniform light at one place, estimated by `sample()` with and without the facets' sampling, three pixel sizes | equal within the errors (0.4 % or better) in both builds |

The averages are unbiased estimates, so they are judged by their standard error (4 standard errors plus 0.3 %); the energy loss is an
intended approximation and is judged by a bound (2 %).
