# LVC-BPT (bidirectional path tracing with a light vertex cache): a prototype

A prototype that renders light shafts and caustics in fog (kiri), on diffuse surfaces and in glass (with dispersion) faster than path tracing (in the same time, the variance of the noise is 10 to 60 times smaller) (2026-11-08). The implementation is `GoLD_rayt/include/Render/VolumeBPT.hpp` (header only). The design learned from the talk at Ray Tracing Camp 10 (listed in the [references](REFERENCES.md)) was rebuilt in this CPU engine.

## 0. Summary

- **N light paths are traced first and their vertices are stored** (Light Vertex Cache). The camera paths are connected to those vertices.
- There are three connection strategies: **NEE, connection to a light vertex, and light tracing**. Every path is weighted by the **balance heuristic MIS** from the probability densities of all strategies.
- **This is not a technique for making the GPU faster.** It is a sampling method that combines the strengths of the light side and the camera side. In the talk (GPU), it appears as a form that keeps the state small and reduces branching (the light vertices are stored in an array). On the CPU here, the advantage of reusing light paths for many camera paths works as it is.
- Glass is put into the sub-paths as a "vertex that cannot be connected". For dispersion, **the wavelength is fixed to one per render pass**.
- Verification is `--scene lvc` (agreement with PT and with the light tracing of photons).

![The same render of about 9 seconds. Top: PT + light tracing of photons, bottom: this prototype (dispersive prism)](img/lvc_equal_time_prism.png)

## 1. What problem this method solves

- **Path tracing (PT)** extends a path from the camera, and it contributes only if it happens to hit a light. It almost never hits point lights, thin lights, or light that passed through glass (caustics). A thin beam in fog has a low probability that a scattering point enters the beam, which gives speckle noise.
- **Light tracing (LT)** extends a path from the light and connects each vertex to the camera. It is strong for caustics and thin beams, but weak for the camera-side appearance (such as multiple reflections on diffuse surfaces).
- **Bidirectional path tracing (BPT)** makes both and connects the vertices to each other. The same path can be made in several ways, so it is weighted by MIS.

## 2. Strategies and MIS

A path with n vertices `x_1 … x_n` (counted from the light side) can be made in `n + 1` ways, with `s` vertices from the light-side sub-path and `n − s` from the camera-side sub-path.

| s | Name | Content |
|---|---|---|
| 0 | NEE | Connect the camera-side sub-path directly to a light |
| 0 < s < n | Connection | Connect a light-side vertex and a camera-side vertex with a shadow ray |
| n | Light tracing | Connect the light-side sub-path to the camera (splatting) |

Let the probability density of each strategy be `p_s(x)` and its number of samples `n_s`. A path made by any strategy then contributes `f(x) / D(x)`, with `D = Σ_s n_s p_s(x)` (the balance heuristic). `n_s` is the number of samples per pixel `spp` for the strategies that start from the camera side, and the number of light paths `N` for light tracing. Recursive MIS (a number per vertex) is not used; the densities of all strategies are computed every time from the positions of the vertices of the path (a prototype in a form that is easy to check for correctness).

## 3. The cache of light vertices

1. `N` light paths are traced and their vertices (position, normal, albedo, transmittance of the edge before it) are stored in an array.
2. One camera sample picks one of the stored light paths at random and **connects all of its own vertices to all vertices of that path**.
3. Sharing the light paths among many camera samples is still unbiased (each pair is an independent sample from the same distribution). The cost on the light side is almost nothing.
4. In the light-tracing-only case, the vertices are not stored (memory ran out at 50 million light paths).

Note: because the light paths are shared, the noise is correlated between pixels (bands and mottling).

## 4. Densities and contributions (fog and diffuse surfaces)

The measure of a vertex is volume in fog and area on a surface. The density of arriving at a vertex `x` from a point `a` with a directional density `p_dir` (per solid angle) is as follows.

| Vertex | Density | Factor of the contribution |
|---|---|---|
| Fog | `p_dir · σ_s · T(a,x) / d²` | `σ_s · phase function` |
| Diffuse surface | `p_dir · ρ · T(a,x) · \|cos\| / d²` | `ρ/π` |

- `σ_s` is the scattering coefficient (`albedo · σ_t`). The choice between scattering and absorption (probability `albedo`) is included in the density. `ρ` is the probability of not being absorbed at a surface (the albedo).
- `p_dir` is the density of the emission for a light (the inverse of the solid angle of the cone; it includes the mixture of importance sampling), the phase function in fog, `cos/π` on a surface, and `1/Ω_pixel` (the inverse of the solid angle of a pixel) for the camera.
- For each edge, `T · |cos_a| · |cos_b| / d²` is multiplied (cos = 1 at a fog vertex).
- The transmittance of an edge is recorded while the sub-path is traced (from the length traveled through the medium). **Only the edge of a connection** is measured with a shadow ray. A surface blocks the shadow ray, and a null interface (the boundary of the fog) lets it through.

## 5. Glass and dispersion

**Glass (a smooth dielectric) is a specular vertex.**
- It cannot be the endpoint of a connection (the directions do not match; probability zero). A strategy whose endpoint would be specular is removed from the sum in `D`. A chain of glass is treated like one long edge.
- The weight is the W returned by `Material::sample` (1 for reflection, 1/η² for refraction). In a chain that goes in and comes out, the η² cancel.
- In the contribution, the **cos on the side that sampled** that vertex is paired with W. In the density, the directional density of a specular vertex is taken as 1 (the delta is omitted), and the geometric conversion is kept.
- **The MIS weights are normalized**, so even if the handling of the delta is somewhat inaccurate, no bias appears (only the variance increases). Bias appears only when the sampling weight (the contribution) is wrong.

For **dispersion**, the shape of a path changes with wavelength, so the wavelengths of the two paths to be connected must be the same. Therefore **the wavelength is fixed to one per render pass** (`--lvc-passes K`).
- In pass k, one wavelength is chosen from a stratified random number, and both the light paths and the camera paths are traced at that wavelength. Each pass builds its own cache (N/K light paths, camera spp/K), and the K images are averaged.
- The density of the wavelength is common to all strategies and cancels in the MIS ratio.
- K is the number that removes the bands (gaps between wavelengths). As measured, it saturates at K = 16, and **the color grain does not decrease with K. It decreases with the number of light paths N (1/√N)**. A scene without dispersion uses K = 1 (the usual 4 wavelengths).

## 6. Importance sampling of the emission

To pick up caustics efficiently from a light with a wide cone that also lights the floor and walls, half of the light paths are sent into a small cone that encloses the object (`PointLight::setEmissionTarget`). The emitted amount (intensity) does not change, and the mixed density is put into `emissionPdf`, so the MIS weights are consistent.

## 7. Verification (`--scene lvc`)

| Configuration | Reference | What is compared |
|---|---|---|
| A ball of fog | PT (with NEE in the medium) | Each strategy alone (NEE, LT), combinations, all. The mean, and every 6×3 block |
| A ball of fog + a diffuse floor | PT | Same |
| A glass ball in fog | Ray tracing of the photon map (an independent implementation) | The image of L S S M… (the first vertex of the light is the glass) |
| Same, the whole picture | PT + light tracing of photons | Connection + LT, all strategies (blocks that do not show the ball) |
| A dispersive glass ball (Abbe 10, K = 32) | Same | Same (wider tolerance. A path of one wavelength has a large color variance) |

It passes in both RGB and spectral. **Comparing each strategy alone with PT** found the following bugs.
- When leaving the boundary of the fog into the air, the medium of the ray was not updated (on the floor outside the fog, scattering happened although it was in air).
- When the number of stored photons was 0, the map was disabled and counted twice.
- For fog behind the glass, the camera side discarded the path, which made the image 2% darker.

The checks of the dispersive configuration use small images, so the ratio fluctuates by about ±5% (there is no sign bias, and the mean RGB of the real scene agrees with the reference within 1%).

## 8. Results (variance of the noise at the same time. PT / this prototype)

| Configuration (640 px, about 3 to 9 seconds) | Beam | Haze | Caustics |
|---|---|---|---|
| A spot in fog | 43 times | 280 times | |
| Fog and a diffuse floor | 10 times | 14 times | The circle of light on the floor (direct light, which PT is good at) 1.8 times |
| A glass ball in fog | 31 times | 58 times | 9 times |
| A dispersive prism in fog | 15 times | 61 times | The fan of the rainbow 12 times |

The same image was rendered twice independently, the noise was estimated from the difference, and the values are the ratio of the values multiplied by the time. It works in the regions where PT is weak (thin beams, indirect light). For direct light, which PT is good at, the difference is small.

## 8b. Vertex merging (specular-diffuse-specular paths)

A path with only one connectable vertex, such as the floor under a glass ball that the camera sees through the ball (`E S S D S S L`), cannot be made by any connection: both neighbors of the floor vertex are specular. **Vertex merging** (Georgiev et al. 2012; on the light vertex cache as in Davidovič et al. 2014) makes it: a camera vertex is paired with every cached light vertex within a radius `r`, as if they were the same vertex.

- **Strategy "merge at vertex k"**: the light sub-path gives `x_1..x_k`, the camera sub-path `x_n..x_k`. Its density is `pdfF[1..k] · pdfC[k..n] · V`, with `V = π r²` on a surface and `4/3 π r³` in the fog (the measure of the kernel), and every camera sample is paired with all `N` light paths (count `spp · N`). It joins the same balance-heuristic sum `D` as NEE, connections and light tracing, for every path (a path with a specular vertex next to the merge vertex is what makes it useful).
- **Bias and consistency**: for a fixed radius the estimate is blurred (biased, `O(r²)` inside a smooth region, `O(r)` near the boundary of the fog or of a surface). Each pass has its own light paths, and the radius shrinks as `r_i = r_1 √(i^(α−1))`, `α = 0.75` (`--lvc-passes`), so the bias vanishes with the number of passes. A single pass keeps its radius `r_1`.
- **Implementation**: a hash grid over the cached non-specular vertices (cell `2r`, rebuilt every pass), a range query at every non-specular camera vertex, and the path rebuilt as `z_1..z_{k−1}, y_t, y_{t−1}..y_1` (the light vertex is replaced by the camera vertex, its edge keeps the recorded transmittance). The same `pathContribution` evaluates it with `sActual = k`.
- **Switch**: `--lvc-merge R` (scene units; 0 = off, the default) and bit 8 of `--lvc-mask` (default mask 15).
- **Checks** (`--scene lvc`): the fog and the fog with a floor (no glass) with merging, the glass ball in fog (light-side glass, merging at fog vertices, also dispersive), a directly lit floor seen through a glass ball (camera-side glass, no caustic: next event estimation, merging and path tracing agree), and the caustic of a ball on a floor seen through the ball (the one case no connection can make; the reference is path tracing with the photon map).
- **Cost**: a camera vertex inside a dense caustic merges thousands of light vertices, each one rebuilding the path and its weights (`O(n)` per pair): about 15 to 20 times the time of the connections alone in the ball-on-a-floor check. The recursive weights of VCM (stored per vertex) would remove that.

## 8c. Manifold connections (two refractions through a flat-faced body)

A connection needs the two vertices to see each other. When a glass body is in the way, the path between them is bent by two refractions, and that chain is fixed by Snell's law once the two ends are fixed: the "generalized half vector" at each refraction must be parallel to the facet normal. **Manifold exploration** (Jakob and Marschner 2012; Hanika et al. 2015 for the light; Speierer et al. 2018 for any two vertices of a bidirectional path tracer) solves it with Newton's method, so the connection is made exactly instead of by chance.

- **Strategies.** *Manifold NEE*: a camera vertex and the light. *Manifold connection*: a camera vertex and a vertex of the cached light path. The chain is the two refractions through ONE body of `SpecularNee` (flat facets, as in the manifold NEE; no reflections, no total internal reflection, no curved glass). `SpecularNee::chains(x, y)` returns every chain with its Fresnel transmittance `(1-F_A)(1-F_B)` and `dω/dA`, the change of the solid angle at `y` per area at `x`, from two re-solved neighbours.
- **Contribution.** The three edges and the two specular vertices of the chain are replaced by `T_out T_out' · (dω/dA) · (1-F_A)(1-F_B)` (the cosine at the camera-side end is inside `dω/dA`). The rest of the path is evaluated as before.
- **Weights.** The strategy joins the balance-heuristic sum `D` of every path that has two refractions through the same body between two connectable positions, with density `pdfF[1..a] · pdfC[b..n] · Γ` (`a`, `b = a + 3` the connected positions). `Γ` stands in for the two solved vertices, which are not sampled (the light-side conversions of the two refractions); the contribution carries the same factor, so the estimate stays unbiased whatever `Γ` is, and only the variance depends on it. A manifold term is counted in `D` only if the solver, started from the two ends, finds exactly the chain of the path (`manifoldReproduces`): the rule that decides which chains the strategy makes must be the rule that weighs them, or energy is lost.
- **Where it applies.** With `manifoldOnlyBlocked` (the default, as in the paper) only where the straight segment between the two vertices is blocked; the tests switch it off to compare with the photon map. The medium of the two outside segments is taken from the camera-side vertex (the body and the vertices share a medium). Dispersive glass needs `--lvc-passes > 1` (the chain is solved at the wavelength of the pass). Lights must be point lights.
- **Switch.** Bits 16 (manifold NEE) and 32 (manifold connections) of `--lvc-mask` (default 15: off); the scene must list its flat-faced glass bodies (the prism of `--scene fogspot --spot-prism`).
- **Checks** (`--scene lvcmanifold`, 2.5 minutes at `--spp 32`):
  1. `chains()` against `SpecularNee::sample()` on 197 floor points: identical (the glue: frame, Fresnel, `dω/dA`).
  2. Manifold NEE alone, with a camera path of one vertex, against the photon map of the two-refraction chains: the sum is within 7 % (the tolerance of the manifold check itself is 10 %), and every chain that is generated is also found again by the weights.
  3. Adding the manifold strategies to the other strategies changes the picture by 0.4 % (noise): adding a strategy does not change the expectation. With merging added, the picture is within 2 % of path tracing plus the photon map.
  4. The camera looks through the prism (a chain on both ends): manifold and merging add the same energy to the strategies without them.
  5. The noise of the picture with and without the manifold strategies (two independent renders).
- **Not shown.** The solver does not find every chain that exists (the same coverage as the manifold NEE: a few percent of the energy of the photon map is missed, and none of the paths that end on a facet's plane within a millimetre).

## 9. Usage

```bash
# A prism in fog (about 9 minutes)
GoLD_rayt.exe --scene fogspot --spot-prism --abbe 6 --spot-deg 1.4 --width 1280 --lvc --lvc-passes 96 --spp 192 --lvc-paths 600000000 --fog-light 150

# A gem (about 11 minutes)
GoLD_rayt.exe --scene caustics --gal-object gem --abbe 9 --width 1000 --lvc --lvc-passes 64 --lvc-depth 16 --spp 512 --lvc-paths 400000000 --gal-dist 0.8 --gal-scale 1.4 --gal-light 8 --gal-albedo 0.08 --gal-spread 8 --gal-look 0.02,0.02,-0.015 --gal-fov 9 --denoise

# Self-check
GoLD_rayt.exe --scene lvc
```

- `--lvc`: render with this method (without it, PT). `--lvc-paths N`: total number of light paths. `--lvc-depth D`: number of vertices per sub-path. `--lvc-passes K`: number of wavelength passes (64 to 128 for dispersive glass). `--lvc-mask M`: strategies (1 = NEE, 2 = connection, 4 = light tracing, 8 = vertex merging, which needs `--lvc-merge R`).
- The render time is decided almost by "the number of light paths" (about 0.3 to 1 microsecond per path, 8 threads).
