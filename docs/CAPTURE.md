# Importing captured objects (scans, 3D Gaussian Splatting) (2026-11-13)

**Status:** Stage 1 (bring the object in as a mesh and replace its material appearance (shitsukan) with a physical material) and stage 2 (place the Gaussians themselves around the scene as emitting, occluding particles) are implemented. 3DGS training is not part of this repository: this document starts from its outputs (a mesh, `point_cloud.ply`), which are produced on a machine with a GPU.

## 1. Aim

3DGS reproduces "how the scene looked under that lighting" from multi-view photographs. It has no physical quantities such as reflectance or scattering coefficients, so mirrors and translucent objects break down and the lighting cannot be changed. This renderer is the opposite: it has physical quantities (spectra, subsurface scattering, layers), but the shapes are made by hand. The sensible approach is to **borrow the shape from 3DGS and replace the material appearance here**, which uses the strengths of both as they are.

| Stage | Content | Status |
|---|---|---|
| 1 | Photograph → 3DGS → mesh → clean up → load → relight with a physical material | Done (this document) |
| 2 | Read a 3DGS `.ply` (Gaussians) directly and place it around the scene as emitting particles (the same idea as NVIDIA's 3DGRT) | Done (section 4) |

## 2. Procedure (stage 1)

**(1) Photograph → 3DGS → mesh** (on a PC with a GPU). Use an implementation that can output a mesh: 2DGS (PLY by `mesh extraction`), Gaussian Opacity Fields, SuGaR (textured OBJ). When the GPU is weak, shrink the images to 1000-1600 px on the long side and keep the iterations to 15k-30k. 2DGS is lighter than the original 3DGS and gives a clean mesh. The subject should be a still object (stone, metal, glass, translucent objects). Human faces and skin move, so they are not suitable.

**(2) Clean up** (`tools/data/mesh_clean.py`; reads and writes OBJ and PLY):

```bat
python tools\data\mesh_clean.py mesh.ply out.ply --keep-uv --up=-y --drop-small 0.01 --height 0.3 --ground
```

- `--up=-y`: rotates the COLMAP / 3DGS coordinates (y down) to the renderer's y up. `--up=z` is for z up (Blender, Open3D).
- `--drop-small 0.01`: discards detached fragments (floaters) that hold less than 1% of all triangles.
- `--height 0.3 --ground`: scales to a height of 0.3 m, puts the bottom at y = 0 and the center of the footprint at x = z = 0 (units are meters).
- `--keep-uv`: keeps the uv when the photographed texture is to be applied (otherwise they are dropped to make the file smaller). `--cap`: closes holes when the inside is filled with a scattering medium or `randomWalk`.
- Normals are always discarded and rebuilt from the geometry (the orientation of scanned normals cannot be trusted). If the signed volume is negative, use `--flip`.

**(3) Place it in the scene** (`docs/SCENE_FILE.md`). A `.ply` can be given directly as the `file` of a `mesh`; an STL scan (binary or ascii) is read the same way, with coincident vertices welded and smooth normals rebuilt from the winding (`--scene mesh`, c3), but it has no uv or color, and its unit is usually mm with Z up, so set `scale` and `rotate`. The appearance as photographed is the `texture` of a `diffuse` material; the replacement can be any material.

```jsonc
"materials": { "captured": { "type": "diffuse", "albedo": 1.0, "texture": "assets/x/texture.jpg" },
               "gold":     { "type": "conductor", "metal": "Au", "roughness": 0.2 } },
"objects":   [ { "shape": "mesh", "file": "assets/x/out.ply", "material": "gold" } ]
```

**Example:** `scenes/capture_erato.json`. The Erato scan (`assets/Meshes/erato/`) is converted to PLY with `--keep-uv`; the left copy keeps the photographed texture, the right copy is gold, both drawn with the same 3 lights. The picture is `out\studies\capture\capture_erato.png`.

## 3. What was checked

- `--scene mesh`, c2: the round trip of PLY writing and reading (ascii: per-vertex uv + skipping the color + a `TextureFile` comment; binary: a per-face `texcoord` list). Triangle count, normals and volume agree, and the v of the uv comes in flipped (texture convention: in `Core/Texture.hpp` v points down; the OBJ reader was made to flip in the same way). PASS in the RGB build.
- Round trip of `mesh_clean.py`: a torus (`--up=-y --height 0.1 --ground` → PLY → OBJ), and Erato (OBJ with 1,237,495 face corners → 227,503 vertices welded with uv kept, 12 seconds, PLY 12.6 MB; 1/8 of the 107 MB OBJ).
- Rendering: the texture is applied in the correct orientation (face, lyre) and is replaced by gold under the same light. 900 × 600, 128 spp, about 3 minutes (development build).

## 4. Stage 2: drawing the Gaussians as particles (2026-11-13)

**What it does:** It reads the `point_cloud.ply` that 3DGS writes during training (position, opacity, 3-axis size, orientation quaternion, color as spherical harmonics up to degree 3) as it is, and path-traces each particle as an **emitting, occluding ellipsoidal Gaussian**. The response of a particle to a ray is the value at t*, the point closest to the center in the covariance distance: α = opacity × exp(−½ |x(t*) − μ|²_Σ) (the "ray-Gaussian response" of 3D Gaussian Ray Tracing, Moenne-Loccoz et al. 2024). The particles a ray meets are sorted by t* and composited front to back, L = Σ T_i α_i c_i, T = Π(1 − α_i). There is no scattering or reflectance: c_i is the radiance the capture saw, and the particles are **background and light source** (visible when a path passes through them, and reflected in metals), while objects behind them are occluded by T. A BVH over the particles' boxes (the range where α reaches alphaMin) collects the particles of a ray. Code: `Geometry/GaussianSplats.hpp`, `Scene::setSplats`, `Integrator.cpp` (compositing per ray) and `Render/Transmittance.hpp` (a shadow ray is only multiplied by T).

**Usage:** the `splats` object of the scene file (`docs/SCENE_FILE.md`). The template is `scenes/capture_splats.json` (`rotate` turns COLMAP's y down upright, `keep` discards distant floaters, `remove` takes out the particles where an object is to be placed, `intensity` matches the brightness). A sample without data is `--scene splats` (a gold sphere and a frosted glass sphere inside a synthetic dome of 40,000 particles; the picture is `out\studies\capture\splats_demo.png`).

**What was checked** (checks a-g of `--scene splats`; PASS in both RGB and spectral): α of one particle (center, 1σ, cut by a surface in front); agreement of the closed form for a rotated anisotropic particle with the numerical maximum along the ray (786 rays, difference 1e-6 or less); the order of compositing 5 particles (insertion order scrambled) and cutting off at a surface; agreement of the BVH with brute force (3000 particles, 5000 rays); the round trip of a 3DGS-format PLY (logit, log, quaternion normalization, the channel-major order of f_rest, the degree-1 color against a hand calculation); the scene transform (the center, the box, and the harmonics rotate together with the cloud); the integrator (a black plate with α = ½ placed between the light and the floor is crossed once by the camera ray and once by the ray to the light, so the floor radiance becomes ¼: estimate 0.09192 ± 0.00011, theory 0.09185). With a synthetic cloud in 3DGS format (20,000 particles + 200 at the object's place + 100 far away), it was also confirmed that `rotate`, `keep` and `remove` work.
