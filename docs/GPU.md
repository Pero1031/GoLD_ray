# GPU

The renderer runs on the CPU. This document is about a second way to run it: **kernel code** that the CPU and the GPU both compile, with checks that the two agree. The CPU renderer stays the reference; features move into the kernel one at a time, each with a check.

**State (2026-11-27): stage 2, conductors and glass.** A scene file made of quads and triangle meshes, diffuse surfaces, isotropic rough conductors, smooth and rough glass, emitting quads and a uniform environment renders on the GPU (`--device gpu`) and, with the same code, on the CPU (`--device cpu-kernel`). The two pictures are the same bits, and both converge to the engine's. The GGX distribution and the conductor's Fresnel term are one implementation that the engine calls too.

## 1. Design, and where each idea comes from

Before writing code I read how three renderers that run on both devices do it: Mitsuba 3 (with its JIT compiler Dr.Jit), pbrt-v4 and Blender's Cycles. Their sources are cloned (git-ignored) under `references/src/`.

| Decision | GoLD_rayt | Taken from |
|---|---|---|
| One source for both devices | Header-only code under `include/Kernel`, every function marked `RAYT_HD` (`__host__ __device__` under nvcc, nothing under MSVC). No virtual functions, exceptions, heap, std containers or `thread_local` in it. | pbrt-v4 (`PBRT_CPU_GPU`), Cycles (`ccl_device`). Mitsuba reaches the same with a JIT compiler, which is far more machinery than this project needs. |
| The GPU computes in float only | Kernel functions are templates on the scalar type; their double versions exist on the CPU only, and in device code the double overloads are deleted. `build_gpu.ps1` makes ptxas reject any double instruction in the CUDA code and searches the OptiX PTX for one, which catches implicit promotions such as `x * 0.1` with a float `x`. | Cycles builds its CPU kernel with `-Werror=double-promotion`; pbrt-v4 defaults to `Float = float`; Mitsuba always traces in float. |
| The same float bits on both devices | IEEE division and square root and **no implicit fused multiply-adds** (`--fmad=false`; an FMA only where kernel code writes `kfma`). Sine and cosine are our own float polynomials (`ksincos`, the CEPHES algorithm as in Dr.Jit), because the library versions differ in the last bits, and CUDA's falls back to double for huge angles (the guard above found it). | Dr.Jit: explicit FMAs and its own transcendental functions, so all of Mitsuba's backends compute the same thing. |
| Fast math as an option | `build_gpu.ps1 -FastMath`. Accurate enough, but no longer bit-identical to the CPU, and only 5 % faster on the test scene (section 3). | pbrt-v4 and Cycles build their GPU code with `--use_fast_math`. |
| CUDA is optional | The GPU code is a separate DLL, `GoLD_rayt_gpu.dll`, loaded at run time (`Gpu/GpuDevice.hpp`). The CLI builds and runs without CUDA; GPU checks print SKIP. A module of another interface version is refused. OptiX comes with the driver; only its headers are needed to build. | Mitsuba and Cycles load the CUDA / OptiX driver at run time; Mitsuba's tests skip missing devices. |
| Random numbers | A PCG32 stream per camera path, keyed by (run seed, pixel, sample) only (`Kernel/Rng.hpp`). The engine's `thread_local` Mersenne Twister cannot run on a GPU. | Cycles (`path_rng_1D`), pbrt-v4 (stateless samplers). |
| The scene | Flat arrays of plain numbers (`Kernel/KernelScene.hpp`), built from the scene file's description (`Cli/KernelScenePack`), not from the engine's objects. Materials are plain data with a type, chosen by a branch, not virtual classes. | Cycles (`KernelData` and flat arrays, closure structs with a `switch`), pbrt-v4 (`TaggedPointer`). |
| Integrator | A megakernel: one OptiX ray-generation program runs a whole camera path (`kernel::tracePath`, the same function the CPU runs); the closest-hit program only reports the hit; one launch per sample adds that sample of every pixel to a sum, so launches stay short and the sums are added in the CPU's order. | Mitsuba (megakernel by default, no shading in hit programs); Cycles runs the same per-path state machine on the CPU. |
| Hit points | Only the triangle is taken from the tracer (Embree on the CPU, OptiX on the GPU); the kernel recomputes the hit point with its own formula, and moves ray origins off surfaces in integer steps of the float representation (Wächter and Binder, Ray Tracing Gems ch. 6). | Cycles re-intersects the triangle to get a hit point it controls, and uses the same integer offset. |
| Tests | A ladder, section 2. | Mitsuba (`test_renders.py`: per-pixel z-test against a CPU reference), and this project's Mitsuba comparison (`MITSUBA_CHECK.md`). |

What the kernel renderer models: quads, triangle meshes (OBJ / PLY / STL, with their transform) and lathes; diffuse surfaces with a constant albedo (smooth vertex normals as the engine uses them); isotropic rough conductors with the Kulla-Conty lobe; dielectrics, smooth (delta reflection and refraction) or rough (GGX reflection and transmission), also against water (`Kernel/Materials.hpp`); emitting quads, one- or two-sided, sampled by next-event estimation (chosen by power, a uniform point on the quad) with the power heuristic against BSDF sampling; a uniform environment; a pinhole camera; the engine's path depth and Russian roulette; the RGB build. Anything else in a scene file or on the command line is refused with a list of what is missing (for example an environment map, dispersion or a medium inside glass, an anisotropic or glinting conductor, `--denoise`); `--device cpu` renders it.

How a material moves into the kernel (stage 2): the physics that both need goes into a kernel header as templates on the scalar type, in the engine's order of operations, and the engine calls it (`Kernel/Microfacet.hpp`: GGX D, Smith Lambda, sampling of the visible normals, the conductor's Fresnel term; `Microfacet/GGX.cpp` and `fresnel::fresnelConductor` now call these). The engine's renders of conductor scenes stayed the same bits (gold, anisotropic, multiple-bounce and glint spheres; RGB and spectral). The kernel's material (`Kernel/Materials.hpp`) repeats the engine class's steps for the case it models and takes its parameters from the engine's object (`RoughConductor::kernelParams`: the RGB index the engine derived from the metal's table, alpha, the Kulla-Conty table row), so nothing is derived twice. The dielectric works the same way (`Dielectric::kernelParams`: the relative index, alpha, smooth or not); for it the path tracer learned delta events (no light sample at a smooth surface; the throughput takes the event's weight, as `PathIntegrator::updateThroughput`) and refraction (a refracted ray leaves from the other side of the surface; a light sample may lie behind a rough surface).

## 2. Checks

`--scene gpu` (both builds; the GPU parts print SKIP without the module or a device):

1. **Kernel = engine, on the CPU in double.** Each kernel sampling function mirrors an engine function (same formula, same order of operations), so their doubles must be the same bits: cosine-weighted hemisphere, uniform sphere, orthonormal basis, dielectric Fresnel; 4096 cases with edges. All bit-identical.
2. **White furnace, kernel renderer on the CPU**: `scenes/gpu_box.json` with every albedo 1 and no lamp in a uniform environment of radiance 1, where every pixel converges to 1. Mean 1.00036 ± 0.00074.
3. **Random numbers, GPU = CPU, bit for bit.**
4. **Floats, GPU = CPU, bit for bit**, per function: 100 % identical on a Quadro M4000. (With `-FastMath`: the GPU's floats must be as close to the double reference as the CPU's; they are, 61 to 84 % bit-identical.)
5. **White furnace on the GPU**: the same 1.00036 ± 0.00074 as on the CPU.
6. **The test scenes on the CPU kernel and on the GPU with the same seed: the same bits** in every pixel, `scenes/gpu_box.json`, `scenes/gpu_metals.json` and `scenes/gpu_glass.json` (`-FastMath`: the same mean within the noise).
7. **The kernel's conductor in double = the engine's RoughConductor**, eval, pdf and sample on 4096 direction pairs (shading normals tilted from the geometric ones, half of the pairs near the mirror direction), for rough gold with the Kulla-Conty lobe, smooth gold, and rough gold without the lobe: 100 % bit-identical.
8. **The kernel's dielectric in double = the engine's Dielectric**, the same way, for smooth glass, rough glass (alpha 0.09) and rough glass in water (relative index 1.125), entering and leaving, with a third of the pairs near the refracted direction: 100 % bit-identical. A rough surface's choice between reflection and refraction is drawn inside `Dielectric::sample` from the engine's own stream; the check seeds that stream, reads the number the engine will draw, seeds it again, and hands the same number to the kernel.
9. **Each material in float, GPU = CPU, bit for bit** (diffuse, the three conductors, the three dielectrics; eval, pdf and sample): 100 %.
10. **The tracers agree:** every ray a CPU path trace of the three test scenes casts (about 60 000 each, closest-hit and shadow) is traced again by OptiX: the same triangle, the same blocked shadow rays, everywhere (at most 1e-4 of the rays may differ). Since the kernel takes only the triangle from the tracer, this is the only place where the devices' paths could part.

The tracer check found the one way they did part: the first glass scene had its vase standing exactly in the floor's plane. Rays inside the glass then hit two coincident surfaces at the same distance, and Embree and OptiX settled the tie differently (both rightly): 83 of 38 000 rays, 2.5 % of the pixels not bit-identical, and a mean 0.1 % off at 8 spp. With the floor 1 mm lower, every pixel is the same bits. Coincident surfaces are ambiguous for any renderer; the engine settles them a third way.

A fixed tolerance on |GPU − CPU| in float would be the wrong test: just past the critical angle (glass to air, cos θ 0.7458 against a critical 0.7454) the Fresnel term takes the square root of a small number, which turns a rounding difference in its input into 1.4e-5 in its output. With the same operations on both devices the difference is zero instead.

**Against the engine** (`tools/test/gpu_compare.py`, the method of `MITSUBA_CHECK.md`: each renderer draws the scene twice with independent seeds, the noise of each pixel comes from the difference of its two renders). `scenes/gpu_box.json`, 400 × 300, 1024 spp, blocks of 16 × 16 pixels, two sets of seeds:

| Seeds | Image mean kernel / engine (R, G, B) | z of the means | Blocks with \|z\| > 3 | Median relative difference / noise |
|---|---|---|---|---|
| 1–4 | 1.00014, 1.00014, 1.00013 | 2.1, 2.0, 1.7 | 0.39 % | 0.049 % / 0.087 % |
| 11–14 | 1.00002, 1.00002, 1.00004 | 0.3, 0.3, 0.4 | 0.77 % (1 of 173) | 0.060 % / 0.088 % |

The difference changes sign between the two sets (mean block z +0.17 and −0.15) and stays within the noise: the float kernel converges to the double engine's picture. The background, which sees only the environment, differs by 1e-5: the rounding of summing 1024 float samples.

`scenes/gpu_glass.json` (the box with a smooth glass vase and a frosted glass cube; the caustics are found by BSDF sampling only), seeds 1–4: image mean kernel / engine 1.00000, 0.99996, 1.00002 (z 0.0, −0.3, 0.2), 0.36 % of the blocks beyond |z| = 3, median relative difference 0.19 % against a noise of 0.28 %. Per render of 1024 spp: engine 77 s, GPU 20 s.

`scenes/gpu_metals.json` (the same box with a rough gold vase, a copper panel and a nearly smooth silver plate), seeds 1–4: image mean kernel / engine 1.00007, 1.00011, 1.00005 (z 0.7, 1.2, 0.5), 0.65 % of the blocks beyond |z| = 3 (1 of 154), median relative difference 0.098 % against a noise of 0.154 %. Per render of 1024 spp: engine 80 s, GPU 18 s.

## 3. Speed

`scenes/gpu_box.json` (780 triangles), 400 × 300, 256 spp, depth 16, on this PC (i7-6700K, 8 threads; Quadro M4000, Maxwell, no RT cores):

| Renderer | Render time |
|---|---|
| engine, `--device cpu` (double, every feature) | 18.2 s |
| kernel on the CPU, `--device cpu-kernel` (float, Embree) | 4.4 s |
| kernel on the GPU, `--device gpu` | 2.2 s (+ 0.3 s to set up OptiX) |
| kernel on the GPU, built with `-FastMath` | 2.1 s |

The GPU is twice as fast as the same kernel on the CPU here; the larger step from the engine to the kernel comes from float and from the much smaller set of features a path has to check. Fast math hardly helps: the time goes into ray traversal and memory, not arithmetic, so the bit-exact default costs about 5 %. A GPU with RT cores traverses rays in hardware; the ratio will look different there.

## 4. Build and run

```
powershell -ExecutionPolicy Bypass -File tools\build_cli.ps1 -Dev      # the CLI (no CUDA needed)
powershell -ExecutionPolicy Bypass -File tools\build_gpu.ps1           # GoLD_rayt_gpu.dll next to it (about 10 s)
cd GoLD_rayt
..\out\cli-dev\GoLD_rayt.exe --scene gpu
..\out\cli-dev\GoLD_rayt.exe scenes\gpu_box.json --device gpu
```

`build_gpu.ps1 -Arch 75,86,89` compiles machine code for other GPUs; `-FastMath` as above. PTX (the GPU's intermediate code) of the lowest architecture is always embedded, so a newer GPU runs the module too: its driver compiles the PTX on the first launch. The design is not tied to the current card. Recent Mitsuba already requires compute capability 7.5 (Turing) or newer, and CUDA 13 no longer compiles for Maxwell, so a Maxwell card like the M4000 needs CUDA 12.x. The OptiX SDK headers are found in `OPTIX_ROOT` (default: OptiX SDK 8.1.0).

Files: `include/Kernel/` (shared code: `HostDevice`, `KernelMath`, `Rng`, `Sampling`, `KernelScene`, `PathTracer`, `SelfTest`), `include/Gpu/` (the module's interface and loader), `GoLD_rayt/gpu/` (the module: `GpuModule.cu` CUDA, `OptixRender.cpp` OptiX host code, `OptixPrograms.cu` the OptiX programs), `src/Cli/KernelScenePack.cpp` (scene file to kernel scene, the CPU kernel renderer, `--device`), `src/Cli/Checks/GpuChecks.cpp`, `tools/build_gpu.ps1`, `tools/test/gpu_compare.py` (against the engine), `tools/test/pfm_compare.py` (two PFM pictures pixel by pixel).

Writing kernel code: float constants as `T(0.1)`, never bare double literals (the build stops otherwise); `kfma` where a fused multiply-add is wanted; `ksin` / `kcos` / `ksincos` instead of the library's.

## 5. Next stages

1. **More of the engine, the same way:** the coated diffuse; environment maps; textures; media inside glass.
2. Spectral mode; timing on larger meshes (the scanned statues); a wavefront integrator only if the timings call for it.
