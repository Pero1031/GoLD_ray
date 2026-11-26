# GPU

The renderer runs on the CPU. This document is about a second way to run it: **kernel code** that the CPU and the GPU both compile, with checks that the two agree. The CPU renderer stays the reference; features move into the kernel one at a time, each with a check.

**State (2026-11-26): stage 0.** The kernel headers, the GPU module and `--scene gpu` exist. The GPU computes in float only; double is for the CPU. No image is rendered on the GPU yet.

## 1. Design, and where each idea comes from

Before writing code I read how three renderers that run on both devices do it: Mitsuba 3 (with its JIT compiler Dr.Jit), pbrt-v4 and Blender's Cycles. Their sources are cloned (git-ignored) under `references/src/`.

| Decision | GoLD_rayt | Taken from |
|---|---|---|
| One source for both devices | Header-only code under `include/Kernel`, every function marked `RAYT_HD` (`__host__ __device__` under nvcc, nothing under MSVC). No virtual functions, exceptions, heap, std containers or `thread_local` in it. | pbrt-v4 (`PBRT_CPU_GPU`), Cycles (`ccl_device`). Mitsuba reaches the same with a JIT compiler, which is far more machinery than this project needs. |
| The GPU computes in float only | Kernel functions are templates on the scalar type. Their double versions (math overloads, double random numbers, the checks' references) exist on the CPU only; in device code the double overloads are deleted, and `build_gpu.ps1` makes ptxas reject any double instruction, which catches implicit promotions such as `x * 0.1` with a float `x`. | Cycles builds its CPU kernel with `-Werror=double-promotion`; pbrt-v4 defaults to `Float = float` and its double build is fragile on the GPU; Mitsuba always traces in float. |
| The same float bits on both devices | The module is built with IEEE division and square root and **without implicit fused multiply-adds** (`--fmad=false`): an FMA happens only where kernel code writes `kfma`, on both devices. Sine and cosine are our own float polynomials (`ksincos`, the CEPHES algorithm as in Dr.Jit), because the C library's and CUDA's `sinf` / `cosf` differ in the last bits, and CUDA's falls back to double for huge angles (the guard above found it). | Dr.Jit: explicit FMAs and its own transcendental functions, so all of Mitsuba's backends compute the same thing. |
| Fast math as an option | `build_gpu.ps1 -FastMath` (approximate division and square root, implicit FMAs, denormals flushed). Accurate enough (section 2), but no longer bit-identical to the CPU. Its speed will be measured in stage 1. | pbrt-v4 and Cycles build their GPU code with `--use_fast_math`. |
| Precision on the CPU | The CPU renderer stays double. The kernel can run on the CPU in float (to pair with the GPU) or in double (to compare with the engine). | Mitsuba's `*_double` variants. |
| CUDA is optional | The GPU code is a separate DLL, `GoLD_rayt_gpu.dll`, loaded at run time (`Gpu/GpuDevice.hpp`). The CLI builds and runs without CUDA; GPU checks print SKIP. A module of another interface version is refused. | Mitsuba and Cycles load the CUDA / OptiX driver at run time; Mitsuba's tests skip missing devices. |
| Random numbers | A PCG32 stream per camera path, keyed by (run seed, pixel, sample) only (`Kernel/Rng.hpp`), so both devices draw the same numbers for the same sample. The engine's `thread_local` Mersenne Twister cannot run on a GPU. | Cycles (`path_rng_1D`), pbrt-v4 (stateless samplers). |
| Materials (planned) | Plain data plus a `switch` on the type, instead of the virtual `Material`. The existing material classes become thin wrappers that call the kernel function, so each material keeps one implementation. | Cycles (closure structs and `switch`), pbrt-v4 (`TaggedPointer`). |
| Integrator (planned) | A megakernel: one OptiX ray-generation program runs the whole path loop, and the hit program only reports the hit; shading runs in the ray-generation program. A wavefront integrator only if profiling asks for it. | Mitsuba (megakernel by default, no hit programs); Cycles runs the same per-path state machine on the CPU. |
| Tests | A ladder, see section 2. | Mitsuba (`test_renders.py`: per-pixel z-test against a CPU reference, pass when 99.75 % of the pixels pass), and this project's Mitsuba comparison (`MITSUBA_CHECK.md`). |

## 2. Checks

`--scene gpu` (both builds; 4096 cases with edges: u = 0, the poles of the basis, n.z = -0, total internal reflection, cos θ < 0):

1. **Kernel = engine, on the CPU in double.** Each kernel function mirrors an engine function (same formula, same order of operations), so their doubles must be the same bits: cosine-weighted hemisphere, uniform sphere, orthonormal basis, dielectric Fresnel. All bit-identical.
2. **Random numbers, GPU = CPU, bit for bit.** Identical.
3. **Floats, GPU = CPU, bit for bit**, per function. Measured on a Quadro M4000: 100 % identical. The float error against the double reference is then the same on both devices (cosine hemisphere 3.7e-7, uniform sphere 4.1e-7, basis 1.4e-7, Fresnel 1.0e-5).

With `-FastMath` check 3 changes to "the GPU's floats are as close to the double reference as the CPU's (at most twice the CPU's error + 1e-7)". It passes, with 61 to 84 % of the values bit-identical; the basis error grows from 1.4e-7 to 2.5e-7, and the Fresnel error even shrinks (3.9e-6, fused multiply-adds round once instead of twice).

Why a fixed tolerance on |GPU − CPU| is the wrong test: just past the critical angle (glass to air, cos θ 0.7458 against a critical 0.7454) the Fresnel term takes the square root of a small number, which turns a rounding difference in its input into 1.4e-5 in its output. That is the problem's conditioning, not a device error. With the same operations on both devices the difference is zero instead.

Planned image checks, each printing SKIP without a device:

- the same scene and seed on the CPU kernel (float) and on the GPU, pixel by pixel (with the same float bits, paths can only diverge where Embree and OptiX disagree on a hit);
- the GPU against the CPU renderer within the Monte Carlo noise, by the method of `MITSUBA_CHECK.md` (two independent renders give the noise; mean ratio ± standard error and the share of 16 × 16 blocks beyond 3 standard errors) and Mitsuba's per-pixel z-test;
- white furnace scenes on the GPU.

## 3. Build and run

```
powershell -ExecutionPolicy Bypass -File tools\build_cli.ps1 -Dev      # the CLI (no CUDA needed)
powershell -ExecutionPolicy Bypass -File tools\build_gpu.ps1           # GoLD_rayt_gpu.dll next to it (8 s)
cd GoLD_rayt
..\out\cli-dev\GoLD_rayt.exe --scene gpu
```

`build_gpu.ps1 -Arch 75,86,89` compiles machine code for other GPUs; `-FastMath` as above. PTX (the GPU's intermediate code) of the lowest architecture is always embedded, so a newer GPU runs the module too: its driver compiles the PTX on the first launch. The design is not tied to the current card. Recent Mitsuba already requires compute capability 7.5 (Turing) or newer, and CUDA 13 no longer compiles for Maxwell, so a Maxwell card like the M4000 needs CUDA 12.x.

Writing kernel code: float constants as `T(0.1)`, never bare double literals (ptxas stops the build otherwise); `kfma` where a fused multiply-add is wanted; `ksin` / `kcos` / `ksincos` instead of the library's.

## 4. Next stages

1. **A minimal path tracer in the kernel:** triangle meshes, diffuse surfaces, area lights and a uniform environment, RGB. It runs on the CPU (OpenMP, Embree in float) and on the GPU (OptiX megakernel), with the image checks above, and the first timings (default against fast math).
2. **Materials one at a time** (rough conductor, dielectric, coated diffuse, ...): the engine class calls the kernel function, the engine's renders stay the same, and eval / pdf / sample agree between the devices on fixed inputs (plus χ² of sample against pdf on the GPU).
3. Spectral mode, textures, timing against the CPU renderer; a wavefront integrator only if the timings call for it.
