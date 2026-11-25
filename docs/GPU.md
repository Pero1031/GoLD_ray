# GPU

The renderer runs on the CPU. This document is about a second way to run it: **kernel code** that the CPU and the GPU both compile, with checks that the two agree. The CPU renderer stays the reference; features move into the kernel one at a time, each with a check.

**State (2026-11-26): stage 0.** The kernel headers, the GPU module and `--scene gpu` exist. No image is rendered on the GPU yet.

## 1. Design, and where each idea comes from

Before writing code I read how three renderers that run on both devices do it: Mitsuba 3 (with its JIT compiler Dr.Jit), pbrt-v4 and Blender's Cycles. Their sources are cloned (git-ignored) under `references/src/`.

| Decision | GoLD_rayt | Taken from |
|---|---|---|
| One source for both devices | Header-only code under `include/Kernel`, every function marked `RAYT_HD` (`__host__ __device__` under nvcc, nothing under MSVC). No virtual functions, exceptions, heap, std containers or `thread_local` in it. | pbrt-v4 (`PBRT_CPU_GPU`), Cycles (`ccl_device`). Mitsuba reaches the same with a JIT compiler, which is far more machinery than this project needs. |
| Precision | The scalar type is a template parameter: float on the GPU, float or double on the CPU. Comparing CPU-float with the GPU separates device bugs from precision; CPU-double with CPU-float measures the cost of float. | Mitsuba's `*_double` variants; all three trace rays in float. |
| CUDA is optional | The GPU code is a separate DLL, `GoLD_rayt_gpu.dll`, loaded at run time (`Gpu/GpuDevice.hpp`). The CLI builds and runs without CUDA; GPU checks print SKIP. | Mitsuba and Cycles load the CUDA / OptiX driver at run time; Mitsuba's tests skip missing devices. |
| Random numbers | A PCG32 stream per camera path, keyed by (run seed, pixel, sample) only (`Kernel/Rng.hpp`), so both devices draw the same numbers for the same sample. The engine's `thread_local` Mersenne Twister cannot run on a GPU. | Cycles (`path_rng_1D`), pbrt-v4 (stateless samplers). |
| Materials (planned) | Plain data plus a `switch` on the type, instead of the virtual `Material`. The existing material classes become thin wrappers that call the kernel function, so each material keeps one implementation. | Cycles (closure structs and `switch`), pbrt-v4 (`TaggedPointer`). |
| Integrator (planned) | A megakernel: one OptiX ray-generation program runs the whole path loop, and the hit program only reports the hit; shading runs in the ray-generation program. A wavefront integrator only if profiling asks for it. | Mitsuba (megakernel by default, no hit programs); Cycles runs the same per-path state machine on the CPU. |
| Tests | A ladder, see section 2. | Mitsuba (`test_renders.py`: per-pixel z-test against a CPU reference, pass when 99.75 % of the pixels pass), and this project's Mitsuba comparison (`MITSUBA_CHECK.md`). |

## 2. Checks

`--scene gpu` (both builds):

1. **Kernel = engine, on the CPU in double.** Each kernel function mirrors an engine function (same formula, same order of operations), so their doubles must be the same bits. 4096 cases, edges included (u = 0, the poles of the basis, n.z = -0, total internal reflection, cos θ < 0): bit-identical.
2. **Random numbers, GPU = CPU, bit for bit.** Integer arithmetic. 4096 streams: identical.
3. **Floats, GPU as accurate as CPU.** The largest |float − double| on the GPU must be at most twice that on the CPU. Measured on a Quadro M4000: GPU 3.9e-6, CPU 1.0e-5. The GPU is the more accurate one, because CUDA fuses multiply-adds (one rounding instead of two).
4. **Doubles, GPU = CPU to rounding** (tolerance 1e-12): largest difference 3.2e-15, 92 % of the values bit-identical.

A fixed tolerance on |GPU − CPU| in float is the wrong test. Just past the critical angle (glass to air, cos θ 0.7458 against a critical 0.7454) the Fresnel term takes the square root of a small number, which turns a rounding difference in its input into 1.4e-5 in its output. The same case is the largest difference in double too. That is the problem's conditioning, not a GPU error, hence check 3.

Planned image checks, each printing SKIP without a device:

- the same scene and seed on the CPU kernel and on the GPU, pixel by pixel (most paths take the same branches, so most pixels agree closely);
- the GPU against the CPU renderer within the Monte Carlo noise, by the method of `MITSUBA_CHECK.md` (two independent renders give the noise; mean ratio ± standard error and the share of 16 × 16 blocks beyond 3 standard errors) and Mitsuba's per-pixel z-test;
- white furnace scenes on the GPU.

## 3. Build and run

```
powershell -ExecutionPolicy Bypass -File tools\build_cli.ps1 -Dev      # the CLI (no CUDA needed)
powershell -ExecutionPolicy Bypass -File tools\build_gpu.ps1           # GoLD_rayt_gpu.dll next to it (8 s)
cd GoLD_rayt
..\out\cli-dev\GoLD_rayt.exe --scene gpu
```

`build_gpu.ps1 -Arch 75,86,89` compiles machine code for other GPUs. PTX (the GPU's intermediate code) of the lowest architecture is always embedded, so a newer GPU runs the module too: its driver compiles the PTX on the first launch. The design is not tied to the current card. Recent Mitsuba already requires compute capability 7.5 (Turing) or newer, and CUDA 13 no longer compiles for Maxwell, so a Maxwell card like the M4000 needs CUDA 12.x.

## 4. Next stages

1. **A minimal path tracer in the kernel:** triangle meshes, diffuse surfaces, area lights and a uniform environment, RGB. It runs on the CPU (OpenMP, Embree in float) and on the GPU (OptiX megakernel), with the image checks above.
2. **Materials one at a time** (rough conductor, dielectric, coated diffuse, ...): the engine class calls the kernel function, the engine's renders stay the same, and eval / pdf / sample agree between the devices on fixed inputs (plus χ² of sample against pdf on the GPU).
3. Spectral mode, textures, timing against the CPU renderer; a wavefront integrator only if the timings call for it.
