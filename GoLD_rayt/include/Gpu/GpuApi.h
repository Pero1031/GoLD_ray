#pragma once
/**
 * @file Gpu/GpuApi.h
 * @brief The functions of the GPU module GoLD_rayt_gpu.dll (docs/GPU.md).
 *
 * The CPU program never links against CUDA: it loads the module at run time (Gpu/GpuDevice.hpp), so it starts on a
 * computer without an NVIDIA GPU or driver, and the GPU parts report "not available" there. The module is built by
 * tools/build_gpu.ps1 from gpu/*.cu. Both sides include this header; the structs hold plain numbers only.
 * Raise RAYT_GPU_API_VERSION whenever a struct or a function here changes: a module of another version is refused.
 */

#include <stdint.h>

#include "Kernel/SelfTest.hpp"

#define RAYT_GPU_API_VERSION 2
#define RAYT_GPU_MODULE_NAME "GoLD_rayt_gpu.dll"

extern "C" {

    struct RaytGpuInfo {
        int deviceCount;           ///< CUDA devices found (the module uses device 0)
        int computeMajor, computeMinor;
        int multiprocessors;
        uint64_t totalMemory;      ///< bytes
        int driverVersion;         ///< CUDA version the driver supports, e.g. 12090
        int runtimeVersion;        ///< CUDA version the module was built with
        int fastMath;              ///< 1 if the module was built with --use_fast_math (tools/build_gpu.ps1 -FastMath)
        char name[256];
        char error[256];           ///< why no device can be used (empty on success)
    };

    /// RAYT_GPU_API_VERSION of the module.
    typedef int (*RaytGpuApiVersionFn)(void);
    /// Fills `info`; returns 0 when device 0 can be used.
    typedef int (*RaytGpuQueryFn)(RaytGpuInfo* info);
    /// Runs kernel::evalMathCase (float) on the GPU for n cases; returns 0 on success, else writes the reason to `error`.
    typedef int (*RaytGpuMathTestFn)(const rayt::kernel::MathCase* cases, rayt::kernel::MathOut* out, int n, char* error, int errorSize);

}
