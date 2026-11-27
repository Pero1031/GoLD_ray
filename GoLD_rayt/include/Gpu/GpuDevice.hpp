#pragma once
/**
 * @file Gpu/GpuDevice.hpp
 * @brief Loads the GPU module GoLD_rayt_gpu.dll at run time (docs/GPU.md).
 *
 * The module is looked for next to the program. Without it, or without a usable CUDA device, the GPU parts are simply
 * not available: the checks print SKIP and the renderer stays on the CPU.
 */

#include <string>

#include "Gpu/GpuApi.h"

namespace rayt::gpu {

    /// The functions of a loaded module.
    struct Module {
        RaytGpuApiVersionFn apiVersion = nullptr;
        RaytGpuQueryFn query = nullptr;
        RaytGpuMathTestFn mathTest = nullptr;
        RaytGpuMaterialTestFn materialTest = nullptr;
        RaytGpuRenderFn render = nullptr;
        RaytGpuTraceTestFn traceTest = nullptr;
        std::string path;   ///< the file that was loaded
    };

    /// The module (loaded on the first call, kept for the run), or null with the reason in `why`: the file is missing,
    /// cannot be loaded, lacks a function or was built for another RAYT_GPU_API_VERSION.
    const Module* LoadModule(std::string& why);

} // namespace rayt::gpu
