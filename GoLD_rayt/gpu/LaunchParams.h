#pragma once
// The parameters of one OptiX launch of the kernel renderer (docs/GPU.md): the scene with its arrays on the GPU, the geometry,
// the accumulation buffer and the sample index of this launch. Shared by the host code (GpuModule.cu) and the programs
// (OptixPrograms.cu).

#include <optix_types.h>

#include "Kernel/KernelScene.hpp"

namespace rayt::gpu {

    struct LaunchParams {
        kernel::KScene scene;              ///< pointers on the GPU
        OptixTraversableHandle handle;     ///< the triangles' GAS, or an instance AS over it and the spheres' GAS
        float* accum;                      ///< width x height x 3: the sum of the samples so far
        float* accumAlbedo;                ///< the same for the denoiser's features (null: not wanted)
        float* accumNormal;
        uint32_t sample;                   ///< the sample index every pixel renders in this launch
        const kernel::KRay* rays;          ///< the tracer check: rays to trace (one per launch index) ...
        kernel::KHit* hits;                ///< ... and their closest hits
    };

} // namespace rayt::gpu
