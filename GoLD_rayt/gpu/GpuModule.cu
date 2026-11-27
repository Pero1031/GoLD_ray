// The GPU module GoLD_rayt_gpu.dll (docs/GPU.md): the CUDA part (device query, the math test kernel of --scene gpu) and the C
// functions of Gpu/GpuApi.h; the OptiX renderer is in OptixRender.cpp. Built by tools/build_gpu.ps1; the CPU program loads the
// module at run time. The device code computes in float only: tools/build_gpu.ps1 makes ptxas reject double instructions.

#include <cuda_runtime.h>

#include <stdio.h>
#include <string.h>

#include "Gpu/GpuApi.h"
#include "GpuCommon.h"

using rayt::gpu::DeviceBuffer;
using rayt::gpu::upload;
using rayt::gpu::writeError;

namespace {

    __global__ void mathTestKernel(const rayt::kernel::MathCase* cases, rayt::kernel::MathOut* out, int n) {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < n) rayt::kernel::evalMathCase(cases[i], out[i]);
    }

    __global__ void materialTestKernel(const rayt::kernel::KMaterial* material, const rayt::kernel::MaterialCase* cases, rayt::kernel::MaterialOut* out, int n) {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < n) rayt::kernel::evalMaterialCase(*material, cases[i], out[i]);
    }

} // namespace

RAYT_GPU_EXPORT int raytGpuMaterialTest(const rayt::kernel::KMaterial* material, const rayt::kernel::MaterialCase* cases, rayt::kernel::MaterialOut* out, int n,
    char* error, int errorSize) {
    if (n <= 0) return 0;
    DeviceBuffer<rayt::kernel::KMaterial> dMaterial;
    DeviceBuffer<rayt::kernel::MaterialCase> dCases;
    DeviceBuffer<rayt::kernel::MaterialOut> dOut;
    RAYT_CUDA(upload(dMaterial, material, 1), "copying the material");
    RAYT_CUDA(upload(dCases, cases, size_t(n)), "copying the cases");
    RAYT_CUDA(dOut.allocate(size_t(n)), "cudaMalloc");
    const int block = 128;
    materialTestKernel<<<(n + block - 1) / block, block>>>(dMaterial.ptr, dCases.ptr, dOut.ptr, n);
    RAYT_CUDA(cudaGetLastError(), "the material test kernel");
    RAYT_CUDA(cudaDeviceSynchronize(), "the material test kernel");
    RAYT_CUDA(cudaMemcpy(out, dOut.ptr, size_t(n) * sizeof(*out), cudaMemcpyDeviceToHost), "cudaMemcpy to the host");
    return 0;
}

RAYT_GPU_EXPORT int raytGpuApiVersion(void) { return RAYT_GPU_API_VERSION; }

RAYT_GPU_EXPORT int raytGpuQuery(RaytGpuInfo* info) {
    memset(info, 0, sizeof(*info));
    cudaRuntimeGetVersion(&info->runtimeVersion);
    cudaError_t e = cudaDriverGetVersion(&info->driverVersion);
    if (e == cudaSuccess && info->driverVersion == 0) {
        snprintf(info->error, sizeof(info->error), "no NVIDIA driver is installed");
        return 1;
    }
    e = cudaGetDeviceCount(&info->deviceCount);
    if (e != cudaSuccess) { writeError(info->error, int(sizeof(info->error)), "cudaGetDeviceCount", e); return 1; }
    if (info->deviceCount == 0) { snprintf(info->error, sizeof(info->error), "no CUDA device found"); return 1; }
    cudaDeviceProp prop;
    e = cudaGetDeviceProperties(&prop, 0);
    if (e != cudaSuccess) { writeError(info->error, int(sizeof(info->error)), "cudaGetDeviceProperties", e); return 1; }
    snprintf(info->name, sizeof(info->name), "%s", prop.name);
    info->computeMajor = prop.major;
    info->computeMinor = prop.minor;
    info->multiprocessors = prop.multiProcessorCount;
    info->totalMemory = uint64_t(prop.totalGlobalMem);
#if defined(RAYT_GPU_FAST_MATH)
    info->fastMath = 1;
#endif
    return 0;
}

RAYT_GPU_EXPORT int raytGpuMathTest(const rayt::kernel::MathCase* cases, rayt::kernel::MathOut* out, int n, char* error, int errorSize) {
    if (n <= 0) return 0;
    DeviceBuffer<rayt::kernel::MathCase> dCases;
    DeviceBuffer<rayt::kernel::MathOut> dOut;
    RAYT_CUDA(upload(dCases, cases, size_t(n)), "copying the cases");
    RAYT_CUDA(dOut.allocate(size_t(n)), "cudaMalloc");
    const int block = 128;
    mathTestKernel<<<(n + block - 1) / block, block>>>(dCases.ptr, dOut.ptr, n);
    RAYT_CUDA(cudaGetLastError(), "the math test kernel");
    RAYT_CUDA(cudaDeviceSynchronize(), "the math test kernel");
    RAYT_CUDA(cudaMemcpy(out, dOut.ptr, size_t(n) * sizeof(*out), cudaMemcpyDeviceToHost), "cudaMemcpy to the host");
    return 0;
}

