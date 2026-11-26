// The GPU module GoLD_rayt_gpu.dll (docs/GPU.md): CUDA device code for the kernel headers (include/Kernel) and the
// C functions of Gpu/GpuApi.h. Built by tools/build_gpu.ps1; the CPU program loads it at run time.
// The device code computes in float only: tools/build_gpu.ps1 makes ptxas reject double instructions.

#include <cuda_runtime.h>

#include <stdio.h>
#include <string.h>

#include "Gpu/GpuApi.h"

#define RAYT_GPU_EXPORT extern "C" __declspec(dllexport)

namespace {

    void writeError(char* error, int errorSize, const char* what, cudaError_t e) {
        if (error && errorSize > 0) snprintf(error, size_t(errorSize), "%s: %s (%s)", what, cudaGetErrorString(e), cudaGetErrorName(e));
    }

    __global__ void mathTestKernel(const rayt::kernel::MathCase* cases, rayt::kernel::MathOut* out, int n) {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < n) rayt::kernel::evalMathCase(cases[i], out[i]);
    }

    /// Frees device memory when the function returns, also on an early error return.
    template <class T> struct DeviceBuffer {
        T* ptr = nullptr;
        cudaError_t allocate(size_t count) { return cudaMalloc(reinterpret_cast<void**>(&ptr), count * sizeof(T)); }
        ~DeviceBuffer() { if (ptr) cudaFree(ptr); }
    };

} // namespace

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
    cudaError_t e = dCases.allocate(size_t(n));
    if (e == cudaSuccess) e = dOut.allocate(size_t(n));
    if (e != cudaSuccess) { writeError(error, errorSize, "cudaMalloc", e); return 1; }
    e = cudaMemcpy(dCases.ptr, cases, size_t(n) * sizeof(*cases), cudaMemcpyHostToDevice);
    if (e != cudaSuccess) { writeError(error, errorSize, "cudaMemcpy to the device", e); return 1; }
    const int block = 128;
    mathTestKernel<<<(n + block - 1) / block, block>>>(dCases.ptr, dOut.ptr, n);
    e = cudaGetLastError();
    if (e == cudaSuccess) e = cudaDeviceSynchronize();
    if (e != cudaSuccess) { writeError(error, errorSize, "the math test kernel", e); return 1; }
    e = cudaMemcpy(out, dOut.ptr, size_t(n) * sizeof(*out), cudaMemcpyDeviceToHost);
    if (e != cudaSuccess) { writeError(error, errorSize, "cudaMemcpy to the host", e); return 1; }
    return 0;
}
