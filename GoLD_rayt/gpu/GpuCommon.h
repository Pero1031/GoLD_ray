#pragma once
// Helpers of the GPU module's host code (docs/GPU.md): error messages and device buffers that free themselves.

#include <cuda_runtime.h>

#include <stdio.h>

#define RAYT_GPU_EXPORT extern "C" __declspec(dllexport)

namespace rayt::gpu {

    inline void writeError(char* error, int errorSize, const char* what, cudaError_t e) {
        if (error && errorSize > 0) snprintf(error, size_t(errorSize), "%s: %s (%s)", what, cudaGetErrorString(e), cudaGetErrorName(e));
    }

    /// Device memory that is freed when it goes out of scope, also on an early error return.
    template <class T> struct DeviceBuffer {
        T* ptr = nullptr;
        cudaError_t allocate(size_t count) { return cudaMalloc(reinterpret_cast<void**>(&ptr), (count ? count : 1) * sizeof(T)); }
        unsigned long long device() const { return reinterpret_cast<unsigned long long>(ptr); }
        DeviceBuffer() = default;
        DeviceBuffer(const DeviceBuffer&) = delete;
        DeviceBuffer& operator=(const DeviceBuffer&) = delete;
        ~DeviceBuffer() { if (ptr) cudaFree(ptr); }
    };

    /// A copy of `count` elements on the GPU.
    template <class T> cudaError_t upload(DeviceBuffer<T>& b, const T* src, size_t count) {
        cudaError_t e = b.allocate(count);
        if (e == cudaSuccess && count) e = cudaMemcpy(b.ptr, src, count * sizeof(T), cudaMemcpyHostToDevice);
        return e;
    }

} // namespace rayt::gpu

#define RAYT_CUDA(call, what) do { const cudaError_t e_ = (call); if (e_ != cudaSuccess) { rayt::gpu::writeError(error, errorSize, what, e_); return 1; } } while (0)
