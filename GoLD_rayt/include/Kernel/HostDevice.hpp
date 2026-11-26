#pragma once
/**
 * @file Kernel/HostDevice.hpp
 * @brief Marks code that runs on the CPU and on the GPU (docs/GPU.md).
 *
 * The headers under include/Kernel are "kernel" code: MSVC compiles them into the CPU program and nvcc compiles the
 * same text into the GPU module (gpu/*.cu). Every function there is marked RAYT_HD. Code in these headers must not use
 * virtual functions, exceptions, the heap, std containers, std::optional or thread_local, and takes its random numbers
 * as arguments or from a kernel::Pcg32 that the caller owns.
 *
 * Precision: the GPU computes in float only; double is for the CPU. Kernel functions are templates on the scalar type T,
 * and everything that only makes sense for double (the double overloads of the math functions, double random numbers,
 * double references in the checks) sits in `#if !defined(RAYT_GPU_CODE)` sections: device code that reaches a double
 * fails to compile, and tools/build_gpu.ps1 makes ptxas reject any double instruction that slips through (an implicit
 * promotion such as `x * 0.5` with a float x). Write float constants as T(0.5), never as bare double literals.
 */

#if defined(__CUDACC__)
#define RAYT_HD __host__ __device__
#define RAYT_HD_INLINE __host__ __device__ __forceinline__
#else
#define RAYT_HD
#define RAYT_HD_INLINE inline
#endif

// Defined while nvcc compiles the device side of a .cu file.
#if defined(__CUDA_ARCH__)
#define RAYT_GPU_CODE 1
#endif
