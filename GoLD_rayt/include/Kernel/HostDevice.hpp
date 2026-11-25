#pragma once
/**
 * @file Kernel/HostDevice.hpp
 * @brief Marks code that runs on the CPU and on the GPU (docs/GPU.md).
 *
 * The headers under include/Kernel are "kernel" code: MSVC compiles them into the CPU program and nvcc compiles the
 * same text into the GPU module (gpu/*.cu). Every function there is marked RAYT_HD. Code in these headers must not use
 * virtual functions, exceptions, the heap, std containers, std::optional or thread_local, and takes its random numbers
 * as arguments or from a kernel::Pcg32 that the caller owns. The scalar type is a template parameter T (float on the
 * GPU; float or double on the CPU), so the same function can be compared across devices and precisions.
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
