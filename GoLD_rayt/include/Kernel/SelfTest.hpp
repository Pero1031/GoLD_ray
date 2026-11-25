#pragma once
/**
 * @file Kernel/SelfTest.hpp
 * @brief One function that --scene gpu runs on the CPU and on the GPU with the same inputs (docs/GPU.md).
 *
 * evalMathCase draws random numbers from the case's stream and evaluates the kernel's sampling functions, in float and
 * in double. The CPU calls it directly, the GPU module calls it from a CUDA kernel, and the check compares the results:
 * the random numbers must be the same bits, the values agree to rounding.
 */

#include <stdint.h>

#include "Kernel/Rng.hpp"
#include "Kernel/Sampling.hpp"

namespace rayt::kernel {

    /// Inputs of one case: plain floats and integers, the same layout under MSVC and nvcc.
    struct MathCase {
        uint64_t seed;
        uint32_t pixel, sample;
        float u1, u2;              ///< the sampling functions' random numbers
        float nx, ny, nz;          ///< a normal (normalized in the precision of the evaluation)
        float cosTheta, etaI, etaT;
    };

    /// The values of one case: the same list in float and in double.
    enum MathValue : int {
        kRngFloat0, kRngFloat1,                    ///< the stream's 5th and 6th draws as floats
        kHemiX, kHemiY, kHemiZ,                    ///< sampleCosineHemisphere(u1, u2)
        kSphereX, kSphereY, kSphereZ,              ///< sampleUniformSphere(u1, u2)
        kBasisSX, kBasisSY, kBasisSZ,              ///< orthonormalBasis(n): s
        kBasisTX, kBasisTY, kBasisTZ,              ///< orthonormalBasis(n): t
        kFresnel,                                  ///< fresnelDielectric(cosTheta, etaI, etaT)
        kMathValueCount
    };

    struct MathOut {
        uint32_t rng[4];                           ///< the stream's first four 32-bit draws
        float f[kMathValueCount];
        double d[kMathValueCount];
    };

    template <class T> RAYT_HD_INLINE void evalMathValues(const MathCase& c, Pcg32 rng, T* v) {
        v[kRngFloat0] = T(rng.nextFloat());
        v[kRngFloat1] = T(rng.nextFloat());
        const T u1 = T(c.u1), u2 = T(c.u2);
        const Vec3<T> h = sampleCosineHemisphere(u1, u2);
        v[kHemiX] = h.x; v[kHemiY] = h.y; v[kHemiZ] = h.z;
        const Vec3<T> sp = sampleUniformSphere(u1, u2);
        v[kSphereX] = sp.x; v[kSphereY] = sp.y; v[kSphereZ] = sp.z;
        Vec3<T> n(T(c.nx), T(c.ny), T(c.nz));
        n = n / ksqrt(n.x * n.x + n.y * n.y + n.z * n.z);
        Vec3<T> s, t;
        orthonormalBasis(n, s, t);
        v[kBasisSX] = s.x; v[kBasisSY] = s.y; v[kBasisSZ] = s.z;
        v[kBasisTX] = t.x; v[kBasisTY] = t.y; v[kBasisTZ] = t.z;
        v[kFresnel] = fresnelDielectric(T(c.cosTheta), T(c.etaI), T(c.etaT));
    }

    RAYT_HD_INLINE void evalMathCase(const MathCase& c, MathOut& out) {
        Pcg32 rng = Pcg32::forSample(c.seed, c.pixel, c.sample);
        for (int k = 0; k < 4; ++k) out.rng[k] = rng.nextUint();
        evalMathValues<float>(c, rng, out.f);
        evalMathValues<double>(c, rng, out.d);
    }

} // namespace rayt::kernel
