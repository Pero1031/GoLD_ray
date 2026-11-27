#pragma once
/**
 * @file Kernel/SelfTest.hpp
 * @brief One function that --scene gpu runs on the CPU and on the GPU with the same inputs (docs/GPU.md).
 *
 * evalMathCase draws random numbers from the case's stream and evaluates the kernel's sampling functions in float. The
 * CPU calls it directly and the GPU module calls it from a CUDA kernel; evalMathCaseDouble (CPU only, the GPU computes in
 * float) gives the reference. The check compares: the random numbers must be the same bits, and the GPU's float values
 * must be as close to the double reference as the CPU's.
 */

#include <stdint.h>

#include "Kernel/Bsdf.hpp"
#include "Kernel/Environment.hpp"
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

    /// The values of one case (the same list in float and in double).
    enum MathValue : int {
        kRngFloat0, kRngFloat1,                    ///< the stream's 5th and 6th draws as floats
        kHemiX, kHemiY, kHemiZ,                    ///< sampleCosineHemisphere(u1, u2)
        kSphereX, kSphereY, kSphereZ,              ///< sampleUniformSphere(u1, u2)
        kBasisSX, kBasisSY, kBasisSZ,              ///< orthonormalBasis(n): s
        kBasisTX, kBasisTY, kBasisTZ,              ///< orthonormalBasis(n): t
        kFresnel,                                  ///< fresnelDielectric(cosTheta, etaI, etaT)
        kAtan2, kAcos,                             ///< katan2(ny, nx), kacos(nz) (the float versions of Kernel/KernelMath.hpp)
        kExp,                                      ///< kexp(-12 u1)
        kMathValueCount
    };

    /// What one device computes for a case.
    struct MathOut {
        uint32_t rng[4];                           ///< the stream's first four 32-bit draws
        float f[kMathValueCount];                  ///< the values in float
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
        v[kAtan2] = katan2(T(c.ny), T(c.nx));
        v[kAcos] = kacos(T(c.nz));
        v[kExp] = kexp(T(-12) * T(c.u1));
    }

    /// The random numbers and the values in float: what the CPU and the GPU both run.
    RAYT_HD_INLINE void evalMathCase(const MathCase& c, MathOut& out) {
        Pcg32 rng = Pcg32::forSample(c.seed, c.pixel, c.sample);
        for (int k = 0; k < 4; ++k) out.rng[k] = rng.nextUint();
        evalMathValues<float>(c, rng, out.f);
    }

#if !defined(RAYT_GPU_CODE)
    /// The same values in double: the reference, computed on the CPU only.
    inline void evalMathCaseDouble(const MathCase& c, double* v) {
        Pcg32 rng = Pcg32::forSample(c.seed, c.pixel, c.sample);
        for (int k = 0; k < 4; ++k) rng.nextUint();
        evalMathValues<double>(c, rng, v);
    }
#endif

    // ---- materials: eval, pdf and sample of one KMaterial on given directions (--scene gpu compares the CPU's and the GPU's) ----

    struct MaterialCase {
        KFloat3 wo, wi, ns, ng;   ///< world space; ns and ng turned toward wo as at a hit
        float u1, u2, u3;         ///< the sample's random numbers
        uint32_t front;           ///< the ray came from the outside
        uint32_t stream;          ///< the random number stream of a layered material's walks
    };

    struct MaterialOut {
        float f[3], pdf;          ///< bsdfEval, bsdfPdf at (wo, wi)
        uint32_t sampled, delta;  ///< bsdfSample succeeded; its event is specular
        float wi[3], sf[3], spdf; ///< its direction, value and density
    };

    RAYT_HD_INLINE void evalMaterialCase(const KMaterial& m, const MaterialCase& c, MaterialOut& out) {
        const Vec3<float> wo(c.wo.x, c.wo.y, c.wo.z), wi(c.wi.x, c.wi.y, c.wi.z), ns(c.ns.x, c.ns.y, c.ns.z), ng(c.ng.x, c.ng.y, c.ng.z);
        const bool front = c.front != 0;
        Pcg32 rng = Pcg32::forSample(0x1a7e5eedull, c.stream, 0);   // a layered material's walks
        const Vec3<float> f = bsdfEval(m, front, wo, wi, ns, ng, rng);
        out.f[0] = f.x; out.f[1] = f.y; out.f[2] = f.z;
        out.pdf = bsdfPdf(m, front, wo, wi, ns, ng);
        BsdfSample<float> s{};
        out.sampled = bsdfSample(m, front, wo, ns, ng, c.u1, c.u2, c.u3, rng, s) ? 1u : 0u;
        out.delta = s.delta ? 1u : 0u;
        out.wi[0] = s.wi.x; out.wi[1] = s.wi.y; out.wi[2] = s.wi.z;
        out.sf[0] = s.f.x; out.sf[1] = s.f.y; out.sf[2] = s.f.z;
        out.spdf = s.pdf;
    }

} // namespace rayt::kernel
