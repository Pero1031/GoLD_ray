#pragma once
/**
 * @file Kernel/Rng.hpp
 * @brief Random numbers that are the same on every device: PCG32 streams keyed by (run seed, pixel, sample).
 *
 * The engine's Random() draws from a thread_local Mersenne Twister, which a GPU thread cannot have, and its numbers
 * depend on which thread renders which row. Kernel code instead owns one Pcg32 per camera path, created from the run
 * seed, the pixel and the sample index alone, so the CPU and the GPU draw the same numbers for the same sample (the
 * integer arithmetic is exact on both). Two paths only drift apart where the devices' floating point makes a path take
 * another branch (a hit against a miss, a Russian roulette decision).
 */

#include <stdint.h>

#include "Kernel/HostDevice.hpp"

namespace rayt::kernel {

    /// The finalizer of SplitMix64: neighbouring inputs give unrelated outputs.
    RAYT_HD_INLINE uint64_t mix64(uint64_t x) {
        x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ull;
        x ^= x >> 27; x *= 0x94d049bb133111ebull;
        x ^= x >> 31;
        return x;
    }

    /// PCG32 (O'Neill 2014, the XSH RR output function): 64 bits of state, 32-bit outputs, period 2^64 per stream.
    struct Pcg32 {
        uint64_t state = 0x853c49e6748fea9bull;
        uint64_t inc = 0xda3e39cb94b95bdbull;

        RAYT_HD_INLINE void seed(uint64_t initState, uint64_t initSeq) {
            state = 0;
            inc = (initSeq << 1u) | 1u;
            nextUint();
            state += initState;
            nextUint();
        }

        /// The stream of one camera path: a function of the run seed, the pixel index and the sample index only.
        RAYT_HD_INLINE static Pcg32 forSample(uint64_t runSeed, uint32_t pixel, uint32_t sample) {
            Pcg32 r;
            r.seed(mix64(runSeed ^ mix64((uint64_t(pixel) << 32) | sample)), mix64(runSeed + 0x9e3779b97f4a7c15ull));
            return r;
        }

        RAYT_HD_INLINE uint32_t nextUint() {
            const uint64_t old = state;
            state = old * 6364136223846793005ull + inc;
            const uint32_t xorshifted = uint32_t(((old >> 18u) ^ old) >> 27u);
            const uint32_t rot = uint32_t(old >> 59u);
            return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31u));
        }

        /// Uniform in [0, 1): the top 24 bits, exact in float.
        RAYT_HD_INLINE float nextFloat() { return float(nextUint() >> 8) * 0x1p-24f; }

#if !defined(RAYT_GPU_CODE)
        /// Uniform in [0, 1) with 53 bits (two draws). CPU only.
        double nextDouble() {
            const uint64_t a = nextUint() >> 5, b = nextUint() >> 6;
            return double(a * 67108864ull + b) * 0x1p-53;
        }
#endif
    };

    /// A uniform number in [0, 1) of the scalar type T (float: one draw; double, CPU only: two).
    template <class T> RAYT_HD_INLINE T uniform(Pcg32& rng);
    template <> RAYT_HD_INLINE float uniform<float>(Pcg32& rng) { return rng.nextFloat(); }
#if !defined(RAYT_GPU_CODE)
    template <> inline double uniform<double>(Pcg32& rng) { return rng.nextDouble(); }
#endif

} // namespace rayt::kernel
