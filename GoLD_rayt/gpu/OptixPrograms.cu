// The OptiX programs of the kernel renderer (docs/GPU.md), compiled to PTX by tools/build_gpu.ps1 and embedded in the GPU module.
//
// A megakernel, as Mitsuba's GPU variants: the ray-generation program runs a whole camera path (kernel::tracePath, the same code
// the CPU runs) and casts its rays with optixTrace; the closest-hit programs only report the hit (distance, barycentrics,
// primitive), and shading happens back in the ray-generation program. Spheres are custom primitives whose intersection program
// is kernel::intersectSphere, the function Embree calls on the CPU (as Mitsuba binds one routine to both tracers). One launch renders one sample of every pixel and adds it
// to the accumulation buffer, so launches stay short (no display watchdog) and the sums are added in the same order as on the CPU.

#include <optix.h>

#include "LaunchParams.h"
#include "Kernel/PathTracer.hpp"

extern "C" {
    __constant__ rayt::gpu::LaunchParams params;
}

namespace {

    struct OptixTracer {
        __device__ rayt::kernel::KHit closest(const rayt::kernel::V3& o, const rayt::kernel::V3& d, float tmin, float tmax) const {
            unsigned int t = 0, b1 = 0, b2 = 0, prim = rayt::kernel::kNoHit;
            optixTrace(params.handle, make_float3(o.x, o.y, o.z), make_float3(d.x, d.y, d.z), tmin, tmax, 0.0f, OptixVisibilityMask(255),
                OPTIX_RAY_FLAG_DISABLE_ANYHIT, 0, 1, 0, t, b1, b2, prim);
            return { __uint_as_float(t), __uint_as_float(b1), __uint_as_float(b2), prim };
        }
        __device__ bool occluded(const rayt::kernel::V3& o, const rayt::kernel::V3& d, float tmin, float tmax) const {
            unsigned int blocked = 1;   // the shadow miss program clears it
            optixTrace(params.handle, make_float3(o.x, o.y, o.z), make_float3(d.x, d.y, d.z), tmin, tmax, 0.0f, OptixVisibilityMask(255),
                OPTIX_RAY_FLAG_DISABLE_ANYHIT | OPTIX_RAY_FLAG_TERMINATE_ON_FIRST_HIT | OPTIX_RAY_FLAG_DISABLE_CLOSESTHIT, 0, 1, 1, blocked);
            return blocked != 0;
        }
    };

} // namespace

extern "C" __global__ void __raygen__render() {
    const uint3 idx = optixGetLaunchIndex();
    const rayt::kernel::V3 L = rayt::kernel::tracePath(params.scene, OptixTracer{}, idx.x, idx.y, params.sample);
    float* a = params.accum + 3 * (size_t(idx.y) * params.scene.width + idx.x);
    a[0] += L.x;
    a[1] += L.y;
    a[2] += L.z;
}

/// The tracer check: the closest hit of one given ray (no shading), to compare OptiX's answers with Embree's.
extern "C" __global__ void __raygen__trace() {
    const uint32_t i = optixGetLaunchIndex().x;
    const rayt::kernel::KRay r = params.rays[i];
    params.hits[i] = OptixTracer{}.closest(rayt::kernel::V3(r.o.x, r.o.y, r.o.z), rayt::kernel::V3(r.d.x, r.d.y, r.d.z), r.tmin, r.tmax);
}

extern "C" __global__ void __miss__radiance() { optixSetPayload_3(rayt::kernel::kNoHit); }

extern "C" __global__ void __miss__shadow() { optixSetPayload_0(0); }

/// A sphere of the scene (custom primitive i of the sphere GAS is scene sphere i): the kernel's test on the ray as it was cast
/// (the instance transform is the identity).
extern "C" __global__ void __intersection__sphere() {
    const rayt::kernel::KSphere& sp = params.scene.spheres[optixGetPrimitiveIndex()];
    const float3 o = optixGetWorldRayOrigin(), d = optixGetWorldRayDirection();
    float t;
    if (rayt::kernel::intersectSphere(rayt::kernel::V3(o.x, o.y, o.z), rayt::kernel::V3(d.x, d.y, d.z), rayt::kernel::V3(sp.center.x, sp.center.y, sp.center.z),
            sp.radius, optixGetRayTmin(), optixGetRayTmax(), t))
        optixReportIntersection(t, 0);
}

extern "C" __global__ void __closesthit__sphere() {
    optixSetPayload_0(__float_as_uint(optixGetRayTmax()));
    optixSetPayload_1(0);
    optixSetPayload_2(0);
    optixSetPayload_3(params.scene.numTriangles + optixGetPrimitiveIndex());   // as the kernel numbers primitives
}

extern "C" __global__ void __closesthit__radiance() {
    const float2 b = optixGetTriangleBarycentrics();
    optixSetPayload_0(__float_as_uint(optixGetRayTmax()));
    optixSetPayload_1(__float_as_uint(b.x));
    optixSetPayload_2(__float_as_uint(b.y));
    optixSetPayload_3(optixGetPrimitiveIndex());
}
