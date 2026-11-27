#pragma once
/**
 * @file Kernel/KernelScene.hpp
 * @brief The scene as the kernel renderer reads it: flat arrays of plain numbers (docs/GPU.md).
 *
 * The CPU builds it from a scene file (Cli/KernelScenePack); the CPU kernel renderer reads it in place, and the GPU module
 * copies the arrays to the GPU and points a copy of this struct at them. Everything is float (the GPU computes in float);
 * KFloat3 is three floats with no padding, the same 12 bytes under MSVC and nvcc (glm's vec3 may be padded on the CPU).
 */

#include <stdint.h>

#include "Kernel/HostDevice.hpp"

namespace rayt::kernel {

    struct KFloat3 { float x, y, z; };

    constexpr uint32_t kNoHit = 0xffffffffu;

    /// A ray's closest hit as a tracer reports it: distance, barycentrics of vertices 1 and 2, triangle index (kNoHit: none).
    struct KHit { float t, b1, b2; uint32_t prim; };

    /// A ray for the tracer checks (--scene gpu): origin, direction, interval.
    struct KRay { KFloat3 o, d; float tmin, tmax; };

    enum : uint32_t { kMatDiffuse = 0, kMatEmission = 1, kMatConductor = 2, kMatDielectric = 3 };

    /// Rows of the GGX directional-albedo table (Microfacet/GGXAlbedo.hpp, GGXAlbedoTable::kMu).
    constexpr int kAlbedoMu = 32;

    /**
     * @brief A rough conductor as the engine's RoughConductor holds it (isotropic GGX, RGB): complex index per channel, the GGX
     *        alpha (already clamped as GGXDistribution clamps it) and the Kulla-Conty multiple-scattering lobe.
     *
     * T = float in a scene; T = double where the checks compare the kernel's formulas with the engine's.
     */
    template <class T> struct ConductorParams {
        T eta[3], k[3];
        T alpha;
        uint32_t multipleScattering;   ///< 1: add the lobe f_ms = F_ms (1 - E(mu_o)) (1 - E(mu_i)) / (pi (1 - E_avg))
        T eAvg;                        ///< E_avg, the cosine-weighted mean of E
        T fms[3];                      ///< F_ms per channel
        T eRow[kAlbedoMu];             ///< E(mu) at mu = i / (kAlbedoMu - 1) for this alpha
    };

    /**
     * @brief A dielectric interface as the engine's Dielectric holds it (isotropic, no dispersion, nothing inside): the relative index
     *        (inside / outside), whether it is smooth (a delta lobe: alpha below 0.001), and the GGX alpha of a rough one (clamped as
     *        GGXDistribution clamps it).
     */
    template <class T> struct DielectricParams {
        T ior;
        T alpha;
        uint32_t smooth;
    };

    struct KMaterial {
        uint32_t type;             ///< kMatDiffuse | kMatEmission | kMatConductor | kMatDielectric
        uint32_t twoSided;         ///< emission: emits from both faces
        float albedo[3];           ///< diffuse: reflectance (linear RGB)
        float emission[3];         ///< emission: radiance (linear RGB)
        ConductorParams<float> conductor;   ///< conductor
        DielectricParams<float> dielectric; ///< dielectric
    };

    enum : uint32_t { kLightQuad = 0, kLightSphere = 1 };

    /// A light that next-event estimation samples: an emitting parallelogram (corner + a edge1 + b edge2, a uniform point on it)
    /// or an emitting sphere (the cone it subtends, as the engine's Sphere::sampleFrom).
    struct KLight {
        uint32_t type;             ///< kLightQuad | kLightSphere
        KFloat3 corner, edge1, edge2;
        KFloat3 normal;            ///< quad: unit cross(edge1, edge2), the front face
        KFloat3 center;            ///< sphere
        float radius;              ///< sphere (as the engine's: a negative radius turns the normals inward)
        float area;
        float radiance[3];
        uint32_t twoSided;
    };

    /// A sphere of the scene (an analytic shape, as the engine's Sphere; the tracers intersect it with kernel::intersectSphere).
    struct KSphere {
        KFloat3 center;
        float radius;              ///< negative: the normals point inward (Geometry/Sphere.hpp)
        uint32_t material;
        int32_t light;             ///< index of its light, -1 if none
    };

    /// A pinhole: film point (s, t) in [0, 1)^2 looks from origin toward lowerLeft + s horizontal + t vertical (t = 0 at the bottom).
    struct KCamera { KFloat3 origin, lowerLeft, horizontal, vertical; };

    /// Everything a render reads. The pointers are valid on the device that runs the render.
    struct KScene {
        const KFloat3* positions;     ///< vertices, world space
        const KFloat3* normals;       ///< per vertex (read for triangles with triSmooth = 1)
        const uint32_t* indices;      ///< 3 per triangle
        const KFloat3* faceNormals;   ///< per triangle: unit geometric normal, turned toward the vertex normals as the engine does
        const uint32_t* triMaterial;  ///< per triangle: index into materials
        const int32_t* triLight;      ///< per triangle: index of the light it belongs to, -1 if none
        const uint8_t* triSmooth;     ///< per triangle: 1 = shade with the interpolated vertex normals
        const KMaterial* materials;
        const KSphere* spheres;       ///< primitive numTriangles + i of the tracers is sphere i
        const KLight* lights;
        const float* lightCdf;        ///< numLights + 1 entries: cdf[0] = 0, cdf[numLights] = 1 (lights chosen by emitted power)
        uint32_t numVertices, numTriangles, numSpheres, numMaterials, numLights;
        // The environment: an equirectangular map as the engine's EnvMap holds it (IO/EnvMap.hpp; a uniform environment is a white
        // map), times envScale, and its importance-sampling distribution (Core/Distribution2D.hpp: luminance x sin(theta) per texel).
        const float* envPixels;       ///< envWidth x envHeight x 3, row 0 = the top (the +y pole)
        const float* envMarginalFunc; ///< envHeight: the integral of each row
        const float* envMarginalCdf;  ///< envHeight + 1
        const float* envCondFunc;     ///< envHeight x envWidth: the weights of each row
        const float* envCondCdf;      ///< envHeight x (envWidth + 1)
        const float* envCondInt;      ///< envHeight: the integral of each row's weights
        float envMarginalInt;
        uint32_t envWidth, envHeight; ///< 0 x 0: no environment
        float envScale;
        float envSelect;              ///< the probability that a light sample picks the environment (the other lights share the rest)
        KCamera camera;
        uint32_t width, height;       ///< pixels
        uint32_t spp, maxDepth;       ///< samples per pixel; path length as the engine counts it (surface vertices)
        uint64_t seed;                ///< the run seed of the random number streams
    };

} // namespace rayt::kernel
