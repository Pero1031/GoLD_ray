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

    enum : uint32_t { kMatDiffuse = 0, kMatEmission = 1 };

    struct KMaterial {
        uint32_t type;             ///< kMatDiffuse | kMatEmission
        uint32_t twoSided;         ///< emission: emits from both faces
        float albedo[3];           ///< diffuse: reflectance (linear RGB)
        float emission[3];         ///< emission: radiance (linear RGB)
    };

    /// An emitting parallelogram that next-event estimation samples: corner + a edge1 + b edge2, a, b in [0, 1).
    struct KQuadLight {
        KFloat3 corner, edge1, edge2;
        KFloat3 normal;            ///< unit cross(edge1, edge2): the front face
        float area;
        float radiance[3];
        uint32_t twoSided;
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
        const int32_t* triLight;      ///< per triangle: index of the quad light it belongs to, -1 if none
        const uint8_t* triSmooth;     ///< per triangle: 1 = shade with the interpolated vertex normals
        const KMaterial* materials;
        const KQuadLight* lights;
        const float* lightCdf;        ///< numLights + 1 entries: cdf[0] = 0, cdf[numLights] = 1 (lights chosen by emitted power)
        uint32_t numVertices, numTriangles, numMaterials, numLights;
        float envRadiance[3];         ///< a uniform environment (0 = none)
        KCamera camera;
        uint32_t width, height;       ///< pixels
        uint32_t spp, maxDepth;       ///< samples per pixel; path length as the engine counts it (surface vertices)
        uint64_t seed;                ///< the run seed of the random number streams
    };

} // namespace rayt::kernel
