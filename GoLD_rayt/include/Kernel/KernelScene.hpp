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

    enum : uint32_t { kMatDiffuse = 0, kMatEmission = 1, kMatConductor = 2, kMatDielectric = 3, kMatLayered = 4 };

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
        uint32_t nullInterface;    ///< Dielectric::isNullInterface: smooth, the same index on both sides (--media-nee crosses it)
    };

    /**
     * @brief A clear coat over a base as the engine's LayeredCoat holds it (RGB): the coat's index and GGX alpha (the perceptual
     *        roughness clamped to 0.02, squared), its optical depth along the normal per channel (-log of the tint), the random
     *        walks per evaluation and their length, and the kind of base (its parameters are KMaterial::albedo or ::conductor).
     */
    template <class T> struct LayeredParams {
        T ior;
        T alpha;
        T sigma[3];
        uint32_t samples, maxDepth;
        uint32_t conductorBase;    ///< 0: a Lambertian base (KMaterial::albedo), 1: a rough conductor (KMaterial::conductor)
        uint32_t absorbing;        ///< 1 if a sigma is not 0 (a tinted coat)
    };

    enum : uint32_t { kPhaseHG = 0, kPhaseTable = 1 };

    /**
     * @brief A scattering medium inside a dielectric (Materials/Medium.hpp: HomogeneousMedium or ChromaticMedium), RGB: the extinction
     *        per channel (1/m; for a density field the majorant, where q = 1), the scattering coefficient sigma_t albedo, the phase
     *        function (two Henyey-Greenstein lobes, or a table such as Mie scattering in KScene::phaseF / phaseCdf), and a density field
     *        (NoiseDensity: fractal value noise mapped by a smooth step, optionally faded out toward a sphere).
     */
    template <class T> struct MediumParams {
        T sigmaT[3], sigmaS[3];
        uint32_t constant;         ///< the same coefficients in every channel: no hero channel, no weights
        uint32_t phaseKind;        ///< kPhaseHG | kPhaseTable
        T g1, g2, w;               ///< p = w HG(g1) + (1 - w) HG(g2)
        uint32_t tableOffset, tableCount;
        T tableStep;               ///< the table's nodes are tableStep radians apart, from 0 to pi
        uint32_t noise;            ///< 1: a density field
        T frequency;
        uint32_t octaves;
        T lo, hi;
        uint32_t seed;
        T maskCenter[3];
        T maskRadius, maskStart;   ///< maskRadius 0: no mask
    };
    using KMedium = MediumParams<float>;

    struct KMaterial {
        uint32_t type;             ///< kMatDiffuse | kMatEmission | kMatConductor | kMatDielectric | kMatLayered
        uint32_t twoSided;         ///< emission: emits from both faces
        float albedo[3];           ///< diffuse, and the diffuse base of a layered material: reflectance (linear RGB)
        float emission[3];         ///< emission: radiance (linear RGB)
        ConductorParams<float> conductor;   ///< conductor, and the conductor base of a layered material (in the coat's medium)
        DielectricParams<float> dielectric; ///< dielectric
        LayeredParams<float> layered;       ///< layered
        int32_t medium;            ///< dielectric: the medium inside it (index into KScene::media), -1 if none
        // The denoiser's features (--aov, --denoise), as the engine's Material reports them: a camera path records them at its first
        // vertex whose material is not specular (Render/Integrator.cpp).
        float aovAlbedo[3];        ///< Material::aovAlbedo: the lighting-free color in [0, 1]
        uint32_t aovSpecular;      ///< Material::isSpecular: 1 = the path looks through (smooth glass) and records the next vertex
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

    /// The denoiser's features of one camera path (the engine's PathAOV): albedo and shading normal (facing the ray) of its first
    /// vertex that is not specular, or the environment's color clamped to [0, 1] (normal 0) if it escapes before one; zero if neither.
    struct KAov { float albedo[3], normal[3]; uint32_t valid; };

    enum : uint32_t { kCamPinhole = 0, kCamThinLens = 1, kCamLens = 2 };
    enum : uint32_t { kLensSphere = 0, kLensFlat = 1, kLensStop = 2 };

    /// A surface of a real lens as the engine's LensSystem holds it (Render/LensSystem.hpp), in the lens frame: metres, z toward the
    /// object, the first vertex at z = 0, the film behind the lens at negative z.
    template <class T> struct LensSurface {
        T radius;                  ///< a sphere's radius (> 0: its centre on the film side)
        T z;                       ///< the vertex on the axis
        T apertureRadius;          ///< the clear semi-diameter (the stop's open radius is this times CameraParams::stopScale)
        T nBehind;                 ///< index of refraction of the medium behind the surface (toward the film); 1 = air
        uint32_t kind;             ///< kLensSphere | kLensFlat | kLensStop
    };

    /**
     * @brief The camera as the engine's Camera holds it (Render/Camera.hpp): a pinhole, a thin lens, or a real lens traced surface by
     *        surface. Film point (s, t) in [0, 1)^2 (t = 0 at the bottom) looks from the lens toward lowerLeft + s horizontal + t vertical,
     *        a point of the plane of focus. A real lens has its surfaces, its exit pupil's bounds per band of film radius and its stop's
     *        edge normals in KScene's lens arrays.
     */
    template <class T> struct CameraParams {
        uint32_t type;             ///< kCamPinhole | kCamThinLens | kCamLens
        T origin[3], lowerLeft[3], horizontal[3], vertical[3];
        T u[3], v[3], w[3];        ///< right, up and backward (the view is -w)
        // the thin lens: a round aperture of radius lensRadius, or a regular polygon of `blades` corners and that circumradius
        T lensRadius;
        uint32_t blades;
        T bladeRotation, bladeRadius;
        // the real lens
        uint32_t numSurfaces, numBands;   ///< surfaces from the front; bands of film radius of the pupil table
        uint32_t stopBlades;       ///< the stop's blades (0: round); their edge normals are in KScene::lensBladeNormals
        T filmZ, filmW, filmH, filmRadius, centreIrradiance, frontPrincipalZ, stopScale, apothemPerRadius;
    };
    using KCamera = CameraParams<float>;

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
        const LensSurface<float>* lensSurfaces;   ///< camera.numSurfaces (a real lens)
        const float* lensPupil;       ///< camera.numBands x 5: the exit pupil's bounds minX, minY, maxX, maxY and 1 if the band is empty
        const float* lensBladeNormals;///< camera.stopBlades x 2
        const KMedium* media;         ///< the media inside dielectrics (KMaterial::medium)
        const float* phaseF;          ///< the phase functions' tables (KMedium::tableOffset, tableCount): density in theta at the nodes ...
        const float* phaseCdf;        ///< ... and the cumulative distribution there
        uint32_t numMedia, numPhaseNodes;
        uint32_t mediaNee;            ///< --media-nee: lights are sampled at every scattering event in a medium, and an index-matched smooth
                                      ///  boundary (relative index 1) is crossed without counting as a bounce, by paths and by shadow rays
        uint32_t volumeCap;           ///< collisions in media per path (--volume-cap)
        uint32_t width, height;       ///< pixels
        uint32_t spp, maxDepth;       ///< samples per pixel; path length as the engine counts it (surface vertices)
        uint64_t seed;                ///< the run seed of the random number streams
    };

} // namespace rayt::kernel
