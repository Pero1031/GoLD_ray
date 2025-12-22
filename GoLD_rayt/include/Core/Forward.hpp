#pragma once

namespace rayt {
    // ---------------------------------------------------------------------
    // Forward Declarations
    // ---------------------------------------------------------------------
    // Declaring classes here allows us to use pointers/references in header files
    // without including the full definition, reducing compilation time and circular dependencies.

    // Scene & Geometry
    class Scene;
    class Shape;
    class Primitive;
    class GeometricPrimitive;
    class Aggregate; // BVH, etc.

    // Main Pipeline
    class Camera;
    class Sampler;
    class Integrator;
    class Film;
    class Filter;

    // Rays & Interactions
    struct Ray;
    class RayDifferential; 
    struct Interaction;
    struct SurfaceInteraction;

    // Lighting
    class Light;
    class AreaLight;
    class VisibilityTester;

    // Shading & Materials
    class Material;
    class BSDF;
    class BxDF;
    class Texture;          // Note: If Texture is a template class later, remove this.
    class Fresnel;
    class MicrofacetDistribution;

    // Volumetrics
    class Medium;
    struct MediumInterface;
    class PhaseFunction;
}