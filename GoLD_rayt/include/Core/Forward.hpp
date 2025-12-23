#pragma once

/**
 * @file forward.hpp
 * @brief Forward declarations for the rayt engine classes and structs.
 * * Purpose:
 * - Reduces compilation time by minimizing header interdependencies.
 * - Prevents circular dependencies between core components (e.g., Scene and Integrator).
 * - Allows the use of pointers and references in headers without full definitions.
 */

namespace rayt {
    
    // ---------------------------------------------------------------------
    // Scene & Geometry
    // ---------------------------------------------------------------------
    // Fundamental building blocks of the 3D environment.
    class Scene;
    class Shape;                 // Pure geometry (Sphere, Triangle, etc.)
    class Primitive;             // Ties Shape and Material together
    class GeometricPrimitive;
    class Aggregate;             // Acceleration structures (BVH, Grid, etc.)

    // ---------------------------------------------------------------------
    // Main Pipeline
    // ---------------------------------------------------------------------
    // Components that drive the rendering process.
    class Camera;        // Generates primary rays
    class Sampler;       // Generates sample patterns (RNG)
    class Integrator;    // Solves the rendering equation
    class Film;          // Stores and processes pixel data
    class Filter;        // Reconstruction filter for antialiasing

    // ---------------------------------------------------------------------
    // Rays & Interactions
    // ---------------------------------------------------------------------
    // Data structures representing light paths and hit information.
    struct Ray;
    struct RayDifferential;   // Ray with derivatives for texture filtering
    struct Interaction;       // Generic hit information
    struct SurfaceInteraction; // Detailed hit data (UV, normal, tangents)

    // ---------------------------------------------------------------------
    // Lighting
    // ---------------------------------------------------------------------
    // Light sources and visibility calculations.
    class Light;
    class AreaLight;          // Light attached to a geometric shape
    class VisibilityTester;   // Connects shading point to light source

    // ---------------------------------------------------------------------
    // Shading & Materials
    // ---------------------------------------------------------------------
    // Physics of light scattering at a surface.
    class Material;                // Describes surface properties
    class BSDF;                    // Bidirectional Scattering Distribution Function
    class BxDF;                    // Base for BRDF and BTDF
    class Texture;                 // Note: Update to template if using Texture<T> later
    class Fresnel;                 // Reflection/refraction ratios
    class MicrofacetDistribution;  // Roughness models (GGX, Beckmann, etc.)

    // ---------------------------------------------------------------------
    // Volumetrics
    // ---------------------------------------------------------------------
    // Handling of light interacting with participating media (fog, smoke).
    class Medium;               // Volume properties
    struct MediumInterface;     // Transition between two media
    class PhaseFunction;        // Angular distribution of volumetric scattering

} // namespace rayt