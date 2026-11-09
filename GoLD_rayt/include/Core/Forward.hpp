/**
 * @file Core/Forward.hpp
 * @brief Forward declarations for the rayt engine classes and structs.
 * * Purpose:
 * - Reduces compilation time by minimizing header interdependencies.
 * - Prevents circular dependencies between core components (e.g., Scene and Integrator).
 * - Allows the use of pointers and references in headers without full definitions.
 */

#pragma once

namespace rayt {

    // ---------------------------------------------------------------------
    // Scene & Geometry
    // ---------------------------------------------------------------------
    
    // Fundamental building blocks of the 3D environment.
    class Scene;
    class Shape;                 // Pure geometry (Sphere, Triangle, etc.)
    class Primitive;             // Ties Shape and Material together
    class Aggregate;             // A plain list of objects (the scene builders fill one; rays go through a BVH over it)

    // ---------------------------------------------------------------------
    // Main Pipeline
    // ---------------------------------------------------------------------
    
    // Components that drive the rendering process.
    class Camera;        // Generates primary rays
    class Integrator;    // Solves the rendering equation
    class Film;          // Stores and processes pixel data

    // ---------------------------------------------------------------------
    // Rays & Interactions
    // ---------------------------------------------------------------------
    
    // Data structures representing light paths and hit information.
    struct Ray;
    struct RayDifferential;    // Ray with derivatives for texture filtering
    struct Interaction;        // Generic hit information
    struct SurfaceInteraction; // Detailed hit data (UV, normal, tangents)

    // ---------------------------------------------------------------------
    // Lighting
    // ---------------------------------------------------------------------
    
    // Light sources and visibility calculations.
    class Light;
    class AreaLight;          // Light attached to a geometric shape

    // ---------------------------------------------------------------------
    // Shading & Materials
    // ---------------------------------------------------------------------
    
    // Physics of light scattering at a surface.
    class Material;                // Describes surface properties (it evaluates and samples the scattering itself)
    class MicrofacetDistribution;  // Roughness models (GGX)

    // ---------------------------------------------------------------------
    // Volumetrics
    // ---------------------------------------------------------------------
    
    // Handling of light interacting with participating media (fog, smoke).
    class Medium;               // Volume properties
    struct PhaseFunction;       // Angular distribution of volumetric scattering

} // namespace rayt