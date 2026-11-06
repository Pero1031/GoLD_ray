#pragma once
// Procedural triangle meshes used by the demo scenes: surfaces of revolution (lathe), a faceted lathe, a rippled pool surface,
// a Gerstner-wave sea and a prism.

#include <vector>

#include "Core/Types.hpp"
#include "Core/Precision.hpp"
#include "Geometry/TriangleMesh.hpp"

namespace rayt {

// Sea surface: a sum of Gerstner (trochoidal) waves, evaluated on a polar grid around the origin whose rings grow geometrically
// (fine near the camera, coarse toward the horizon). Wavelengths 0.5-64 m, four random components per octave, directions spread
// around the wind direction (normal distribution), the longest waves damped (Phillips-like cutoff for a 7 m/s wind). `steep` is
// k * A of every component. A component that the local ring spacing cannot resolve (fewer than ~2-5 samples per wavelength) fades
// out, so the distance does not alias. Normals are analytic; the mesh is one open sheet facing up (t = 0, a still frame).
struct SeaParams {
    Real steep = Real(0.03);
    unsigned seed = 1;
    Real windDeg = Real(0);       // direction the waves travel to, as an angle in the x-z plane (atan2(z, x))
    Real rInner = Real(0.8), rOuter = Real(3000);
    int rings = 400, sectors = 512;
};

// Surface of revolution about the y axis. `profile` is a polyline of (radius, height) points ordered counter-clockwise
// with the solid on the left, so the outward normal of a segment is (dy, -dr). The segment from the last point back to the
// first is not generated: both must lie on the axis (r = 0), otherwise the mesh is open. Every segment has its own
// vertices; the normal is flat along the profile and smooth around the axis (analytic, so a coarse ring count is enough).
TriangleMesh MakeLatheMesh(const std::vector<Point2>& profile, int segments);

// Surface of revolution with FLAT facets: the profile (same convention as MakeLatheMesh: (radius, height) points, counter-clockwise,
// solid on the left) is revolved into `sides` planar facets, each with its own vertices and its own face normal, so the faces
// refract light in one direction and throw sharp-edged caustics. sides = 6 turns a straight profile into a hexagonal column,
// and a profile that ends on the axis into a pyramid top.
TriangleMesh MakeFacetedLathe(const std::vector<Point2>& profile, int sides);

// Water surface of a pool: gentle ripples on the rectangle |x| <= w/2, |z| <= d/2 at y = 0, a sum of 14 plane waves (wavelengths
// 0.12-1 m, any direction, k * A = steep each) on a regular grid; analytic normals, one open sheet facing up.
TriangleMesh MakePoolSurface(Real w, Real d, Real steep, unsigned seed, Real spacing);

TriangleMesh MakeSeaMesh(const SeaParams& P);

// -----------------------------------------------------------------------------
// Main Entry Point
// -----------------------------------------------------------------------------
// A closed prism: the polygon `xz` (counter-clockwise or not: the winding is fixed per face) extruded from y = 0 to `height`, flat faces
// with their own normals (outward). The body of the specular manifold NEE prototype.
TriangleMesh MakePrismMesh(const std::vector<Point2>& xz, Real height);

} // namespace rayt
