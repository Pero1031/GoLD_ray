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

// Ocean surface from a wave spectrum (spectral synthesis; Tessendorf 2001 evaluates the same sum on a grid with an FFT).
//
// Wind speed U (m/s at 10 m) sets the Pierson-Moskowitz spectrum S(omega) = alpha g^2 omega^-5 exp(-beta (g / (U omega))^4), optionally
// sharpened at the peak (JONSWAP, `peakEnhance` = gamma). Waves travel with the wind, spread over directions as cos^(2s)(theta / 2).
// Deep-water dispersion omega^2 = g k turns the spectrum into S(k) and gives every component its phase speed. The surface is the sum of
// kBins x dirBins plane waves, one random wavenumber, direction and phase per (log k, direction) bin, each with the amplitude
// a = sqrt(2 S(k) dk D(theta) dtheta): the energy of its bin. The sum is not periodic (no tiling), and its height and slope
// statistics follow the spectrum (checked by --scene ocean). `chop` scales the horizontal (Gerstner) displacement that sharpens crests.
// The mesh is the polar grid of the sea (rings growing geometrically); a component the local spacing cannot resolve fades out.
struct OceanParams {
    Real windSpeed = Real(8);       // m/s
    Real windDeg = Real(0);         // direction the waves travel to, atan2(z, x) in degrees
    Real spread = Real(4);          // s of the directional spreading cos^(2s)(theta / 2)
    Real peakEnhance = Real(1);     // JONSWAP gamma (1 = Pierson-Moskowitz, 3.3 = a growing sea)
    Real chop = Real(0.8);          // 0 = pure heights, 1 = full Gerstner displacement
    Real time = Real(0);            // seconds (waves move at their phase speed)
    Real minWavelength = Real(0.05);
    unsigned seed = 1;
    int kBins = 56, dirBins = 24;
    Real rInner = Real(0.8), rOuter = Real(3000);
    int rings = 400, sectors = 512;
};

struct OceanWave {
    Real kx, kz, k, a, phase, omega;
};

/// The plane waves of the spectrum (random phases, jittered wavenumbers and directions): the same set the mesh is built from.
std::vector<OceanWave> MakeOceanWaves(const OceanParams& P);

/// Height variance and mean-square slope of the set (sum of a^2 / 2, sum of k^2 a^2 / 2) next to the values of the continuous spectrum
/// between the same wavenumber limits (numerical integral), for the self-check.
struct OceanStats {
    double heightVarianceSet = 0, heightVarianceSpectrum = 0;
    double slopeVarianceSet = 0, slopeVarianceSpectrum = 0;
    double peakWavenumber = 0, kMin = 0, kMax = 0;
};
OceanStats OceanSpectrumStats(const OceanParams& P, const std::vector<OceanWave>& waves);

/// Height of the surface (no displacement, no fading) at (x, z): the sum of the components.
Real OceanHeightAt(const std::vector<OceanWave>& waves, Real time, Real x, Real z);

TriangleMesh MakeOceanMesh(const OceanParams& P);

// -----------------------------------------------------------------------------
// Main Entry Point
// -----------------------------------------------------------------------------
// A closed prism: the polygon `xz` (counter-clockwise or not: the winding is fixed per face) extruded from y = 0 to `height`, flat faces
// with their own normals (outward). The body of the specular manifold NEE prototype.
TriangleMesh MakePrismMesh(const std::vector<Point2>& xz, Real height);

} // namespace rayt
