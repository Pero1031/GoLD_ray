#pragma once
// Building blocks shared by the demo scenes and the self-checks: meshes as hittables, lights, assets.

#include <memory>
#include <string>
#include "Core/Texture.hpp"
#include "Color/EmissionSpectrum.hpp"
#include "Color/SpectralCurve.hpp"
#include "Geometry/TriangleMesh.hpp"
#include "Materials/Material.hpp"
#include "Materials/Dielectric.hpp"
#include "Scene/Aggregate.hpp"
#include "Scene/Hittable.hpp"
#include "Scene/Scene.hpp"
#include "Cli/Options.hpp"

namespace rayt::cli {

// Triangle mesh -> Hittable with the chosen intersection backend (own SAH BVH or Embree).
std::shared_ptr<Hittable> MakeMeshHittable(const std::shared_ptr<const TriangleMesh>& mesh, const std::shared_ptr<Material>& material, const std::string& accel);

// Procedural wood floor (size x size texels): `planks` boards side by side with their own tint, growth rings bent by low-frequency
// noise, fibers along the boards and dark seams between them. Dark brown on purpose (about 0.2 reflectance on average).
std::shared_ptr<Texture2D> MakeWoodTexture(int size, int planks, unsigned seed);

// Absorption coefficient of pure water (1/m) from the table in assets/Water (the file is in 1/cm).
bool LoadWaterAbsorption(color::SpectralCurve& sigmaA);

// Spherical emitter. registerAsLight = false leaves it out of the light list (no NEE: pure BSDF sampling).
void AddSphereLight(Aggregate& world, Scene& scene, const Point3& center, Real radius, const EmissionSpectrum& color, Real intensity, bool registerAsLight);

// Rectangular emitter (corner + two edges; emits toward cross(edge1, edge2)).
std::shared_ptr<Light> AddQuadLight(Aggregate& world, Scene* scene, const Point3& corner, const Vector3& e1, const Vector3& e2, const EmissionSpectrum& color, Real intensity, bool registerAsLight, bool twoSided = false);

// The prototype scene of the specular manifold NEE: a dispersive triangular prism (IOR 1.55, --abbe) standing on a wooden-gray disc, one
// point light beside it. Shared by --scene manifold (image) and --scene manifoldtest (comparison with the photon map).
struct PrismWorld {
    std::shared_ptr<Aggregate> world;
    std::shared_ptr<const TriangleMesh> mesh;
    std::shared_ptr<Dielectric> glass;
    // The light shines on the left face at about 45 degrees from its normal (near the minimum deviation of a 60 degree prism): at 30 degrees
    // the ray would hit the second face beyond the critical angle (IOR 1.75) and no two-refraction path would exist at all.
    Point3 lamp = Point3(-0.19, 0.11, -0.064);
};

PrismWorld BuildPrismWorld(const CliOptions& opt, const std::string& accel);
} // namespace rayt::cli
