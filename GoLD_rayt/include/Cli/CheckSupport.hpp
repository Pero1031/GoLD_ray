#pragma once
// Helpers of the self-checks: estimating the radiance of a ray, small test scenes, closed-form references.

#include <functional>
#include <memory>
#include <string>
#include "Core/Types.hpp"
#include "Core/Ray.hpp"
#include "Color/EmissionSpectrum.hpp"
#include "Scene/Aggregate.hpp"
#include "Scene/Scene.hpp"

namespace rayt {
    class EnvMap;
    class CausticPhotonMap;
    class SpecularNee;
}

namespace rayt::cli {

// -----------------------------------------------------------------------------
// Light tests (--scene lights)
// -----------------------------------------------------------------------------
// Mean linear RGB (and its standard error) of Li() along one fixed ray.
void EstimateRadiance(const Scene& scene, const Ray& ray, int n, int maxDepth, Vector3& mean, Vector3& se, int& badSamples, const CausticPhotonMap* caustics = nullptr, const SpecularNee* nee = nullptr, bool skipDirectCaustics = false, bool mediaNee = false);

// Passes when every channel is within 4 standard errors (+0.5 % of the value) of the expectation.
bool ReportLightCase(const std::string& name, const Vector3& expected, const Vector3& mean, const Vector3& se, int bad);

// The variance (per RGB channel) of the photon map's radiance at the first hit of `ray`, from the photons' own randomness. A map is built once, so
// this noise is the same offset in every camera sample and does not appear in the standard error that EstimateRadiance returns; a comparison
// against the map adds it. Zero when the ray hits nothing.
Vector3 PhotonVarianceAt(const Scene& scene, const CausticPhotonMap& map, const Ray& ray);

// Lambertian ground (a huge sphere, locally flat) with whatever `addLights` puts above it.
std::shared_ptr<Scene> MakeGroundScene(Real albedo, const std::function<void(Aggregate&, Scene&)>& addLights, const char* lightSampler = "spatial");

// Independent reference for a quad light: radiance of the ground point at the origin (normal +y)
//   Lo = albedo / pi * L * Integral cos(receiver) cos(light) / r^2 dA     (midpoint rule on the quad)
Real QuadGroundRadiance(const Point3& corner, const Vector3& e1, const Vector3& e2, bool twoSided, Real L, Real albedo);

// Environment map: dim bluish sky that brightens upward + a small bright sun (elevation 50 deg).
std::shared_ptr<EnvMap> MakeSunSkyEnv(Vector3& sunDirOut);

// Independent reference: radiance leaving a white-ish Lambertian ground point (normal +y) under
// an environment map = albedo / pi * Integral_{upper hemisphere} Le(w) cos(theta) dw, evaluated
// on a dense (theta, phi) grid with EnvMap::eval (no importance sampling involved).
Vector3 GroundRadianceFromEnv(const EnvMap& env, Real albedo);

// Consistency of EnvMap::sample()/pdf(): the pdf integrates to 1, sample() and pdf() agree, and
// Le/pdf estimates the same integral as the grid.
bool CheckEnvSampling(const EnvMap& env, const std::string& label, int n);
} // namespace rayt::cli
