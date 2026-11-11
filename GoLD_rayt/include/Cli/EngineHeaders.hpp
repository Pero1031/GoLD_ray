#pragma once
// The engine headers the CLI sources work with, in one place (the CLI is an application on top of the header-only engine).
// A source file of the CLI includes this after "pch.h" instead of repeating the list.

#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

// Core
#include "Core/Core.hpp"
#include "Core/Fresnel.hpp"
#include "Core/Image.hpp"
#include "Core/Interaction.hpp"
#include "Core/RenderSettings.hpp"
#include "Core/Sampling.hpp"

// Scene, geometry, acceleration
#include "Accelerators/BVH.hpp"
#include "Accelerators/EmbreeMesh.hpp"
#include "Geometry/Frame.hpp"
#include "Geometry/HairCurves.hpp"
#include "Geometry/HairGenerator.hpp"
#include "Geometry/Quad.hpp"
#include "Geometry/Sphere.hpp"
#include "Geometry/TriangleMesh.hpp"
#include "Scene/Aggregate.hpp"
#include "Scene/Hittable.hpp"
#include "Scene/Primitive.hpp"
#include "Scene/Scene.hpp"

// Rendering
#include "Render/Camera.hpp"
#include "Render/Denoiser.hpp"
#include "Render/Film.hpp"
#include "Render/Integrator.hpp"
#include "Render/VolumeBPT.hpp"

// Materials and lights
#include "Materials/CoatedDiffuse.hpp"
#include "Materials/Hair.hpp"
#include "Materials/SkinRelief.hpp"
#include "Materials/WrinkleField.hpp"
#include "Materials/SkinPigment.hpp"
#include "Materials/LayeredCoat.hpp"
#include "Materials/DiffuseLight.hpp"
#include "Materials/Dielectric.hpp"
#include "Materials/Lambertian.hpp"
#include "Materials/LayeredSkin.hpp"
#include "Materials/Material.hpp"
#include "Materials/RoughConductor.hpp"
#include "Materials/ThinFilm.hpp"
#include "Lights/AreaLight.hpp"
#include "Lights/PointLight.hpp"

// Color and input / output
#include "Color/ConductorFromReflectance.hpp"
#include "Color/EmissionSpectrum.hpp"
#include "Color/IorRgbApprox.hpp"
#include "IO/ComplexIorLoader.hpp"
#include "IO/EnvMap.hpp"
#include "IO/ImageLoader.hpp"
#include "IO/ObjLoader.hpp"
#include "IO/PmxLoader.hpp"
#include "IO/PmxMorph.hpp"
#include "IO/PmxPose.hpp"
#include "IO/SpectralCurveLoader.hpp"
#include "IO/FadingSeries.hpp"

#include "DebugTools/FrameDebug.hpp"
