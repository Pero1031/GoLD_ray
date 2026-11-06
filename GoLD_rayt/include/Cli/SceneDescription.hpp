#pragma once
// A scene written as data (the "objects" form of a scene file, see SceneFile.hpp): plain structs, no engine objects.
// SceneFile.cpp fills them from JSON and checks the document; Scenes/FileScene.cpp builds the scene from them.

#include <map>
#include <string>
#include <vector>

#include "Core/Types.hpp"
#include "Core/Precision.hpp"

namespace rayt::cli {

struct CameraSpec {
    Point3 from = Point3(0, 0.5, 2.5);
    Point3 at = Point3(0, 0, -1);
    Vector3 up = Vector3(0, 1, 0);
    Real fovDeg = 35.0;
};

// A light spectrum: a preset (d65 | tungsten | candle | led), a blackbody temperature, an RGB color (x D65) or a table from a CSV file.
struct SpectrumSpec {
    enum class Kind { Preset, Blackbody, Rgb, Curve };
    Kind kind = Kind::Rgb;
    std::string preset;
    double kelvin = 0;
    Vector3 rgb = Vector3(1.0);
    std::string file;
};

struct MaterialSpec {
    enum class Type { Diffuse, Conductor, Dielectric, Emission };
    Type type = Type::Diffuse;
    Vector3 albedo = Vector3(0.5);                 // diffuse
    std::string metal;                             // conductor: Au | Ag | Cu (the Johnson & Christy tables in assets/Metal_data) ...
    std::string table;                             // ... or a complex IOR table (RefractiveIndex.info CSV)
    double roughness = 0;                          // conductor, dielectric
    bool multipleScattering = true;                // conductor: Kulla-Conty compensation
    double ior = 1.5;                              // dielectric
    double abbe = 0;                               // dielectric: Abbe number (0 = no dispersion)
    SpectrumSpec spectrum;                         // emission
    double radiance = 1.0;                         // emission
    bool twoSided = false;                         // emission
};

struct ObjectSpec {
    enum class Shape { Sphere, Quad, Mesh, Lathe };
    Shape shape = Shape::Sphere;
    std::string material;
    int light = -1;                                // emission only: 1 = register as a light for next-event estimation, 0 = do not (-1: yes)
    Point3 center = Point3(0.0);                   // sphere
    Real radius = 1;
    Point3 corner = Point3(0.0);                   // quad (corner and two edges; faces toward cross(edge1, edge2))
    Vector3 edge1 = Vector3(1, 0, 0), edge2 = Vector3(0, 1, 0);
    std::string file;                              // mesh: an OBJ file
    double scale = 1;                              // mesh: uniform scale, then rotation, then translation
    Vector3 rotateAxis = Vector3(1, 0, 0);
    double rotateDegrees = 0;
    Vector3 translate = Vector3(0.0);
    std::vector<Point2> profile;                   // lathe: (radius, height) pairs of the profile that is turned around the y axis
    int segments = 96;                             // lathe
};

struct PointLightSpec {
    Point3 position = Point3(0.0);
    Vector3 color = Vector3(1.0);
    double intensity = 1.0;
};

struct EnvironmentSpec {
    enum class Type { None, Hdr, Uniform };
    Type type = Type::None;
    double scale = 1.0;
};

struct SceneDescription {
    CameraSpec camera;
    EnvironmentSpec environment;
    std::map<std::string, MaterialSpec> materials;
    std::vector<ObjectSpec> objects;
    std::vector<PointLightSpec> lights;
};

} // namespace rayt::cli
