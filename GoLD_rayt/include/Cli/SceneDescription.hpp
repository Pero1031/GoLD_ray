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
    enum class Type { Diffuse, Conductor, Dielectric, Emission, ThinFilm, Layered };
    Type type = Type::Diffuse;
    Vector3 albedo = Vector3(0.5);                 // diffuse
    std::string metal;                             // conductor: Au | Ag | Cu (the Johnson & Christy tables in assets/Metal_data) ...
    std::string table;                             // ... or a complex IOR table (RefractiveIndex.info CSV)
    std::string reflectance;                       // ... or a measured normal-incidence reflectance R(lambda) (CSV; Gulbrandsen's method)
    double roughness = 0;                          // conductor, dielectric
    bool multipleScattering = true;                // conductor: Kulla-Conty compensation
    double anisotropy = 0;                         // conductor: -1..1; > 0 rougher along the tangent, < 0 smoother (brushed along it)
    double rotation = 0;                           // conductor: degrees, rotates the tangent about the normal
    // layered: a clear coat over a base (Materials/LayeredCoat.hpp). The base reuses albedo (diffuse) or metal / table / reflectance + roughness (conductor).
    double coatIor = 1.5;
    double coatRoughness = 0.1;
    Vector3 coatTint = Vector3(1.0);               //   transmittance of the coat along the normal (1 = clear)
    int coatSamples = 1;                           //   random walks per evaluation
    bool baseIsConductor = false;
    double farRoughness = -1;                      // dielectric: >= 0: the roughness grows to this (perceptual) between farStart and farEnd (m from the y axis)
    double farStart = 50, farEnd = 500;
    double ior = 1.5;                              // dielectric
    double abbe = 0;                               // dielectric: Abbe number (0 = no dispersion)
    std::string sellmeierPreset;                   // dielectric: Sellmeier dispersion of a named glass (diamond, sapphire, bk7, silica) ...
    bool hasSellmeier = false;                     //   ... or its coefficients (b[3] and c[3] in um^2)
    double sellmeierB[3] = { 0, 0, 0 }, sellmeierC[3] = { 0, 0, 0 };
    bool hasMedium = false;                        // dielectric: the (closed) inside is filled with a scattering medium
    Vector3 mediumSigmaT = Vector3(100.0);         //   extinction (1/m), per channel
    Vector3 mediumAlbedo = Vector3(0.9);           //   single-scattering albedo, per channel
    double mediumG = 0;                            //   Henyey-Greenstein anisotropy (the first lobe)
    double mediumG2 = 0;                           //   second lobe (used when mediumWeight < 1)
    double mediumWeight = 1;                       //   weight of the first lobe
    std::string mediumPhaseTable;                  //   a CSV "theta_deg, p" (Mie scattering...): replaces g, g2 and weight
    bool mediumNoise = false;                      //   a heterogeneous medium: sigmaT is the extinction at the densest point, a noise field gives q(x)
    double noiseFrequency = 20;                    //     features per metre (the first octave)
    int noiseOctaves = 4;
    double noiseLo = 0.45, noiseHi = 0.65;         //     the noise below lo gives density 0, above hi density 1
    int noiseSeed = 1;
    double maskRadius = 0;                         //     > 0: the density fades to 0 at this distance from maskCenter (a soft-edged cloud in a sphere)
    Vector3 maskCenter = Vector3(0.0);
    double maskStart = 0.5;                        //     the fade starts at this fraction of the radius
    // thinFilm: a smooth surface under one thin transparent layer (Materials/ThinFilm.hpp, docs/THIN_FILM.md)
    double filmIor = 1.33;                         //   index of the film (soap water 1.33, oil 1.45)
    double filmThickness = 400;                    //   nm; the start of the gradient when filmGradient is set
    bool filmGradient = false;                     //   thickness changes linearly along filmAxis ...
    Vector3 filmAxis = Vector3(0, 1, 0);
    double filmRange0 = 0, filmRange1 = 1;         //   ... from filmThickness at dot(p, axis) = range0 to filmTop at range1
    double filmTop = 400;
    double filmExponent = 1;                       //   gradient profile t^exponent (> 1: the top stays thin for longer)
    double filmWarpAmp = 0;                        //   the position is displaced by up to +- this many metres (bands wobble and swirl)
    double filmWarpFreq = 3;
    int filmWarpSeed = 5;
    Vector3 filmNoiseStretch = Vector3(1.0);       //   the noise coordinates are multiplied by this ((1, 3, 1): horizontal streaks)
    double filmNoiseAmp = 0;                       //   +- nm of fractal noise (swirls)
    double filmNoiseFreq = 10;                     //   features per metre
    int filmNoiseOctaves = 3;
    int filmNoiseSeed = 1;
    bool filmCoverage = false;                     //   the film covers only patches (a noise field; elsewhere the surface is not there)
    double covFrequency = 0.5, covLo = 0.5, covHi = 0.58;
    int covOctaves = 4, covSeed = 3;
    std::string substrateKind = "free";            //   free (a bubble) | dielectric (oil on water) | conductor (coated metal)
    double substrateIor = 1.33;                    //   dielectric
    std::string substrateMetal;                    //   conductor: Au | Ag | Cu
    std::string substrateTable;                    //   ... or a complex IOR table
    SpectrumSpec spectrum;                         // emission
    double radiance = 1.0;                         // emission
    bool twoSided = false;                         // emission
};

struct ObjectSpec {
    enum class Shape { Sphere, Quad, Mesh, Lathe, Ocean };
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
    // ocean: a sea surface from a wave spectrum (Geometry/ProceduralMesh.hpp, docs/WATER.md), a polar grid around the origin
    double windSpeed = 8, windDirection = 0, spread = 4, peakEnhance = 1, chop = 0.8, time = 0, minWavelength = 0.05;
    int oceanSeed = 1, rings = 400, sectors = 512;
    double rInner = 0.8, rOuter = 3000;
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
