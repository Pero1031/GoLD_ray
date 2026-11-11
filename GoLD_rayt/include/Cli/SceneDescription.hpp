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

// A light spectrum: a preset (d65 | tungsten | candle | led), a blackbody temperature, an RGB color (x D65), a table from a CSV file,
// or a weighted sum of tables (a mix of LEDs).
struct SpectrumSpec {
    enum class Kind { Preset, Blackbody, Rgb, Curve, Mix };
    Kind kind = Kind::Rgb;
    std::string preset;
    double kelvin = 0;
    Vector3 rgb = Vector3(1.0);
    std::string file;
    std::vector<std::pair<std::string, double>> mix;   // Mix: (CSV file, weight)
};

// One level of the wrinkle simulation (Materials/WrinkleField.hpp); the fields are those of WrinkleField::Level.
struct WrinkleLevelSpec {
    double spacingMm = 1.5, angleDeg = 0, base = 0.2, distance = 0.7, continuation = 0.3, similarity = 0.3, crossing = 0.5, deposit = 0.25;
    double width = 0.15, poreWidth = 0.1, poreMin = 0.1, poreBlend = 0.5, poreSkew = 0, wobble = 0.15, wobbleFreq = 2, depthUm = 30;
    bool bulge = true, darkPores = true;
};

struct MaterialSpec {
    enum class Type { Diffuse, Conductor, Dielectric, Emission, ThinFilm, Layered, Hair, Skin };
    Type type = Type::Diffuse;
    Vector3 albedo = Vector3(0.5);                 // diffuse
    std::string metal;                             // conductor: Au | Ag | Cu (the Johnson & Christy tables in assets/Metal_data) ...
    std::string table;                             // ... or a complex IOR table (RefractiveIndex.info CSV)
    std::string reflectance;                       // ... or a measured normal-incidence reflectance R(lambda) (CSV; Gulbrandsen's method)
    std::string reflectanceTexture;                // diffuse: a measured reflectance per texel R(u, v, lambda) (.sptex, Core/SpectralTexture.hpp)
    std::string fadingSeries, fadingColor;         // diffuse: a measured reflectance that changes with the time of exposure (IO/FadingSeries.hpp):
    double fadingHours = 0;                        //   the series folder, the patch in it, and the hours of exposure ("fading")
    double fadingUnexposedHours = 0;               //   with an exposure field: the hours where nothing reached the surface (hours = where everything did)
    std::string fadingMap;                         //   exposure field: a gray image over the surface's uv (white = exposed)
    bool fadingMat = false;                        //   ... a covered border in uv (a mat board): left, right, bottom, top, and the softness of its edge
    double fadingMatInset[4] = { 0, 0, 0, 0 }, fadingMatSoft = 0.002;
    bool fadingGradient = false;                   //   ... exposure 1 at fadingBright and 0 at fadingDark, linear between (world space)
    Vector3 fadingBright = Vector3(0.0), fadingDark = Vector3(0, 1, 0);
    double fadingGradientExponent = 1;
    double fadingNoiseAmount = 0, fadingNoiseFrequency = 10;   //   ... uneven fading: 1 - amount * noise (world space)
    int fadingNoiseOctaves = 3, fadingNoiseSeed = 1;
    bool hasExposureField() const { return !fadingMap.empty() || fadingMat || fadingGradient || fadingNoiseAmount > 0; }
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
    // hair: a dielectric fiber with pigment (Materials/Hair.hpp, docs/HAIR.md); the absorption is given as melanin, as the apparent color, or directly
    Vector3 hairSigmaA = Vector3(0.0);             //   absorption per hair radius (R, G, B); filled from melanin or color when those are used
    double hairEta = 1.55;
    double hairBetaM = 0.3, hairBetaN = 0.3;       //   longitudinal / azimuthal roughness, 0..1
    double hairTiltDeg = 2;                        //   cuticle scale tilt
    unsigned hairLobes = 0xF;                      //   lobes switched on: bit 0 R, 1 TT, 2 TRT, 3 the rest
    // skin: a coated diffuse (albedo, roughness, specularWeight), optionally with subsurface scattering (Burley profile) and a micro-relief
    double skinSpecular = 0.4;                     //   weight of the oily coat (0 = matte)
    bool skinSss = false;                          //   subsurface scattering on
    Vector3 skinScatterMm = Vector3(3.67, 1.37, 0.68);   //   scatter distance per channel (mm), x skinSssScale
    double skinSssScale = 1;
    bool skinWalk = false;                         //   the subsurface light walks through the inside of the object ("randomWalk") instead of the BSSRDF
    double skinWalkG = 0;                          //   mean cosine of its phase function
    bool skinRelief = false;                       //   furrows and pores (Materials/SkinRelief.hpp)
    double reliefCellMm = 0.6, reliefFurrowDepthMm = 0.04, reliefFurrowWidth = 0.12, reliefPoreDepthMm = 0.08, reliefPoreRadiusMm = 0.07, reliefPoreFraction = 0.35;
    double reliefStretch = 1, reliefCoarseScale = 3, reliefCoarseWeight = 0.8;
    double reliefPoreSizeVariation = 0.5, reliefPoreDensityVariation = 0.6, reliefPoreDensityFreq = 60, reliefPoreDarkness = 0;
    double skinGlossVariation = 0, skinGlossFreq = 60;   //   oily and dry patches (CoatedDiffuse::setGlossVariation)
    int skinGlossSeed = 1;
    bool skinFilmOn = false;                       //   a liquid film that pools in the furrows (Materials/SkinFilm.hpp, docs/SKIN_MOISTURE.md); needs the relief
    double skinFilmThicknessUm = 5, skinFilmIor = 1.33, skinFilmRoughness = 0.02, skinFilmSoftnessUm = 1.5;
    bool skinFoundationOn = false;                 //   a layer of foundation make-up on top (Materials/CosmeticLayer.hpp): "foundation"
    double fndThickness = 0.5, fndDiffuserFraction = 0.9, fndG1 = 0.5, fndG2 = -0.2, fndLobeWeight = 0.7;   //   its optical thickness, share of diffusers, their two lobes
    double fndPlateletRoughness = 0.15, fndPlateletTiltDeg = 0;                                               //   the platelets
    double fndLeveling = 1.0;                                                                                 //   how far the layer levels the skin relief (0: follows it, 1: smooth)
    Vector3 fndDiffuserAlbedo = Vector3(0.97, 0.94, 0.90), fndPlateletAlbedo = Vector3(0.85);                 //   the colours of the two particles
    bool skinScatterOn = false;                    //   the white scatter of flaking skin (Materials/SkinScatter.hpp)
    double scatRho0 = 0.08, scatFrequency = 300, scatThreshold = 0.52, scatSoftness = 0.08, scatHydration = 0, scatKappaH = 0.3, scatGapUm = 1.0;
    int scatSeed = 1;
    double skinEmulsionTau = 0;                    //   transport optical depth of the oil droplets of an emulsion across the film (Materials/SkinLotion.hpp; set by "lotion" with "oil")
    std::string skinLotionNote;                   //   a lotion at one time sets the film and the hydration above (Materials/SkinLotion.hpp); this describes the state it chose
    bool skinPigment = false;                    //   uneven color (Materials/SkinPigment.hpp)
    double pigFreckleCellMm = 2.5, pigFreckleFraction = 0.5, pigFreckleRadiusMm = 0.8, pigFreckleStrength = 0.8, pigClusterFreq = 40;
    double pigBlotchFreq = 150, pigBlotchStrength = 0.15, pigRednessFreq = 35, pigRednessStrength = 0.5;
    int pigSeed = 1;
    int reliefSeed = 1;
    bool reliefWrinkles = false;                   //   furrows and pores from a graph simulation (replaces the Worley furrows and pores)
    double wrTileMm = 24, wrDepthScale = 1;
    int wrResolution = 2048, wrSteps = 300, wrSeed = 1;
    std::vector<WrinkleLevelSpec> wrLevels;        //   empty: the default skin (three levels)
    SpectrumSpec spectrum;                         // emission
    double radiance = 1.0;                         // emission
    bool twoSided = false;                         // emission
};

// A set of procedural strands (Geometry/HairGenerator.hpp): a hanging "lock", or hair falling over a spherical "head".
struct HairSpec {
    std::string style = "lock";
    int strands = 3000, segments = 24;
    double length = 0.3;
    double rootRadius = 8e-5, tipRadius = 5e-5;    // fiber radius (m)
    int seed = 1;
    // lock: hangs down from `root`
    Point3 root = Point3(0.0);
    double width = 0.012, tipWidth = 1.4, waveAmp = 0.018, waveLength = 0.16, twist = 0.5, frizz = 0.0006, flyaway = 0.03;
    // head: roots on a sphere
    Point3 center = Point3(0.0);
    double headRadius = 0.09, droop = 9, outward = 0.55, sway = 5, swayFreq = 28, headFrizz = 1.2, maxPolarDeg = 105, hairline = 0.55, clumpStrength = 0.65;
    int clumps = -1;                               // sub-locks (lock) or guide paths (head); -1: 9 for a lock, 500 for a head
};

struct ObjectSpec {
    enum class Shape { Sphere, Quad, Mesh, Lathe, Ocean, Hair };
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
    HairSpec hair;                                 // hair
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
