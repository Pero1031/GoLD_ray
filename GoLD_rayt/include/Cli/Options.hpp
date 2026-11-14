#pragma once
// Command line options of the CLI renderer: the defaults, the option struct and the parser.

#include <memory>
#include <optional>
#include <string>

#include "Core/Precision.hpp"
#include "Cli/Lens.hpp"

namespace rayt::cli {

struct SceneDescription;   // a scene written as data (SceneDescription.hpp)


// -----------------------------------------------------------------------------
// Scene Configuration
// -----------------------------------------------------------------------------
const int IMAGE_WIDTH = 800;

const int IMAGE_HEIGHT = 450;      // 16:9 Aspect Ratio

const int SAMPLES_PER_PIXEL = 256; // Higher = less noise, slower

const int MAX_DEPTH = 16;          // Max recursion depth for rays


// env path
const std::string ENV_HDR_PATH = "assets/env/bryanston_park_sunrise_2k.hdr";


const std::string AU_IOR_CSV_PATH = "assets/Metal_data/Johnson_Au.csv";


// -----------------------------------------------------------------------------
// Command line
//   main_cli [--scene gold|glass|furnace] [--width N] [--spp N] [--depth N]
//            [--env file.hdr] [--abbe V] [--out file.png] [--trace x,y]
//
//   --trace x,y : debug aid. Prints a few light paths of the pixel (x,y) (image coordinates,
//                 origin at the top-left) with the object hit at every bounce, then exits.
//
//   gold    : three gold spheres (default, as before)
//   glass   : clear / dispersive / frosted glass spheres on a floor
//   furnace : objects in a uniform white environment. The mean color of every object must match
//             its expected value (1.0 for lossless materials; a value below/above means energy is
//             lost/created, i.e. a bug). Prints OK/FAIL per object and exits with code 2 on FAIL.
//             Debug builds also assert that all ray directions have unit length.
//   spd     : colored diffuse spheres + gold under ONE sphere light whose spectrum is chosen with
//             --light (d65 | tungsten = illuminant A | candle = 1850 K blackbody | led = synthetic
//             SPD loaded from assets/Illuminants/*.csv); --wb white-balances the film to that light.
//             No environment light (black background). Spectral mode shows metamerism best.
//   goldbowl: a gold bowl with a gold ball in it, a white reference ball and a gray floor under one sphere light (--light as in spd):
//             light reflected between gold surfaces under a non-D65 light is where RGB and spectral rendering differ most.
//   pmx     : an MMD model from MMD_model/ (--pmx-list, --pmx-index N, --pmx-info, --pmx-yaw deg). Diffuse textures (uv, repeat, bilinear)
//             and alpha cut-outs; "+" overlay materials are skipped unless --pmx-all. HDR environment lighting.
//   goldrough: gold spheres of roughness 0.4 / 0.7 / 1.0 under the HDR environment; --no-ms = without the
//             multiple-scattering compensation (to compare).
//   studio  : gold / dispersive glass / red diffuse sphere lit by a tungsten rectangular softbox
//             and a small cool-white LED sphere light (SPDs from EmissionSpectrum / CSV). No environment.
//   mesh    : triangle meshes: runs mesh checks first (geometry, OBJ, BVH vs brute force, furnace,
//             emissive triangles), then renders a tilted gold torus (assets/Meshes/torus.obj) and a
//             flat- vs smooth-shaded low-poly sphere pair under the studio lights.
//   bench   : BVH benchmark (recursive reference vs SAH flat BVH) on icospheres of 20k / 82k / 328k triangles;
//             rays = --spp x 2000. Use a Release build for meaningful speeds.
//   metals  : metal spectral data tests (no image): file layouts / units, reflectance input,
//             edge-tint accuracy against the Johnson tables, rendered gold. Samples = --spp x 1024.
//   metalprobe: a smooth metal sphere (--metal Au|Ag|Cu) inside a uniform emissive shell with the spectrum of
//             --light. One specular bounce, so the mean film color per angle of incidence (printed) is
//             F(theta) x light. Build it in RGB and in spectral mode and compare with tools/metal_accuracy.py.
//   skindemo: a skin-colored sphere and ring (centimeter scale) under a key light, a fill and a back light, to see
//             subsurface scattering. Plain diffuse-like skin without options; --sss = real skin scatter distances
//             (3.7 / 1.4 / 0.7 mm for R / G / B), --sss-scale S multiplies them (3 = anime-like glow).
//   skinlayers: Monte Carlo radial profile of a layered skin (epidermis with melanin over a dermis with blood), written to a
//             table (--skin-out) for skindemo --skin-table; prints the profile per RGB channel, the Burley profile with the
//             same albedo and mean radius, and a check against diffusion theory. --melanin --blood --oxygen --epidermis-mm
//             --photons --single-layer. skindemo --skin-table f.csv renders it, --skin-matched the matched Burley profile.
//   volume  : volumetric random walk tests (no image): Beer-Lambert transmittance of an absorbing slab, white furnaces of
//             scattering media (slab, glass sphere), and the reflectance of a semi-infinite slab against a standalone Monte Carlo.
//   blob    : a closed mesh (--mesh file.obj, default an icosphere of 5 cm) of glass filled with a homogeneous scattering
//             medium (--sigma-t 1/m, --albedo, --phase-g, --phase-g2/--phase-w for a two-lobe phase function, --glass-roughness; --sigma-t-rgb r,g,b and --albedo-rgb r,g,b make the medium chromatic), lit by the environment (--env); --bg-gray v
//             shows a uniform gray backdrop to camera rays that touch nothing. --volume-cap limits collisions per path.
//   water   : a glass tumbler of water (metres) in front of a striped wall, lit by softboxes. Water = smooth dielectric (IOR 1.333)
//             with the measured absorption spectrum of pure water (Pope & Fry 1997, assets/Water) applied as a per-wavelength
//             transmittance (AbsorbingMedium: red is absorbed first, deep water turns cyan). --absorb-scale S multiplies the
//             absorption (default 40: a 6 cm glass is colorless at 1, so the tint is exaggerated); 0 = no absorption.
//             --water-sky: lit by the environment map (--env) instead of softboxes, no wall, camera looking down at --cam-elev
//             degrees (default 38) from --cam-yaw (default: the side of the map's brightest spot), to see the caustic on the table.
//             --caustics N: caustic photon map with N emitted photons (Render/CausticPhotonMap.hpp) for clean caustics even
//             from small / concentrated lights; --caustic-radius m = gather radius (default 0.0015, in the test 0.004).
//             --caustic-passes P: P passes of spp / P samples, each with its own independent photon set, the images averaged: the
//             grain of a photon map does not shrink with spp (the photons are fixed), but it does shrink as 1/sqrt(P), at the
//             memory of one photon set (each pass pays the photon tracing time once).
//             --water-lamp: one small point light (left, about 33 degrees up) instead of the softboxes / wall, camera looking down on
//             the caustic behind the glass (same --cam-elev / --cam-yaw). A point light can only be seen with --caustics.
//   sea     : prototype of a shallow tropical sea at sunrise: a Gerstner-wave surface (polar mesh, ~410k triangles, fine waves
//             fade out with distance) of smooth water (IOR 1.333, physical absorption of pure water) over sand and rocks
//             --sea-depth m above the seabed, lit by a procedural sky (blue gradient, a 1.2 degree sun at --sea-sun-elev degrees
//             that glitters on the waves; --sea-hdr = the --env map instead), camera 2 m above the water looking toward the
//             sun (--cam-yaw to turn). --sea-steep s = wave steepness (rough sea), --sea-seed n = another wave set.
//             No caustics on the seabed yet (that needs the photon map in absorbing media).
//   pool    : a shallow tiled swimming pool (3.2 x 5.2 m, --pool-depth m) in the --env map (the park HDR by default), seen from
//             the deck. Water = a rippled open sheet (IOR 1.333, physical absorption, --pool-steep s) over a tiled basin; a
//             sun (--pool-sun-elev / --pool-sun-az, 1.2 degree radius) is added to the map's own light and its caustic network
//             on the tiles comes from the caustic photon map (--caustics N, default 8 M; 0 = off; --caustic-radius, default
//             0.015 m). The map's sky has no photons: its light through the water is found by BSDF sampling (Light::photonCaustics).
//   crystal : prototype of faceted glass on a dark stage: a hexagonal crystal column, a cut-glass tumbler (8 flat facets), an
//             octagonal block and a pyramid on a wooden-colored disc in a black space, lit by one small point light (the
//             flat facets throw sharp-edged polygonal caustics). Dispersive glass (--abbe V, default 20; spectral mode only) gives
//             the rainbow edges. Caustics come from the photon map (--caustics N default 16 M, --caustic-radius default 0.0006 m).
//             --crystal-inside: the lamp sits INSIDE the empty tumbler and is the only light (all of its light leaves through the
//             glass facets as beams; default 24 M photons, radius 0.0012 m).
//   manifold: PROTOTYPE of a deterministic caustic: a dispersive triangular prism (--abbe) and one point light on a gray disc. The
//             caustic of the point light through the two refractions is solved exactly (specular manifold NEE, Render/SpecularNee.hpp:
//             Newton's method on Snell's law, no photons, no blur) with --caustic-method manifold (default); photons = the caustic photon
//             map for comparison (--caustics N); none = direct light only (the prism's shadow).
//             Also for --scene crystal --crystal-inside --caustic-method manifold: the lamp inside the empty tumbler, every ray through
//             the wall is a T T chain solved exactly (the wall chord must stay in the glass: the chord across the cavity is rejected); the
//             photon map (--caustics N, default 8 M) then only carries the other chains (it skips the T T chains of the registered bodies).
//   manifoldtest: no image: the prism caustic at a few lit floor points from the manifold estimator against the photon map (they
//             must agree), plus solver statistics.
//   manifoldcup: the same comparison for a hollow faceted tumbler with a lamp inside (T T through the wall).
//   caustic : photon map test (no image): floor radiance beside / under a glass sphere in a uniform sky, plain path tracing
//             (BSDF-sampled caustics, unbiased) against path tracing + caustic photons. Photons emitted = --spp x 8000.
//   lights  : light source tests (no image). A Lambertian ground is lit by one light and the
//             radiance of one ray is compared with a closed-form / independent numerical value
//             (point light, sphere area light with and without NEE, colored emission, environment
//             maps). Samples per test = --spp x 4096. Exits with code 2 on FAIL.
//   aov     : denoiser feature tests (no image): the albedo / normal recorded for a camera path on
//             a diffuse surface, a coated sphere, a light, through smooth glass and for escaping rays.
//   sss     : subsurface scattering tests (no image): the radial profile (normalization, inverse CDF,
//             mean radius), a white furnace over a skin plane and the profile across a shadow edge
//             against a numerically integrated reference. Samples per test = --spp x 4096.
// Without options the behavior is the same as before (gold scene, constants above).
// -----------------------------------------------------------------------------
struct CliOptions {
    std::string scene = "gold";
    int width = IMAGE_WIDTH;           // height follows the 16:9 aspect ratio of the constants
    int height = 0;                    // --height N: set the height explicitly (0 = from the format of the scene)
    std::string sceneFile;             // the scene file the options came from (empty: none)
    std::shared_ptr<const SceneDescription> sceneDescription;   // its "objects" form (scene "file"), read once while parsing
    int spp = SAMPLES_PER_PIXEL;
    long long seed = -1;               // --seed N: the random numbers of the path tracer (-1: a new seed every run, printed at the start)
    int maxDepth = MAX_DEPTH;
    std::string env = ENV_HDR_PATH;
    Real abbe = 20.0;                  // Abbe number of the dispersive glass (smaller = stronger)
    std::string out;                   // default: result_<scene>.png ("gold" keeps result_gold_pbr.png)
    int pmxIndex = -1;                     // --pmx-index N (see --pmx-list); -1 = the largest .pmx
    bool pmxList = false;                  // --pmx-list : print the .pmx files found in MMD_model/ and exit
    Real pmxZoom = 1;                      // --pmx-zoom Z : camera distance / Z
    Real pmxFocus = 0.5;                   // --pmx-focus f : look at this fraction of the model height (0.9 = face)
    int poseIndex = -1;                    // --pose-index N (see --pose-list): apply a .vpd pose from MMD_model/pose
    bool poseList = false;                 // --pose-list : print the .vpd files found and exit
    bool morphList = false;                // --morph-list : print the model's morphs and exit
    std::string morphSpec;                 // --morph "12=1.0,15=0.5" : morph weights by index
    std::string exprFile;                  // --expr file.txt : "name = weight" lines (UTF-8) from MMD_model/expression
    bool pmxCheck = false;                 // --pmx-check : skeleton summary + pose self-tests, then exit
    bool pmxBench = false;                 // --pmx-bench : time the stages of one pose update (morphs, pose + skinning, mesh, BVH), then exit
    std::string studioLook = "default";    // --studio-look default|warm|cool|dramatic|soft : lighting of --pmx-studio
    bool aov = false;                      // --aov : also write <out>_albedo.png and <out>_normal.png (denoiser features)
    bool denoise = false;                  // --denoise : also write <out>_denoised.* (Intel OIDN, uses albedo + normal)
    std::string gamut = "srgb";            // --gamut srgb|p3|rec2020 : primaries of the written image / .pfm / .hdr (the film stays linear sRGB)
    Real fogPrismDeg = 21.0;               // --fog-prism-deg : rotation of the prism in the fogbeam scene (degrees)
    Real fogFill = 0.4;                    // --fog-fill : intensity of the faint fill light in the fogbeam scene (0 = none)
    Real fogSigma = 1.5;                   // --fog-sigma : extinction of the fog in the fogbeam scene (1/m)
    Real fogLight = 60.0;                  // --fog-light : intensity of its lamp
    bool lvc = false;                      // --lvc : bidirectional path tracing with a cache of light paths (VolumeBPT) in the fogspot scene
    long long lvcPaths = 2000000;          // --lvc-paths N : light paths in the cache
    int lvcPasses = 1;                     // --lvc-passes K : K passes with one wavelength each (needed with dispersive glass)
    int lvcDepth = 10;                      // --lvc-depth N : scattering vertices per sub-path
    int lvcMask = 7;                       // --lvc-mask M : strategies, 1 = next event estimation, 2 = connections, 4 = light tracing
    std::string galObject = "ball";        // --gal-object ball|prism|column|tumbler|gem : the glass object of the caustics scene
    Real galLight = 3.0;                   // --gal-light : intensity of its lamp
    Real galFog = 0.2;                     // --gal-fog : extinction of its haze (1/m)
    Real galFov = 11.0;
    Vector3 galLook = Vector3(0.0, -999.0, 0.0);   // --gal-look x,y,z : where the camera looks (y = -999: the object's default)
    int galFloor = 1;                      // --gal-floor 0 : no floor and no wall: the glass in the haze only (use --gal-fog 1 or more)
    double sssrwAlbedo = 0.99;             // --sssrw-albedo : single-scattering albedo of the medium in --scene sssrw
    double sssrwG = 0.0;                   // --sssrw-g : its Henyey-Greenstein anisotropy
    double sssrwEll = 0.002;               // --sssrw-ell : the reduced mean free path (m) = the scatter distance of the BSSRDF
    double sssrwA = 0.0;                   // --sssrw-A : the diffuse color of the BSSRDF (0: Monte Carlo of the semi-infinite medium)
    double sssrwLampDeg = 0.0;             // --sssrw-lamp-deg : the lamp's angle from the axis of the discs (0: on the axis, in line with the camera)
    double sssrwSphere = 5.0;              // --sssrw-sphere : the radius of the sphere in scatter distances
    Real galScale = 1.0;                   // --gal-scale : size of the object (the wall moves with it)
    Real galDist = 0.0;                    // --gal-dist : distance of the lamp and the camera from the object (0: the scene's default)
    Real galAlbedo = 0.3;                  // --gal-albedo : albedo of the floor and the wall
    Real galSpread = 11.0;                 // --gal-spread : half angle of the lamp's cone (degrees)                    // --gal-fov : vertical field of view (degrees)
    Real galAz = 999, galEl = 999;         // --gal-az, --gal-el : the lamp's direction (degrees; 999 = the object's default)
    bool spotPrism = false;                // --spot-prism : a dispersive prism in the fogspot scene (the lamp aims at it); use --lvc-passes
    bool spotGlass = false;                // --spot-glass : a glass ball in the fogspot scene (the lamp aims at it)
    bool spotFloor = false;                // --spot-floor : a diffuse floor in the fogspot scene (the lamp aims at it)
    Real spotDeg = 3.0;                    // --spot-deg : half angle of the lamp's cone in the fogspot scene
    bool causticTraceOnly = false;         // --caustic-trace-only : write only the light-traced image (to look at it alone)
    bool causticTrace = false;             // --caustic-trace : light tracing for the caustics on surfaces the camera sees directly (unbiased, sharp edges); the photons carry the rest
    int pngBits = 8;                       // --bits 8|16 : bits per channel of a written .png (16: no banding in a wide gamut)
    bool cleanAux = false;                // --clean-aux : tell OIDN the albedo / normal are noise-free (skips its prefilter)
    bool pmxStudio = false;                // --pmx-studio : studio lighting (key/fill/rim soft boxes, dim backdrop), framing by pose width
    bool pmxLight = false;                 // --pmx-light : add a white softbox in front-left (highlights show the gloss)
    bool pmxMatte = false;                 // --pmx-matte : every material plain diffuse (no gloss presets)
    bool sss = false;                      // --sss : subsurface scattering for the PMX skin materials
    Real sssScale = 1.0;                   // --sss-scale S : multiplies the skin scatter distance (3.7 / 1.4 / 0.7 mm); artistic control
    bool pmxAll = false;                   // --pmx-all : also draw the "+" overlay materials
    bool pmxInfoOnly = false;              // --pmx-info : load, print the model summary, exit (no render)
    Real pmxYaw = 0;                       // --pmx-yaw deg : turn the model about +Y
    bool noMs = false;                     // --no-ms : single-scattering GGX only (goldrough scene)
    std::string accel = "own";             // own | embree (triangle mesh intersection)
    std::string lightSampler = "spatial";  // uniform | power | spatial (light selection for NEE)
    std::string light = "tungsten";    // spd scene: d65 | tungsten | candle | led
    std::string metal = "Au";          // metalprobe scene: Au | Ag | Cu (Johnson tables)
    Real probeRoughness = 0.01;        // metalprobe scene: --roughness r of the metal sphere
    Real melanin = 0.04;               // skinlayers: melanosome volume fraction of the epidermis
    Real bloodFraction = 0.02;         // skinlayers: blood volume fraction of the dermis
    Real oxygen = 0.75;                // skinlayers: oxygen saturation of the blood
    Real epidermisMm = 0.1;            // skinlayers: epidermis thickness (mm)
    long long photons = 2000000;       // skinlayers: photons per wavelength
    bool singleLayer = false;          // skinlayers: dermis only (homogeneous), for comparison
    std::string skinOut = "assets/Skin/skin_layers.csv";   // skinlayers: output table
    std::string skinTable;             // skindemo: layered profile table (--skin-table file.csv)
    Real sigmaT = 128;                 // blob: extinction coefficient of the medium (1/m)
    Real mediumAlbedo = 0.99;          // blob: single-scattering albedo of the medium
    Real phaseG = 0;                   // blob: Henyey-Greenstein anisotropy (0 = isotropic)
    Real phaseG2 = 0;                  // blob: second (backward) lobe of the phase function (--phase-g2); used with --phase-w
    Real phaseW = 1;                   // blob: weight of the first lobe (--phase-w; 1 = a single lobe)
    std::optional<Vector3> sigmaTRgb;  // blob: chromatic extinction (1/m) per channel (--sigma-t-rgb); with --albedo-rgb, a ChromaticMedium
    std::optional<Vector3> albedoRgb;  // blob: chromatic single-scattering albedo (--albedo-rgb)
    Real glassRoughness = 0.224;       // blob: perceptual roughness of the glass boundary
    std::string mesh;                  // blob: closed OBJ mesh with outward normals (default: an icosphere)
    Real bgGray = -1;                  // blob: gray backdrop for camera rays that touch nothing (< 0: the environment)
    int volumeCap = 100000;            // blob: limit of collisions inside the medium per path
    bool mediaNee = false;             // --media-nee: sample the lights at every collision in a medium, and cross index-matched boundaries (clouds, fog; scene files)
    Real absorbScale = -1;             // water / sea: multiplies the absorption of pure water (1 = physical); < 0: the scene's default (water 40, sea 1)
    Real seaDepth = 2.5;               // sea: depth of the water above the sand (m)
    Real seaSteep = 0.03;              // sea: steepness k*A of every wave component (about 0.1 rms slope in total)
    unsigned seaSeed = 1;              // sea: random seed of the wave set and the rocks
    Real poolDepth = 0.9;              // pool: water depth (m)
    Real poolSteep = 0.025;            // pool: steepness k*A of every ripple (about 0.07 rms slope in total)
    Real poolSunElev = 50;             // pool: elevation of the sun (degrees)
    Real poolSunAz = 25;               // pool: azimuth of the sun about +y, 0 = +z (the camera stands at +z looking toward -z)
    bool seaHdr = false;               // sea: light with the --env map instead of the procedural sky
    Real seaSunElev = 14;              // sea: elevation of the procedural sky's sun (degrees)
    bool waterSky = false;             // water: light the tumbler with the environment map (--env) instead of softboxes; the wall is removed
    // manifold / crystal scenes: how the caustic of the point light is computed: manifold (specular manifold NEE for the T T chains, plus the photon map
    // for the other chains when --caustics N > 0) | photons | none. Empty: the scene's default (manifold: manifold; crystal: photons).
    std::string causticMethod;
    bool crystalInside = false;       // crystal: the lamp sits inside the tumbler (no water, no other light): beams radiate through the facets
    bool waterLamp = false;           // water: one small point light at about 33 degrees elevation (a sharp caustic; needs --caustics), no wall, camera looks down
    long long causticPhotons = -1;    // water / pool: photons emitted for the caustic photon map (water: off unless > 0; pool: 8 M unless 0 = off)
    int causticPasses = 1;             // water / pool / crystal: render in P passes with independent photon sets (spp / P each), images averaged
    Real causticRadius = -1;          // gather radius of the caustic photon map (m); < 0: the scene's default (water 1.5 mm, caustic test 4 mm)
    Real camElev = -1;                 // water --water-sky: camera elevation above the table in degrees (< 0: 38)
    Real camYaw = 1e9;                 // water --water-sky: camera azimuth in degrees about +y, 0 = on the +z side (1e9: from the brightest spot of the map)
    bool skinMatched = false;          // skindemo: Burley profile with the same albedo and mean radius as the table
    bool whiteBalance = false;         // --wb: white-balance the film to the light (spd scene; scene files: the first emission material)
    std::optional<Real> fadeHours;     // --fade-hours H: the hours of exposure of every "fading" material of a scene file (overrides the file)
    int traceX = -1, traceY = -1;      // --trace x,y : print a few light paths of one pixel and exit
    // The lens: --focal-length MM (replaces the field of view of the scene), --fstop N (defocus blur), --focus-dist D | --focus-at x,y,z,
    // --blades N (polygonal aperture), --blade-rotation DEG. Whatever is given overrides the scene's own lens (see Lens.hpp).
    LensSpec lens;
};

// Parses the command line, and the scene file if the first argument is one (the command line wins over the file).
// Returns an exit code when the program should stop (0 after --help / --list-scenes, 1 after an error), std::nullopt to carry on.
std::optional<int> ParseArgs(int argc, char** argv, CliOptions& o);

} // namespace rayt::cli
