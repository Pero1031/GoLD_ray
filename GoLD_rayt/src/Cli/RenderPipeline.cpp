// The common steps of rendering a scene: standard materials and environment, the scene builder, lights, caustics, camera, integrator,
// render, save, denoise. What differs between the scenes is in src/Cli/Scenes.
#include "pch.h"
#include "Cli/RenderPipeline.hpp"
#include "Cli/Scene.hpp"
#include "Cli/Tools.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Denoiser.hpp"
#include "Render/Integrator.hpp"
#include "Render/MaskDepth.hpp"
#include "Render/SpecularNee.hpp"
#include "Render/VolumeBPT.hpp"

#include <algorithm>
#include <sstream>

namespace rayt::cli {

namespace {

// The lens as the Camera takes it. The scene gives a field of view (or a lens); the command line's lens replaces it field by field.
struct ResolvedLens {
    Real orthoHeight = 0;     // > 0: the orthographic view, this high (scene units)
    Real fovDeg = 35;
    Real focalLengthMm = 0;
    Real fStop = 0;           // 0: a pinhole
    Real aperture = 0;        // diameter in scene units (metres)
    Real focusDist = 1;       // along the view axis
    int blades = 0;
    Real bladeRotationDeg = 0;
    Real k1 = 0, k2 = 0;      // radial distortion
    Real lateralCA = 0;       // chromatic aberrations (per 100 nm)
    Real axialCA = 0;
    Real barrelRadius = 0;    // vignetting: the opening of the lens barrel (metres; 0 = none), its distance in front of the stop (metres), the cos^4 share
    Real barrelDist = 0;
    Real naturalVignetting = 0;
    Real shutter = 0;         // motion blur: seconds the shutter is open (0 = off), the velocity of the camera (m/s) and its pan rates (degrees per second)
    Vector3 velocity = Vector3(0.0);
    Real panYaw = 0, panPitch = 0;
};

ResolvedLens ResolveLens(const CameraSpec& cam, const LensSpec& fromOptions, Real aspect) {
    LensSpec s = cam.lens;
    s.overrideWith(fromOptions);

    ResolvedLens r;
    // The 35 mm convention: the long side of the picture is the 36 mm side of the film. The vertical field of view follows from it.
    const Real sensorHeightMm = (aspect >= Real(1)) ? Real(36) / aspect : Real(36);
    if (s.focalLengthMm) {
        r.focalLengthMm = *s.focalLengthMm;
        r.fovDeg = Real(2) * std::atan(sensorHeightMm / (Real(2) * r.focalLengthMm)) * constants::RAD_TO_DEG;
    }
    else {
        r.fovDeg = cam.fovDeg;
        r.focalLengthMm = sensorHeightMm / (Real(2) * std::tan(r.fovDeg * constants::DEG_TO_RAD * Real(0.5)));
    }

    // The plane of focus: the point or the distance given ("focusAt", "focusDistance", --focus-at, --focus-dist), the scene's "focus",
    // or the point the camera looks at.
    const Vector3 viewDir = glm::normalize(cam.at - cam.from);
    r.focusDist = cam.focus > 0 ? cam.focus : glm::length(cam.from - cam.at);
    if (s.focusDist) r.focusDist = *s.focusDist;
    else if (s.focusAt) r.focusDist = glm::dot(*s.focusAt - cam.from, viewDir);
    if (!(r.focusDist > 0)) {
        std::cout << "[Camera] the point to focus on is behind the camera: focusing on the point it looks at\n";
        r.focusDist = glm::length(cam.from - cam.at);
    }

    // The aperture: an f-number ("fStop", --fstop; 0 = a pinhole) wins over the scene's lens diameter ("aperture"). Scenes are in metres.
    if (s.fStop) r.aperture = (*s.fStop > 0) ? r.focalLengthMm * Real(0.001) / *s.fStop : Real(0);
    else r.aperture = cam.aperture;
    if (r.aperture > 0) r.fStop = r.focalLengthMm * Real(0.001) / r.aperture;
    r.blades = s.blades.value_or(0);
    r.bladeRotationDeg = s.bladeRotationDeg.value_or(Real(0));
    r.k1 = s.distortionK1.value_or(Real(0));
    r.k2 = s.distortionK2.value_or(Real(0));
    r.lateralCA = s.lateralCA.value_or(Real(0));
    r.axialCA = s.axialCA.value_or(Real(0));
    // The lens barrel: an opening of diameter focal length / N, by default half a focal length in front of the stop
    if (s.vignetteFStop && *s.vignetteFStop > 0) {
        r.barrelRadius = r.focalLengthMm * Real(0.001) / (Real(2) * *s.vignetteFStop);
        r.barrelDist = s.vignetteDistMm.value_or(r.focalLengthMm * Real(0.5)) * Real(0.001);
    }
    r.naturalVignetting = s.naturalVignetting.value_or(Real(0));
    r.orthoHeight = s.orthoHeight.value_or(Real(0));
    r.shutter = s.shutter.value_or(Real(0));
    r.velocity = s.velocity.value_or(Vector3(0.0));
    r.panYaw = s.panYawRate.value_or(Real(0));
    r.panPitch = s.panPitchRate.value_or(Real(0));
    return r;
}

// --mask: the materials named in opt.mask (the names of the scene file). False, with a message, for a built-in scene or an unknown name.
bool ResolveMaskMaterials(const CliOptions& opt, const SceneSetup& setup, std::vector<const Material*>& out) {
    if (opt.mask.empty()) return true;
    if (setup.namedMaterials.empty()) {
        std::cerr << "[Mask] --mask names materials of a scene file; the built-in scene \"" << opt.scene << "\" has no named materials\n";
        return false;
    }
    std::stringstream list(opt.mask);
    for (std::string name; std::getline(list, name, ',');) {
        if (name.empty()) continue;
        const auto it = setup.namedMaterials.find(name);
        if (it == setup.namedMaterials.end()) {
            std::cerr << "[Mask] --mask: the scene file has no material \"" << name << "\"; it has:";
            for (const auto& [n, m] : setup.namedMaterials) std::cerr << " " << n;
            std::cerr << "\n";
            return false;
        }
        out.push_back(it->second);
    }
    return true;
}

// The camera of the scene with the lens of the command line, and the notes on what the chosen render mode does not show of it.
// lightTracing: light paths are connected to a pinhole camera (--caustic-trace, --lvc); lvcCamera: --lvc renders the picture.
std::shared_ptr<Camera> BuildCamera(const CliOptions& opt, const SceneSetup& setup, int imageWidth, int imageHeight, bool lightTracing, bool lvcCamera) {
    const Point3 lookFrom = setup.camera.from;
    const Point3 lookAt = setup.camera.at;
    const Vector3 vUp = setup.camera.up;

    const ResolvedLens lens = ResolveLens(setup.camera, opt.lens, double(imageWidth) / imageHeight);
    if (lens.orthoHeight > 0) {
        if (lightTracing) {
            std::cerr << "[Camera] the orthographic view cannot be used with light tracing (--caustic-trace, --lvc, or a scene that light-traces its caustics):"
                " it weighs its paths for a pinhole camera\n";
            return nullptr;
        }
        std::cout << "[Camera] orthographic: a view " << lens.orthoHeight * imageWidth / imageHeight << " x " << lens.orthoHeight
            << " (scene units) wide and high, rays parallel to the view axis";
        if (lens.aperture > 0) std::cout << "; f/" << lens.fStop << " = aperture " << lens.aperture * 1000.0 << " mm, focus at " << lens.focusDist << " m";
        std::cout << "\n";
        if (lens.barrelRadius > 0 || lens.naturalVignetting > 0) std::cout << "[Camera] note: the vignetting has no effect in the orthographic view (every chief ray is on the axis)\n";
    }
    else if (lens.aperture > 0 || opt.lens.focalLengthMm || setup.camera.lens.focalLengthMm) {
        std::cout << "[Camera] focal length " << lens.focalLengthMm << " mm (vertical field of view " << lens.fovDeg << " deg)";
        if (lens.aperture > 0) {
            std::cout << ", f/" << lens.fStop << " = aperture " << lens.aperture * 1000.0 << " mm, focus at " << lens.focusDist << " m, "
                << (lens.blades >= 3 ? std::to_string(lens.blades) + " blades" : std::string("round aperture"));
        }
        std::cout << std::endl;
    }

    auto camera = std::make_shared<Camera>(
        lookFrom, lookAt, vUp,
        lens.fovDeg,
        double(imageWidth) / imageHeight,
        lens.aperture,
        lens.focusDist
    );
    if (lens.blades >= 3) camera->setApertureShape(lens.blades, lens.bladeRotationDeg * constants::DEG_TO_RAD);
    camera->setAberrations(lens.k1, lens.k2, lens.lateralCA, lens.axialCA, lens.focalLengthMm * Real(0.001));
    if (lens.k1 != 0 || lens.k2 != 0) {
        std::cout << "[Camera] radial distortion k1 = " << lens.k1 << ", k2 = " << lens.k2 << (lens.k1 > 0 ? " (barrel)" : lens.k1 < 0 ? " (pincushion)" : "") << "\n";
        if (camera->distortionFolds()) std::cout << "[Camera] warning: this distortion folds the picture over itself near the edge; use a smaller k1 / k2\n";
    }
    if (lens.lateralCA != 0 || lens.axialCA != 0) {
#if defined(RAYT_SPECTRAL_MODE)
        std::cout << "[Camera] chromatic aberration: lateral " << lens.lateralCA << ", axial " << lens.axialCA
            << " per 100 nm (one wavelength per path: the color noise is higher)\n";
        if (lens.axialCA != 0 && !(lens.aperture > 0)) std::cout << "[Camera] note: axial chromatic aberration shows as colored defocus: it needs --fstop\n";
        if (lvcCamera) std::cout << "[Camera] note: --lvc connects to a pinhole camera of one wavelength: no chromatic aberration there\n";
#else
        std::cout << "[Camera] warning: chromatic aberration needs the wavelength of the path: it is ignored in the RGB build (use the spectral build)\n";
#endif
    }
    camera->setOrthographic(lens.orthoHeight);
    camera->setShutter(lens.shutter);
    camera->setMotion(lens.velocity, lens.panYaw, lens.panPitch);
    if (camera->hasShutter()) {
        std::cout << "[Camera] shutter " << lens.shutter << " s";
        if (lens.shutter < 1.0) std::cout << " (1/" << std::lround(1.0 / lens.shutter) << ")";
        if (camera->moving()) std::cout << "; the camera moves " << lens.velocity.x << ", " << lens.velocity.y << ", " << lens.velocity.z << " m/s, pans " << lens.panYaw << " / " << lens.panPitch << " deg/s";
        std::cout << "\n";
        if (lvcCamera) std::cout << "[Camera] note: --lvc connects to a camera at time 0: no motion blur there\n";
    }
    else if (camera->moving()) std::cout << "[Camera] note: the camera has a velocity or a pan but the shutter is closed (--shutter, \"shutter\"): it does not move\n";
    camera->setVignetting(lens.barrelRadius, lens.barrelDist, lens.naturalVignetting);
    if (lens.barrelRadius > 0 && !camera->orthographic()) {
        if (lens.aperture > 0) {
            std::cout << "[Camera] vignetting: lens barrel opening f/" << lens.focalLengthMm * 0.001 / (2 * lens.barrelRadius) << " (diameter " << 2000.0 * lens.barrelRadius
                << " mm) " << 1000.0 * lens.barrelDist << " mm in front of the stop; of the round aperture, " << std::lround(100.0 * camera->apertureFraction(1, 1))
                << " % is left open at the corner (about)\n";
            if (2 * lens.barrelRadius <= lens.aperture) std::cout << "[Camera] note: the barrel opening is not larger than the aperture: it acts as the stop. A smaller --vignette N opens it\n";
        }
        else std::cout << "[Camera] note: the barrel vignetting needs an aperture (--fstop): a pinhole is not clipped\n";
    }
    if (lens.naturalVignetting > 0 && !camera->orthographic()) std::cout << "[Camera] natural vignetting: the film is lit in proportion to cos^" << 4 * lens.naturalVignetting << " of the angle off the axis\n";
    if ((lens.k1 != 0 || lens.k2 != 0 || lens.lateralCA != 0) && lightTracing)
        std::cout << "[Camera] note: light tracing (--caustic-trace, --lvc) connects photons to the undistorted pinhole picture: its caustics are not distorted\n";
    if (lens.aperture > 0 && lightTracing)
        std::cout << "[Camera] note: light tracing (--caustic-trace, --lvc) connects photons to a pinhole camera: its caustics stay sharp while the rest is defocused\n";
    return camera;
}

// --spectral-image: the bands of the film's spectral image, and how many wavelength samples the poorest band gets per pixel (its noise
// goes as one over the square root of that): spp x the wavelengths of a path x the chance that one falls into the band.
void EnableSpectralImage(const CliOptions& opt, Film& film) {
    const double first = opt.spectralFirst, step = opt.spectralStep;
    const int bands = static_cast<int>(std::lround((opt.spectralLast - first) / step));
    film.enableSpectral(Real(first), Real(step), bands);
    double fewest = 1e300, fewestAt = first;
    for (int k = 0; k < bands; ++k) {
        double chance = 0;
        for (int s = 0; s < 16; ++s) chance += visibleWavelengthPdf(Real(first + (k + (s + 0.5) / 16) * step)) * step / 16;
        const double samples = double(opt.spp) * kSpectralSamples * chance;
        if (samples < fewest) { fewest = samples; fewestAt = first + (k + 0.5) * step; }
    }
    std::cout << "[Spectral] spectral image: " << bands << " bands of " << step << " nm from " << first << " to " << opt.spectralLast
        << " nm; the poorest band (at " << fewestAt << " nm) gets " << fewest << " wavelength samples per pixel\n";
}

// Writes the picture and everything asked for beside it (mask, depth, spectral image, AOVs, Stokes images, polarizer views, the denoised copy).
// photonsAdded: caustic photons or light paths were added to the picture (they carry no polarization).
void SaveOutputs(const CliOptions& opt, Film& film, bool photonsAdded) {
    std::cout << "[Output] Saving images..." << std::endl;
    color::Gamut gamut = color::Gamut::SRGB;
    color::parseGamut(opt.gamut, gamut);   // already validated by the option parser
    // Every written picture (the image, the polarizer images, the denoised copy) gets the same output settings.
    auto setOutput = [&](Film& f) {
        f.setOutputGamut(gamut);
        f.setPng16(opt.pngBits == 16);
        f.setDisplay(opt.exposure, opt.tonemap == "clip" ? Film::ToneMap::Clip : Film::ToneMap::Reinhard);
    };
    setOutput(film);
    film.save(opt.out);
    const std::string outStem = std::filesystem::path(opt.out).replace_extension("").string();
    const std::string outExt = std::filesystem::path(opt.out).extension().string();
    std::string outExtLower = outExt;
    std::transform(outExtLower.begin(), outExtLower.end(), outExtLower.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    const bool floatOut = (outExtLower == ".pfm" || outExtLower == ".hdr");
    if (floatOut && (opt.exposure != 0 || opt.tonemap != "reinhard"))
        std::cout << "[Output] note: --exposure and --tonemap change only .png / .bmp / .jpg; " << outExtLower << " holds the radiance as rendered\n";
    if (film.hasMask()) film.saveMask(outStem + (floatOut ? "_mask.pfm" : "_mask.png"));
    if (film.hasDepth()) film.saveDepth(outStem + "_depth.pfm");
    if (film.hasSpectral()) film.saveSpectral(outStem);
    if (opt.aov) {
        const std::string aovExt = (outExt == ".hdr") ? ".hdr" : ".png";
        film.saveAlbedo(outStem + "_albedo" + aovExt);
        film.saveNormal(outStem + "_normal" + aovExt);
    }
    if (opt.stokes) {
        film.saveStokes(outStem);
        film.saveLinearPolarization(outStem);
        for (const Real angle : opt.polarizers) {
            Film filtered(film.width(), film.height());
            setOutput(filtered);
            film.polarizerImage(angle, filtered);
            std::ostringstream name;
            name << outStem << "_pol" << angle;
            std::string stem = name.str();
            std::replace(stem.begin() + static_cast<std::ptrdiff_t>(outStem.size()), stem.end(), '.', 'p');   // 22.5 -> 22p5: keeps the extension clear
            filtered.save(stem + outExt);
        }
        if (photonsAdded)
            std::cout << "[Polarization] note: the caustic photons and the light-traced image are added without polarization" << std::endl;
    }
    if (opt.denoise) {
        Film denoised(film.width(), film.height());
        setOutput(denoised);
        if (DenoiseFilm(film, denoised, opt.cleanAux, gamut)) denoised.save(outStem + "_denoised" + outExt);
    }
}

} // namespace

ImageSize ComputeImageSize(const CliOptions& opt) {
    if (opt.height > 0) return { opt.width, opt.height };
    const SceneDef* def = FindScene(opt.scene);
    const double heightPerWidth = def ? def->heightPerWidth : double(IMAGE_HEIGHT) / IMAGE_WIDTH;
    return { opt.width, std::max(1, static_cast<int>(std::lround(double(opt.width) * heightPerWidth))) };
}

int RenderScene(CliOptions& opt) {
    const SceneDef* def = FindScene(opt.scene);
    if (!def) { std::cerr << "[Scene] unknown scene: " << opt.scene << "\n"; return 1; }
    const ImageSize size = ComputeImageSize(opt);
    const int imageWidth = size.width;
    const int imageHeight = size.height;

    // ---- the environment and the shared materials, then the scene itself ----
    std::cout << "CWD = " << std::filesystem::current_path() << std::endl;
    Environment env = MakeEnvironment(opt, def->envKind(opt));
    const StdMaterials mats = MakeStandardMaterials(opt);
    SceneSetup setup;
    SceneContext ctx{ opt, mats, env };
    if (const std::optional<int> exitCode = def->build(ctx, setup)) return *exitCode;
    const CausticSetup& caustics = setup.caustics;       // what the scene asks of the photon map and the manifold NEE

    // --mask: the materials whose objects the mask covers; checked before the render starts
    std::vector<const Material*> maskMaterials;
    if (!ResolveMaskMaterials(opt, setup, maskMaterials)) return 1;

    if (caustics.lightTraceOnly) opt.causticTrace = true;

    // ---- the scene: lights, caustics ----
    // The builders collect the objects in a plain list (spheres, quads, whole meshes with their own BVH, hair). Rays go through
    // one BVH over all of them, so a scene of many small objects does not test every object for every ray.
    const auto worldBvh = std::make_shared<BVH>(setup.world->objects, 4, /*recordInPlace=*/true);
    std::cout << "[Scene] " << worldBvh->numPrimitives() << " objects in a BVH of " << worldBvh->numNodes() << " nodes.\n";
    auto scene = std::make_shared<Scene>(worldBvh);

    if (setup.background) scene->setCameraBackground(*setup.background);
    if (setup.splats) {
        scene->setSplats(setup.splats);
        std::cout << "[Scene] Gaussian splats: " << setup.splats->count() << " particles in a BVH of " << setup.splats->numNodes() << " nodes.\n";
    }

    // the environment light, if there is one
    if (env.light) {
        scene->setEnvLight(env.light);
        std::cout << "[Scene] Environment light set.\n";
    }

    for (const auto& l : setup.lights) scene->addLight(l);

    // IMPORTANT: finalize the scene (builds the light sampler)
    const std::string& lightSampler = setup.lightSampler.empty() ? opt.lightSampler : setup.lightSampler;
    scene->finalize(lightSampler);
    std::cout << "[Scene] Finalized (light sampler: " << lightSampler << ") with "
        << scene->numLights() << " scene lights";
    if (scene->hasEnvLight()) {
        std::cout << " + environment light";
    }
    std::cout << ".\n";

    // Caustic photon map (water / pool scene, --caustics N): photons from the lights through the water onto the table / the basin.
    // Caustic method of the point-light scenes: manifold = the T T chains of the flat-faced glass bodies are solved exactly (specular manifold NEE),
    // the photon map (if --caustics N > 0) carries the other chains only; photons = the photon map alone; none = direct light only.
    // The scene's CausticSetup holds its defaults (the prism scene uses manifold alone, the crystal scene the photon map) and says
    // when the manifold NEE cannot be used.
    std::string causticMethod = opt.causticMethod.empty() ? caustics.defaultMethod : opt.causticMethod;
    if (causticMethod == "manifold" && !caustics.manifoldUnavailable.empty()) {
        std::cout << "[Caustics] " << caustics.manifoldUnavailable << "; using the photon map\n";
        causticMethod = "photons";
    }
    const bool useManifold = (causticMethod == "manifold") && !setup.manifoldBodies.empty();   // only scenes with flat-faced glass bodies list them
    std::unique_ptr<SpecularNee> specularNee;
    if (useManifold) {
        specularNee = std::make_unique<SpecularNee>(setup.manifoldBodies, *scene);
        std::cout << "[Manifold] specular manifold NEE: " << setup.manifoldBodies.size() << " bodies, " << specularNee->numFacets() << " facets, "
            << specularNee->numLights() << " point light(s)\n";
    }

    // ---- the camera ----
    const bool lvcCamera = opt.lvc && setup.bptLight;   // --lvc connects to a pinhole camera of its own
    const auto camera = BuildCamera(opt, setup, imageWidth, imageHeight, opt.causticTrace || caustics.lightTraceOnly || lvcCamera, lvcCamera);
    if (!camera) return 1;

    // --spectral-image comes from the path tracer's own samples: refused where other parts of the picture are added in RGB (checked
    // before any photon is traced)
    if (opt.spectralStep > 0) {
        const char* why = (opt.lvc && setup.bptLight && setup.bptMedium) ? "--lvc renders the picture with a film of its own"
            : opt.causticTrace ? "the light-traced caustics (--caustic-trace, or a scene that light-traces them) are added in RGB"
            : opt.causticPasses > 1 ? "--caustic-passes averages RGB pictures" : nullptr;
        if (why) { std::cerr << "[Spectral] --spectral-image cannot be used here: " << why << "\n"; return 1; }
    }

    CausticPhotonMap causticMap;
    CausticPhotonParams causticParams;   // kept for --caustic-passes: every pass traces a fresh photon set with the same parameters
    bool causticOn = false;
    bool traceOn = false;                  // --caustic-trace: light tracing of the caustics on directly visible surfaces
    std::vector<Vector3> traceSum;         // the light-traced images of all photon sets (summed)
    // The number of photons: --caustics N, or the scene's default (a smaller one while the manifold NEE solves the T T chains).
    long long nCausticPhotons = opt.causticPhotons;
    const long long defaultPhotons = (useManifold && caustics.defaultPhotonsWithManifold >= 0) ? caustics.defaultPhotonsWithManifold : caustics.defaultPhotons;
    if (nCausticPhotons < 0 || (nCausticPhotons == 0 && causticMethod == "photons" && caustics.photonMethodNeedsPhotons)) nCausticPhotons = defaultPhotons;
    if (causticMethod == "none") nCausticPhotons = 0;
    if (caustics.photonMap && nCausticPhotons > 0) {
        causticParams.numPhotons = nCausticPhotons;
        causticParams.radius = opt.causticRadius > 0 ? opt.causticRadius : caustics.defaultRadius;
        causticParams.targetCenter = setup.causticTarget.center;
        causticParams.targetRadius = setup.causticTarget.radius;
        if (specularNee) causticParams.skipTransmissionPairs = specularNee->materials();   // their T T chains are solved: do not store them twice
        if (opt.causticTrace) {   // also connect every caustic photon to the camera (pinhole): unbiased, sharp, on directly visible surfaces
            causticParams.lightTraceCamera = camera.get();
            causticParams.traceWidth = imageWidth;
            causticParams.traceHeight = imageHeight;
            traceOn = true;
            causticParams.traceMedia = caustics.throughMedia;   // the photons go on through the fog and are connected to the camera there
        }
        causticMap.build(*scene, causticParams);
        causticOn = true;
        if (traceOn) traceSum = causticMap.splat();
    }

    if (opt.traceX >= 0) {
        TracePixelDebug(*scene, *camera, opt.traceX, opt.traceY, imageWidth, imageHeight, 12, 24,
            { { "clear", mats.glass.get() }, { "disp", mats.glassDisp.get() }, { "frosted", mats.frosted.get() },
              { "white", mats.white.get() }, { "floor", mats.floor.get() },
              { "goldS", mats.goldSmooth.get() }, { "goldM", mats.goldMedium.get() }, { "goldR", mats.goldRough.get() } });
        return 0;
    }

    // ---- the integrator ----
    Film film(imageWidth, imageHeight);
    if (opt.aov || opt.denoise) film.enableAOVs();   // the denoiser uses them even if they are not written out
    if (opt.stokes) film.enableStokes();             // --stokes / --polarizer: the Stokes vector of every pixel (docs/POLARIZATION.md)
    if (opt.spectralStep > 0) EnableSpectralImage(opt, film);

    // the scene's default path depth, unless --depth was given
    if (def->defaultMaxDepth > 0 && opt.maxDepth == MAX_DEPTH) opt.maxDepth = def->defaultMaxDepth;
    // --caustic-passes P: P independent photon sets, spp / P samples per pass, the images averaged (see the render below)
    const int causticPasses = (causticOn && opt.causticPasses > 1) ? opt.causticPasses : 1;
    const int sppPerPass = (causticPasses > 1) ? std::max(1, opt.spp / causticPasses) : opt.spp;
    auto integrator = std::make_unique<PathIntegrator>(
        camera,
        opt.maxDepth,
        sppPerPass
    );

    integrator->setMaxVolumeEvents(opt.volumeCap);
    integrator->setCausticMap(&causticMap);   // ignored while the map is empty
    if (traceOn) {
        integrator->setSkipDirectCaustics(true);
        std::cout << "[Caustics] light tracing: the caustics on surfaces the camera sees directly are splatted from the photon paths (no gather there)\n";
    }
    else if (opt.causticTrace) std::cout << "[Caustics] --caustic-trace needs the photon map (--caustics N > 0, caustic method photons): ignored\n";
    integrator->setSpecularNee(specularNee.get());
    if (setup.mediaNee || opt.mediaNee) integrator->setMediaNee(true);   // the fog is lit by next event estimation, its boundary does not block shadow rays
    std::cout << "[Integrator] PathIntegrator created (maxDepth="
        << opt.maxDepth << ", spp=" << opt.spp << ").\n";
    if (causticPasses > 1)
        std::cout << "[Caustics] " << causticPasses << " passes of " << sppPerPass << " spp, each with its own " << causticParams.numPhotons
            << " photons (" << double(causticPasses) * double(causticParams.numPhotons) / 1e6 << " M in total, memory of one set)\n";

    // ---- render and save ----
    std::cout << "[Render] Start PBR rendering..." << std::endl;
    if (causticPasses > 1) {
        // The noise of a photon map is frozen: more samples per pixel re-read the same photons and do not remove it. Independent
        // photon sets per pass decorrelate it, so the average of the passes improves by 1/sqrt(passes) at the memory of one set.
        // The last pass renders into `film` (so the denoiser features come from it); the earlier ones go to a scratch film.
        std::vector<Vector3> sum(static_cast<size_t>(imageWidth) * static_cast<size_t>(imageHeight), Vector3(0.0));
        for (int pass = 0; pass < causticPasses; ++pass) {
            if (pass > 0) {
                causticMap.build(*scene, causticParams);   // a fresh photon set (the random number streams continue)
                if (traceOn) for (size_t i = 0; i < traceSum.size(); ++i) traceSum[i] += causticMap.splat()[i];
            }
            std::cout << "[Caustics] pass " << pass + 1 << " / " << causticPasses << "\n";
            if (pass + 1 < causticPasses) {
                Film scratch(imageWidth, imageHeight);
                integrator->render(*scene, scratch);
                for (size_t i = 0; i < sum.size(); ++i) sum[i] += scratch.getData()[i];
            }
            else integrator->render(*scene, film);
        }
        for (int y = 0; y < imageHeight; ++y)
            for (int x = 0; x < imageWidth; ++x) {
                const size_t i = static_cast<size_t>(y) * static_cast<size_t>(imageWidth) + static_cast<size_t>(x);
                film.setPixel(x, y, (sum[i] + film.getData()[i]) / Real(causticPasses));
            }
    }
    else if (opt.lvc && setup.bptLight && setup.bptMedium) {
        VolumeBptParams bp;
        bp.width = imageWidth; bp.height = imageHeight; bp.spp = opt.spp; bp.lightPaths = opt.lvcPaths; bp.maxVertices = opt.lvcDepth; bp.passes = opt.lvcPasses; bp.surfaceAlbedo = setup.bptAlbedo;
        bp.useNee = (opt.lvcMask & 1) != 0; bp.useConnections = (opt.lvcMask & 2) != 0; bp.useLightTracing = (opt.lvcMask & 4) != 0;
        bp.mergeRadius = Real(opt.lvcMerge); bp.useMerging = (opt.lvcMask & 8) != 0;
        bp.useManifoldNee = (opt.lvcMask & 16) != 0; bp.useManifoldConn = (opt.lvcMask & 32) != 0;
        VolumeBPT vbpt(*scene, *camera, *setup.bptLight, *setup.bptMedium, bp);
        std::unique_ptr<SpecularNee> lvcManifold;   // the flat-faced glass bodies of the scene (the prism): manifold connections
        if ((opt.lvcMask & 48) != 0 && !setup.manifoldBodies.empty()) {
            lvcManifold = std::make_unique<SpecularNee>(setup.manifoldBodies, *scene);
            vbpt.setManifold(lvcManifold.get());
        }
        else if ((opt.lvcMask & 48) != 0) std::cout << "[VolumeBPT] note: this scene lists no flat-faced glass body: no manifold connections\n";
        vbpt.render(film);
    }
    else integrator->render(*scene, film);
    if (traceOn) {   // add the light-traced caustic image (already an average radiance per pixel)
        const Real inv = Real(1) / Real(causticPasses);
        for (int y = 0; y < imageHeight; ++y)
            for (int x = 0; x < imageWidth; ++x) {
                const size_t i = static_cast<size_t>(y) * static_cast<size_t>(imageWidth) + static_cast<size_t>(x);
                film.setPixel(x, y, (opt.causticTraceOnly ? Vector3(0.0) : film.getData()[i]) + traceSum[i] * inv);
            }
    }
    // --mask, --depth-map: a pass of camera rays alone, with its own random numbers (Render/MaskDepth.hpp); as many mask samples as the picture's
    if (!maskMaterials.empty()) film.enableMask();
    if (opt.depthMap) film.enableDepth();
    if (film.hasMask() || film.hasDepth()) {
        const auto t0 = std::chrono::steady_clock::now();
        const int side = std::clamp(static_cast<int>(std::lround(std::sqrt(double(opt.spp)))), 1, 16);
        RenderMaskAndDepth(*scene, *camera, film, maskMaterials, side);
        std::cout << "[MaskDepth] " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s\n";
    }
    if (setup.afterRender) setup.afterRender({ *scene, *camera, film, *integrator });

    if (specularNee) {
        const auto& c = specularNee->counters();
        std::cout << "[Manifold] " << specularNee->queries() << " receiver evaluations, " << specularNee->solutions() << " solutions; pairs tried "
            << c.pairs << " -> Newton converged " << c.converged << " -> inside " << c.insideFacets << " -> refraction " << c.refraction
            << " -> unblocked " << c.unblocked << " -> focusing ok " << c.focusOk << "\n";
    }

    SaveOutputs(opt, film, opt.lvc || causticOn);

    bool hooksPassed = true;
    if (setup.afterSave) hooksPassed = setup.afterSave({ *scene, *camera, film, *integrator });

    std::cout << "[System] Finished." << std::endl;

    // exit code 2 when a check failed (furnace test, mesh checks), so a script can see it
    return (hooksPassed && setup.checksPassed) ? 0 : 2;
}

} // namespace rayt::cli
