// The common steps of rendering a scene: standard materials and environment, the scene builder, lights, caustics, camera, integrator,
// render, save, denoise. What differs between the scenes is in src/Cli/Scenes.
#include "pch.h"
#include "Cli/KernelScenePack.hpp"
#include "Cli/RenderPipeline.hpp"
#include "Cli/Scene.hpp"
#include "Cli/Tools.hpp"
#include "Color/CameraResponse.hpp"
#include "Color/SpectralConversion.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Denoiser.hpp"
#include "Render/Diffraction.hpp"
#include "Render/Exposure.hpp"
#include "Render/LensSystem.hpp"
#include "Render/LensFlare.hpp"
#include "Render/Integrator.hpp"
#include "Render/ToonIntegrator.hpp"
#include "Render/MaskDepth.hpp"
#include "Render/Sensor.hpp"
#include "IO/DngWriter.hpp"
#include "Render/SpecularNee.hpp"
#include "Render/VolumeBPT.hpp"

#include <algorithm>
#include <sstream>

namespace rayt::cli {

namespace {

// The lens as the Camera takes it. The scene gives a field of view (or a lens); the command line's lens replaces it field by field.
struct ResolvedLens {
    Real orthoHeight = 0;     // > 0: the orthographic view, this high (scene units)
    Real diffraction = 0;     // > 0: diffraction at the aperture, the pattern this many times its real size
    std::string ghosts;       // the ghosts of a real lens: "coated", "uncoated" ("" or "off": none)
    std::string lensFile;     // a real lens instead of the thin one ("" = none)
    bool fStopGiven = false;  // an f-number was asked for ("fStop", --fstop), not only a lens diameter
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
    Real iso = 0;             // > 0: physical exposure (Render/Exposure.hpp) with this ISO speed
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
    r.diffraction = s.diffraction.value_or(Real(0));
    r.ghosts = s.ghosts.value_or(std::string());
    r.lensFile = s.lensFile.value_or(std::string());
    r.fStopGiven = s.fStop.has_value() && *s.fStop > 0;
    r.shutter = s.shutter.value_or(Real(0));
    r.iso = s.iso.value_or(Real(0));
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

// --camera-response: the film as the camera of this file sees (Color/CameraResponse.hpp), and what the set-up found.
bool ApplyCameraResponse(const std::string& path) {
    std::string error;
    if (!color::setCameraResponse(path, error)) { std::cerr << "[Camera] --camera-response: " << error << "\n"; return false; }
    const color::CameraResponse& cam = color::cameraResponse();
    const Vector3 white = cam.d65Raw / cam.d65Raw.y;
    std::cout << "[Camera] response " << cam.name << ": raw of D65 R " << white.x << " G 1 B " << white.z
        << " (balanced to 1 1 1); colour matrix fitted on " << cam.trainingCount << " reflectances under D65, RMS error " << cam.fitRms << ":\n";
    for (int r = 0; r < 3; ++r)
        std::cout << "[Camera]   " << cam.matrix[0][r] << "  " << cam.matrix[1][r] << "  " << cam.matrix[2][r] << "\n";
    return true;
}

// Diffraction at an aperture of f-number N (--diffraction): the pattern's setup, and a line on its size.
DiffractionSetup MakeDiffraction(const ResolvedLens& lens, Real fNumber, int imageWidth, int imageHeight, Real filmLongSide = Real(0.036)) {
    DiffractionSetup d;
    d.fNumber = fNumber;
    d.blades = lens.blades;
    d.bladeRotation = lens.bladeRotationDeg * constants::DEG_TO_RAD;
    d.pixelPitch = filmLongSide / Real(std::max(imageWidth, imageHeight));   // the film's long side (the 35 mm convention: 36 mm)
    d.scale = lens.diffraction;
    d.radius = std::min(255, std::max(imageWidth, imageHeight));
    std::cout << "[Camera] diffraction at the aperture (f/" << fNumber << ", " << (lens.blades >= 3 ? std::to_string(lens.blades) + " blades" : std::string("round"))
        << "): Airy radius 1.22 lambda N = " << 1.22 * 0.55 * fNumber << " um at 550 nm, pixel " << d.pixelPitch * 1e6 << " um";
    if (lens.diffraction > 1) std::cout << ", the pattern enlarged " << lens.diffraction << " times";
    else std::cout << " (real size)";
    std::cout << "\n";
    return d;
}

// Physical exposure (--iso, camera "iso"): the gain from the scene's radiance to the camera's linear value, in stops, for the irradiance
// factor G (sr) at the film centre, the shutter and --nits. Prints the settings; false when something it needs is missing (the run stops).
bool PhysicalExposureStops(const CliOptions& opt, const ResolvedLens& lens, Real fNumber, Real irradianceFactor, Real& stops) {
    stops = 0;
    if (!(lens.iso > 0)) return true;
    if (!(lens.shutter > 0)) { std::cerr << "[Camera] --iso needs an exposure time (--shutter, camera \"shutter\")\n"; return false; }
    if (!(opt.nits > 0)) {
        std::cerr << "[Camera] --iso needs the scene's luminance scale: --nits K, the luminance (cd/m^2) of a radiance whose Y is 1\n";
        return false;
    }
    if (!(irradianceFactor > 0)) { std::cerr << "[Camera] --iso: no light reaches the film centre through this lens\n"; return false; }
    const Real gain = PhysicalExposureGain(lens.iso, lens.shutter, opt.nits, irradianceFactor);
    const Real working = std::sqrt(Real(3.14159265358979323846) / (Real(4) * irradianceFactor));
    const Real ev100 = std::log2(fNumber * fNumber / lens.shutter) - std::log2(lens.iso / Real(100));
    stops = std::log2(gain);
    std::cout << "[Camera] physical exposure: f/" << fNumber << " (working f/" << working << " at the focus distance), " << lens.shutter << " s";
    if (lens.shutter < 1.0) std::cout << " (1/" << std::lround(1.0 / lens.shutter) << ")";
    std::cout << ", ISO " << lens.iso << " (EV100 " << ev100 << "); a radiance of Y 1 (" << opt.nits << " cd/m^2) is recorded as " << gain
        << " (" << (stops >= 0 ? "+" : "") << stops << " EV on the .png / .bmp / .jpg; the .pfm keeps the radiance); mid grey (0.18) is "
        << MidGreyLuminance(lens.iso, lens.shutter, irradianceFactor) << " cd/m^2\n";
    if (opt.tonemap != "clip") std::cout << "[Camera] note: a camera clips at its full well: --tonemap clip shows over-exposure as a camera does (Reinhard compresses it)\n";
    return true;
}

// A real lens (--lens, camera "lens"): its prescription, the stop closed to the f-number asked for, focused on the focus distance (from the
// front principal plane, which sits at the camera position), the film of the 35 mm convention, or the picture's rectangle on the lens's own
// image circle when its file gives one (image-height: a phone lens and its sensor). Its field of view comes from the lens and the film.
std::shared_ptr<Camera> BuildLensCamera(const CliOptions& opt, const SceneSetup& setup, const ResolvedLens& lens, int imageWidth, int imageHeight,
    bool lightTracing, std::optional<DiffractionSetup>& diffraction, Real& exposureStops) {
    if (lens.orthoHeight > 0) { std::cerr << "[Camera] --lens and --ortho exclude each other\n"; return nullptr; }
    if (lightTracing) {
        std::cerr << "[Camera] a real lens (--lens) cannot be used with light tracing (--caustic-trace, --lvc, or a scene that light-traces its caustics):"
            " it weighs its paths for a pinhole camera\n";
        return nullptr;
    }
    std::string error;
    const auto system = LensSystem::load(lens.lensFile, error);
    if (!system) { std::cerr << "[Lens] " << error << "\n"; return nullptr; }
    const double aspect = double(imageWidth) / imageHeight;
    Real filmW = aspect >= 1 ? Real(0.036) : Real(0.036 * aspect), filmH = aspect >= 1 ? Real(0.036 / aspect) : Real(0.036);
    if (system->imageHeight() > 0) {   // the rectangle whose diagonal is the image circle's diameter
        filmW = Real(2) * system->imageHeight() * Real(aspect / std::sqrt(1 + aspect * aspect));
        filmH = filmW / Real(aspect);
    }
    const Real wideOpen = system->wideOpenFNumber();
    const Real fNumber = lens.fStopGiven ? system->setFNumber(lens.fStop) : wideOpen;
    system->setStopShape(lens.blades, lens.bladeRotationDeg * constants::DEG_TO_RAD);
    const bool focused = system->focus(lens.focusDist);
    system->prepare(filmW, filmH);

    auto camera = std::make_shared<Camera>(setup.camera.from, setup.camera.at, setup.camera.up, lens.fovDeg, aspect, Real(0), lens.focusDist);
    camera->setLensSystem(system);
    camera->setShutter(lens.shutter);
    camera->setMotion(lens.velocity, lens.panYaw, lens.panPitch);

    const Real f = system->focalLength();
    std::cout << "[Lens] " << lens.lensFile << ": " << system->surfaces().size() << " surfaces, focal length " << f * 1000 << " mm, f/" << fNumber
        << " (wide open f/" << wideOpen << ")" << (lens.blades >= 3 ? ", " + std::to_string(lens.blades) + " blades" : std::string())
        << "; film " << filmW * 1000 << " x " << filmH * 1000 << " mm, paraxial field of view " << 2 * std::atan(filmW / (2 * f)) * constants::RAD_TO_DEG
        << " x " << 2 * std::atan(filmH / (2 * f)) * constants::RAD_TO_DEG << " deg; focused at " << lens.focusDist << " m, the film "
        << (system->surfaces().back().z - system->filmZ()) * 1000 << " mm behind the rear vertex\n";
    if (lens.fStopGiven && lens.fStop < wideOpen * Real(0.999)) std::cout << "[Lens] note: the lens opens no wider than f/" << wideOpen << "\n";
    if (!focused) std::cout << "[Lens] warning: the focus distance is not beyond the focal length: focused at infinity\n";
    if (lens.k1 != 0 || lens.k2 != 0 || lens.lateralCA != 0 || lens.axialCA != 0 || lens.barrelRadius > 0 || lens.naturalVignetting > 0)
        std::cout << "[Lens] note: a real lens has its own distortion, chromatic aberration and vignetting: those settings are not used\n";
    if (opt.lens.focalLengthMm || setup.camera.lens.focalLengthMm) std::cout << "[Lens] note: the focal length is the lens's own; \"focalLength\" / --focal-length is not used\n";
    const bool hasAbbe = std::any_of(system->surfaces().begin(), system->surfaces().end(), [](const LensSystem::Surface& s) { return s.abbe > 0 && s.nd > 1; });
#if defined(RAYT_SPECTRAL_MODE)
    if (system->dispersive()) std::cout << "[Lens] dispersive glass: every wavelength takes its own path (chromatic aberration; one wavelength per path)\n";
    else std::cout << "[Lens] note: the file gives no Abbe numbers: the glass is not dispersive\n";
#else
    if (hasAbbe) std::cout << "[Lens] note: the glass is dispersive, but the RGB build traces n_d only (use the spectral build for chromatic aberration)\n";
#endif
    (void)hasAbbe;
    if (lens.diffraction > 0) diffraction = MakeDiffraction(lens, fNumber, imageWidth, imageHeight, std::max(filmW, filmH));
    if (!PhysicalExposureStops(opt, lens, fNumber, system->centreIrradiance(), exposureStops)) return nullptr;
    return camera;
}

// The camera of the scene with the lens of the command line, and the notes on what the chosen render mode does not show of it.
// lightTracing: light paths are connected to a pinhole camera (--caustic-trace, --lvc); lvcCamera: --lvc renders the picture.
// diffraction: set when the lens asks for diffraction at the aperture (applied to the pictures when they are written).
std::shared_ptr<Camera> BuildCamera(const CliOptions& opt, const SceneSetup& setup, int imageWidth, int imageHeight, bool lightTracing, bool lvcCamera,
    std::optional<DiffractionSetup>& diffraction, Real& exposureStops) {
    exposureStops = 0;
    const Point3 lookFrom = setup.camera.from;
    const Point3 lookAt = setup.camera.at;
    const Vector3 vUp = setup.camera.up;

    const ResolvedLens lens = ResolveLens(setup.camera, opt.lens, double(imageWidth) / imageHeight);
    if (!lens.lensFile.empty()) return BuildLensCamera(opt, setup, lens, imageWidth, imageHeight, lightTracing, diffraction, exposureStops);
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
    if (lens.diffraction > 0) {
        if (!(lens.aperture > 0)) {
            std::cerr << "[Camera] --diffraction needs an aperture (--fstop, \"fStop\"): a pinhole has no diffraction pattern of finite size\n";
            return nullptr;
        }
        diffraction = MakeDiffraction(lens, lens.fStop, imageWidth, imageHeight);
    }
    if (lens.iso > 0) {
        if (lens.orthoHeight > 0) { std::cerr << "[Camera] --iso: the orthographic view has no film irradiance of a photographic lens (use the perspective view)\n"; return nullptr; }
        if (!(lens.aperture > 0)) { std::cerr << "[Camera] --iso needs an f-number (--fstop, \"fStop\"): a pinhole lets no measurable light through\n"; return nullptr; }
        const Real g = ThinLensIrradianceFactor(lens.fStop, lens.focalLengthMm * Real(0.001), lens.focusDist);
        if (!PhysicalExposureStops(opt, lens, lens.fStop, g, exposureStops)) return nullptr;
        if (!(lens.naturalVignetting > 0)) std::cout << "[Camera] note: the exposure is the centre's; the thin lens darkens the corners only with --natural-vignette 1 (cos^4)\n";
    }
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

// --sensor: the model and its numbers at this ISO in the log.
void LogSensor(const SensorModel& m, const SensorSetting& s, bool physicalExposure) {
    std::cout << "[Sensor] full well " << m.fullWell << " e-, base ISO " << m.baseIso << ", read noise " << m.readNoise << " e- before the amplifier and "
        << m.adcNoise << " e- after it (at the base ISO), " << m.bits << "-bit ADC\n";
    const Real grey = Real(0.18) * s.whiteElectrons;
    const Real snr = grey / std::sqrt(grey + s.readElectrons * s.readElectrons);
    std::cout << "[Sensor] ISO " << s.iso << (physicalExposure ? "" : " (the base ISO: no physical exposure, --exposure sets the light)") << ": the picture's white is "
        << s.whiteElectrons << " e- (" << s.whiteElectrons / m.fullWell << " of the full well), " << s.dnPerElectron << " DN per e- above a black level of "
        << s.black << " DN, read noise " << s.readElectrons << " e-; mid grey " << grey << " e-, signal-to-noise ratio " << snr << " (shot noise alone "
        << std::sqrt(grey) << ")\n";
    if (s.whiteElectrons > m.fullWell) std::cout << "[Sensor] below the base ISO the well is full at " << m.fullWell / s.whiteElectrons << " of the white\n";
    std::cout << "[Sensor] " << (m.mosaic ? std::string("Bayer mosaic RGGB, ") + (m.demosaic == DemosaicMethod::Malvar ? "Malvar" : "bilinear")
        + " demosaicing (a raw file <out>.dng beside the picture)" : std::string("no mosaic: every pixel records all three channels")) << "\n";
}

// The DNG tags of the sensor's mosaic: its levels, the camera's colours under illuminant A and D65 (XYZ -> raw, each scaled so that its
// light's white has a largest raw channel of 1) and the white it is balanced to (without --camera-response the raw channels are the
// film's linear sRGB).
io::DngInfo SensorDngInfo(const SensorModel& m, const SensorSetting& s) {
    io::DngInfo info;
    info.black = s.black;
    info.white = s.maxDn;
    auto calibration = [](bool illuminantA, double (&out)[9]) {
        const color::Matrix3 a = color::cameraRawFromXYZ(illuminantA);
        const Vector3 lightWhite = a * (illuminantA ? Vector3(1.09850, 1, 0.35585) : Vector3(0.95047, 1, 1.08883));   // XYZ of A / D65, Y 1
        const Real top = std::max(lightWhite.x, std::max(lightWhite.y, lightWhite.z));
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) out[3 * i + j] = double(a[j][i] / top);   // row i, column j (glm keeps columns)
    };
    calibration(true, info.colorMatrixA);
    calibration(false, info.colorMatrix);
    // ForwardMatrix: white-balanced raw -> the film's RGB (the calibration's colour matrix) -> XYZ (D65 of Y 1) -> Bradford to D50, the
    // DNG's connection space; a reader goes back to D65 by Bradford too, so it gets the renderer's colours.
    auto fromRows = [](std::initializer_list<Real> v) { color::Matrix3 m(Real(0)); auto it = v.begin(); for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) m[c][r] = *it++; return m; };
    const color::Matrix3 bradford = fromRows({ 0.8951, 0.2664, -0.1614, -0.7502, 1.7135, 0.0367, 0.0389, -0.0685, 1.0296 });
    const color::Matrix3 toXYZ = color::filmToXYZ() / (color::filmToXYZ() * Vector3(1.0)).y;   // the film's RGB -> XYZ, D65 of Y 1
    const Vector3 coneD65 = bradford * (toXYZ * Vector3(1.0)), coneD50 = bradford * Vector3(0.96422, 1.0, 0.82521);
    color::Matrix3 scale(Real(0));
    for (int i = 0; i < 3; ++i) scale[i][i] = coneD50[i] / coneD65[i];
    const color::Matrix3 toD50 = glm::inverse(bradford) * scale * bradford * toXYZ;
    auto forward = [&](bool illuminantA, double (&out)[9]) {
        const color::Matrix3 f = toD50 * color::cameraCalibrationMatrix(illuminantA);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) out[3 * i + j] = double(f[j][i]);
    };
    forward(true, info.forwardMatrixA);
    forward(false, info.forwardMatrix);
    const Vector3 white = Vector3(1.0) / color::cameraWhiteBalanceGains();   // the raw of the balanced white
    const Real top = std::max(white.x, std::max(white.y, white.z));
    for (int i = 0; i < 3; ++i) info.asShotNeutral[i] = double(white[i] / top);
    const color::CameraResponse& cam = color::cameraResponse();
    info.model = "GoLD_rayt " + (cam.active ? cam.name : std::string("ideal sRGB camera"));
    std::ostringstream d;
    d << "GoLD_rayt sensor at ISO " << s.iso << ": full well " << m.fullWell << " e- at ISO " << m.baseIso << ", read noise " << m.readNoise
      << " + " << m.adcNoise << " e-, " << m.bits << "-bit ADC, " << s.whiteElectrons << " e- at white";
    info.description = d.str();
    return info;
}

// --awb: a camera's auto white balance. The light is guessed from the picture by the grey world assumption (on average a scene
// reflects grey, so the picture's mean is the light's colour), and the picture is balanced to it: with a camera's response by gains on
// its raw channels and a matrix interpolated between its A and D65 fits by the white's colour temperature; for the CIE film by CAT16
// from that white to D65. The denoised copy gets the same balance.
void LogAutoWhiteBalance(Film& film, Film* denoised) {
    const AutoBalance b = AutoWhiteBalance(film, denoised);
    if (!b.ok) { std::cout << "[AWB] the picture has no usable white to balance to; left as it is\n"; return; }
    const bool camera = color::cameraResponse().active;
    std::cout << "[AWB] grey world: the picture's mean " << (camera ? "raw " : "RGB ") << b.mean.x / b.mean.y << " 1 " << b.mean.z / b.mean.y
        << ", a light of about " << b.cct << " K; ";
    if (camera) {
        const Vector3 gains = color::cameraWhiteBalanceGains() / color::cameraWhiteBalanceGains().y;
        std::cout << "gains R " << gains.x << " G 1 B " << gains.z << ", the colour matrix interpolated between the A and D65 fits\n";
    }
    else std::cout << "CAT16 from it to D65\n";
}

// Writes the picture and everything asked for beside it (mask, depth, spectral image, AOVs, Stokes images, polarizer views, the denoised copy).
// photonsAdded: caustic photons or light paths were added to the picture (they carry no polarization). iso: the physical exposure's (0 = none).
void SaveOutputs(const CliOptions& opt, Film& film, bool photonsAdded, const std::optional<DiffractionSetup>& diffraction, Real exposureStops, Real iso,
    const std::vector<Vector3>& flare) {
    std::cout << "[Output] Saving images..." << std::endl;
    color::Gamut gamut = color::Gamut::SRGB;
    color::parseGamut(opt.gamut, gamut);   // already validated by the option parser
    // Every written picture (the image, the polarizer images, the denoised copy) gets the same output settings.
    auto setOutput = [&](Film& f) {
        f.setOutputGamut(gamut);
        f.setPng16(opt.pngBits == 16);
        f.setDisplay(opt.exposure + exposureStops, opt.tonemap == "clip" ? Film::ToneMap::Clip : Film::ToneMap::Reinhard);
    };
    setOutput(film);

    // The denoiser works on the picture as rendered (a noisy bright pixel would otherwise turn into a star of its own).
    std::unique_ptr<Film> denoised;
    if (opt.denoise) {
        denoised = std::make_unique<Film>(film.width(), film.height());
        setOutput(*denoised);
        if (!DenoiseFilm(film, *denoised, opt.cleanAux, gamut)) denoised.reset();
    }
    // The ghosts of the lens (--ghosts): optics, added after the denoiser (they are free of noise) and before diffraction
    if (!flare.empty()) {
        for (Film* f : { &film, denoised.get() }) {
            if (!f) continue;
            for (int y = 0; y < f->height(); ++y)
                for (int x = 0; x < f->width(); ++x) {
                    const size_t i = static_cast<size_t>(y) * static_cast<size_t>(f->width()) + static_cast<size_t>(x);
                    f->setPixel(x, y, f->getData()[i] + flare[i]);
                }
        }
    }
    // Diffraction at the aperture: optics, so it reaches every picture, .pfm included (not the AOVs, the mask, the depth or the spectral image)
    std::optional<DiffractionKernels> kernels;
    if (diffraction) {
        const auto t0 = std::chrono::steady_clock::now();
        kernels = ComputeDiffractionKernels(*diffraction);
        ApplyDiffraction(film, *kernels);
        if (denoised) ApplyDiffraction(*denoised, *kernels);
        std::cout << "[Diffraction] kernel of " << 2 * kernels->radius + 1 << " px holding " << 100.0 * kernels->captured << " % of the light, applied in "
            << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s\n";
    }
    if (opt.autoWhiteBalance) LogAutoWhiteBalance(film, denoised.get());   // before the sensor, whose processing uses the balance
    // The sensor: after the optics (diffraction) and after the denoiser, which would otherwise take its noise away. With a physical
    // exposure the sensor gets the exposure's light and --exposure pushes the recording afterwards, noise and all; without one --exposure
    // sets the light and the picture's white is the full well at the base ISO. The film keeps radiance units (.pfm).
    if (!opt.sensor.empty()) {
        SensorModel model;
        std::string error;
        ParseSensorSpec(opt.sensor, model, error);   // validated by the option parser
        const SensorSetting setting = MakeSensorSetting(model, iso);
        const Real gain = std::pow(Real(2), iso > 0 ? exposureStops : opt.exposure);
        LogSensor(model, setting, iso > 0);
        const std::uint64_t seed = static_cast<std::uint64_t>(opt.seed);
        std::vector<std::uint16_t> mosaic, mosaicDenoised;
        const Real clipped = ApplySensor(film, model, setting, gain, seed, &mosaic);
        if (denoised) ApplySensor(*denoised, model, setting, gain, seed ^ 0xd3e01c0ull, &mosaicDenoised);   // its own noise
        std::cout << "[Sensor] " << Real(100) * clipped << " % of the pixels have a channel at the top of the range\n";
        if (!mosaic.empty()) {   // the raw files: the mosaic as recorded, for a raw developer
            const io::DngInfo info = SensorDngInfo(model, setting);
            const std::string stem = std::filesystem::path(opt.out).replace_extension("").string();
            auto writeDng = [&](const std::string& path, const std::vector<std::uint16_t>& dn) {
                if (io::WriteDng(path, film.width(), film.height(), dn, info)) std::cout << "[Sensor] Saved raw file (DNG, the RGGB mosaic as recorded): " << path << "\n";
                else std::cerr << "[Sensor] Error: could not write " << path << "\n";
            };
            writeDng(stem + ".dng", mosaic);
            if (denoised) writeDng(stem + "_denoised.dng", mosaicDenoised);
        }
        if (!color::cameraResponse().active)
            std::cout << "[Sensor] note: its channels are the film's linear sRGB (an ideal camera); --camera-response gives a real camera's (spectral build)\n";
        if (!denoised) std::cout << "[Sensor] note: the render's own noise adds to the sensor's; with --denoise the _denoised picture has the sensor's alone\n";
        if (opt.tonemap != "clip") std::cout << "[Sensor] note: --tonemap clip shows the sensor's clip as it is (Reinhard compresses it)\n";
    }

    film.save(opt.out);
    const std::string outStem = std::filesystem::path(opt.out).replace_extension("").string();
    const std::string outExt = std::filesystem::path(opt.out).extension().string();
    std::string outExtLower = outExt;
    std::transform(outExtLower.begin(), outExtLower.end(), outExtLower.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    const bool floatOut = (outExtLower == ".pfm" || outExtLower == ".hdr");
    if (floatOut && (opt.exposure != 0 || opt.tonemap != "reinhard" || exposureStops != 0))
        std::cout << "[Output] note: --exposure, --iso and --tonemap change only .png / .bmp / .jpg; " << outExtLower << " holds the radiance as rendered\n";
    if (film.hasMask()) film.saveMask(outStem + (floatOut ? "_mask.pfm" : "_mask.png"));
    if (film.hasDepth()) film.saveDepth(outStem + "_depth.pfm");
    if (film.hasSpectral()) film.saveSpectral(outStem);
    if (color::cameraResponse().active) {   // the camera's raw RGB before its white balance and colour matrix (G of D65 = 1)
        Film raw(film.width(), film.height());
        for (int y = 0; y < film.height(); ++y)
            for (int x = 0; x < film.width(); ++x)
                raw.setPixel(x, y, color::filmToCameraRaw(film.getData()[static_cast<size_t>(y) * film.width() + x]));
        raw.save(outStem + "_camraw.pfm");
    }
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
            if (kernels) ApplyDiffraction(filtered, *kernels);
            std::ostringstream name;
            name << outStem << "_pol" << angle;
            std::string stem = name.str();
            std::replace(stem.begin() + static_cast<std::ptrdiff_t>(outStem.size()), stem.end(), '.', 'p');   // 22.5 -> 22p5: keeps the extension clear
            filtered.save(stem + outExt);
        }
        if (photonsAdded)
            std::cout << "[Polarization] note: the caustic photons and the light-traced image are added without polarization" << std::endl;
    }
    if (denoised) denoised->save(outStem + "_denoised" + outExt);
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
#if defined(RAYT_SPECTRAL_MODE)
    // The toon mode paints display colors: it has no use for spectra and lives in the RGB build only.
    if (opt.toon) { std::cerr << "[Toon] --toon is a mode of the RGB build (out\\cli or out\\cli-dev)\n"; return 1; }
#endif
    const ImageSize size = ComputeImageSize(opt);
    const int imageWidth = size.width;
    const int imageHeight = size.height;

    // ---- the environment and the shared materials, then the scene itself ----
    std::cout << "CWD = " << std::filesystem::current_path() << std::endl;
    if (!opt.cameraResponse.empty() && !ApplyCameraResponse(opt.cameraResponse)) return 1;   // before anything turns a spectrum into the film's RGB
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

    if (const color::CameraResponse& cam = color::cameraResponse(); cam.active && cam.balancedTo != "D65") {   // --wb with a camera
        const Vector3 gains = color::cameraWhiteBalanceGains() / color::cameraWhiteBalanceGains().y;
        std::cout << "[Camera] balanced to " << cam.balancedTo << ": gains R " << gains.x << " G 1 B " << gains.z << ", colour matrix fitted under it (RMS "
            << cam.fitRms << " against the CIE film adapted to D65 with CAT16)\n";
    }

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
    std::optional<DiffractionSetup> diffraction;
    Real exposureStops = 0;   // --iso: the physical exposure's gain on the written pictures
    const auto camera = BuildCamera(opt, setup, imageWidth, imageHeight, opt.causticTrace || caustics.lightTraceOnly || lvcCamera, lvcCamera, diffraction,
        exposureStops);
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

    // --device gpu | cpu-kernel: the kernel renderer (docs/GPU.md) draws the picture instead of the PathIntegrator; it is written as usual
    if (opt.device != "cpu") {
        if (const int exitCode = RenderOnKernel(opt, setup.namedMaterials, env.light.get(), *camera, film)) return exitCode;
        const ResolvedLens resolvedLens = ResolveLens(setup.camera, opt.lens, double(imageWidth) / imageHeight);
        SaveOutputs(opt, film, false, diffraction, exposureStops, resolvedLens.iso, {});
        std::cout << "[System] Finished." << std::endl;
        return 0;
    }

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
    if (setup.pixelFootprints) integrator->setPixelFootprints(true);
    if (setup.mediaNee || opt.mediaNee) integrator->setMediaNee(true);   // the fog is lit by next event estimation, its boundary does not block shadow rays
    std::cout << "[Integrator] PathIntegrator created (maxDepth="
        << opt.maxDepth << ", spp=" << opt.spp << ").\n";
    if (causticPasses > 1)
        std::cout << "[Caustics] " << causticPasses << " passes of " << sppPerPass << " spp, each with its own " << causticParams.numPhotons
            << " photons (" << double(causticPasses) * double(causticParams.numPhotons) / 1e6 << " M in total, memory of one set)\n";

    // ---- render and save ----
    std::cout << "[Render] Start PBR rendering..." << std::endl;
    if (opt.toon) {
        ToonSettings ts;
        ts.lineWidth = Real(opt.toonLine);
        ts.looks = setup.toonLooks;
        ToonIntegrator(camera, opt.spp, ts).render(*scene, film);
        opt.tonemap = "clip";   // a display-referred picture in [0, 1]
    }
    else if (causticPasses > 1) {
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
    if (traceOn && !opt.toon) {   // add the light-traced caustic image (already an average radiance per pixel)
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

    const ResolvedLens resolved = ResolveLens(setup.camera, opt.lens, double(imageWidth) / imageHeight);
    const Real iso = resolved.iso;   // the physical exposure's, for the sensor
    // The lens flare (--ghosts): the ghosts of the scene's lights through the real lens, traced from the lights
    std::vector<Vector3> flare;
    if (!resolved.ghosts.empty() && resolved.ghosts != "off") {
        if (!camera->lensSystem()) std::cout << "[Ghosts] --ghosts needs a real lens (--lens): ignored\n";
        else {
            LensFlareSettings fs;
            fs.coating = resolved.ghosts == "uncoated" ? LensSystem::Coating::Uncoated : LensSystem::Coating::SingleLayer;
            LensFlareStats st;
            flare = RenderLensFlare(*scene, *camera, imageWidth, imageHeight, fs, &st);
            std::cout << "[Ghosts] " << st.ghosts << " pairs of surfaces x " << st.lights << " light point(s), "
                << (resolved.ghosts == "uncoated" ? "bare glass" : "a single MgF2 coating") << ": " << st.rays / 1e6 << " M rays in " << st.seconds
                << " s; the brightest pixel " << st.peak << "\n";
            if (st.lights == 0) std::cout << "[Ghosts] no light of the scene reaches the lens (environment lights cast no ghosts here)\n";
        }
    }
    SaveOutputs(opt, film, opt.lvc || causticOn, diffraction, exposureStops, iso, flare);

    bool hooksPassed = true;
    if (setup.afterSave) hooksPassed = setup.afterSave({ *scene, *camera, film, *integrator });

    std::cout << "[System] Finished." << std::endl;

    // exit code 2 when a check failed (furnace test, mesh checks), so a script can see it
    return (hooksPassed && setup.checksPassed) ? 0 : 2;
}

} // namespace rayt::cli
