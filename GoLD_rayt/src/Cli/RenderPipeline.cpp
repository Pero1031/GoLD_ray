// The common steps of rendering a scene: standard materials and environment, the scene builder, lights, caustics, camera, integrator,
// render, save, denoise. What differs between the scenes is in src/Cli/Scenes.
#include "pch.h"
#include "Cli/RenderPipeline.hpp"
#include "Cli/Scene.hpp"
#include "Cli/Tools.hpp"
#include "Render/CausticPhotonMap.hpp"
#include "Render/Denoiser.hpp"
#include "Render/Integrator.hpp"
#include "Render/SpecularNee.hpp"
#include "Render/VolumeBPT.hpp"

namespace rayt::cli {

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
    const Point3 lookFrom = setup.camera.from;
    const Point3 lookAt = setup.camera.at;
    const Vector3 vUp = setup.camera.up;
    const Real fovDeg = setup.camera.fovDeg;

    Real distToFocus = glm::length(lookFrom - lookAt);
    Real aperture = 0.0;   // pinhole camera (no depth of field)

    auto camera = std::make_shared<Camera>(
        lookFrom, lookAt, vUp,
        fovDeg, // FOV
        double(imageWidth) / imageHeight,
        aperture,
        distToFocus
    );

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
        VolumeBPT(*scene, *camera, *setup.bptLight, *setup.bptMedium, bp).render(film);
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
    if (setup.afterRender) setup.afterRender({ *scene, *camera, film, *integrator });

    if (specularNee) {
        const auto& c = specularNee->counters();
        std::cout << "[Manifold] " << specularNee->queries() << " receiver evaluations, " << specularNee->solutions() << " solutions; pairs tried "
            << c.pairs << " -> Newton converged " << c.converged << " -> inside " << c.insideFacets << " -> refraction " << c.refraction
            << " -> unblocked " << c.unblocked << " -> focusing ok " << c.focusOk << "\n";
    }

    std::cout << "[Output] Saving images..." << std::endl;
    color::Gamut gamut = color::Gamut::SRGB;
    color::parseGamut(opt.gamut, gamut);   // already validated by the option parser
    film.setOutputGamut(gamut);
    film.setPng16(opt.pngBits == 16);
    film.save(opt.out);
    const std::string outStem = std::filesystem::path(opt.out).replace_extension("").string();
    const std::string outExt = std::filesystem::path(opt.out).extension().string();
    if (opt.aov) {
        const std::string aovExt = (outExt == ".hdr") ? ".hdr" : ".png";
        film.saveAlbedo(outStem + "_albedo" + aovExt);
        film.saveNormal(outStem + "_normal" + aovExt);
    }
    if (opt.denoise) {
        Film denoised(imageWidth, imageHeight);
        denoised.setOutputGamut(gamut);
        denoised.setPng16(opt.pngBits == 16);
        if (DenoiseFilm(film, denoised, opt.cleanAux, gamut)) denoised.save(outStem + "_denoised" + outExt);
    }

    bool hooksPassed = true;
    if (setup.afterSave) hooksPassed = setup.afterSave({ *scene, *camera, film, *integrator });

    std::cout << "[System] Finished." << std::endl;

    // exit code 2 when a check failed (furnace test, mesh checks), so a script can see it
    return (hooksPassed && setup.checksPassed) ? 0 : 2;
}

} // namespace rayt::cli
