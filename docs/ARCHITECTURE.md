# Code structure (2026-11-16)

A map of what is where. The theory and results of each feature are in the respective `docs/*.md`. The whole is about 37 thousand lines (engine about 24 thousand, CLI about 6 thousand, self-checks about 6 thousand).

## 1. Two layers, and the flow of drawing one image

The **engine** (mainly the headers in `GoLD_rayt/include/`; where the body is large, a .cpp of the same name in `src/`: `Materials/` (CoatedDiffuse, Dielectric, LayeredCoat, RoughConductor, ThinFilm, Hair, WrinkleField), `Render/` (Integrator, Film, CausticPhotonMap), `Geometry/` (ProceduralMesh, GaussianSplats), `IO/` (PmxLoader, PlyLoader, StlLoader, ImageIO), `Accelerators/` (BVH)) and the **application** on top of it (`include/Cli` and `src/Cli`; the command line `main_cli.cpp`) form two layers. The GUI (`Application`, `UI`, `main_ui.cpp`) is a separate thin layer that calls the same engine. **The GUI is frozen** (the CLI is the main line; real-time use will pay off once a GPU implementation exists. It may later be repurposed as a "control panel for deciding what to draw", so `Cli/SceneDescription` (data independent of the engine) is kept as its foundation).

**A new .cpp must be registered in `GoLD_rayt.vcxproj` (and `.filters`), or it is not linked.** It begins with `#include "pch.h"`.

The flow of drawing one image:

1. `main_cli` → `Cli/Options` (parsing the arguments)
2. `Cli/Scenes/SceneRegistry` (chooses the builder function from the `--scene` name; JSON scenes go through `Scenes/FileScene`)
3. `Cli/RenderPipeline` (the common procedure: standard materials, environment, the scene-wide BVH, registering the lights, camera (`BuildCamera`), choosing the integrator, drawing, the mask / depth pass, saving (`SaveOutputs`))
4. `Render/Integrator` (the core of path tracing, `PathIntegrator::Li`)
5. `Render/Film` (accumulating pixels; gamut, exposure and tone map of the written values; the files are encoded by `IO/ImageWriter`) → `Render/Denoiser`

## 2. Role of each directory

| Directory | Contents |
|---|---|
| `Core` | Types, random numbers and sampling, AABB, `Ray`, `SurfaceInteraction` (information about the hit point; many fields were added by the skin features), `Spectrum` (switches between RGB and wavelength samples with the macro `RAYT_SPECTRAL_MODE`; this is why there are two kinds of build), textures (`SpectralTexture` is a spectral reflectance image), polarization (`Polarization`: Stokes vectors, Mueller matrices, Fresnel amplitudes) |
| `Color` | Conversion between spectrum and RGB, color matching functions, gamuts (sRGB / P3 / Rec.2020), CAT16 white balance, illumination spectra, tables of the complex refractive index of metals, Sellmeier dispersion |
| `Materials` | `Material` is the common BSDF interface (eval / sample / pdf, plus hooks for emission, opacity, medium boundaries, subsurface scattering, normal perturbation and so on). Implementations: diffuse, rough conductor, dielectric, thin film, coat layer, hair, skin (`CoatedDiffuse` is the skin itself and combines `SkinRelief`, `WrinkleField`, `SkinPigment`, `SkinFilm`, `SkinScatter`, `CosmeticLayer` and `KubelkaMunkLayer`; `SkinLotion` gives the time course of a lotion with which the film and the scatter are built; `LayeredSkin` is the Monte Carlo radial profile of a layered skin), subsurface scattering (the Burley type of `Subsurface` and the volumetric walk of `SubsurfaceWalk`), media (`Medium`, `DensityField`, `TabulatedPhase`), fading (`FadingReflectance`, `FadingDye`), measured BRDF (`MeasuredHpBRDF`: a table of hpBRDF compressed by PCA) |
| `Microfacet` | GGX and its directional albedo (multiple-scattering compensation) |
| `Lights` | Area light, point light, environment light, how lights are chosen (`LightSampler`: uniform, emitted power, spatial importance) |
| `Geometry` | Sphere, quad, triangle mesh, hair curves and the generator, meshes of revolution and of the ocean, Gaussian particles (3DGS) |
| `Scene` | `Hittable` (the common intersection interface), `Primitive` (bundles a shape, a material and a light), `Aggregate` (a list), `MovingObject` (motion blur), `Scene` (the whole: BVH, lights, environment, particles) |
| `Accelerators` | BVH (own SAH) and Embree |
| `Render` | Path tracing (`Integrator.cpp`; the hooks for media, particles and caustics are also here), transmittance of shadow rays, exit search for subsurface scattering, caustic photon map, manifold NEE, a prototype of LVC-BPT, camera, film (the values of the written pictures), the mask and depth pass (`MaskDepth`), denoiser |
| `IO` | OBJ, PLY, STL, PMX (model, morphs, pose), environment maps, CSV spectra, measured series of fading; `ImageWriter` encodes the image files (PNG 8 / 16 bit with cICP, BMP, JPG, PFM, HDR) |

## 3. The CLI side

| Where | Contents |
|---|---|
| `Cli/Options` | Default values of the arguments, structures, parsing |
| `Cli/Commands` | Commands that do not draw (checks, measurements, listing models) |
| `Cli/Scenes/` | Building the built-in scenes: metal, glass and media, studio, water, PMX, particles. `SceneRegistry` is the table |
| `Cli/SceneParts` | Building blocks shared by the demo scenes and the self-checks (meshes as hittables, lights, assets) |
| `Cli/SceneFile.cpp`, `SceneFileMaterials.cpp`, `SceneFileReaders.hpp`, `SceneDescription.hpp`, `Scenes/FileScene.cpp` | JSON scenes. Three stages: read → plain structures → build the engine objects |
| `Cli/RenderPipeline` | The common drawing procedure |
| `Cli/Checks/` | Self-checks by topic (`--scene furnace \| lights \| aov \| sss \| volume \| caustic \| manifoldtest \| metals \| mesh \| splats \| cosmetic \| relief \| camera \| film \| ...`) |
| `Cli/Tools/` | Measurement and diagnostics (`--scene skinlayers`, `sssrw`, `bench`, `--trace`) |

## 4. Features and where they are

| Feature | Where | Document |
|---|---|---|
| Spectral rendering, gamut, white balance | `Core/Spectrum`, `Color/*` | `METAL_ACCURACY.md`, `SPECTRUM_FROM_RGB.md` |
| Metal accuracy, n and k from measured reflectance, measured BRDF | `Materials/RoughConductor`, `MeasuredHpBRDF`, `Color/ComplexIorTable`, `ConductorFromReflectance` | `METAL_ACCURACY.md`, `MEASURED_GOLD.md` |
| Thin film, dispersion, coat | `Materials/ThinFilm`, `Color/Sellmeier`, `Materials/LayeredCoat` | `THIN_FILM.md`, `GLOSS.md` |
| Subsurface scattering | `Materials/Subsurface*`, `Render/SubsurfaceProbe` | `SSS.md` |
| Skin (relief, pigment, moisture, makeup) | `Materials/CoatedDiffuse` + `Skin*`, `WrinkleField`, `CosmeticLayer` | `SKIN.md`, `SKIN_MOISTURE.md`, `COSMETICS.md`, `WETNESS.md` |
| Hair | `Materials/Hair`, `Geometry/HairCurves`, `HairGenerator` | `HAIR.md` |
| Water, ocean, cloud and fog | `Geometry/ProceduralMesh`, `Materials/Medium`, `DensityField`, `TabulatedPhase` | `WATER.md`, `CLOUD.md` |
| Caustics, manifold NEE, LVC-BPT | `Render/CausticPhotonMap`, `SpecularNee`, `VolumeBPT` | `LVC_BPT.md` |
| Importing captured objects, 3DGS | `IO/PlyLoader`, `IO/StlLoader`, `Geometry/GaussianSplats`, `Cli/Scenes/FileScene` | `CAPTURE.md` |
| Polarization (Stokes, `--stokes`, `--polarizer`) | `Core/Polarization`, `Render/Integrator` | `POLARIZATION.md` |
| Camera and film output (lens, orthographic view, aberrations, diffraction, exposure and tone map, mask, depth, spectral image) | `Render/Camera`, `Render/Diffraction`, `Render/Film`, `Render/MaskDepth`, `IO/ImageWriter`, `Cli/RenderPipeline` | `CAMERA.md` |
| MMD | `IO/Pmx*`, `Cli/PmxSupport`, `Scenes/PmxScene` | — |
| Scene file | `Cli/SceneFile*`, `SceneDescription` | `SCENE_FILE.md` |

## 5. Records of tidy-ups

**Done (2026-11-13, splitting):** things whose body is long, which are edited often and whose interface is stable were split into .hpp and .cpp: `WrinkleField` (construction in the .cpp, only the per-ray lookup in the header), `Hair` (the lobe formulas `hairmath` in the header, the material in the .cpp), `CausticPhotonMap`, `PmxLoader`, `PlyLoader`, `GaussianSplats`. Not split: `VolumeBPT` and `SpecularNee` (prototypes that are left alone), `CosmeticLayer`, `SkinRelief` and the like (many small functions called per ray, which can be slow in the development build (no LTCG)), templates (`Spectrum`). Regression was checked with `tools/test/regress.py` (the mean color of small renderings + 15 self-checks). **Build time (development build):** rebuilding after editing the inside (.cpp) of a split item went from the full build (about 80 seconds) to about 4 seconds. When a header is touched, the whole is rebuilt both before and after the split (because `Cli/EngineHeaders.hpp` distributed all the headers to every CLI file).

**Done (2026-11-13, tidying includes):** `Cli/EngineHeaders.hpp` (a bundle of 63 headers, read by all 41 CLI files) was abolished, and each CLI file now includes only the headers it uses (a per-file list was made mechanically from the declared names in the headers, and the rest was fixed by building). Types used only by pointer or reference became forward declarations: `Render/Integrator.hpp` (forward declarations for the photon map, manifold NEE, medium, camera and film; `setCausticMap` and `mediaTraced` moved to the .cpp), `Cli/Scene.hpp`, `CheckSupport.hpp`, `Checks.hpp`, `Tools.hpp`, `Materials/Material.hpp` (`SubsurfaceTable.hpp` removed). `SpecularNee::Body` was moved out to `Render/SpecularBody.hpp` (only the shared_ptr of the mesh and glass), so the scene side does not read the prototype's header. **Effect (of 44 CLI files, the number that the header reaches):** the skin layers `CoatedDiffuse`, `WrinkleField`, `SkinRelief` 41 → 11, `Lambertian` 41 → 24, `Dielectric` 41 → 22, `Integrator` 41 → 21, `RoughConductor` 41 → 15, `PmxLoader` 41 → 8, `GaussianSplats` 41 → 5, `Hair` 41 → 3, `SpecularNee` 41 → 3, `VolumeBPT` 41 → 2. `Core/*`, `Color/*`, `Material.hpp` and `Scene.hpp` are foundations and still reach almost all (touching them rebuilds everything).
