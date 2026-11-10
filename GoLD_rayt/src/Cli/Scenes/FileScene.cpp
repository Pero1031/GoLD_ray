// The scene "file": a scene written as data in a JSON file (see Cli/SceneFile.hpp and Cli/SceneDescription.hpp).
#include "pch.h"
#include "Cli/Scenes.hpp"
#include "Cli/LightPresets.hpp"
#include "Cli/SceneParts.hpp"
#include "Geometry/ProceduralMesh.hpp"
#include "Cli/EngineHeaders.hpp"
#include "stb_image_write.h"

namespace rayt::cli {

namespace {

// Debugging aid: with the environment variable RAYT_WRINKLE_DUMP = file.png the height map of the wrinkle simulation is written as a picture
// (dark = a valley), and RAYT_WRINKLE_DUMP_PORES = file.png writes the pore mask.
void DumpWrinkleMap(const WrinkleField& field) {
    auto write = [&](const char* var, const std::vector<float>& map, bool heights) {
        char* env = nullptr;
        size_t envLen = 0;
        if (_dupenv_s(&env, &envLen, var) != 0 || !env || !*env) { std::free(env); return; }
        const std::string pathStr(env);
        std::free(env);
        const char* path = pathStr.c_str();
        const int R = field.resolution();
        float lo = 0, hi = 0;
        for (float v : map) { lo = std::min(lo, v); hi = std::max(hi, v); }
        std::vector<unsigned char> img(map.size());
        for (size_t i = 0; i < map.size(); ++i) {
            const float t = heights ? (map[i] - lo) / std::max(hi - lo, 1e-12f) : map[i];
            img[i] = static_cast<unsigned char>(std::clamp(t, 0.f, 1.f) * 255.f + 0.5f);
        }
        const int ok = stbi_write_png(path, R, R, 1, img.data(), R);
        std::cout << "[Wrinkles] " << (ok ? "wrote " : "could not write ") << path << " (" << lo * 1e6f << " .. " << hi * 1e6f << " um)\n";
    };
    write("RAYT_WRINKLE_DUMP", field.heightMap(), true);
    write("RAYT_WRINKLE_DUMP_PORES", field.poreMap(), false);
}

// The spectrum of an emission material.
bool MakeSpectrum(const SpectrumSpec& s, EmissionSpectrum& out) {
    switch (s.kind) {
    case SpectrumSpec::Kind::Preset:
        return MakeLightEmission(s.preset, out);
    case SpectrumSpec::Kind::Blackbody:
        out = EmissionSpectrum::Blackbody(Real(s.kelvin));
        return true;
    case SpectrumSpec::Kind::Rgb:
        out = EmissionSpectrum(RgbColor(s.rgb));
        return true;
    case SpectrumSpec::Kind::Curve: {
        color::SpectralCurve curve;
        std::string err;
        if (!io::loadSpectralCurveCsv(s.file, curve, &err)) { std::cerr << "[SceneFile] spectrum " << s.file << ": " << err << "\n"; return false; }
        out = EmissionSpectrum::FromCurve(curve);
        return true;
    }
    }
    return false;
}

} // namespace

std::optional<int> BuildFromDescription(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    if (!opt.sceneDescription) {
        std::cerr << "[SceneFile] the scene \"file\" is a scene written as data: GoLD_rayt scene.json (a file with \"objects\")\n";
        return 1;
    }
    const SceneDescription& d = *opt.sceneDescription;
    const auto& world = out.world;
    out.camera = d.camera;
    if (ctx.env.light && d.environment.scale != 1.0) ctx.env.light->setScale(Real(d.environment.scale));
    std::cout << "[SceneFile] " << opt.sceneFile << ": " << d.materials.size() << " materials, " << d.objects.size() << " objects, "
        << d.lights.size() << " point lights\n";

    // ---- materials ----
    std::map<std::string, std::shared_ptr<Material>> materials;
    std::map<std::string, std::shared_ptr<color::ComplexIorTable>> tables;   // one table per metal / file, shared by its materials
    for (const auto& [name, m] : d.materials) {
        switch (m.type) {
        case MaterialSpec::Type::Diffuse:
            if (m.reflectance.empty()) {
                materials[name] = std::make_shared<Lambertian>(RgbColor(m.albedo));
                break;
            }
            else {   // a measured reflectance R(lambda): used as is by the spectral build, as its D65 color by the RGB build
                auto curve = std::make_shared<color::SpectralCurve>();
                std::string err;
                if (!io::loadReflectanceCsv(m.reflectance, *curve, &err)) {
                    std::cerr << "[SceneFile] material \"" << name << "\": " << m.reflectance << ": " << err << "\n";
                    return 1;
                }
                Real lo = Real(1e9), hi = Real(0);
                for (const auto& s : curve->samples()) { lo = std::min(lo, s.value); hi = std::max(hi, s.value); }
                std::cout << "[SceneFile] material \"" << name << "\": reflectance " << m.reflectance << " (" << curve->firstNm() << "-"
                    << curve->lastNm() << " nm, " << lo << " to " << hi << ")\n";
                if (hi > Real(1)) std::cout << "[SceneFile]   warning: R above 1 (a diffuse surface cannot reflect more than it receives)\n";
                materials[name] = std::make_shared<Lambertian>(std::shared_ptr<const color::SpectralCurve>(curve));
            }
            break;
        case MaterialSpec::Type::Conductor: {
            const std::string path = !m.reflectance.empty() ? m.reflectance
                : m.table.empty() ? "assets/Metal_data/Johnson_" + m.metal + ".csv" : m.table;
            auto& table = tables[path];
            if (!table) {
                table = std::make_shared<color::ComplexIorTable>();
                std::string err;
                if (!m.reflectance.empty()) {   // a measured reflectance R(lambda): the (n, k) that give this color at normal incidence
                    color::SpectralCurve reflectance;
                    if (!io::loadReflectanceCsv(path, reflectance, &err)) {
                        std::cerr << "[SceneFile] material \"" << name << "\": " << path << ": " << err << "\n";
                        return 1;
                    }
                    *table = color::complexIorTableFromReflectance(reflectance);
                }
                else if (!io::loadComplexIor_RefractiveIndexInfoCsv(path, *table, &err)) {
                    std::cerr << "[SceneFile] material \"" << name << "\": " << path << ": " << err << "\n";
                    return 1;
                }
            }
            auto metal = std::make_shared<RoughConductor>(table, Real(m.roughness), Real(m.anisotropy));
            metal->setTangentRotation(Real(m.rotation));
            metal->setMultipleScattering(m.multipleScattering);
            materials[name] = metal;
            break;
        }
        case MaterialSpec::Type::Dielectric: {
            auto glass = std::make_shared<Dielectric>(Real(m.ior), Real(m.roughness), Real(0), Real(m.abbe));
            if (!m.sellmeierPreset.empty() || m.hasSellmeier) {
                Sellmeier sm;
                if (!m.sellmeierPreset.empty()) sellmeierPreset(m.sellmeierPreset, sm);
                else for (int i = 0; i < 3; ++i) { sm.b[i] = Real(m.sellmeierB[i]); sm.c[i] = Real(m.sellmeierC[i]); }
                glass->setSellmeier(sm);
            }
            if (m.hasMedium) {   // a scattering medium inside the closed surface; equal channels make the cheaper gray medium
                PhaseFunction phase(Real(m.mediumG), Real(m.mediumG2), Real(m.mediumWeight));
                if (!m.mediumPhaseTable.empty()) {   // a tabulated phase function (Mie scattering of droplets...)
                    std::string err;
                    auto table = TabulatedPhase::fromCsv(m.mediumPhaseTable, &err);
                    if (!table) {
                        std::cerr << "[SceneFile] material \"" << name << "\": " << m.mediumPhaseTable << ": " << err << "\n";
                        return 1;
                    }
                    phase = PhaseFunction(table);
                }
                if (m.mediumNoise) {   // a density field: always the walk with null collisions
                    auto field = std::make_shared<NoiseDensity>(Real(m.noiseFrequency), m.noiseOctaves, Real(m.noiseLo), Real(m.noiseHi), std::uint32_t(m.noiseSeed));
                    if (m.maskRadius > 0) field->setSphereMask(m.maskCenter, Real(m.maskRadius), Real(m.maskStart));
                    glass->setMedium(std::make_shared<ChromaticMedium>(m.mediumSigmaT, m.mediumAlbedo, phase, field));
                }
                else if (m.mediumSigmaT.x == m.mediumSigmaT.y && m.mediumSigmaT.y == m.mediumSigmaT.z
                    && m.mediumAlbedo.x == m.mediumAlbedo.y && m.mediumAlbedo.y == m.mediumAlbedo.z)
                    glass->setMedium(std::make_shared<HomogeneousMedium>(Real(m.mediumSigmaT.x), Real(m.mediumAlbedo.x), phase));
                else
                    glass->setMedium(std::make_shared<ChromaticMedium>(m.mediumSigmaT, m.mediumAlbedo, phase));
            }
            if (m.farRoughness >= 0) glass->setDistanceRoughness(Real(m.farRoughness), Real(m.farStart), Real(m.farEnd));
            materials[name] = glass;
            break;
        }
        case MaterialSpec::Type::Layered: {   // a clear coat over a diffuse or a metal base
            LayeredCoat::Params lp;
            lp.ior = Real(m.coatIor);
            lp.roughness = Real(m.coatRoughness);
            lp.tint = m.coatTint;
            lp.samples = m.coatSamples;
            if (!m.baseIsConductor) {
                materials[name] = std::make_shared<LayeredCoat>(lp, RgbColor(m.albedo));
                break;
            }
            const std::string path = !m.reflectance.empty() ? m.reflectance
                : m.table.empty() ? "assets/Metal_data/Johnson_" + m.metal + ".csv" : m.table;
            // The base metal sits under the coat: its Fresnel reflectance uses the index relative to the coat, so it gets its own table entry.
            auto& table = tables[path];
            if (!table) {
                table = std::make_shared<color::ComplexIorTable>();
                std::string err;
                if (!m.reflectance.empty()) {
                    color::SpectralCurve reflectance;
                    if (!io::loadReflectanceCsv(path, reflectance, &err)) {
                        std::cerr << "[SceneFile] material \"" << name << "\": " << path << ": " << err << "\n";
                        return 1;
                    }
                    *table = color::complexIorTableFromReflectance(reflectance);
                }
                else if (!io::loadComplexIor_RefractiveIndexInfoCsv(path, *table, &err)) {
                    std::cerr << "[SceneFile] material \"" << name << "\": " << path << ": " << err << "\n";
                    return 1;
                }
            }
            auto base = std::make_shared<RoughConductor>(table, Real(m.roughness), Real(m.anisotropy), Real(m.coatIor));
            base->setTangentRotation(Real(m.rotation));
            base->setMultipleScattering(m.multipleScattering);
            materials[name] = std::make_shared<LayeredCoat>(lp, base);
            break;
        }
        case MaterialSpec::Type::ThinFilm: {
            FilmThickness th;
            th.base = Real(m.filmThickness);
            if (m.filmGradient) {
                th.gradient = true;
                th.axis = m.filmAxis;
                th.s0 = Real(m.filmRange0); th.s1 = Real(m.filmRange1);
                th.top = Real(m.filmTop);
                th.exponent = Real(m.filmExponent);
            }
            if (m.filmWarpAmp != 0) {
                th.warp = Real(m.filmWarpAmp);
                th.warpX = std::make_shared<NoiseDensity>(Real(m.filmWarpFreq), 3, Real(0.25), Real(0.75), std::uint32_t(m.filmWarpSeed));
                th.warpY = std::make_shared<NoiseDensity>(Real(m.filmWarpFreq), 3, Real(0.25), Real(0.75), std::uint32_t(m.filmWarpSeed + 101));
                th.warpZ = std::make_shared<NoiseDensity>(Real(m.filmWarpFreq), 3, Real(0.25), Real(0.75), std::uint32_t(m.filmWarpSeed + 202));
            }
            th.noiseStretch = m.filmNoiseStretch;
            if (m.filmNoiseAmp != 0) {
                th.noiseAmplitude = Real(m.filmNoiseAmp);
                th.noise = std::make_shared<NoiseDensity>(Real(m.filmNoiseFreq), m.filmNoiseOctaves, Real(0.3), Real(0.7), std::uint32_t(m.filmNoiseSeed));
            }
            std::shared_ptr<const DensityField> coverage;
            if (m.filmCoverage) coverage = std::make_shared<NoiseDensity>(Real(m.covFrequency), m.covOctaves, Real(m.covLo), Real(m.covHi), std::uint32_t(m.covSeed));
            if (m.substrateKind == "dielectric") materials[name] = std::make_shared<ThinFilm>(Real(m.filmIor), th, Real(m.substrateIor));
            else if (m.substrateKind == "conductor") {
                const std::string path = m.substrateTable.empty() ? "assets/Metal_data/Johnson_" + m.substrateMetal + ".csv" : m.substrateTable;
                auto& table = tables[path];
                if (!table) {
                    table = std::make_shared<color::ComplexIorTable>();
                    std::string err;
                    if (!io::loadComplexIor_RefractiveIndexInfoCsv(path, *table, &err)) {
                        std::cerr << "[SceneFile] material \"" << name << "\": " << path << ": " << err << "\n";
                        return 1;
                    }
                }
                materials[name] = std::make_shared<ThinFilm>(Real(m.filmIor), th, std::shared_ptr<const color::ComplexIorTable>(table));
            }
            else materials[name] = std::make_shared<ThinFilm>(Real(m.filmIor), th);
            if (coverage) std::static_pointer_cast<ThinFilm>(materials[name])->setCoverage(coverage);
            break;
        }
        case MaterialSpec::Type::Skin: {
            auto skin = std::make_shared<CoatedDiffuse>(RgbColor(m.albedo), Real(m.roughness), Real(m.skinSpecular), Real(0));
            if (!m.skinLotionNote.empty()) std::cout << "[Lotion] material \"" << name << "\": " << m.skinLotionNote << "\n";
            if (m.skinSss)   // the exit point is searched on every surface that shares this material
                skin->setSubsurface(std::make_shared<SubsurfaceProfile>(m.skinScatterMm * Real(0.001) * Real(m.skinSssScale), Real(1.4)));
            if (m.skinRelief) {
                SkinRelief::Params rp;
                rp.cell = Real(m.reliefCellMm * 1e-3); rp.furrowDepth = Real(m.reliefFurrowDepthMm * 1e-3); rp.furrowWidth = Real(m.reliefFurrowWidth);
                rp.poreDepth = Real(m.reliefPoreDepthMm * 1e-3); rp.poreRadius = Real(m.reliefPoreRadiusMm * 1e-3);
                rp.stretch = Real(m.reliefStretch); rp.coarseScale = Real(m.reliefCoarseScale); rp.coarseWeight = Real(m.reliefCoarseWeight);
                rp.poreSizeVariation = Real(m.reliefPoreSizeVariation); rp.poreDensityVariation = Real(m.reliefPoreDensityVariation); rp.poreDensityFreq = Real(m.reliefPoreDensityFreq);
                rp.poreFraction = Real(m.reliefPoreFraction); rp.seed = unsigned(m.reliefSeed);
                if (m.reliefWrinkles) {   // the furrows and pores come from the graph simulation: the Worley ones are off
                    rp.furrowDepth = 0; rp.coarseWeight = 0; rp.poreDepth = 0;
                }
                auto relief = std::make_shared<SkinRelief>(rp);
                if (m.reliefWrinkles) {
                    WrinkleField::Params wp;
                    wp.tileMm = Real(m.wrTileMm); wp.resolution = m.wrResolution; wp.steps = m.wrSteps; wp.depthScale = Real(m.wrDepthScale); wp.seed = unsigned(m.wrSeed);
                    for (const WrinkleLevelSpec& ls : m.wrLevels) {
                        WrinkleField::Level l;
                        l.spacingMm = Real(ls.spacingMm); l.angleDeg = Real(ls.angleDeg); l.base = Real(ls.base); l.distance = Real(ls.distance);
                        l.continuation = Real(ls.continuation); l.similarity = Real(ls.similarity); l.crossing = Real(ls.crossing); l.deposit = Real(ls.deposit);
                        l.width = Real(ls.width); l.poreWidth = Real(ls.poreWidth); l.poreMin = Real(ls.poreMin); l.poreBlend = Real(ls.poreBlend); l.poreSkew = Real(ls.poreSkew);
                        l.wobble = Real(ls.wobble); l.wobbleFreq = Real(ls.wobbleFreq); l.depthUm = Real(ls.depthUm); l.bulge = ls.bulge; l.darkPores = ls.darkPores;
                        wp.levels.push_back(l);
                    }
                    auto field = std::make_shared<WrinkleField>(wp);
                    DumpWrinkleMap(*field);
                    relief->setWrinkles(field);
                }
                skin->setRelief(relief);
                if (m.skinFilmOn) {
                    SkinFilm::Params fp;
                    fp.thickness = Real(m.skinFilmThicknessUm * 1e-6); fp.ior = Real(m.skinFilmIor); fp.roughness = Real(m.skinFilmRoughness);
                    fp.softness = Real(m.skinFilmSoftnessUm * 1e-6);
                    skin->setFilm(std::make_shared<SkinFilm>(fp, relief));
                    if (m.skinEmulsionTau > 0) skin->setEmulsion(Real(m.skinEmulsionTau));
                }
            }
            if (m.skinScatterOn) {
                SkinScatter::Params sp;
                sp.rho0 = Real(m.scatRho0); sp.frequency = Real(m.scatFrequency); sp.threshold = Real(m.scatThreshold); sp.softness = Real(m.scatSoftness);
                sp.hydration = Real(m.scatHydration); sp.kappaH = Real(m.scatKappaH); sp.gap = Real(m.scatGapUm * 1e-6); sp.seed = std::uint32_t(m.scatSeed);
                skin->setScatter(std::make_shared<SkinScatter>(sp));
            }
            if (m.skinFoundationOn) {
                CosmeticLayer::Params cp;
                cp.thickness = Real(m.fndThickness); cp.diffuserFraction = Real(m.fndDiffuserFraction); cp.g1 = Real(m.fndG1); cp.g2 = Real(m.fndG2);
                cp.lobeWeight = Real(m.fndLobeWeight); cp.plateletRoughness = Real(m.fndPlateletRoughness); cp.plateletTiltDeg = Real(m.fndPlateletTiltDeg);
                cp.diffuserAlbedo = m.fndDiffuserAlbedo; cp.plateletAlbedo = m.fndPlateletAlbedo;
                const auto coat = std::make_shared<CosmeticCoat>(cp);
                skin->setCosmetic(coat, Real(m.fndLeveling));
                const Vector3 r = coat->reflected(Real(1)), t = coat->through(Real(1));
                std::cout << "[Foundation] a layer of thickness " << cp.thickness << ", " << 100 * cp.diffuserFraction << " % diffusers, platelet roughness " << cp.plateletRoughness
                    << ": seen from the normal it reflects (" << r.x << ", " << r.y << ", " << r.z << "), lets through (" << t.x << ", " << t.y << ", " << t.z
                    << "), without a collision " << coat->straight(Real(1)) << "\n";
            }
            if (m.skinGlossVariation > 0) skin->setGlossVariation(Real(m.skinGlossVariation), Real(m.skinGlossFreq), std::uint32_t(m.skinGlossSeed));
            if (m.skinPigment || (m.skinRelief && m.reliefPoreDarkness > 0)) {
                std::shared_ptr<const SkinPigment> pigment;
                if (m.skinPigment) {
                    SkinPigment::Params gp;
                    gp.freckleCell = Real(m.pigFreckleCellMm * 1e-3); gp.freckleFraction = Real(m.pigFreckleFraction); gp.freckleRadius = Real(m.pigFreckleRadiusMm * 1e-3);
                    gp.freckleStrength = Real(m.pigFreckleStrength); gp.clusterFreq = Real(m.pigClusterFreq); gp.blotchFreq = Real(m.pigBlotchFreq);
                    gp.blotchStrength = Real(m.pigBlotchStrength); gp.rednessFreq = Real(m.pigRednessFreq); gp.rednessStrength = Real(m.pigRednessStrength);
                    gp.seed = unsigned(m.pigSeed);
                    pigment = std::make_shared<SkinPigment>(gp);
                }
                skin->setPigment(pigment, Real(m.reliefPoreDarkness));
            }
            if (m.skinWalk) {   // after setSubsurface(): the medium is made from the albedo and the profile's scatter distance
                skin->setSubsurfaceWalk(Real(m.skinWalkG));
                if (const SubsurfaceWalkMedium* w = skin->subsurfaceWalk()) {
                    const Vector3 sT = w->sigmaTRgb(), sS = w->sigmaSRgb();
                    std::cout << "[Skin] material \"" << name << "\": random walk inside the object, extinction (" << sT.x << ", " << sT.y << ", " << sT.z
                        << ") 1/m, single-scattering albedo (" << sS.x / sT.x << ", " << sS.y / sT.y << ", " << sS.z / sT.z << "), g " << w->anisotropy() << "\n";
                }
            }
            materials[name] = skin;
            break;
        }
        case MaterialSpec::Type::Hair: {
            Hair::Params hp;
            hp.eta = Real(m.hairEta); hp.sigmaA = m.hairSigmaA; hp.betaM = Real(m.hairBetaM); hp.betaN = Real(m.hairBetaN); hp.alphaDeg = Real(m.hairTiltDeg); hp.lobeMask = m.hairLobes;
            materials[name] = std::make_shared<Hair>(hp);
            break;
        }
        case MaterialSpec::Type::Emission: {
            EmissionSpectrum spectrum;
            if (!MakeSpectrum(m.spectrum, spectrum)) return 1;
            materials[name] = std::make_shared<DiffuseLight>(spectrum, Real(m.radiance), m.twoSided);
            break;
        }
        }
    }

    // ---- objects ----
    for (const ObjectSpec& o : d.objects) {
        const MaterialSpec& spec = d.materials.at(o.material);
        const std::shared_ptr<Material>& mat = materials.at(o.material);
        const bool emission = spec.type == MaterialSpec::Type::Emission;
        const bool asLight = emission && o.light != 0;   // an emitter is a light for next-event estimation unless it says "light": false
        switch (o.shape) {
        case ObjectSpec::Shape::Sphere: {
            auto sphere = std::make_shared<Sphere>(o.center, o.radius, mat);
            if (!emission) world->add(sphere);
            else if (asLight) {
                auto light = std::make_shared<AreaLight>(sphere, mat);
                world->add(std::make_shared<Primitive>(sphere, mat, light));
                out.lights.push_back(light);
            }
            else world->add(std::make_shared<Primitive>(sphere, mat));
            break;
        }
        case ObjectSpec::Shape::Quad:
            if (!emission) world->add(std::make_shared<Quad>(o.corner, o.edge1, o.edge2, mat));
            else {
                EmissionSpectrum spectrum;
                if (!MakeSpectrum(spec.spectrum, spectrum)) return 1;
                const auto light = AddQuadLight(*world, nullptr, o.corner, o.edge1, o.edge2, spectrum, Real(spec.radiance), asLight, spec.twoSided);
                if (asLight && light) out.lights.push_back(light);
            }
            break;
        case ObjectSpec::Shape::Mesh: {
            TriangleMesh mesh;
            std::string err;
            if (!io::loadObj(o.file, mesh, &err)) { std::cerr << "[SceneFile] mesh " << o.file << ": " << err << "\n"; return 1; }
            const glm::dmat3 rotation(glm::rotate(glm::dmat4(1.0), glm::radians(Real(o.rotateDegrees)), Vector3(o.rotateAxis)));
            mesh.transform(rotation * Real(o.scale), Vector3(o.translate));
            world->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(mesh)), mat, opt.accel));
            break;
        }
        case ObjectSpec::Shape::Lathe: {
            TriangleMesh turned = MakeLatheMesh(o.profile, o.segments);
            turned.transform(glm::dmat3(1.0), Vector3(o.translate));
            world->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(turned)), mat, opt.accel));
            break;
        }
        case ObjectSpec::Shape::Hair: {
            const HairSpec& hs = o.hair;
            auto curves = std::make_shared<HairCurves>(mat);
            if (hs.style == "lock") {
                hairgen::LockParams lp;
                lp.strands = hs.strands; lp.length = Real(hs.length); lp.segments = hs.segments;
                lp.rootRadius = Real(hs.rootRadius); lp.tipRadius = Real(hs.tipRadius);
                lp.width = Real(hs.width); lp.tipWidth = Real(hs.tipWidth); lp.waveAmp = Real(hs.waveAmp); lp.waveLength = Real(hs.waveLength);
                lp.twist = Real(hs.twist); lp.frizz = Real(hs.frizz); lp.flyaway = Real(hs.flyaway); lp.root = hs.root; lp.seed = unsigned(hs.seed);
                if (hs.clumps > 0) lp.clumps = hs.clumps;
                hairgen::MakeLock(*curves, lp);
            }
            else {
                hairgen::HeadHairParams hp;
                hp.headCenter = hs.center; hp.headRadius = Real(hs.headRadius);
                hp.strands = hs.strands; hp.length = Real(hs.length); hp.segments = hs.segments;
                hp.rootRadius = Real(hs.rootRadius); hp.tipRadius = Real(hs.tipRadius);
                hp.droop = Real(hs.droop); hp.outward = Real(hs.outward); hp.waveAmp = Real(hs.sway); hp.waveFreq = Real(hs.swayFreq); hp.frizz = Real(hs.headFrizz);
                hp.maxPolarDeg = Real(hs.maxPolarDeg); hp.hairlineHeight = Real(hs.hairline); hp.seed = unsigned(hs.seed);
                hp.clumpStrength = Real(hs.clumpStrength);
                if (hs.clumps > 0) hp.clumps = hs.clumps;
                hairgen::MakeHeadHair(*curves, hp);
            }
            curves->build();
            std::cout << "[Hair] " << curves->numStrands() << " strands, " << curves->numSegments() << " segments (" << hs.style << ")\n";
            world->add(curves);
            break;
        }
        case ObjectSpec::Shape::Ocean: {
            OceanParams op;
            op.windSpeed = Real(o.windSpeed); op.windDeg = Real(o.windDirection); op.spread = Real(o.spread);
            op.peakEnhance = Real(o.peakEnhance); op.chop = Real(o.chop); op.time = Real(o.time); op.minWavelength = Real(o.minWavelength);
            op.seed = static_cast<unsigned>(o.oceanSeed); op.rings = o.rings; op.sectors = o.sectors;
            op.rInner = Real(o.rInner); op.rOuter = Real(o.rOuter);
            TriangleMesh sea = MakeOceanMesh(op);
            std::cout << "[Ocean] " << sea.numTriangles() << " triangles, wind " << o.windSpeed << " m/s toward " << o.windDirection << " deg\n";
            world->add(MakeMeshHittable(std::make_shared<const TriangleMesh>(std::move(sea)), mat, opt.accel));
            break;
        }
        }
    }

    // ---- point lights ----
    for (const PointLightSpec& l : d.lights)
        out.lights.push_back(std::make_shared<PointLight>(l.position, RgbColor(l.color), Real(l.intensity)));
    return std::nullopt;
}

} // namespace rayt::cli
