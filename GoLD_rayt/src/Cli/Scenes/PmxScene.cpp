// An MMD (PMX) character: model, pose, expressions, materials, studio lights.
#include "pch.h"
#include "Cli/Scenes.hpp"
#include "Cli/Checks.hpp"
#include "Cli/LightPresets.hpp"
#include "Cli/PmxSupport.hpp"
#include "Cli/SceneParts.hpp"
#include "Cli/Tools.hpp"
#include "Geometry/ProceduralMesh.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

// --scene pmx
std::optional<int> BuildPmx(SceneContext& ctx, SceneSetup& out) {
    const CliOptions& opt = ctx.opt;
    auto& worldList = out.world;
    const auto& envLight = ctx.env.light;
    const auto& matFloor = ctx.mats.floor;
    auto& studioLights = out.lights;
    Point3 pmxCamFrom(0.0), pmxCamTarget(0.0);                    // camera framing the model

    const auto files = FindPmxFiles();
    if (files.empty()) { std::cerr << "[PMX] no .pmx files found in MMD_model/\n"; return 1; }
    size_t pick = 0;
    if (opt.pmxIndex >= 0) pick = std::min(static_cast<size_t>(opt.pmxIndex), files.size() - 1);
    else for (size_t i = 1; i < files.size(); ++i) if (std::filesystem::file_size(files[i]) > std::filesystem::file_size(files[pick])) pick = i;

    io::PmxModel model;
    std::string pmxErr;
    if (!io::loadPmx(files[pick], model, &pmxErr)) { std::cerr << "[PMX] load failed: " << pmxErr << "\n"; return 1; }
    PrintPmxInfo(model, files[pick]);
    if (opt.pmxCheck) return RunPmxChecks(model) ? 0 : 2;
    if (opt.morphList) {
        static const char* kTypes[] = { "group", "vertex", "bone", "uv", "uv1", "uv2", "uv3", "uv4", "material", "flip", "impulse" };
        static const char* kPanels[] = { "-", "brow", "eye", "mouth", "other" };
        for (size_t i = 0; i < model.morphs.size(); ++i) {
            const auto& m = model.morphs[i];
            std::cout << "[Morph] " << i << ": " << m.name << "  [" << kPanels[std::clamp(m.panel, 0, 4)] << ", " << kTypes[std::clamp(m.type, 0, 10)]
                << ", " << m.offsetCount() << " offsets]\n";
        }
        return 0;
    }

    // ---- pose file (bones, and expression entries) ----
    io::VpdPose vpd;
    bool havePose = false;
    std::string poseLabel;
    if (opt.poseIndex >= 0) {
        const auto poses = FindMmdFiles("pose", ".vpd");
        if (poses.empty()) { std::cerr << "[Pose] no .vpd files found in MMD_model/pose\n"; return 1; }
        const auto& poseFile = poses[std::min(static_cast<size_t>(opt.poseIndex), poses.size() - 1)];
        std::string poseErr;
        if (!io::loadVpd(poseFile, vpd, &poseErr)) { std::cerr << "[Pose] load failed: " << poseErr << "\n"; return 1; }
        havePose = true;
        poseLabel = PathToUtf8(poseFile.parent_path().filename()) + "/" + PathToUtf8(poseFile.filename());
    }

    // ---- morph weights: --morph i=w,..  /  --expr name  /  Morph blocks of the pose file ----
    std::map<int, Real> morphWeights;
    auto morphByName = [&](const std::string& name) {
        for (size_t i = 0; i < model.morphs.size(); ++i) if (model.morphs[i].name == name) return static_cast<int>(i);
        return -1;
    };
    auto setMorph = [&](const std::string& key, Real w, const char* from) {
        int idx = -1;
        const bool numeric = !key.empty() && std::all_of(key.begin(), key.end(), [](unsigned char c) { return std::isdigit(c) != 0; });
        if (numeric) idx = std::stoi(key);
        else idx = morphByName(key);
        if (idx < 0 || idx >= static_cast<int>(model.morphs.size())) { std::cout << "[Morph] " << from << ": no morph " << key << " in this model\n"; return; }
        morphWeights[idx] = w;
    };
    if (!opt.morphSpec.empty()) {
        std::stringstream ss(opt.morphSpec);
        std::string item;
        while (std::getline(ss, item, ',')) {
            const size_t eq = item.find('=');
            if (eq != std::string::npos) setMorph(item.substr(0, eq), std::stod(item.substr(eq + 1)), "--morph");
        }
    }
    if (!opt.exprFile.empty()) {
        const auto files = FindMmdFiles("expression", ".txt");
        const auto it = std::find_if(files.begin(), files.end(), [&](const std::filesystem::path& f) { return PathToUtf8(f.stem()) == opt.exprFile; });
        std::vector<std::pair<std::string, Real>> items;
        std::string exprErr;
        if (it == files.end() || !io::loadExpression(*it, items, &exprErr)) { std::cerr << "[Morph] expression " << opt.exprFile << " not found in MMD_model/expression\n"; return 1; }
        for (const auto& [name, w] : items) setMorph(name, w, "--expr");
    }
    if (havePose) for (const auto& [name, w] : vpd.morphs) setMorph(name, w, "pose file");

    const io::PmxModel restModel = opt.pmxBench ? model : io::PmxModel();   // the rest pose, kept for --pmx-bench
    std::vector<io::VpdBone> morphBones;
    if (!morphWeights.empty()) {
        model = io::applyMorphs(model, morphWeights, &morphBones);
        std::cout << "[Morph] applied:";
        for (const auto& [idx, w] : morphWeights) std::cout << " " << model.morphs[static_cast<size_t>(idx)].name << "=" << w;
        std::cout << "\n";
    }
    if (havePose || !morphBones.empty()) {
        io::mergePoseBones(vpd, morphBones);
        io::PoseReport rep;
        model = io::applyPose(model, vpd, &rep);
        std::cout << "[Pose] " << (havePose ? poseLabel : std::string("(bone morphs only)")) << "  made for " << vpd.parentModel
            << ": " << vpd.bones.size() << " bones, " << rep.applied << " applied";
        if (!rep.missing.empty()) {
            std::cout << ", " << rep.missing.size() << " not in this model:";
            for (const auto& n : rep.missing) std::cout << " " << n;
        }
        std::cout << "\n";
    }
    if (opt.pmxInfoOnly) return 0;

    Real height = 1;
    Vector3 modelSize(1.0);
    const std::shared_ptr<TriangleMesh> editableMesh = PmxToMesh(model, opt.pmxYaw, &height, &modelSize);
    if (opt.pmxBench) {
        // The cost of one pose update, stage by stage: what an interactive posing tool pays for every new pose,
        // before a single ray is traced.
        using clock = std::chrono::steady_clock;
        const int reps = 10;
        auto timeIt = [&](const char* what, const std::function<void()>& fn) {
            fn();   // warm-up
            const auto t0 = clock::now();
            for (int i = 0; i < reps; ++i) fn();
            const double ms = std::chrono::duration<double, std::milli>(clock::now() - t0).count() / reps;
            std::cout << "[PmxBench]   " << what << ": " << ms << " ms\n";
            return ms;
        };
        std::cout << "[PmxBench] " << model.vertices.size() << " vertices, " << model.indices.size() / 3 << " triangles, "
            << model.bones.size() << " bones; mean of " << reps << " runs, one thread\n";
        double common = 0;
        io::PmxModel morphed = restModel;
        if (!morphWeights.empty()) {
            common += timeIt("morphs (vertex / UV / material / bone)", [&]() { std::vector<io::VpdBone> b; morphed = io::applyMorphs(restModel, morphWeights, &b); });
        }
        common += timeIt("pose: bones (FK + grant) and linear blend skinning", [&]() { io::PoseReport r; const io::PmxModel posed = io::applyPose(morphed, vpd, &r); (void)posed; });
        common += timeIt("convert to the renderer mesh (mirror, rotate, normals)", [&]() { Real h = 1; Vector3 sz(1.0); const auto m = PmxToMesh(model, opt.pmxYaw, &h, &sz); (void)m; });
        const std::shared_ptr<const TriangleMesh> constMesh = editableMesh;
        auto white = std::make_shared<Lambertian>(RgbColor(1.0));
        const double ownBvh = timeIt("own SAH BVH build (incl. one Triangle object per face)", [&]() { auto b = BuildMeshBVH(constMesh, white); (void)b; });
        const double embreeBuild = timeIt("Embree scene build (incl. one Triangle object per face)", [&]() { EmbreeMesh e(constMesh, white); });
        std::cout << "[PmxBench] one pose update: " << common + ownBvh << " ms with the own BVH, " << common + embreeBuild << " ms with Embree\n";
        return 0;
    }
    {
        // A translucent layer (blush ...) is a copy of the surface it covers, lying exactly on it. A ray that
        // passes through the layer must still find the surface below, so the layer is lifted 0.5 mm along
        // its normals (far more than the 1e-5 ray offset, invisible in the picture).
        std::vector<char> lifted(editableMesh->p.size(), 0);
        for (const auto& pm : model.materials) {
            const bool overlayName = !pm.name.empty() && pm.name.back() == '+';
            if (pm.diffuse.w <= Real(0.001) || !(pm.alphaFromMorph && (overlayName || pm.diffuse.w < Real(0.999)))) continue;
            for (int i = pm.firstIndex; i < pm.firstIndex + pm.indexCount; ++i) {
                const size_t v = static_cast<size_t>(editableMesh->indices[static_cast<size_t>(i)]);
                if (!lifted[v]) { editableMesh->p[v] += editableMesh->n[v] * Real(5e-4); lifted[v] = 1; }
            }
        }
    }
    const std::shared_ptr<const TriangleMesh> mesh = editableMesh;
    std::vector<std::shared_ptr<Hittable>> prims;
    std::map<int, std::shared_ptr<Texture2D>> textureCache;   // PMX texture index -> loaded texture (nullptr = failed)
    // Skin scatter distance (mean free path) in meters: 3.7 / 1.4 / 0.7 mm for R / G / B, times --sss-scale.
    std::shared_ptr<const SubsurfaceProfile> skinProfile;
    if (opt.sss && !opt.pmxMatte) {
        skinProfile = std::make_shared<SubsurfaceProfile>(Vector3(3.67, 1.37, 0.68) * Real(0.001) * opt.sssScale, Real(1.4));
        std::cout << "[PMX] subsurface scattering on the skin materials (scale " << opt.sssScale << ")\n";
    }
    int skippedOverlay = 0, skippedHidden = 0, alphaMaterials = 0, missingTextures = 0, translucentLayers = 0;
    for (const auto& pm : model.materials) {
        if (pm.diffuse.w <= Real(0.001)) { ++skippedHidden; continue; }   // fully transparent material
        // MMD models stack extra layers on top of the real surface and switch them with morphs; names
        // ending in '+' mark them (sphere-map highlights, blush ...). A layer a morph has switched on
        // (diffuse alpha changed by a morph) is drawn as a translucent layer; the others would cover the
        // real surface as opaque diffuse, so they are left out unless --pmx-all.
        const bool overlayName = !pm.name.empty() && pm.name.back() == '+';
        bool translucent = pm.alphaFromMorph && (overlayName || pm.diffuse.w < Real(0.999));
        if (overlayName && !pm.alphaFromMorph && !opt.pmxAll) { ++skippedOverlay; continue; }

        std::shared_ptr<Texture2D> tex;
        if (pm.texture >= 0) {
            auto it = textureCache.find(pm.texture);
            if (it == textureCache.end()) {
                std::string rel = model.textures[static_cast<size_t>(pm.texture)];
                std::replace(rel.begin(), rel.end(), '\\', '/');
                it = textureCache.emplace(pm.texture, Texture2D::load(files[pick].parent_path() / PathFromUtf8(rel))).first;
                if (!it->second) { ++missingTextures; std::cout << "[PMX] could not load texture " << rel << "\n"; }
            }
            tex = it->second;
        }

        const RgbColor baseColor(SrgbToLinear(pm.diffuse.x), SrgbToLinear(pm.diffuse.y), SrgbToLinear(pm.diffuse.z));
        std::shared_ptr<Material> mat;
        if (translucent) {
            auto layer = std::make_shared<Lambertian>(baseColor, tex);
            layer->setTranslucent(std::clamp(pm.diffuse.w, Real(0), Real(1)));
            mat = layer;
            ++translucentLayers;
        }
        else {
            mat = MakePmxMaterial(pm, baseColor, tex, opt.pmxMatte, skinProfile);
        }
        const bool cutout = tex && tex->hasAlpha() && !translucent;   // a translucent layer uses the texture alpha as opacity instead
        if (cutout) ++alphaMaterials;
        for (int t = pm.firstIndex / 3; t < (pm.firstIndex + pm.indexCount) / 3; ++t) {
            const Vector3 e1 = mesh->p[static_cast<size_t>(mesh->indices[3 * t + 1])] - mesh->p[static_cast<size_t>(mesh->indices[3 * t])];
            const Vector3 e2 = mesh->p[static_cast<size_t>(mesh->indices[3 * t + 2])] - mesh->p[static_cast<size_t>(mesh->indices[3 * t])];
            if (glm::length(glm::cross(e1, e2)) < 1e-12) continue;   // degenerate triangle
            auto tri = std::make_shared<Triangle>(mesh, static_cast<size_t>(t), mat);
            if (cutout) tri->setAlphaTexture(tex, Real(0.5));
            if (translucent) tri->setCastsShadow(false);
            prims.push_back(tri);
        }
    }
    size_t texBytes = 0;
    for (const auto& [idx, t] : textureCache) if (t) texBytes += static_cast<size_t>(t->width()) * static_cast<size_t>(t->height()) * 4;
    std::cout << "[PMX] " << prims.size() << " triangles in the scene, height " << height << " m;  " << textureCache.size() - missingTextures
        << " textures loaded (" << texBytes / (1024 * 1024) << " MB), " << alphaMaterials << " materials with alpha cut-outs, skipped "
        << skippedOverlay << " overlay ('+') and " << skippedHidden << " fully transparent materials, " << translucentLayers << " translucent layers\n";
    worldList->add(std::make_shared<BVH>(std::move(prims)));
    worldList->add(std::make_shared<Sphere>(Point3(0, -100, 0), 100.0, matFloor));
    // Framing: the picture is 3:4 with a 35 degree vertical field of view; a wide pose (a raised leg) is
    // fitted by its width as well as its height.
    const Real zoom = std::max(opt.pmxZoom, Real(0.1));
    const Real extent = std::max(height, modelSize.x / Real(0.75));
    pmxCamTarget = Point3(0, height * opt.pmxFocus, 0);
    pmxCamFrom = Point3(0, height * (opt.pmxFocus + 0.05 / zoom), extent * 1.22 / (2 * std::tan(glm::radians(17.5))) / zoom);
    if (opt.pmxStudio) {
        // Studio: a dim backdrop, key / fill soft boxes and two rim lights; the HDR only adds a faint ambient.
        // Looks (--studio-look): color temperature (K; 0 = neutral white) and radiance of each light,
        // the ambient HDR scale and the backdrop reflectance. Positions and sizes are shared.
        struct Lamp { Real kelvin, radiance; };
        struct StudioLook { Lamp key, fill, rimLeft, rimRight; Real env; RgbColor wall; };
        static const std::map<std::string, StudioLook> looks = {
            { "default",  { { 5200, 13 },   { 0, 3.2 },     { 0, 28 },      { 7500, 22 },    0.10, RgbColor(0.16, 0.17, 0.21) } },
            { "warm",     { { 3200, 14 },   { 4500, 2.0 },  { 2500, 30 },   { 3800, 22 },    0.08, RgbColor(0.20, 0.16, 0.13) } },
            { "cool",     { { 8500, 8 },    { 10000, 3.0 }, { 12000, 28 },  { 9000, 22 },    0.12, RgbColor(0.13, 0.16, 0.22) } },
            { "dramatic", { { 4300, 20 },   { 0, 0.4 },     { 6000, 35 },   { 6500, 28 },    0.02, RgbColor(0.07, 0.07, 0.09) } },
            { "soft",     { { 5600, 7 },    { 5600, 6.0 },  { 5600, 10 },   { 5600, 8 },     0.35, RgbColor(0.45, 0.46, 0.50) } },
        };
        const StudioLook& look = looks.at(opt.studioLook);
        auto spectrum = [](const Lamp& l) {
            return l.kelvin > 0 ? EmissionSpectrum::Blackbody(l.kelvin) : EmissionSpectrum(RgbColor(1.0));
        };
        if (envLight) envLight->setScale(look.env);
        const Point3 target(0, height * 0.55, 0);
        auto softbox = [&](const Point3& center, Real w, Real h, const EmissionSpectrum& e, Real radiance) {
            const Vector3 n = glm::normalize(target - center);
            const Vector3 e1 = glm::normalize(glm::cross(Vector3(0, 1, 0), n)) * w;
            const Vector3 e2 = glm::normalize(glm::cross(n, e1)) * h;
            studioLights.push_back(AddQuadLight(*worldList, nullptr, center - e1 * 0.5 - e2 * 0.5, e1, e2, e, radiance, true));
        };
        const Real H = height;
        softbox(Point3(-1.7, H * 1.05, 2.3), 1.3, 1.3, spectrum(look.key), look.key.radiance);               // key: front-left, above
        softbox(Point3(2.3, H * 0.60, 2.4), 1.8, 1.2, spectrum(look.fill), look.fill.radiance);              // fill: front-right, soft
        softbox(Point3(-1.5, H * 1.10, -1.7), 0.5, 1.8, spectrum(look.rimLeft), look.rimLeft.radiance);      // rim: back-left strip
        softbox(Point3(1.8, H * 1.00, -1.6), 0.5, 1.8, spectrum(look.rimRight), look.rimRight.radiance);     // rim: back-right strip
        const Real wallZ = -(modelSize.z * 0.5 + 1.8);
        worldList->add(std::make_shared<Quad>(Point3(-14, 0, wallZ), Vector3(28, 0, 0), Vector3(0, 12, 0),
            std::make_shared<Lambertian>(look.wall)));
        // Dark glossy studio floor (normal +y: cross(z, x)) lying 1 mm above the plain floor: soft reflections of the lights and the figure.
        worldList->add(std::make_shared<Quad>(Point3(-14, 0.001, wallZ), Vector3(0, 0, 16), Vector3(28, 0, 0),
            std::make_shared<CoatedDiffuse>(RgbColor(0.09, 0.09, 0.11), Real(0.20), Real(1.0), Real(0))));
    }
    if (opt.pmxLight) {
        // Softbox 1.0 x 0.7 m, 2.2 m from the focus point, front-left and a little above, facing it.
        const Point3 target = pmxCamTarget;
        const Point3 center = target + Vector3(-1.3, 0.9, 1.6);
        const Vector3 n = glm::normalize(target - center);
        const Vector3 e1 = glm::normalize(glm::cross(Vector3(0, 1, 0), n)) * 1.0;
        const Vector3 e2 = glm::normalize(glm::cross(n, e1)) * 0.7;
        studioLights.push_back(AddQuadLight(*worldList, nullptr, center - e1 * 0.5 - e2 * 0.5, e1, e2,
            EmissionSpectrum(RgbColor(1.0)), 8.0, true));
    }
    out.camera.from = pmxCamFrom;
    out.camera.at = pmxCamTarget;
    return std::nullopt;
}

} // namespace rayt::cli
