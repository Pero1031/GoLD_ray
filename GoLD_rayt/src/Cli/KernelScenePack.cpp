// The kernel renderer's scene from a scene file, and the --device cpu-kernel | gpu render (Cli/KernelScenePack.hpp, docs/GPU.md).
#include "pch.h"
#include "Cli/KernelScenePack.hpp"

#include "Accelerators/EmbreeMesh.hpp"   // detail::embreeDevice
#include "Geometry/ProceduralMesh.hpp"
#include "Geometry/TriangleMesh.hpp"
#include "Gpu/GpuDevice.hpp"
#include "IO/PlyLoader.hpp"
#include "Kernel/PathTracer.hpp"
#include "Render/Camera.hpp"
#include "Render/Film.hpp"

#include <embree4/rtcore.h>
#include <omp.h>

#include <chrono>
#include <cmath>
#include <glm/gtc/matrix_transform.hpp>

namespace rayt::cli {

namespace {

    kernel::KFloat3 F3(const Vector3& v) { return { float(v.x), float(v.y), float(v.z) }; }

    /// Appends a triangle mesh (world space) with one material; the face normals follow the engine's Triangle (Geometry/TriangleMesh.hpp).
    void AddMesh(const TriangleMesh& mesh, uint32_t material, KernelSceneHost& out) {
        const uint32_t base = uint32_t(out.positions.size());
        const bool smooth = mesh.hasNormals();
        for (size_t i = 0; i < mesh.p.size(); ++i) {
            out.positions.push_back(F3(mesh.p[i]));
            out.normals.push_back(smooth ? F3(mesh.n[i]) : kernel::KFloat3{ 0, 0, 0 });
        }
        for (size_t t = 0; t < mesh.numTriangles(); ++t) {
            const int i0 = mesh.indices[3 * t], i1 = mesh.indices[3 * t + 1], i2 = mesh.indices[3 * t + 2];
            const Vector3 e1 = mesh.p[size_t(i1)] - mesh.p[size_t(i0)], e2 = mesh.p[size_t(i2)] - mesh.p[size_t(i0)];
            const Vector3 c = glm::cross(e1, e2);
            const Real area = Real(0.5) * glm::length(c);
            Vector3 ng = area > Real(0) ? c / (Real(2) * area) : Vector3(0, 1, 0);
            if (smooth && glm::dot(ng, mesh.n[size_t(i0)] + mesh.n[size_t(i1)] + mesh.n[size_t(i2)]) < Real(0)) ng = -ng;
            out.indices.insert(out.indices.end(), { base + uint32_t(i0), base + uint32_t(i1), base + uint32_t(i2) });
            out.faceNormals.push_back(F3(ng));
            out.triMaterial.push_back(material);
            out.triLight.push_back(-1);
            out.triSmooth.push_back(smooth ? 1 : 0);
        }
    }

    /// Appends a quad (two triangles that face cross(edge1, edge2), as the engine's Quad); light >= 0: the quad light they belong to.
    void AddQuad(const Point3& c, const Vector3& e1, const Vector3& e2, uint32_t material, int32_t light, KernelSceneHost& out) {
        const uint32_t base = uint32_t(out.positions.size());
        for (const Point3& p : { c, c + e1, c + e1 + e2, c + e2 }) {
            out.positions.push_back(F3(p));
            out.normals.push_back({ 0, 0, 0 });
        }
        const kernel::KFloat3 n = F3(glm::normalize(glm::cross(e1, e2)));
        const uint32_t tris[2][3] = { { 0, 1, 2 }, { 0, 2, 3 } };
        for (const auto& t : tris) {
            out.indices.insert(out.indices.end(), { base + t[0], base + t[1], base + t[2] });
            out.faceNormals.push_back(n);
            out.triMaterial.push_back(material);
            out.triLight.push_back(light);
            out.triSmooth.push_back(0);
        }
    }

    /// Reasons a diffuse material cannot be modeled (empty: it can).
    std::string DiffuseLimits(const MaterialSpec& m) {
        if (!m.texture.empty()) return "an albedo texture";
        if (!m.reflectanceTexture.empty()) return "a reflectance texture";
        if (!m.reflectance.empty()) return "a measured reflectance";
        if (m.fadingDye || !m.fadingSeries.empty()) return "fading";
        return "";
    }

    // ---- the CPU kernel renderer: Embree in float, OpenMP over rows, the same tracePath as the GPU ----

    struct EmbreeTracer {
        RTCScene scene;
        kernel::KHit closest(const kernel::V3& o, const kernel::V3& d, float tmin, float tmax) const {
            RTCRayHit rh;
            rh.ray.org_x = o.x; rh.ray.org_y = o.y; rh.ray.org_z = o.z;
            rh.ray.dir_x = d.x; rh.ray.dir_y = d.y; rh.ray.dir_z = d.z;
            rh.ray.tnear = tmin; rh.ray.tfar = tmax;
            rh.ray.time = 0.0f; rh.ray.mask = 0xffffffffu; rh.ray.id = 0; rh.ray.flags = 0;
            rh.hit.geomID = RTC_INVALID_GEOMETRY_ID;
            rh.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
            rtcIntersect1(scene, &rh);
            if (rh.hit.geomID == RTC_INVALID_GEOMETRY_ID) return { 0.0f, 0.0f, 0.0f, kernel::kNoHit };
            return { rh.ray.tfar, rh.hit.u, rh.hit.v, rh.hit.primID };
        }
        bool occluded(const kernel::V3& o, const kernel::V3& d, float tmin, float tmax) const {
            RTCRay r;
            r.org_x = o.x; r.org_y = o.y; r.org_z = o.z;
            r.dir_x = d.x; r.dir_y = d.y; r.dir_z = d.z;
            r.tnear = tmin; r.tfar = tmax;
            r.time = 0.0f; r.mask = 0xffffffffu; r.id = 0; r.flags = 0;
            rtcOccluded1(scene, &r);
            return r.tfar < 0.0f;   // Embree sets tfar = -inf on a hit
        }
    };

    void GpuProgress(uint32_t done, uint32_t total) {
        if (done % 8 == 0 || done == total) std::cout << "\rProgress: " << done << "/" << total << " samples per pixel" << std::flush;
    }

} // namespace

bool RenderKernelCpu(const kernel::KScene& s, std::vector<float>& mean, std::string& why, bool progress) {
    RTCDevice device = detail::embreeDevice();
    if (!device) { why = "Embree: rtcNewDevice failed"; return false; }
    RTCScene scene = rtcNewScene(device);
    rtcSetSceneFlags(scene, RTC_SCENE_FLAG_ROBUST);
    rtcSetSceneBuildQuality(scene, RTC_BUILD_QUALITY_HIGH);
    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
    // Copied into buffers Embree allocates (it reads vertices 16 bytes at a time and pads its own buffers for that).
    float* vb = static_cast<float*>(rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float), s.numVertices));
    memcpy(vb, s.positions, size_t(s.numVertices) * sizeof(kernel::KFloat3));
    unsigned* ib = static_cast<unsigned*>(rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned), s.numTriangles));
    memcpy(ib, s.indices, size_t(s.numTriangles) * 3 * sizeof(unsigned));
    rtcCommitGeometry(geom);
    rtcAttachGeometry(scene, geom);
    rtcReleaseGeometry(geom);
    rtcCommitScene(scene);

    const EmbreeTracer tracer{ scene };
    mean.assign(size_t(s.width) * s.height * 3, 0.0f);
    std::atomic<int> rows{ 0 };
#pragma omp parallel for schedule(dynamic, 1)
    for (int y = 0; y < int(s.height); ++y) {
        for (uint32_t x = 0; x < s.width; ++x) {
            // Summed sample by sample in float, in the same order as the GPU accumulates its launches.
            kernel::V3 sum(0.0f);
            for (uint32_t k = 0; k < s.spp; ++k) sum += kernel::tracePath(s, tracer, x, uint32_t(y), k);
            const kernel::V3 m = sum / float(s.spp);
            float* px = &mean[(size_t(y) * s.width + x) * 3];
            px[0] = m.x; px[1] = m.y; px[2] = m.z;
        }
        const int done = ++rows;
        if (progress && (done % 16 == 0 || done == int(s.height))) {
#pragma omp critical
            std::cout << "\rProgress: " << done << "/" << s.height << " rows" << std::flush;
        }
    }
    if (progress) std::cout << std::endl;
    rtcReleaseScene(scene);
    return true;
}

bool RenderKernelGpu(const kernel::KScene& s, std::vector<float>& mean, std::string& why, bool progress, bool* unavailable, std::string* device) {
    if (unavailable) *unavailable = false;
    const gpu::Module* mod = gpu::LoadModule(why);
    if (!mod) { if (unavailable) *unavailable = true; why = "no GPU module: " + why; return false; }
    RaytGpuInfo info{};
    if (mod->query(&info) != 0) { if (unavailable) *unavailable = true; why = std::string("no usable GPU: ") + info.error; return false; }
    if (device) *device = std::string(info.name) + " (compute capability " + std::to_string(info.computeMajor) + "." + std::to_string(info.computeMinor) +
        ", OptiX, " + (info.fastMath ? "fast math" : "IEEE float, no implicit FMA") + ")";
    mean.assign(size_t(s.width) * s.height * 3, 0.0f);
    char err[512] = {};
    if (mod->render(&s, mean.data(), progress ? &GpuProgress : nullptr, err, int(sizeof(err))) != 0) { why = err; return false; }
    if (progress) std::cout << std::endl;
    return true;
}

kernel::KScene KernelSceneHost::view() const {
    kernel::KScene s = header;
    s.positions = positions.data();
    s.normals = normals.data();
    s.indices = indices.data();
    s.faceNormals = faceNormals.data();
    s.triMaterial = triMaterial.data();
    s.triLight = triLight.data();
    s.triSmooth = triSmooth.data();
    s.materials = materials.data();
    s.lights = lights.data();
    s.lightCdf = lightCdf.data();
    s.numVertices = uint32_t(positions.size());
    s.numTriangles = uint32_t(triMaterial.size());
    s.numMaterials = uint32_t(materials.size());
    s.numLights = uint32_t(lights.size());
    return s;
}

bool PackKernelScene(const SceneDescription& d, const Camera& camera, const CliOptions& opt, int width, int height, KernelSceneHost& out, std::string& why) {
    std::vector<std::string> missing;
    out = KernelSceneHost{};

    // ---- camera ----
    Point3 origin, lowerLeft;
    Vector3 horizontal, vertical;
    if (!camera.pinholeFrame(origin, lowerLeft, horizontal, vertical))
        missing.push_back("the camera is not a plain pinhole (a lens, an aperture, distortion, vignetting, motion or an orthographic view)");
    out.header.camera = { F3(origin), F3(lowerLeft), F3(horizontal), F3(vertical) };

    // ---- environment and point lights ----
    if (d.environment.type == EnvironmentSpec::Type::Hdr) missing.push_back("an environment map (only a uniform environment)");
    const float env = d.environment.type == EnvironmentSpec::Type::Uniform ? float(d.environment.scale) : 0.0f;
    out.header.envRadiance[0] = out.header.envRadiance[1] = out.header.envRadiance[2] = env;
    if (!d.lights.empty()) missing.push_back("point lights");

    // ---- materials ----
    std::map<std::string, uint32_t> index;
    std::map<std::string, const MaterialSpec*> specs;
    for (const auto& [name, m] : d.materials) {
        kernel::KMaterial km{};
        if (m.type == MaterialSpec::Type::Diffuse) {
            const std::string limit = DiffuseLimits(m);
            if (!limit.empty()) { missing.push_back("material \"" + name + "\": " + limit); continue; }
            km.type = kernel::kMatDiffuse;
            km.albedo[0] = float(m.albedo.x); km.albedo[1] = float(m.albedo.y); km.albedo[2] = float(m.albedo.z);
        }
        else if (m.type == MaterialSpec::Type::Emission) {
            EmissionSpectrum spectrum;
            if (!MakeEmissionSpectrum(m.spectrum, spectrum)) { missing.push_back("material \"" + name + "\": its spectrum cannot be read"); continue; }
            const Vector3 rgb = spectrum.rgb() * Real(m.radiance);   // what the RGB build's light emits
            km.type = kernel::kMatEmission;
            km.twoSided = m.twoSided ? 1u : 0u;
            km.emission[0] = float(rgb.x); km.emission[1] = float(rgb.y); km.emission[2] = float(rgb.z);
        }
        else { missing.push_back("material \"" + name + "\": only diffuse and emission materials"); continue; }
        index[name] = uint32_t(out.materials.size());
        specs[name] = &m;
        out.materials.push_back(km);
    }

    // ---- objects ----
    std::vector<float> power;
    for (size_t k = 0; k < d.objects.size(); ++k) {
        const ObjectSpec& o = d.objects[k];
        const std::string what = "object " + std::to_string(k + 1);
        const auto mat = index.find(o.material);
        if (mat == index.end()) continue;   // its material was refused above (or the scene reader rejected it before)
        const bool emission = out.materials[mat->second].type == kernel::kMatEmission;
        if (glm::dot(o.velocity, o.velocity) > 0) { missing.push_back(what + ": motion"); continue; }
        switch (o.shape) {
        case ObjectSpec::Shape::Quad: {
            int32_t light = -1;
            if (emission && o.light != 0) {
                if (!std::isnan(o.polarizer)) { missing.push_back(what + ": a polarized light"); continue; }
                const kernel::KMaterial& km = out.materials[mat->second];
                kernel::KQuadLight q{};
                q.corner = F3(o.corner); q.edge1 = F3(o.edge1); q.edge2 = F3(o.edge2);
                q.normal = F3(glm::normalize(glm::cross(o.edge1, o.edge2)));
                q.area = float(glm::length(glm::cross(o.edge1, o.edge2)));
                for (int c = 0; c < 3; ++c) q.radiance[c] = km.emission[c];
                q.twoSided = km.twoSided;
                light = int32_t(out.lights.size());
                out.lights.push_back(q);
                power.push_back(q.area * (0.2126f * q.radiance[0] + 0.7152f * q.radiance[1] + 0.0722f * q.radiance[2]) * (q.twoSided ? 2.0f : 1.0f));
            }
            AddQuad(o.corner, o.edge1, o.edge2, mat->second, light, out);
            break;
        }
        case ObjectSpec::Shape::Mesh:
        case ObjectSpec::Shape::Lathe: {
            if (emission) { missing.push_back(what + ": an emitting mesh (only emitting quads)"); continue; }
            TriangleMesh mesh;
            if (o.shape == ObjectSpec::Shape::Mesh) {
                if (!o.cutout.empty()) { missing.push_back(what + ": a cutout"); continue; }
                std::string err;
                if (!io::loadMesh(o.file, mesh, &err)) { missing.push_back(what + ": " + o.file + ": " + err); continue; }
                const glm::dmat3 rotation(glm::rotate(glm::dmat4(1.0), glm::radians(Real(o.rotateDegrees)), Vector3(o.rotateAxis)));
                mesh.transform(rotation * Real(o.scale), Vector3(o.translate));   // as Scenes/FileScene.cpp
            }
            else mesh = MakeLatheMesh(o.profile, o.segments);
            AddMesh(mesh, mat->second, out);
            break;
        }
        default:
            missing.push_back(what + ": only quads, meshes and lathes (no spheres, sea, hair or splats)");
        }
    }
    if (out.triMaterial.empty() && missing.empty()) missing.push_back("no geometry");

    // ---- lights chosen by emitted power ----
    float total = 0;
    for (float p : power) total += p;
    out.lightCdf.assign(1, 0.0f);
    for (size_t i = 0; i < power.size(); ++i) {
        const float p = total > 0 ? power[i] / total : 1.0f / float(power.size());
        out.lightCdf.push_back(i + 1 == power.size() ? 1.0f : out.lightCdf.back() + p);
    }

    // ---- options the kernel renderer does not model ----
    if (opt.aov || opt.denoise) missing.push_back("--aov / --denoise");
    if (opt.stokes) missing.push_back("--stokes");
    if (opt.spectralStep > 0) missing.push_back("--spectral-image");
    if (opt.toon) missing.push_back("--toon");
    if (!opt.mask.empty() || opt.depthMap) missing.push_back("--mask / --depth-map");

    out.header.width = uint32_t(width);
    out.header.height = uint32_t(height);
    out.header.spp = uint32_t(std::max(opt.spp, 1));
    out.header.maxDepth = uint32_t(std::max(opt.maxDepth, 0));
    out.header.seed = uint64_t(opt.seed);

    if (!missing.empty()) {
        why.clear();
        for (const std::string& m : missing) why += "\n    - " + m;
        return false;
    }
    return true;
}

int RenderOnKernel(const CliOptions& opt, const Camera& camera, Film& film) {
    const char* tag = opt.device == "gpu" ? "[GPU]" : "[Kernel]";
#if defined(RAYT_SPECTRAL_MODE)
    std::cerr << tag << " --device " << opt.device << " renders in RGB: use the RGB build (out\\cli or out\\cli-dev)\n";
    return 1;
#endif
    if (!opt.sceneDescription) {
        std::cerr << tag << " --device " << opt.device << " renders scene files only (GoLD_rayt scene.json --device " << opt.device << ")\n";
        return 1;
    }
    const int width = film.width(), height = film.height();
    KernelSceneHost host;
    std::string why;
    if (!PackKernelScene(*opt.sceneDescription, camera, opt, width, height, host, why)) {
        std::cerr << tag << " the kernel renderer (docs/GPU.md, stage 1) does not model:" << why << "\n" << tag << " render it with --device cpu\n";
        return 1;
    }
    const kernel::KScene s = host.view();
    std::cout << tag << " kernel scene: " << s.numTriangles << " triangles, " << s.numMaterials << " materials, " << s.numLights << " quad lights"
        << (s.envRadiance[0] > 0 ? ", uniform environment" : "") << "; " << width << " x " << height << ", " << s.spp << " spp, depth " << s.maxDepth << ", float\n";

    std::vector<float> mean;
    const auto t0 = std::chrono::steady_clock::now();
    if (opt.device == "cpu-kernel") {
        std::cout << tag << " on the CPU (" << omp_get_max_threads() << " threads, Embree in float)" << std::endl;
        if (!RenderKernelCpu(s, mean, why, true)) { std::cerr << tag << " " << why << std::endl; return 1; }
    }
    else {
        std::string device;
        if (!RenderKernelGpu(s, mean, why, true, nullptr, &device)) { std::cerr << std::endl << tag << " " << why << std::endl; return 1; }
        std::cout << tag << " on " << device << std::endl;
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << tag << " rendered in " << seconds << " s (" << double(width) * height * s.spp / seconds / 1e6 << " M samples/s)\n";

    // Into the film: row 0 of the kernel is the bottom (t = 0), the film's is the top, as in PathIntegrator::render.
    size_t nonFinite = 0;
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x) {
            const float* p = &mean[(size_t(y) * size_t(width) + size_t(x)) * 3];
            Vector3 c(p[0], p[1], p[2]);
            if (!std::isfinite(c.x) || !std::isfinite(c.y) || !std::isfinite(c.z)) { c = Vector3(0.0); ++nonFinite; }
            film.setPixel(x, height - 1 - y, c);
        }
    if (nonFinite) std::cout << tag << " " << nonFinite << " pixels were not finite and were set to 0\n";
    return 0;
}

} // namespace rayt::cli
