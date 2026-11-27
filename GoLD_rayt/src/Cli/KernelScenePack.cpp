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
#include "Cli/Scene.hpp"
#include "Cli/SceneFile.hpp"
#include "Lights/EnvLight.hpp"
#include "Materials/Dielectric.hpp"
#include "Materials/RoughConductor.hpp"

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

    /// Geometry ids in the Embree scene: the triangles, then the spheres (user geometry).
    constexpr unsigned kEmbreeTriangles = 0, kEmbreeSpheres = 1;

    struct EmbreeTracer {
        RTCScene scene;
        uint32_t numTriangles;   ///< a sphere i is reported as primitive numTriangles + i, as the kernel numbers them
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
            if (rh.hit.geomID == kEmbreeSpheres) return { rh.ray.tfar, 0.0f, 0.0f, numTriangles + rh.hit.primID };
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

    // ---- spheres as Embree user geometry: the box and the test are the kernel's (Kernel/Lights.hpp), as OptiX's programs ----
    void SphereBounds(const RTCBoundsFunctionArguments* args) {
        const kernel::KSphere& sp = static_cast<const kernel::KSphere*>(args->geometryUserPtr)[args->primID];
        float lo[3], hi[3];
        kernel::sphereBounds(sp, lo, hi);
        RTCBounds& b = *args->bounds_o;
        b.lower_x = lo[0]; b.lower_y = lo[1]; b.lower_z = lo[2];
        b.upper_x = hi[0]; b.upper_y = hi[1]; b.upper_z = hi[2];
    }

    void SphereIntersect(const RTCIntersectFunctionNArguments* args) {
        if (args->N != 1 || !args->valid[0]) return;   // rtcIntersect1 only
        const kernel::KSphere& sp = static_cast<const kernel::KSphere*>(args->geometryUserPtr)[args->primID];
        RTCRayHit& rh = *reinterpret_cast<RTCRayHit*>(args->rayhit);
        float t;
        if (!kernel::intersectSphere(kernel::V3(rh.ray.org_x, rh.ray.org_y, rh.ray.org_z), kernel::V3(rh.ray.dir_x, rh.ray.dir_y, rh.ray.dir_z),
                kernel::V3(sp.center.x, sp.center.y, sp.center.z), sp.radius, rh.ray.tnear, rh.ray.tfar, t)) return;
        rh.ray.tfar = t;
        rh.hit.u = 0.0f; rh.hit.v = 0.0f;
        rh.hit.Ng_x = 0.0f; rh.hit.Ng_y = 0.0f; rh.hit.Ng_z = 0.0f;   // not read: the kernel computes the normal
        rh.hit.primID = args->primID;
        rh.hit.geomID = args->geomID;
        rh.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    }

    void SphereOccluded(const RTCOccludedFunctionNArguments* args) {
        if (args->N != 1 || !args->valid[0]) return;   // rtcOccluded1 only
        const kernel::KSphere& sp = static_cast<const kernel::KSphere*>(args->geometryUserPtr)[args->primID];
        RTCRay& r = *reinterpret_cast<RTCRay*>(args->ray);
        float t;
        if (kernel::intersectSphere(kernel::V3(r.org_x, r.org_y, r.org_z), kernel::V3(r.dir_x, r.dir_y, r.dir_z),
                kernel::V3(sp.center.x, sp.center.y, sp.center.z), sp.radius, r.tnear, r.tfar, t))
            r.tfar = -std::numeric_limits<float>::infinity();
    }

    /// The triangles (robust, high-quality build) and the spheres (user geometry) of a kernel scene in an Embree scene, released
    /// with this object.
    struct EmbreeScene {
        RTCScene scene = nullptr;
        bool build(const kernel::KScene& s, std::string& why) {
            RTCDevice device = detail::embreeDevice();
            if (!device) { why = "Embree: rtcNewDevice failed"; return false; }
            scene = rtcNewScene(device);
            rtcSetSceneFlags(scene, RTC_SCENE_FLAG_ROBUST);
            rtcSetSceneBuildQuality(scene, RTC_BUILD_QUALITY_HIGH);
            if (s.numTriangles > 0) {
                RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
                // Copied into buffers Embree allocates (it reads vertices 16 bytes at a time and pads its own buffers for that).
                float* vb = static_cast<float*>(rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float), s.numVertices));
                memcpy(vb, s.positions, size_t(s.numVertices) * sizeof(kernel::KFloat3));
                unsigned* ib = static_cast<unsigned*>(rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned), s.numTriangles));
                memcpy(ib, s.indices, size_t(s.numTriangles) * 3 * sizeof(unsigned));
                rtcCommitGeometry(geom);
                rtcAttachGeometryByID(scene, geom, kEmbreeTriangles);
                rtcReleaseGeometry(geom);
            }
            if (s.numSpheres > 0) {
                RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
                rtcSetGeometryUserPrimitiveCount(geom, s.numSpheres);
                rtcSetGeometryUserData(geom, const_cast<kernel::KSphere*>(s.spheres));
                rtcSetGeometryBoundsFunction(geom, &SphereBounds, nullptr);
                rtcSetGeometryIntersectFunction(geom, &SphereIntersect);
                rtcSetGeometryOccludedFunction(geom, &SphereOccluded);
                rtcCommitGeometry(geom);
                rtcAttachGeometryByID(scene, geom, kEmbreeSpheres);
                rtcReleaseGeometry(geom);
            }
            rtcCommitScene(scene);
            return true;
        }
        ~EmbreeScene() { if (scene) rtcReleaseScene(scene); }
    };

    /// An Embree tracer that also writes down every query: closest hits as they are, shadow queries as a hit (prim 0) or none.
    struct RecordingTracer {
        EmbreeTracer base;
        std::vector<kernel::KRay>* rays;
        std::vector<kernel::KHit>* hits;
        std::vector<uint8_t>* shadow;
        kernel::KHit closest(const kernel::V3& o, const kernel::V3& d, float tmin, float tmax) const {
            const kernel::KHit h = base.closest(o, d, tmin, tmax);
            rays->push_back({ { o.x, o.y, o.z }, { d.x, d.y, d.z }, tmin, tmax });
            hits->push_back(h);
            shadow->push_back(0);
            return h;
        }
        bool occluded(const kernel::V3& o, const kernel::V3& d, float tmin, float tmax) const {
            const bool blocked = base.occluded(o, d, tmin, tmax);
            rays->push_back({ { o.x, o.y, o.z }, { d.x, d.y, d.z }, tmin, tmax });
            hits->push_back({ 0.0f, 0.0f, 0.0f, blocked ? 0u : kernel::kNoHit });
            shadow->push_back(1);
            return blocked;
        }
    };

} // namespace

bool RecordPathRays(const kernel::KScene& s, std::vector<kernel::KRay>& rays, std::vector<kernel::KHit>& hits, std::vector<uint8_t>& shadow, std::string& why) {
    EmbreeScene es;
    if (!es.build(s, why)) return false;
    rays.clear(); hits.clear(); shadow.clear();
    const RecordingTracer tracer{ EmbreeTracer{ es.scene, s.numTriangles }, &rays, &hits, &shadow };
    for (uint32_t y = 0; y < s.height; ++y)
        for (uint32_t x = 0; x < s.width; ++x)
            for (uint32_t k = 0; k < s.spp; ++k) kernel::tracePath(s, tracer, x, y, k);
    return true;
}

bool RenderKernelCpu(const kernel::KScene& s, std::vector<float>& mean, std::string& why, bool progress) {
    EmbreeScene es;
    if (!es.build(s, why)) return false;
    const EmbreeTracer tracer{ es.scene, s.numTriangles };
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

bool PackEnvironment(const EnvLight* envLight, KernelSceneHost& out) {
    const EnvMap* map = envLight ? envLight->envMap() : nullptr;
    const Distribution2D* dist = map ? map->distribution() : nullptr;
    if (!map || !map->image().isValid() || !dist) return false;
    const Image& img = map->image();
    const int w = img.width(), h = img.height();
    out.envPixels.resize(size_t(w) * size_t(h) * 3);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            const Vector3& c = img.at(x, y);   // floats from the file, held as doubles: converting back is exact
            float* p = &out.envPixels[(size_t(y) * size_t(w) + size_t(x)) * 3];
            p[0] = float(c.x); p[1] = float(c.y); p[2] = float(c.z);
        }
    out.envMarginalFunc = dist->pMarginal->func;
    out.envMarginalCdf = dist->pMarginal->cdf;
    out.header.envMarginalInt = dist->pMarginal->funcInt;
    out.envCondFunc.clear(); out.envCondCdf.clear(); out.envCondInt.clear();
    for (int y = 0; y < h; ++y) {
        const Distribution1D& row = *dist->pConditionalV[size_t(y)];
        out.envCondFunc.insert(out.envCondFunc.end(), row.func.begin(), row.func.end());
        out.envCondCdf.insert(out.envCondCdf.end(), row.cdf.begin(), row.cdf.end());
        out.envCondInt.push_back(row.funcInt);
    }
    out.header.envWidth = uint32_t(w);
    out.header.envHeight = uint32_t(h);
    out.header.envScale = float(envLight->scale());
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
    s.spheres = spheres.data();
    s.lights = lights.data();
    s.lightCdf = lightCdf.data();
    s.envPixels = envPixels.data();
    s.envMarginalFunc = envMarginalFunc.data();
    s.envMarginalCdf = envMarginalCdf.data();
    s.envCondFunc = envCondFunc.data();
    s.envCondCdf = envCondCdf.data();
    s.envCondInt = envCondInt.data();
    s.numVertices = uint32_t(positions.size());
    s.numTriangles = uint32_t(triMaterial.size());
    s.numSpheres = uint32_t(spheres.size());
    s.numMaterials = uint32_t(materials.size());
    s.numLights = uint32_t(lights.size());
    return s;
}

bool PackKernelScene(const SceneDescription& d, const std::map<std::string, const Material*>& engineMaterials, const EnvLight* envLight, const Camera& camera,
    const CliOptions& opt,
    int width, int height, KernelSceneHost& out, std::string& why) {
    std::vector<std::string> missing;
    out = KernelSceneHost{};

    // ---- camera ----
    Point3 origin, lowerLeft;
    Vector3 horizontal, vertical;
    if (!camera.pinholeFrame(origin, lowerLeft, horizontal, vertical))
        missing.push_back("the camera is not a plain pinhole (a lens, an aperture, distortion, vignetting, motion or an orthographic view)");
    out.header.camera = { F3(origin), F3(lowerLeft), F3(horizontal), F3(vertical) };

    // ---- environment (the engine's EnvLight: its map, its sampling distribution, its scale) and point lights ----
    if (d.environment.type != EnvironmentSpec::Type::None && !PackEnvironment(envLight, out)) missing.push_back("the engine's environment is not available");
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
        else if (m.type == MaterialSpec::Type::Conductor) {
            // The engine derives the RGB index from the metal's table (and divides it by a surrounding medium's index): read it from there.
            const auto it = engineMaterials.find(name);
            const RoughConductor* rc = it == engineMaterials.end() ? nullptr : dynamic_cast<const RoughConductor*>(it->second);
            if (!rc) { missing.push_back("material \"" + name + "\": the engine's conductor is not available"); continue; }
#if defined(RAYT_SPECTRAL_MODE)
            missing.push_back("material \"" + name + "\": a conductor in the spectral build");
            continue;
#else
            kernel::ConductorParams<Real> c{};
            std::string limit;
            if (!rc->kernelParams(c, limit)) { missing.push_back("material \"" + name + "\": a conductor with " + limit); continue; }
            km.type = kernel::kMatConductor;
            kernel::ConductorParams<float>& k = km.conductor;
            for (int i = 0; i < 3; ++i) { k.eta[i] = float(c.eta[i]); k.k[i] = float(c.k[i]); k.fms[i] = float(c.fms[i]); }
            k.alpha = float(c.alpha);
            k.multipleScattering = c.multipleScattering;
            k.eAvg = float(c.eAvg);
            for (int i = 0; i < kernel::kAlbedoMu; ++i) k.eRow[i] = float(c.eRow[i]);
#endif
        }
        else if (m.type == MaterialSpec::Type::Dielectric) {
            const auto it = engineMaterials.find(name);
            const Dielectric* de = it == engineMaterials.end() ? nullptr : dynamic_cast<const Dielectric*>(it->second);
            if (!de) { missing.push_back("material \"" + name + "\": the engine's dielectric is not available"); continue; }
#if defined(RAYT_SPECTRAL_MODE)
            missing.push_back("material \"" + name + "\": a dielectric in the spectral build");
            continue;
#else
            kernel::DielectricParams<Real> c{};
            std::string limit;
            if (!de->kernelParams(c, limit)) { missing.push_back("material \"" + name + "\": a dielectric with " + limit); continue; }
            km.type = kernel::kMatDielectric;
            km.dielectric.ior = float(c.ior);
            km.dielectric.alpha = float(c.alpha);
            km.dielectric.smooth = c.smooth;
#endif
        }
        else { missing.push_back("material \"" + name + "\": only diffuse, conductor, dielectric and emission materials"); continue; }
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
                kernel::KLight q{};
                q.type = kernel::kLightQuad;
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
        case ObjectSpec::Shape::Sphere: {
            // As the engine's Sphere (Scenes/FileScene.cpp): an analytic sphere, a light when it emits (unless "light": false).
            kernel::KSphere sp{};
            sp.center = F3(o.center);
            sp.radius = float(o.radius);
            sp.material = mat->second;
            sp.light = -1;
            if (emission && o.light != 0) {
                const kernel::KMaterial& km = out.materials[mat->second];
                kernel::KLight q{};
                q.type = kernel::kLightSphere;
                q.center = sp.center;
                q.radius = sp.radius;
                q.area = float(4.0 * constants::PI * double(o.radius) * double(o.radius));
                for (int c = 0; c < 3; ++c) q.radiance[c] = km.emission[c];
                q.twoSided = km.twoSided;
                sp.light = int32_t(out.lights.size());
                out.lights.push_back(q);
                power.push_back(q.area * (0.2126f * q.radiance[0] + 0.7152f * q.radiance[1] + 0.0722f * q.radiance[2]) * (q.twoSided ? 2.0f : 1.0f));
            }
            out.spheres.push_back(sp);
            break;
        }
        case ObjectSpec::Shape::Mesh:
        case ObjectSpec::Shape::Lathe: {
            if (emission) { missing.push_back(what + ": an emitting mesh (only emitting quads and spheres)"); continue; }
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
            missing.push_back(what + ": only quads, spheres, meshes and lathes (no sea, hair or splats)");
        }
    }
    if (out.triMaterial.empty() && out.spheres.empty() && missing.empty()) missing.push_back("no geometry");

    // ---- lights chosen by emitted power ----
    float total = 0;
    for (float p : power) total += p;
    out.lightCdf.assign(1, 0.0f);
    for (size_t i = 0; i < power.size(); ++i) {
        const float p = total > 0 ? power[i] / total : 1.0f / float(power.size());
        out.lightCdf.push_back(i + 1 == power.size() ? 1.0f : out.lightCdf.back() + p);
    }
    // The environment's share of the light samples, by emitted power as pbrt's PowerLightSampler weighs an infinite light: the
    // power that crosses the scene's bounding sphere, pi^2 r^2 L_avg, against pi A L of each light (once more for two-sided ones).
    // Each side keeps at least a tenth when both emit (any choice is unbiased; the engine's light sampler uses another rule).
    out.header.envSelect = 0.0f;
    if (out.header.envWidth > 0) {
        Vector3 lo(1e30), hi(-1e30);
        for (const kernel::KFloat3& p : out.positions) {
            lo = glm::min(lo, Vector3(p.x, p.y, p.z));
            hi = glm::max(hi, Vector3(p.x, p.y, p.z));
        }
        for (const kernel::KSphere& sp : out.spheres) {
            const Vector3 c(sp.center.x, sp.center.y, sp.center.z), r(std::abs(double(sp.radius)));
            lo = glm::min(lo, c - r);
            hi = glm::max(hi, c + r);
        }
        const double r2 = (out.positions.empty() && out.spheres.empty()) ? 1.0 : 0.25 * glm::dot(hi - lo, hi - lo);
        const double envPower = constants::PI * constants::PI * r2 * double(out.header.envScale) * double(envLight->envMap()->averageLuminance());
        const double lightPower = constants::PI * double(total);
        if (power.empty() || lightPower <= 0) out.header.envSelect = envPower > 0 ? 1.0f : 0.0f;
        else if (envPower > 0) out.header.envSelect = float(std::clamp(envPower / (envPower + lightPower), 0.1, 0.9));
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

bool PackSceneFile(const std::string& path, const CliOptions& base, int width, int height, const std::function<void(SceneDescription&)>& edit,
    KernelSceneHost& out, std::string& why) {
    std::vector<std::string> args;
    std::shared_ptr<const SceneDescription> loaded;
    if (!LoadSceneFile(path, args, loaded, why) || !loaded) { why = path + ": " + why; return false; }
    auto d = std::make_shared<SceneDescription>(*loaded);
    if (edit) edit(*d);
    CliOptions opt = base;
    opt.scene = "file";
    opt.sceneFile = path;
    opt.sceneDescription = d;
    const SceneDef* def = FindScene("file");
    if (!def) { why = "the scene \"file\" is not registered"; return false; }
    SceneSetup setup;
    Environment env;
    {
        // The engine's scene build reports what it loads; a check has no use for that.
        struct Quiet {
            std::ostringstream sink;
            std::streambuf* old = std::cout.rdbuf(sink.rdbuf());
            ~Quiet() { std::cout.rdbuf(old); }
        } quiet;
        env = MakeEnvironment(opt, def->envKind(opt));
        const StdMaterials mats = MakeStandardMaterials(opt);
        SceneContext ctx{ opt, mats, env };
        if (def->build(ctx, setup)) { why = path + ": the engine could not build the scene"; return false; }
    }
    const CameraSpec& c = d->camera;
    const Camera camera(c.from, c.at, c.up, c.fovDeg, Real(width) / Real(height), Real(0), glm::length(c.at - c.from));
    return PackKernelScene(*d, setup.namedMaterials, env.light.get(), camera, opt, width, height, out, why);
}

int RenderOnKernel(const CliOptions& opt, const std::map<std::string, const Material*>& engineMaterials, const EnvLight* envLight, const Camera& camera,
    Film& film) {
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
    if (!PackKernelScene(*opt.sceneDescription, engineMaterials, envLight, camera, opt, width, height, host, why)) {
        std::cerr << tag << " the kernel renderer (docs/GPU.md) does not model:" << why << "\n" << tag << " render it with --device cpu\n";
        return 1;
    }
    const kernel::KScene s = host.view();
    std::cout << tag << " kernel scene: " << s.numTriangles << " triangles, " << s.numSpheres << " spheres, " << s.numMaterials << " materials, " << s.numLights << " lights"
        << (kernel::hasEnvironment(s) ? ", an environment of " + std::to_string(s.envWidth) + " x " + std::to_string(s.envHeight) + " (" +
            std::to_string(int(std::lround(100.0 * s.envSelect))) + " % of the light samples)" : std::string()) << "; " << width << " x " << height << ", " << s.spp << " spp, depth " << s.maxDepth << ", float\n";

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
