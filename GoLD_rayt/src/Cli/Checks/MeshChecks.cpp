// Triangle mesh checks.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Geometry/TriangleMesh.hpp"
#include "Materials/Lambertian.hpp"
#include "Accelerators/EmbreeMesh.hpp"
#include "IO/ObjLoader.hpp"
#include "IO/PlyLoader.hpp"
#include "Lights/AreaLight.hpp"
#include "Materials/DiffuseLight.hpp"
#include "Scene/Primitive.hpp"

namespace rayt::cli {

bool RunMeshChecks(const CliOptions& opt) {
    bool allPassed = true;
    auto check = [&](const std::string& what, bool pass, const std::string& detail) {
        std::cout << "[Mesh]   " << what << ": " << detail << "   " << (pass ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && pass;
    };
    auto white = std::make_shared<Lambertian>(RgbColor(1.0));

    // a. geometry of a procedural icosphere (an inscribed polyhedron: a bit smaller than the sphere)
    auto sphere = std::make_shared<TriangleMesh>(TriangleMesh::Icosphere(3, 1.0, Point3(0.0), true));
    {
        const Real a = sphere->surfaceArea() / (4 * constants::PI);
        const Real v = sphere->signedVolume() / (4.0 / 3.0 * constants::PI);
        check("icosphere (" + std::to_string(sphere->numTriangles()) + " triangles) area / sphere, volume / sphere",
            a > 0.97 && a <= 1.0 && v > 0.95 && v <= 1.0, std::to_string(a) + ", " + std::to_string(v) + "  (volume > 0 = outward winding)");
    }

    // b. one ray: hit distance and the interpolated (smooth) normal
    const std::shared_ptr<const TriangleMesh> sphereC = sphere;
    auto bvh = BuildMeshBVH(sphereC, white);
    {
        SurfaceInteraction rec;
        const bool hit = bvh->hit(Ray(Point3(0, 0, 5), Vector3(0, 0, -1)), rec);
        const bool ok = hit && rec.t > 3.98 && rec.t <= 4.0 + 1e-9 && glm::dot(rec.n, Vector3(0, 0, 1)) > 0.999 && glm::dot(rec.gn, Vector3(0, 0, 1)) > 0.99;
        check("ray at the pole: t, shading normal . z, geometric normal . z", ok,
            hit ? std::to_string(rec.t) + ", " + std::to_string(glm::dot(rec.n, Vector3(0, 0, 1))) + ", " + std::to_string(glm::dot(rec.gn, Vector3(0, 0, 1))) : std::string("no hit"));
    }

    // c. OBJ round trip
    {
        const std::string path = (std::filesystem::temp_directory_path() / "rayt_test_mesh.obj").string();
        {
            std::ofstream f(path);
            f.precision(12);
            for (const auto& p : sphere->p) f << "v " << p.x << " " << p.y << " " << p.z << "\n";
            for (const auto& n : sphere->n) f << "vn " << n.x << " " << n.y << " " << n.z << "\n";
            for (size_t i = 0; i < sphere->numTriangles(); ++i)
                f << "f " << sphere->indices[3 * i] + 1 << "//" << sphere->indices[3 * i] + 1 << " "
                  << sphere->indices[3 * i + 1] + 1 << "//" << sphere->indices[3 * i + 1] + 1 << " "
                  << sphere->indices[3 * i + 2] + 1 << "//" << sphere->indices[3 * i + 2] + 1 << "\n";
        }
        TriangleMesh back;
        std::string err;
        const bool loaded = io::loadObj(path, back, &err);
        const bool ok = loaded && back.numTriangles() == sphere->numTriangles() && back.hasNormals()
            && std::abs(back.signedVolume() - sphere->signedVolume()) < 1e-9;
        check("OBJ write + load (triangles, normals, volume)", ok,
            loaded ? std::to_string(back.numTriangles()) + " triangles" : err);
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }

    // c2. PLY round trip: ascii with per-vertex uv, and binary little endian with MeshLab's per-face "texcoord" list.
    //     Both must give the sphere back (triangles, normals, volume) with the uv's v flipped to the renderer's convention.
    {
        const std::string base = (std::filesystem::temp_directory_path() / "rayt_test_mesh").string();
        std::vector<Point2> uvs;   // a uv per vertex: longitude, latitude
        for (const auto& p : sphere->p) uvs.emplace_back(0.5 + std::atan2(p.z, p.x) / (2 * constants::PI), 0.5 + std::asin(std::clamp(p.y, -1.0, 1.0)) / constants::PI);
        const size_t nt = sphere->numTriangles();
        {   // ascii
            std::ofstream f(base + "_a.ply");
            f.precision(12);
            f << "ply\nformat ascii 1.0\ncomment TextureFile sphere.png\nelement vertex " << sphere->p.size()
              << "\nproperty float x\nproperty float y\nproperty float z\nproperty float nx\nproperty float ny\nproperty float nz\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nproperty float s\nproperty float t\n"
              << "element face " << nt << "\nproperty list uchar int vertex_indices\nproperty int flags\nend_header\n";
            for (size_t i = 0; i < sphere->p.size(); ++i)
                f << sphere->p[i].x << " " << sphere->p[i].y << " " << sphere->p[i].z << " " << sphere->n[i].x << " " << sphere->n[i].y << " " << sphere->n[i].z
                  << " 200 100 50 " << uvs[i].x << " " << uvs[i].y << "\n";
            for (size_t i = 0; i < nt; ++i) f << "3 " << sphere->indices[3 * i] << " " << sphere->indices[3 * i + 1] << " " << sphere->indices[3 * i + 2] << " 0\n";
        }
        {   // binary little endian, uv per face corner (the vertices are not split here: each vertex has one uv)
            std::ofstream f(base + "_b.ply", std::ios::binary);
            f << "ply\nformat binary_little_endian 1.0\nelement vertex " << sphere->p.size()
              << "\nproperty double x\nproperty double y\nproperty double z\nproperty float nx\nproperty float ny\nproperty float nz\n"
              << "element face " << nt << "\nproperty list uchar int vertex_indices\nproperty list uchar float texcoord\nend_header\n";
            auto put = [&](const void* p, size_t n) { f.write(static_cast<const char*>(p), static_cast<std::streamsize>(n)); };
            for (size_t i = 0; i < sphere->p.size(); ++i) {
                const double xyz[3] = { sphere->p[i].x, sphere->p[i].y, sphere->p[i].z };
                const float n3[3] = { float(sphere->n[i].x), float(sphere->n[i].y), float(sphere->n[i].z) };
                put(xyz, 24); put(n3, 12);
            }
            for (size_t i = 0; i < nt; ++i) {
                const unsigned char three = 3, six = 6;
                const std::int32_t idx[3] = { sphere->indices[3 * i], sphere->indices[3 * i + 1], sphere->indices[3 * i + 2] };
                float uv6[6];
                for (int k = 0; k < 3; ++k) { uv6[2 * k] = float(uvs[static_cast<size_t>(idx[k])].x); uv6[2 * k + 1] = float(uvs[static_cast<size_t>(idx[k])].y); }
                put(&three, 1); put(idx, 12); put(&six, 1); put(uv6, 24);
            }
        }
        for (const char* suffix : { "_a.ply", "_b.ply" }) {
            TriangleMesh back;
            std::string err, texture;
            const bool loaded = io::loadPly(base + suffix, back, &err, &texture);
            bool uvOk = loaded && back.uv.size() == back.p.size();
            double maxN = 0;
            if (uvOk) {   // every loaded vertex must carry the uv of the original vertex at its position, with v flipped
                std::map<std::tuple<float, float, float>, size_t> byPos;
                for (size_t i = 0; i < sphere->p.size(); ++i) byPos[std::make_tuple(float(sphere->p[i].x), float(sphere->p[i].y), float(sphere->p[i].z))] = i;
                for (size_t i = 0; i < back.p.size() && uvOk; ++i) {
                    const auto it = byPos.find(std::make_tuple(float(back.p[i].x), float(back.p[i].y), float(back.p[i].z)));
                    if (it == byPos.end()) { uvOk = false; break; }
                    const Point2 want(uvs[it->second].x, 1.0 - uvs[it->second].y);
                    if (std::abs(back.uv[i].x - want.x) > 1e-5 || std::abs(back.uv[i].y - want.y) > 1e-5) uvOk = false;
                    maxN = std::max(maxN, 1.0 - glm::dot(back.n[i], sphere->n[it->second]));
                }
            }
            const bool ok = loaded && back.numTriangles() == nt && back.hasNormals() && maxN < 1e-6
                && std::abs(back.signedVolume() - sphere->signedVolume()) < 1e-6 && uvOk && (suffix[1] != 'a' || texture == "sphere.png");
            check(std::string("PLY write + load, ") + (suffix[1] == 'a' ? "ascii (vertex uv, skipped colors, TextureFile comment)" : "binary (per-face texcoord list)")
                + ": triangles, normals, volume, uv flipped", ok, loaded ? std::to_string(back.numTriangles()) + " triangles, " + std::to_string(back.p.size()) + " vertices" : err);
            std::error_code ec;
            std::filesystem::remove(base + suffix, ec);
        }
    }

    // d. BVH vs brute force on the torus file (2304 triangles): same hit / miss / distance for random rays
    {
        TriangleMesh torus;
        std::string err;
        if (!io::loadObj("assets/Meshes/torus.obj", torus, &err)) {
            check("torus.obj", false, err);
        }
        else {
            auto mesh = std::make_shared<const TriangleMesh>(std::move(torus));
            std::vector<std::shared_ptr<Triangle>> tris;
            auto torusBvh = BuildMeshBVH(mesh, white, &tris);
            Aggregate brute;
            for (const auto& t : tris) brute.add(t);
            int mismatches = 0, hits = 0;
            for (int i = 0; i < 20000; ++i) {
                const Vector3 target(sampling::Random(-0.7, 0.7), sampling::Random(-0.25, 0.25), sampling::Random(-0.7, 0.7));
                const Vector3 origin = glm::normalize(Vector3(sampling::Random(-1, 1), sampling::Random(-1, 1), sampling::Random(-1, 1))) * 3.0;
                const Ray ray(origin, glm::normalize(target - origin));
                SurfaceInteraction a, b;
                const bool ha = torusBvh->hit(ray, a), hb = brute.hit(ray, b);
                hits += hb ? 1 : 0;
                if (ha != hb || (ha && std::abs(a.t - b.t) > 1e-9)) ++mismatches;
                // any-hit with a random shortened range must agree with the closest-hit distance
                Ray shortRay = ray;
                shortRay.tMax = sampling::Random(0.5, 5.0);
                if (torusBvh->occluded(shortRay) != brute.occluded(shortRay) || brute.occluded(shortRay) != (hb && b.t <= shortRay.tMax)) ++mismatches;
            }
            check("torus.obj " + std::to_string(mesh->numTriangles()) + " triangles, volume " + std::to_string(mesh->signedVolume())
                + ": BVH vs brute force on 20000 rays (" + std::to_string(hits) + " hits), mismatches",
                mismatches == 0 && hits > 1000 && mesh->signedVolume() > 0, std::to_string(mismatches));

            // Embree (float traversal, Real hit record) vs brute force on the same kind of rays.
            {
                EmbreeMesh embree(mesh, white);
                int diff = 0, normalDiff = 0, n = 0;
                for (int i = 0; i < 20000; ++i) {
                    const Vector3 target(sampling::Random(-0.7, 0.7), sampling::Random(-0.25, 0.25), sampling::Random(-0.7, 0.7));
                    const Vector3 origin = glm::normalize(Vector3(sampling::Random(-1, 1), sampling::Random(-1, 1), sampling::Random(-1, 1))) * 3.0;
                    const Ray ray(origin, glm::normalize(target - origin));
                    SurfaceInteraction a, b;
                    const bool ha = embree.hit(ray, a), hb = brute.hit(ray, b);
                    ++n;
                    if (ha != hb || (ha && std::abs(a.t - b.t) > 1e-6)) { ++diff; continue; }
                    if (ha && (glm::dot(a.n, b.n) < 0.9999 || a.matPtr != b.matPtr)) ++normalDiff;
                    Ray shortRay = ray;
                    shortRay.tMax = sampling::Random(0.5, 5.0);
                    if (embree.occluded(shortRay) != brute.occluded(shortRay)) ++diff;
                }
                check("Embree vs brute force on " + std::to_string(n) + " rays (hit/miss, distance, any-hit; normals + material)",
                    diff <= 5 && normalDiff == 0, std::to_string(diff) + " ray(s) differ at float/double edge cases, " + std::to_string(normalDiff) + " normal/material differences");
            }
        }
    }

    // e. white (albedo 1) diffuse mesh sphere in a uniform white environment must stay exactly white
    {
        auto envMap = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        for (const bool smooth : { true, false }) {
            auto mesh = std::make_shared<const TriangleMesh>(TriangleMesh::Icosphere(2, 1.0, Point3(0.0), smooth));
            auto world = std::make_shared<Aggregate>();
            world->add(BuildMeshBVH(mesh, white));
            auto scene = std::make_shared<Scene>(world);
            scene->setEnvLight(std::make_shared<EnvLight>(envMap, 1.0));
            scene->finalize();
            double sum = 0;
            const int n = opt.spp * 256;
            Vector3 m, se;
            int bad = 0;
            EstimateRadiance(*scene, Ray(Point3(0.3, 0.2, 5), Vector3(0, 0, -1)), n, 16, m, se, bad);
            const bool ok = std::abs(m.x - 1.0) <= 4 * se.x + 0.01 && bad == 0;
            check(std::string("furnace, white diffuse mesh sphere, ") + (smooth ? "smooth normals" : "flat normals"), ok,
                "mean " + std::to_string(m.x) + " +- " + std::to_string(se.x));
            (void)sum;
        }
    }

    // f. two emissive triangles = the Quad light of --scene lights (numerical reference)
    {
        const Point3 corner(-0.5, 2, -0.5);
        const Vector3 e1(1, 0, 0), e2(0, 0, 1);
        auto mesh = std::make_shared<TriangleMesh>();
        mesh->p = { corner, corner + e1, corner + e1 + e2, corner + e2 };
        mesh->indices = { 0, 1, 2, 0, 2, 3 };
        const std::shared_ptr<const TriangleMesh> meshC = mesh;
        const Real albedo = 0.5;
        auto scene = MakeGroundScene(albedo, [&](Aggregate& w, Scene& s) {
            auto mat = std::make_shared<DiffuseLight>(EmissionSpectrum(RgbColor(1.0)), 10.0);
            for (size_t i = 0; i < meshC->numTriangles(); ++i) {
                auto tri = std::make_shared<Triangle>(meshC, i, mat);
                auto light = std::make_shared<AreaLight>(tri, mat);
                w.add(std::make_shared<Primitive>(tri, mat, light));
                s.addLight(light);
            }
        });
        Vector3 m, se;
        int bad = 0;
        EstimateRadiance(*scene, Ray(Point3(-3, 3, 0), glm::normalize(Vector3(1, -1, 0))), opt.spp * 1024, 8, m, se, bad);
        const Real expected = QuadGroundRadiance(corner, e1, e2, false, 10.0, albedo);
        check("emissive triangle pair (area lights, NEE+MIS) vs quad reference " + std::to_string(expected),
            std::abs(m.x - expected) <= 4 * se.x + 0.005 * expected && bad == 0,
            std::to_string(m.x) + " +- " + std::to_string(se.x));
    }

    std::cout << "[Mesh] checks: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed;
}

} // namespace rayt::cli
