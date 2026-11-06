// Triangle mesh checks.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Cli/CheckSupport.hpp"
#include "Cli/EngineHeaders.hpp"

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
