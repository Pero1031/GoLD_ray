// BVH benchmark (--scene bench).
#include "pch.h"
#include "Cli/Tools.hpp"
#include "Accelerators/EmbreeMesh.hpp"
#include "Materials/Lambertian.hpp"

namespace rayt::cli {

int RunBvhBenchmark(const CliOptions& opt) {
    using clock = std::chrono::steady_clock;
    auto ms = [](clock::time_point a, clock::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    const int nRays = opt.spp * 2000;
    bool ok = true;
    std::cout << "[Bench] " << nRays << " random rays per mesh, single thread (use a Release build)\n";

    auto white = std::make_shared<Lambertian>(RgbColor(1.0));
    for (const int subdiv : { 1, 2, 3, 4, 5, 6, 7 }) {   // 80 ... 328k triangles (small meshes show the per-ray overhead)
        const auto mesh = std::make_shared<const TriangleMesh>(TriangleMesh::Icosphere(subdiv, 1.0, Point3(0.0), true));
        std::vector<std::shared_ptr<Hittable>> prims;
        for (size_t i = 0; i < mesh->numTriangles(); ++i) prims.push_back(std::make_shared<Triangle>(mesh, i, white));

        // Rays from a sphere of radius 3 toward random points of a ball of radius 1.3 (about half hit).
        std::vector<Ray> rays(static_cast<size_t>(nRays));
        for (Ray& r : rays) {
            const Vector3 from = glm::normalize(Vector3(sampling::Random(-1, 1), sampling::Random(-1, 1), sampling::Random(-1, 1))) * 3.0;
            Vector3 to(0.0);
            do { to = Vector3(sampling::Random(-1.3, 1.3), sampling::Random(-1.3, 1.3), sampling::Random(-1.3, 1.3)); } while (glm::dot(to, to) > 1.69);
            r = Ray(from, glm::normalize(to - from));
        }

        auto v1 = prims;
        const auto t0 = clock::now();
        BVHNode legacy(v1, 0, v1.size());
        const auto t1 = clock::now();
        BVH flat(prims);
        const auto t2 = clock::now();

        SurfaceInteraction rec;
        int hitsLegacy = 0, hitsFlat = 0, mismatches = 0;
        std::vector<Real> tLegacy(rays.size()), tFlat(rays.size());
        const auto c0 = clock::now();
        for (size_t i = 0; i < rays.size(); ++i) { const bool h = legacy.hit(rays[i], rec); tLegacy[i] = h ? rec.t : -1; hitsLegacy += h; }
        const auto c1 = clock::now();
        for (size_t i = 0; i < rays.size(); ++i) { const bool h = flat.hit(rays[i], rec); tFlat[i] = h ? rec.t : -1; hitsFlat += h; }
        const auto c2 = clock::now();
        int occ = 0;
        for (size_t i = 0; i < rays.size(); ++i) occ += flat.occluded(rays[i]) ? 1 : 0;
        const auto c3 = clock::now();

        for (size_t i = 0; i < rays.size(); ++i) {
            if (std::abs(tLegacy[i] - tFlat[i]) > 1e-9) ++mismatches;
        }
        ok = ok && mismatches == 0 && occ == hitsFlat;

        // Embree (float traversal; the hit record is built in Real precision by the Triangle).
        const auto e0 = clock::now();
        EmbreeMesh embree(mesh, white);
        const auto e1 = clock::now();
        int hitsEmbree = 0, embreeMismatch = 0;
        std::vector<Real> tEmbree(rays.size());
        const auto e2 = clock::now();
        for (size_t i = 0; i < rays.size(); ++i) { const bool h = embree.hit(rays[i], rec); tEmbree[i] = h ? rec.t : -1; hitsEmbree += h; }
        const auto e3 = clock::now();
        int occEmbree = 0;
        for (size_t i = 0; i < rays.size(); ++i) occEmbree += embree.occluded(rays[i]) ? 1 : 0;
        const auto e4 = clock::now();
        for (size_t i = 0; i < rays.size(); ++i) {
            if (std::abs(tEmbree[i] - tFlat[i]) > 1e-6) ++embreeMismatch;
        }
        const bool embreeOk = embreeMismatch <= std::max(2, nRays / 10000);   // a ray grazing an edge may differ in float vs double
        ok = ok && embreeOk;
        std::cout << "[Bench]     Embree: build " << ms(e0, e1) << " ms (incl. Triangle objects, all threads);  closest hit "
            << double(nRays) / 1e6 / (ms(e2, e3) / 1000.0) << " Mray/s (x" << ms(c1, c2) / ms(e2, e3) << " vs flat BVH);  any-hit "
            << double(nRays) / 1e6 / (ms(e3, e4) / 1000.0) << " Mray/s;  " << hitsEmbree << " hits, " << embreeMismatch
            << " differ from the double BVH (|dt|>1e-6)   " << (embreeOk ? "OK" : "FAIL") << "\n";
        (void)occEmbree;

        // Embree without the robust (watertight) flag: how much of its time is that flag?
        {
            EmbreeMesh fast(mesh, white, /*robust=*/false);
            int h = 0;
            const auto f0 = clock::now();
            for (size_t i = 0; i < rays.size(); ++i) h += fast.hit(rays[i], rec) ? 1 : 0;
            const auto f1 = clock::now();
            std::cout << "[Bench]     Embree (no robust flag): closest hit " << double(nRays) / 1e6 / (ms(f0, f1) / 1000.0)
                << " Mray/s (x" << ms(c1, c2) / ms(f0, f1) << " vs flat BVH, " << h << " hits)\n";
        }

        // All threads at once, as in a render: do both backends scale?
        {
            const int count = static_cast<int>(rays.size());
            int hf = 0, he = 0;
            const auto m0 = clock::now();
#pragma omp parallel for reduction(+:hf) schedule(static)
            for (int i = 0; i < count; ++i) { SurfaceInteraction r2; hf += flat.hit(rays[static_cast<size_t>(i)], r2) ? 1 : 0; }
            const auto m1 = clock::now();
#pragma omp parallel for reduction(+:he) schedule(static)
            for (int i = 0; i < count; ++i) { SurfaceInteraction r2; he += embree.hit(rays[static_cast<size_t>(i)], r2) ? 1 : 0; }
            const auto m2 = clock::now();
            const double flatMt = double(nRays) / 1e6 / (ms(m0, m1) / 1000.0), embreeMt = double(nRays) / 1e6 / (ms(m1, m2) / 1000.0);
            std::cout << "[Bench]     all threads: flat BVH " << flatMt << " Mray/s (x" << ms(c1, c2) / ms(m0, m1) << " of one thread),  Embree "
                << embreeMt << " Mray/s (x" << ms(e2, e3) / ms(m1, m2) << " of one thread),  Embree / flat x" << embreeMt / flatMt
                << "  (" << hf << " / " << he << " hits)\n";
        }

        if (subdiv == 6) {
            // Coherent rays (camera-like scanline order from one point) on the same BVH: shows how much
            // of the random-ray cost is cache misses.
            const int side = static_cast<int>(std::sqrt(double(nRays)));
            std::vector<Ray> cam;
            cam.reserve(static_cast<size_t>(side) * static_cast<size_t>(side));
            for (int y = 0; y < side; ++y)
                for (int x = 0; x < side; ++x) {
                    const Vector3 target(-1.2 + 2.4 * (x + 0.5) / side, -1.2 + 2.4 * (y + 0.5) / side, 0.0);
                    cam.emplace_back(Point3(0, 0, 3), glm::normalize(target - Vector3(0, 0, 3)));
                }
            int h = 0;
            const auto s0 = clock::now();
            for (const Ray& r : cam) h += flat.hit(r, rec) ? 1 : 0;
            const auto s1 = clock::now();
            std::cout << "[Bench]     coherent camera rays (" << cam.size() << "): " << double(cam.size()) / 1e6 / (ms(s0, s1) / 1000.0)
                << " Mray/s (" << h << " hits)\n";

            // Leaf size sensitivity (closest hit, same rays).
            for (const int leaf : { 1, 2, 4, 8 }) {
                BVH b(prims, leaf);
                int h = 0;
                const auto s0 = clock::now();
                for (size_t i = 0; i < rays.size(); ++i) h += b.hit(rays[i], rec) ? 1 : 0;
                const auto s1 = clock::now();
                std::cout << "[Bench]     leaf size " << leaf << ": " << double(nRays) / 1e6 / (ms(s0, s1) / 1000.0) << " Mray/s, "
                    << b.numNodes() << " nodes (" << h << " hits)\n";
            }
        }

        const double mr = double(nRays) / 1e6;
        std::cout << "[Bench] " << mesh->numTriangles() << " triangles:  build  recursive " << ms(t0, t1) << " ms,  SAH flat " << ms(t1, t2)
            << " ms  (" << flat.numNodes() << " nodes, depth " << flat.maxDepth() << ")\n"
            << "[Bench]     closest hit  recursive " << mr / (ms(c0, c1) / 1000.0) << " Mray/s,  flat " << mr / (ms(c1, c2) / 1000.0)
            << " Mray/s  (x" << ms(c0, c1) / ms(c1, c2) << ");  any-hit flat " << mr / (ms(c2, c3) / 1000.0) << " Mray/s\n"
            << "[Bench]     results: " << hitsLegacy << " / " << hitsFlat << " hits, " << mismatches << " distance mismatches, any-hit = "
            << occ << "   " << ((mismatches == 0 && occ == hitsFlat) ? "OK" : "FAIL") << "\n";
    }
    std::cout << "[Bench] RESULT: " << (ok ? "PASS" : "FAIL") << "\n";
    return ok ? 0 : 2;
}

} // namespace rayt::cli
