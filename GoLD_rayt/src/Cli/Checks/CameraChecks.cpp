// Camera checks (--scene camera): the perspective and the orthographic view against closed forms.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Geometry/Sphere.hpp"
#include "Materials/Lambertian.hpp"
#include "Render/Camera.hpp"
#include "Render/Film.hpp"
#include "Render/MaskDepth.hpp"
#include "Scene/Aggregate.hpp"
#include "Scene/Scene.hpp"

namespace rayt::cli {

int RunCameraTests() {
    bool allPassed = true;
    auto report = [&](const std::string& name, bool ok, const std::string& detail = "") {
        std::cout << "[Camera]   " << name << (detail.empty() ? "" : "   (" + detail + ")") << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };
    auto fmt = [](double v) { std::ostringstream s; s << std::setprecision(3) << v; return s.str(); };
    // A camera at (0.3, 0.2, 4) looking at the origin at a slant, so the basis is not the world's axes; 4:3 picture.
    const Point3 from(0.3, 0.2, 4.0), at(0.0, 0.0, 0.0);
    const Vector3 up(0, 1, 0);
    const Real aspect = Real(4) / Real(3), focus = glm::length(from - at);
    const Vector3 axis = glm::normalize(at - from);
    // Where a ray meets the plane of focus (perpendicular to the view axis, `focus` in front of the camera position).
    auto onFocusPlane = [&](const Ray& r, const Point3& planePoint) { return r.o + r.d * (glm::dot(planePoint - r.o, axis) / glm::dot(r.d, axis)); };
    const Point2 films[5] = { Point2(0.5, 0.5), Point2(0.05, 0.1), Point2(0.93, 0.2), Point2(0.15, 0.88), Point2(0.97, 0.96) };

    // (a) Perspective: project() undoes the ray through a pixel centre.
    {
        const Camera cam(from, at, up, 40.0, aspect, 0.0, focus);
        Real worst = 0;
        for (const Point2& f : films) {
            const Ray r = cam.centerRay(f.x, f.y);
            Real s, t, cosCam, dist;
            const bool in = cam.project(r.o + r.d * Real(3.7), s, t, cosCam, dist);
            worst = std::max(worst, in ? std::max(std::fabs(s - f.x), std::fabs(t - f.y)) : Real(1));
        }
        report("a  perspective: project() is the inverse of the pixel's ray", worst < 1e-9, "largest error " + fmt(worst));
    }

    // (b) A lens (f/2.8 of a 50 mm lens): the rays of every lens point of a film position meet on the plane of focus, where the pixel's
    //     centre ray meets it; in the perspective view and in the orthographic one.
    for (const bool ortho : { false, true }) {
        Camera cam(from, at, up, 40.0, aspect, 0.05 / 2.8, focus);
        if (ortho) cam.setOrthographic(1.5);
        Real worst = 0, spread = 0;
        for (const Point2& f : films) {
            const Ray c = cam.centerRay(f.x, f.y);
            const Point3 sharp = c.o + c.d * (focus / glm::dot(c.d, axis));
            for (int k = 0; k < 64; ++k) {
                const Ray r = cam.getRay(f.x, f.y, sampling::Random2D());
                worst = std::max(worst, Real(glm::length(onFocusPlane(r, sharp) - sharp)));
                spread = std::max(spread, Real(glm::length(r.o - c.o)));
            }
        }
        report(std::string(ortho ? "b2 orthographic" : "b1 perspective") + ": with an aperture, the rays of a film position meet on the plane of focus",
            worst < 1e-9 && spread > 0.005, "largest miss " + fmt(worst) + " m, lens points up to " + fmt(spread * 1000) + " mm apart");
    }

    // (c) Orthographic: every ray runs along the view axis, from a rectangle 1.5 high and 2 wide centred on the camera position.
    {
        Camera cam(from, at, up, 40.0, aspect, 0.0, focus);
        cam.setOrthographic(1.5);
        Real worstDir = 0, worstPos = 0;
        for (const Point2& f : films) {
            const Ray r = cam.centerRay(f.x, f.y);
            worstDir = std::max(worstDir, Real(1) - glm::dot(r.d, axis));
            const Vector3 offset = r.o - from;
            const Vector3 right = glm::normalize(glm::cross(axis, up)), upCam = glm::cross(right, axis);
            const Real ex = (f.x - Real(0.5)) * Real(1.5) * aspect, ey = (f.y - Real(0.5)) * Real(1.5);
            worstPos = std::max(worstPos, Real(glm::length(offset - (right * ex + upCam * ey))));
        }
        Real s, t, c, d;
        report("c  orthographic: rays parallel to the axis, from a 2 x 1.5 rectangle; project() refuses",
            worstDir < 1e-12 && worstPos < 1e-12 && !cam.project(at, s, t, c, d), "direction " + fmt(worstDir) + ", position " + fmt(worstPos));
    }

    // (d) Orthographic: a sphere covers the same pixels at 3 m and at 9 m, pi R^2 over the pixel area; its depth follows the sphere.
    {
        const int size = 65;
        const Real height = 2.0, R = 0.5;
        const auto mat = std::make_shared<Lambertian>(RgbColor(Vector3(0.5)));
        const Camera cam(Point3(0.0), Point3(0, 0, -1), up, 40.0, 1.0, 0.0, 1.0);
        Camera ortho = cam;
        ortho.setOrthographic(height);
        const Real pixel = height / Real(size);
        const Real expected = constants::PI * R * R / (pixel * pixel);
        double area[2] = { 0, 0 };
        Real worstDepth = 0;
        const Real dist[2] = { 3.0, 9.0 };
        for (int k = 0; k < 2; ++k) {
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Sphere>(Point3(0, 0, -dist[k]), R, mat));
            auto scene = std::make_shared<Scene>(world);
            scene->finalize("uniform");
            Film f(size, size);
            f.enableMask();
            f.enableDepth();
            RenderMaskAndDepth(*scene, ortho, f, { mat.get() }, 16);
            for (const float c : f.mask()) area[k] += c;
            // depth at a pixel centre that sees the sphere: D - sqrt(R^2 - x^2 - y^2), x and y its offset from the axis
            for (int y = 0; y < size; ++y)
                for (int x = 0; x < size; ++x) {
                    const Real px = (Real(x) + Real(0.5)) * pixel - height / 2, py = height / 2 - (Real(y) + Real(0.5)) * pixel;   // row 0 is the top
                    const Real q = R * R - px * px - py * py;
                    const Real z = f.depth()[static_cast<size_t>(y) * size + x];
                    worstDepth = std::max(worstDepth, q > Real(1e-4) ? std::fabs(z - (dist[k] - std::sqrt(q))) : Real(0));
                }
        }
        report("d  orthographic: a sphere's pixel area is pi R^2 at 3 m and at 9 m (0.5 %), its depth follows the sphere",
            std::fabs(area[0] - expected) < 0.005 * expected && std::fabs(area[1] - expected) < 0.005 * expected && worstDepth < 1e-5,
            "areas " + fmt(area[0]) + " and " + fmt(area[1]) + " px against " + fmt(expected) + ", depth error " + fmt(worstDepth));
    }

    std::cout << "[Camera] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
