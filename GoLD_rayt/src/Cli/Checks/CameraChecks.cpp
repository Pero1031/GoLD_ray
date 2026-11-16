// Camera checks (--scene camera): the perspective and the orthographic view against closed forms.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Color/ChromaticAdaptation.hpp"
#include "Color/ColorSpace.hpp"
#include "Geometry/Sphere.hpp"
#include "Materials/Lambertian.hpp"
#include "Render/Camera.hpp"
#include "Render/Diffraction.hpp"
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

    // ---- Diffraction at the aperture (--diffraction; Render/Diffraction.hpp) ----
    const double pi = constants::PI;

    // (e) The Fourier transform of the aperture against a direct sum over a fine grid of the aperture (disc, hexagon, heptagon).
    {
        double worst = 0;
        for (const int blades : { 0, 6, 7 }) {
            const double rot = 0.3, rc = blades ? std::sqrt(pi / (2.0 * blades * std::sin(2 * pi / blades))) : 0.5;
            auto inside = [&](double x, double y) {
                if (!blades) return x * x + y * y <= 0.25;
                for (int k = 0; k < blades; ++k) {
                    const double a = rot + 2 * pi * (k + 0.5) / blades;
                    if (x * std::cos(a) + y * std::sin(a) > rc * std::cos(pi / blades)) return false;
                }
                return true;
            };
            const int m = 1200;
            const double cell = 1.2 / m;
            for (const Point2 q : { Point2(0.0, 0.0), Point2(7.0, 2.0), Point2(-3.0, 11.0), Point2(20.0, -9.0) }) {
                std::complex<double> sum(0.0);
                for (int i = 0; i < m; ++i)
                    for (int j = 0; j < m; ++j) {
                        const double x = -0.6 + (i + 0.5) * cell, y = -0.6 + (j + 0.5) * cell;
                        if (inside(x, y)) sum += std::complex<double>(std::cos(q.x * x + q.y * y), -std::sin(q.x * x + q.y * y)) * (cell * cell);
                    }
                worst = std::max(worst, std::abs(ApertureTransform(blades, rot, q.x, q.y) - sum) / (pi / 4));
            }
        }
        report("e  aperture transform (disc, 6 and 7 blades) against a sum over the aperture", worst < 2e-3, "largest difference " + fmt(worst) + " of the area");
    }

    // (f) The pattern of one wavelength, in units of lambda N: the Airy disc holds 83.8 % of the light inside its first dark ring (1.2197),
    //     and disc and heptagon hold nearly all of it within 60 units (the tail falls as 1 / radius).
    {
        auto encircled = [&](int blades, double radius) {
            const int nr = 2400, na = 3600;
            double sum = 0;
            for (int i = 0; i < nr; ++i)
                for (int j = 0; j < na; ++j) {
                    const double rr = (i + 0.5) * radius / nr, a = (j + 0.5) * 2 * pi / na;
                    sum += std::norm(ApertureTransform(blades, 0.3, 2 * pi * rr * std::cos(a), 2 * pi * rr * std::sin(a))) / (pi / 4) * rr * (radius / nr) * (2 * pi / na);
                }
            return sum;
        };
        const double airy = encircled(0, 1.2197), disc = encircled(0, 60), hept = encircled(7, 60);
        report("f  Airy disc: 83.8 % inside the first dark ring; disc and heptagon near 1 within 60 lambda N",
            std::fabs(airy - 0.8378) < 0.002 && disc > 0.99 && disc < 1.0005 && hept > 0.99 && hept < 1.0005,
            "first ring " + fmt(airy) + ", within 60: " + fmt(disc) + " and " + fmt(hept));
    }

    // (g) The spikes of an enlarged pattern: 6 for a hexagon, 14 for a heptagon, perpendicular to the edges (Y kernel on a ring of 25 px).
    {
        std::string detail;
        bool ok = true;
        for (const int blades : { 6, 7 }) {
            DiffractionSetup ds;
            ds.fNumber = 8; ds.blades = blades; ds.bladeRotation = 0.3; ds.pixelPitch = 30e-6; ds.scale = 12; ds.radius = 40;
            const DiffractionKernels k = ComputeDiffractionKernels(ds);
            const int size = 2 * k.radius + 1, na = 1440;
            std::vector<double> ring(na);
            for (int j = 0; j < na; ++j) {   // bilinear on the kernel; angle counter-clockwise from the right, the kernel's rows go down
                const double a = 2 * pi * j / na, fx = k.radius + 25 * std::cos(a), fy = k.radius - 25 * std::sin(a);
                const int x0 = int(fx), y0 = int(fy);
                const double tx = fx - x0, ty = fy - y0;
                auto at = [&](int x, int y) { return k.xyz[1][static_cast<size_t>(y) * size + x]; };
                ring[j] = (1 - tx) * (1 - ty) * at(x0, y0) + tx * (1 - ty) * at(x0 + 1, y0) + (1 - tx) * ty * at(x0, y0 + 1) + tx * ty * at(x0 + 1, y0 + 1);
            }
            std::vector<double> sorted = ring;
            std::nth_element(sorted.begin(), sorted.begin() + na / 2, sorted.end());
            const double median = sorted[na / 2];
            int spikes = 0;
            double worstAngle = 0;
            const int expected = (blades % 2) ? 2 * blades : blades;
            for (int j = 0; j < na; ++j) {
                const double v = ring[j];
                if (v < 4 * median || v < ring[(j + na - 1) % na] || v < ring[(j + 1) % na]) continue;
                ++spikes;
                // the nearest edge normal: 0.3 + pi (2k + 1) / n, every pi / n for an odd n (both ends of every spike)
                const double a = 2 * pi * j / na, period = (blades % 2) ? pi / blades : 2 * pi / blades;
                const double off = std::remainder(a - (0.3 + pi / blades), period);
                worstAngle = std::max(worstAngle, std::fabs(off) * 180 / pi);
            }
            ok = ok && spikes == expected && worstAngle < 1.5;
            detail += std::to_string(blades) + " blades: " + std::to_string(spikes) + " spikes, off by up to " + fmt(worstAngle) + " deg; ";
        }
        report("g  spikes: 6 for 6 blades, 14 for 7, across the edges", ok, detail.substr(0, detail.size() - 2));
    }

    // (h) The convolution: a uniform picture stays uniform; one bright white pixel spreads as the kernel (its Y, with the kernel's orientation:
    //     the heptagon is not mirror-symmetric about the horizontal, so a flipped kernel would not match).
    {
        DiffractionSetup ds;
        ds.fNumber = 8; ds.blades = 7; ds.bladeRotation = 0.3; ds.pixelPitch = 30e-6; ds.scale = 12; ds.radius = 30;
        const DiffractionKernels k = ComputeDiffractionKernels(ds);
        const int w = 101, h = 81, size = 2 * k.radius + 1;
        Film flat(w, h), spot(w, h);
        for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x) flat.setPixel(x, y, Vector3(0.3, 0.5, 0.7));
        spot.setPixel(50, 40, Vector3(1000.0));
        ApplyDiffraction(flat, k);
        ApplyDiffraction(spot, k);
        double flatErr = 0, spotErr = 0, mirror = 0;
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x) flatErr = std::max(flatErr, double(glm::length(flat.getData()[static_cast<size_t>(y) * w + x] - Vector3(0.3, 0.5, 0.7))));
        // The expected spread: the white's X, Y, Z, each times its kernel, back to linear sRGB (the matrices of Render/Diffraction.cpp)
        color::Matrix3 toRgb;
        for (int c = 0; c < 3; ++c) { color::CIE_XYZ e; e.X = c == 0; e.Y = c == 1; e.Z = c == 2; toRgb[c] = color::XYZ_to_sRGB(e); }
        const Vector3 white = glm::inverse(toRgb) * Vector3(1000.0);
        auto expected = [&](int dx, int dy) {
            const size_t i = static_cast<size_t>(dy + k.radius) * size + dx + k.radius;
            return toRgb * Vector3(white.x * k.xyz[0][i], white.y * k.xyz[1][i], white.z * k.xyz[2][i]);
        };
        for (int dy = -k.radius; dy <= k.radius; ++dy)
            for (int dx = -k.radius; dx <= k.radius; ++dx) {
                const Vector3 got = spot.getData()[static_cast<size_t>(40 + dy) * w + 50 + dx];
                spotErr = std::max(spotErr, double(glm::length(got - expected(dx, dy))));
                mirror = std::max(mirror, double(glm::length(expected(dx, dy) - expected(dx, -dy))));
            }
        report("h  convolution: a uniform picture stays uniform, a bright pixel spreads as the kernel (not mirrored)",
            flatErr < 1e-9 && spotErr < 1e-6 && mirror > 1e-3,
            "uniform " + fmt(flatErr) + ", spot " + fmt(spotErr) + " against a mirror difference of " + fmt(mirror));
    }

    std::cout << "[Camera] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
