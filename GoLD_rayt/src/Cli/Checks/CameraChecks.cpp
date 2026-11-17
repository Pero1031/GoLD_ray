// Camera checks (--scene camera): the perspective and the orthographic view against closed forms.
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Color/ChromaticAdaptation.hpp"
#include "Color/ColorSpace.hpp"
#include "Geometry/Quad.hpp"
#include "Geometry/Sphere.hpp"
#include "Materials/Lambertian.hpp"
#include "Render/LensSystem.hpp"
#include "Render/Camera.hpp"
#include "Render/Diffraction.hpp"
#include "Render/Film.hpp"
#include "Render/MaskDepth.hpp"
#include "Scene/Aggregate.hpp"
#include "Scene/Scene.hpp"

namespace rayt::cli {

namespace {
    /// A deterministic number in [0, 1) for the checks' lens samples.
    Real Hash01Check(int i, int k) {
        std::uint64_t h = static_cast<std::uint64_t>(i) * 0x9E3779B97F4A7C15ull + static_cast<std::uint64_t>(k) * 0xC2B2AE3D27D4EB4Full;
        h ^= h >> 33; h *= 0xFF51AFD7ED558CCDull; h ^= h >> 33; h *= 0xC4CEB9FE1A85EC53ull; h ^= h >> 33;
        return Real(double(h >> 11) * (1.0 / 9007199254740992.0));
    }
}

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

    // ---- A real lens (--lens; Render/LensSystem.hpp) ----
    // The singlet of assets/Lenses/singlet_bk7_100mm.txt: a 25 mm stop 2 mm in front of an equiconvex N-BK7 lens, R = 103.36 mm, 5 mm thick.
    const Real nBk7 = 1.5168, rLens = 0.10336, thick = 0.005;
    auto singlet = [&]() {
        std::vector<LensSystem::Surface> s(3);
        s[0].radius = 0;       s[0].z = 0;                   s[0].apertureRadius = 0.0125; s[0].nd = 1;
        s[1].radius = rLens;   s[1].z = -0.002;              s[1].apertureRadius = 0.015;  s[1].nd = nBk7; s[1].abbe = 64.17;
        s[2].radius = -rLens;  s[2].z = -0.002 - thick;      s[2].apertureRadius = 0.015;  s[2].nd = 1;
        return std::make_shared<LensSystem>(s);
    };

    // (i) Its focal length against the thick-lens formula 1/f = (n-1) (1/R1 - 1/R2 + (n-1) d / (n R1 R2)); the principal planes of the
    //     symmetric lens lie symmetric to its vertices.
    {
        const auto lens = singlet();
        const Real f = Real(1) / ((nBk7 - 1) * (Real(2) / rLens - (nBk7 - 1) * thick / (nBk7 * rLens * rLens)));
        const Real frontDepth = lens->frontPrincipalZ() - Real(-0.002), rearDepth = Real(-0.002 - thick) - lens->rearPrincipalZ();
        report("i  real lens: a singlet's focal length is the thick-lens formula's, its principal planes symmetric",
            std::fabs(lens->focalLength() - f) < 1e-8 && std::fabs(frontDepth - rearDepth) < 1e-8,
            "f " + fmt(lens->focalLength() * 1000) + " mm against " + fmt(f * 1000) + " mm");
    }

    // (j) Stopped down to f/50 and focused at 1 m: rays from the axis point 1 m in front of the front principal plane, near the axis, cross
    //     it on the film (within 0.1 um). Farther out they cross it nearer the lens, by an amount that grows with the square of the height
    //     in the stop (third-order spherical aberration of an uncorrected positive lens): 0.8 mm misses about 4 times as much as 0.4 mm.
    {
        const auto lens = singlet();
        lens->setFNumber(50);
        lens->focus(1.0);
        const Point3 object(0, 0, lens->frontPrincipalZ() + 1.0);
        auto crossing = [&](Real h) {   // where the ray through the stop at height h crosses the axis, relative to the film (+: toward the lens)
            Point3 o = object;
            Vector3 d = glm::normalize(Point3(h, 0, Real(0)) - object);
            if (!lens->traceFromObject(o, d, 587.56)) return Real(1);
            return o.z + (-o.x / d.x) * d.z - lens->filmZ();
        };
        Real paraxialMiss = 0;
        for (const Real h : { 10e-6, 30e-6, 60e-6 }) paraxialMiss = std::max(paraxialMiss, std::fabs(crossing(h)));
        const Real m4 = crossing(0.4e-3), m8 = crossing(0.8e-3);
        report("j  real lens: near the axis the rays meet on the film; farther out spherical aberration grows as the height squared",
            paraxialMiss < 1e-7 && m4 > 0 && std::fabs(m8 / m4 - 4) < 0.4,
            "paraxial " + fmt(paraxialMiss * 1e6) + " um; at 0.4 / 0.8 mm " + fmt(m4 * 1e6) + " / " + fmt(m8 * 1e6) + " um nearer the lens");
    }

    // (k) The f-number: wide open the stop is the entrance pupil (nothing in front of it), f / 25 mm; closed to f/8 the pupil is f / 8.
    {
        const auto lens = singlet();
        const Real open = lens->fNumber(), want = lens->focalLength() / 0.025;
        const Real set = lens->setFNumber(8);
        report("k  real lens: f-number wide open is f / stop diameter, and closes to f/8",
            std::fabs(open / want - 1) < 1e-6 && std::fabs(set - 8) < 1e-6, "f/" + fmt(open) + " against f/" + fmt(want) + ", then f/" + fmt(set));
    }

    // (l) Dispersion (spectral build): the focal lengths at the F and C lines differ by about f / V_d.
#if defined(RAYT_SPECTRAL_MODE)
    {
        const auto lens = singlet();
        const Real fF = lens->focalLengthAt(486.13), fC = lens->focalLengthAt(656.27);
        const Real ratio = (fC - fF) / lens->focalLength() * 64.17;
        report("l  real lens: the focal lengths at the F and C lines differ by f / V_d (2 %)", lens->dispersive() && std::fabs(ratio - 1) < 0.02,
            "(f_C - f_F) V_d / f = " + fmt(ratio));
    }
#else
    std::cout << "[Camera]   l  real lens dispersion: needs the spectral build (skipped)\n";
#endif

    // (n) Dispersion through a whole lens (spectral build): the double-Gauss with the patent's glasses. The real ray traced a micrometre
    //     off the axis crosses it where a paraxial y-u trace of the prescription (a separate, matrix-like calculation) puts the back focus,
    //     at the F, d and C lines. The lens is partly achromatised: its F and C foci differ by 0.19 % of f, about a ninth of what a single
    //     lens of these glasses would show (1 / V_d, 1.6 to 3.3 %); a fast normal lens wide open shows this as coloured bokeh fringes.
#if defined(RAYT_SPECTRAL_MODE)
    {
        std::string error;
        const auto lens = LensSystem::load("assets/Lenses/dgauss_50mm.txt", error);
        if (!lens) report("n  dispersion through a lens (" + error + ")", false);
        else {
            const auto& S = lens->surfaces();
            auto realFocus = [&](Real nm) {   // where a ray parallel to the axis, 1 um off it, crosses the axis behind the last vertex
                Point3 o(Real(1e-6), 0, Real(0.01));
                Vector3 d(0, 0, -1);
                if (!lens->traceFromObject(o, d, nm) || std::fabs(d.x) < Real(1e-15)) return Real(0);
                return S.back().z - (o.z + (-o.x / d.x) * d.z);
            };
            auto paraxialFocus = [&](Real nm) {   // y-u trace: refraction n'u' = nu - y c (n' - n), transfer y += u t
                const Real lf = 0.48613, ld = 0.58756, lc = 0.65627, l = nm / 1000;
                Real y = 1, u = 0, n = 1;
                for (size_t i = 0; i < S.size(); ++i) {
                    Real n2 = S[i].nd;
                    if (n2 > 1 && S[i].abbe > 0) {
                        const Real B = (S[i].nd - 1) / (S[i].abbe * (1 / (lf * lf) - 1 / (lc * lc)));
                        n2 = S[i].nd - B / (ld * ld) + B / (l * l);
                    }
                    if (n2 < 1) n2 = 1;   // the stop
                    const Real c = S[i].radius == 0 ? Real(0) : 1 / S[i].radius;
                    u = (n * u - y * c * (n2 - n)) / n2;
                    n = n2;
                    if (i + 1 < S.size()) y += u * (S[i].z - S[i + 1].z);
                }
                return -y / u;
            };
            Real worst = 0;
            for (const Real nm : { 486.13, 587.56, 656.27 }) worst = std::max(worst, std::fabs(realFocus(nm) - paraxialFocus(nm)));
            const Real f = lens->focalLength();
            const Real fc = (paraxialFocus(486.13) - paraxialFocus(656.27)) / f;
            report("n  dispersion through a lens: the double-Gauss's foci at F, d, C as a paraxial y-u trace (0.1 um); f = 50 mm",
                lens->dispersive() && worst < 1e-7 && std::fabs(f - 0.050) < 1e-4,
                "largest difference " + fmt(worst * 1e6) + " um, F - C " + fmt(fc * 100) + " % of f, f = " + fmt(f * 1000) + " mm");
        }
    }
#endif

    // (m) The double-Gauss in a camera at the origin looking down -z: the weight averages 1 at the film centre and falls in the corner
    //     (vignetting); the picture is upright (its upper right looks up and to the right); a wall 10 m away has depth 10 in the middle
    //     and the corners.
    {
        std::string error;
        const auto lens = LensSystem::load("assets/Lenses/dgauss_50mm.txt", error);
        if (!lens) report("m  real lens in a camera (" + error + ")", false);
        else {
            lens->focus(10.0);
            lens->prepare(0.036, 0.024);
            auto cam = std::make_shared<Camera>(Point3(0.0), Point3(0, 0, -1), Vector3(0, 1, 0), 40.0, 1.5, 0.0, 10.0);
            cam->setLensSystem(lens);
            auto meanWeight = [&](Real s, Real t) {
                double sum = 0;
                for (int i = 0; i < 20000; ++i) {
                    Real w = 0;
                    cam->getRay(s, t, Point2(Hash01Check(i, 0), Hash01Check(i, 1)), 0, 550, &w);
                    sum += w;
                }
                return sum / 20000;
            };
            const double centre = meanWeight(0.5, 0.5), corner = meanWeight(0.97, 0.97);
            const Ray upRight = cam->centerRay(0.9, 0.9);
            auto mat = std::make_shared<Lambertian>(RgbColor(Vector3(0.5)));
            auto world = std::make_shared<Aggregate>();
            world->add(std::make_shared<Quad>(Point3(-50, -50, -10), Vector3(100, 0, 0), Vector3(0, 100, 0), mat));
            auto scene = std::make_shared<Scene>(world);
            scene->finalize("uniform");
            Film f(30, 20);
            f.enableDepth();
            RenderMaskAndDepth(*scene, *cam, f, {}, 1);
            const float dc = f.depth()[10 * 30 + 15], d0 = f.depth()[0], d1 = f.depth()[30 * 20 - 1];
            const bool ok = std::fabs(centre - 1) < 0.02 && corner < 0.9 * centre && upRight.d.x > 0 && upRight.d.y > 0
                && std::fabs(dc - 10) < 1e-4 && std::fabs(d0 - 10) < 1e-4 && std::fabs(d1 - 10) < 1e-4;
            report("m  real lens in a camera: weight 1 at the centre, less in the corner; upright picture; depth of a wall",
                ok, "weight " + fmt(centre) + " / " + fmt(corner) + ", depth " + fmt(dc) + " " + fmt(d0) + " " + fmt(d1));
        }
    }

    std::cout << "[Camera] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
