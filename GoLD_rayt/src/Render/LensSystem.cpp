// A real lens traced ray by ray (--lens); see Render/LensSystem.hpp.
#include "pch.h"
#include "Render/LensSystem.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

namespace rayt {

    namespace {

        constexpr Real kLineD = Real(587.56), kLineF = Real(486.13), kLineC = Real(656.27);   // Fraunhofer lines (nm)
        constexpr Real kPi = Real(3.14159265358979323846);

        /// Refracts the unit direction d at a surface whose unit normal n faces against d; eta = n_from / n_to. False on total internal reflection.
        bool Refract(Vector3& d, const Vector3& n, Real eta) {
            const Real cosi = -glm::dot(n, d);
            const Real sin2t = eta * eta * (Real(1) - cosi * cosi);
            if (sin2t >= Real(1)) return false;
            d = glm::normalize(eta * d + (eta * cosi - std::sqrt(Real(1) - sin2t)) * n);
            return true;
        }

        void SetCauchy(LensSystem::Surface& s) {
            // the dispersion of Materials/Dielectric: n = A + B / lambda^2 (micrometres) through n_d with n_F - n_C = (n_d - 1) / V_d
            s.cauchyA = s.nd;
            s.cauchyB = 0;
            if (s.abbe > Real(0) && s.nd > Real(1)) {
                const Real d = kLineD * Real(1e-3), f = kLineF * Real(1e-3), c = kLineC * Real(1e-3);
                s.cauchyB = (s.nd - Real(1)) / (s.abbe * (Real(1) / (f * f) - Real(1) / (c * c)));
                s.cauchyA = s.nd - s.cauchyB / (d * d);
            }
        }

    } // namespace

    std::shared_ptr<LensSystem> LensSystem::load(const std::string& path, std::string& error) {
        std::ifstream in(path);
        if (!in) { error = "cannot open " + path; return nullptr; }
        std::vector<Surface> surfaces;
        Real z = 0;
        int lineNo = 0;
        for (std::string line; std::getline(in, line);) {
            ++lineNo;
            if (const size_t hash = line.find('#'); hash != std::string::npos) line.erase(hash);
            std::istringstream ss(line);
            std::vector<double> v;
            for (double x; ss >> x;) v.push_back(x);
            if (v.empty()) continue;
            if (v.size() < 4 || v.size() > 5) { error = path + ", line " + std::to_string(lineNo) + ": expected radius, thickness, n_d, diameter [, Abbe number]"; return nullptr; }
            Surface s;
            s.radius = Real(v[0] * 1e-3);
            s.z = z;
            s.nd = v[2] > 1 ? Real(v[2]) : Real(1);
            s.apertureRadius = Real(v[3] * 0.5e-3);
            s.abbe = v.size() == 5 ? Real(v[4]) : Real(0);
            if (!(s.apertureRadius > 0)) { error = path + ", line " + std::to_string(lineNo) + ": the diameter must be above 0"; return nullptr; }
            if (s.radius != 0 && std::fabs(s.radius) < s.apertureRadius) { error = path + ", line " + std::to_string(lineNo) + ": the radius is smaller than the semi-diameter"; return nullptr; }
            surfaces.push_back(s);
            z -= Real(v[1] * 1e-3);
        }
        if (surfaces.empty()) { error = path + ": no surfaces"; return nullptr; }
        auto lens = std::make_shared<LensSystem>(std::move(surfaces));
        if (!(lens->focalLength() > 0)) { error = path + ": not a converging lens (no positive focal length at 587.6 nm)"; return nullptr; }
        return lens;
    }

    LensSystem::LensSystem(std::vector<Surface> surfaces) : m_surfaces(std::move(surfaces)) {
        for (Surface& s : m_surfaces) SetCauchy(s);
        paraxial();
        m_filmZ = m_rearFocal;
        m_wideOpen = fNumber();
    }

    Real LensSystem::indexAfter(int i, Real nm) const {
        if (i < 0) return Real(1);
        const Surface& s = m_surfaces[static_cast<size_t>(i)];
        if (s.nd <= Real(1)) return Real(1);
#if defined(RAYT_SPECTRAL_MODE)
        if (s.cauchyB != Real(0) && nm > Real(0)) {
            const Real um = nm * Real(1e-3);
            return s.cauchyA + s.cauchyB / (um * um);
        }
#else
        (void)nm;
#endif
        return s.nd;
    }

    bool LensSystem::dispersive() const {
#if defined(RAYT_SPECTRAL_MODE)
        for (const Surface& s : m_surfaces) if (s.cauchyB != Real(0)) return true;
#endif
        return false;
    }

    bool LensSystem::insideStop(const Surface& s, Real x, Real y, const StopShape& stop) const {
        const Real r = s.apertureRadius * stop.scale;
        if (stop.blades < 3) return x * x + y * y <= r * r;
        // a regular polygon of the circle's area (as Camera::setApertureShape): every edge's inner side
        const Real n = Real(stop.blades);
        const Real circumradius = r * std::sqrt(kPi / (Real(0.5) * n * std::sin(Real(2) * kPi / n)));
        const Real apothem = circumradius * std::cos(kPi / n);
        for (int k = 0; k < stop.blades; ++k) {
            const Real a = stop.rotation + Real(2) * kPi * (Real(k) + Real(0.5)) / n;
            if (x * std::cos(a) + y * std::sin(a) > apothem) return false;
        }
        return true;
    }

    bool LensSystem::crossSurface(int i, Point3& o, Vector3& d, Real nFrom, Real nTo, const StopShape& stop) const {
        const Surface& s = m_surfaces[static_cast<size_t>(i)];
        if (s.radius == Real(0)) {   // the aperture stop: a plane, no refraction
            if (std::fabs(d.z) < Real(1e-12)) return false;
            const Real t = (s.z - o.z) / d.z;
            if (t <= Real(0)) return false;
            const Point3 p = o + d * t;
            if (!insideStop(s, p.x, p.y, stop)) return false;
            o = p;
            return true;
        }
        // a sphere through the vertex, its centre on the axis; the intersection on the vertex's side of the sphere
        const Point3 centre(0, 0, s.z - s.radius);
        const Vector3 oc = o - centre;
        const Real b = glm::dot(oc, d), c = glm::dot(oc, oc) - s.radius * s.radius;
        const Real disc = b * b - c;
        if (disc < Real(0)) return false;
        const Real root = std::sqrt(disc);
        const Real t = ((d.z > Real(0)) == (s.radius > Real(0))) ? -b + root : -b - root;
        if (t <= Real(1e-12)) return false;
        const Point3 p = o + d * t;
        if (p.x * p.x + p.y * p.y > s.apertureRadius * s.apertureRadius) return false;
        Vector3 normal = glm::normalize(p - centre);
        if (glm::dot(normal, d) > Real(0)) normal = -normal;
        if (!Refract(d, normal, nFrom / nTo)) return false;
        o = p;
        return true;
    }

    bool LensSystem::traceFromFilm(Point3& o, Vector3& d, Real nm) const {
        return traceFromFilmWith(o, d, nm, StopShape{ m_stopScale, m_blades, m_bladeRotation });
    }

    bool LensSystem::traceFromFilmWith(Point3& o, Vector3& d, Real nm, const StopShape& stop) const {
        for (int i = static_cast<int>(m_surfaces.size()) - 1; i >= 0; --i)
            if (!crossSurface(i, o, d, indexAfter(i, nm), indexAfter(i - 1, nm), stop)) return false;
        return true;
    }

    bool LensSystem::traceFromObject(Point3& o, Vector3& d, Real nm) const {
        const StopShape stop{ m_stopScale, m_blades, m_bladeRotation };
        for (int i = 0; i < static_cast<int>(m_surfaces.size()); ++i)
            if (!crossSurface(i, o, d, indexAfter(i - 1, nm), indexAfter(i, nm), stop)) return false;
        return true;
    }

    void LensSystem::paraxial() {
        // A ray parallel to the axis a micrometre off it, from each side: where it crosses the axis is a focal point, where its extension
        // meets the height it came in at is a principal plane.
        const Real h = Real(1e-6);
        Point3 o(h, 0, Real(0.01));
        Vector3 d(0, 0, -1);
        if (traceFromObject(o, d, kLineD) && std::fabs(d.x) > Real(1e-15)) {
            m_rearFocal = o.z + (-o.x / d.x) * d.z;
            m_rearPrincipal = o.z + ((h - o.x) / d.x) * d.z;
            m_focal = m_rearPrincipal - m_rearFocal;
        }
        o = Point3(h, 0, m_surfaces.back().z - Real(0.01));
        d = Vector3(0, 0, 1);
        if (traceFromFilm(o, d, kLineD) && std::fabs(d.x) > Real(1e-15)) {
            m_frontFocal = o.z + (-o.x / d.x) * d.z;
            m_frontPrincipal = o.z + ((h - o.x) / d.x) * d.z;
        }
    }

    Real LensSystem::focalLengthAt(Real nm) const {
        const Real h = Real(1e-6);
        Point3 o(h, 0, Real(0.01));
        Vector3 d(0, 0, -1);
        if (!traceFromObject(o, d, nm) || std::fabs(d.x) < Real(1e-15)) return Real(0);
        return ((h - o.x) / d.x) * d.z - (-o.x / d.x) * d.z;   // principal minus focal point
    }

    Real LensSystem::entrancePupilRadius() const {
        auto passes = [&](Real h) {
            Point3 o(h, 0, Real(0.01));
            Vector3 d(0, 0, -1);
            return traceFromObject(o, d, kLineD);
        };
        const Real top = m_surfaces.front().apertureRadius;
        if (!passes(top * Real(1e-6))) return Real(0);
        if (passes(top)) return top;
        Real lo = 0, hi = top;
        for (int k = 0; k < 60; ++k) {
            const Real mid = Real(0.5) * (lo + hi);
            (passes(mid) ? lo : hi) = mid;
        }
        return lo;
    }

    Real LensSystem::fNumber() const {
        const Real r = entrancePupilRadius();
        return r > Real(0) ? m_focal / (Real(2) * r) : Real(1e30);
    }

    Real LensSystem::setFNumber(Real n) {
        const bool hasStop = std::any_of(m_surfaces.begin(), m_surfaces.end(), [](const Surface& s) { return s.radius == Real(0); });
        m_stopScale = 1;
        if (!hasStop || !(n > m_wideOpen)) return fNumber();
        Real lo = Real(1e-4), hi = Real(1);   // the f-number falls as the stop opens
        for (int k = 0; k < 60; ++k) {
            m_stopScale = Real(0.5) * (lo + hi);
            (fNumber() > n ? lo : hi) = m_stopScale;
        }
        m_stopScale = Real(0.5) * (lo + hi);
        return fNumber();
    }

    bool LensSystem::focus(Real distance) {
        if (!(distance > m_focal * Real(1.0001))) { m_filmZ = m_rearFocal; return false; }
        m_filmZ = m_rearPrincipal - m_focal * distance / (distance - m_focal);   // 1/D + 1/D' = 1/f, from the principal planes
        return true;
    }

    LensSystem::Bounds LensSystem::pupilBounds(Real r0, Real r1) const {
        // Film points on the +x axis at four radii of the band, rays to a 128 x 128 grid on the plane of the rear vertex: the box of the
        // grid points whose rays get through, widened by one grid step. A polygonal stop is tested as its circumcircle here, so the box
        // holds the pupil of a film point at any angle.
        StopShape circle{ m_stopScale, 0, Real(0) };
        if (m_blades >= 3) {
            const Real n = Real(m_blades);
            circle.scale *= std::sqrt(kPi / (Real(0.5) * n * std::sin(Real(2) * kPi / n)));
        }
        const Surface& rear = m_surfaces.back();
        const Real extent = Real(1.5) * rear.apertureRadius;
        const int grid = 128;
        const Real step = Real(2) * extent / Real(grid);
        Bounds b;
        for (int j = 0; j < 4; ++j) {
            const Real r = r0 + (r1 - r0) * Real(j) / Real(3);
            for (int gy = 0; gy < grid; ++gy)
                for (int gx = 0; gx < grid; ++gx) {
                    const Real x = -extent + (Real(gx) + Real(0.5)) * step, y = -extent + (Real(gy) + Real(0.5)) * step;
                    Point3 o(r, 0, m_filmZ);
                    Vector3 d = glm::normalize(Vector3(x - r, y, rear.z - m_filmZ));
                    if (!traceFromFilmWith(o, d, kLineD, circle)) continue;
                    if (b.empty) { b.minX = b.maxX = x; b.minY = b.maxY = y; b.empty = false; }
                    b.minX = std::min(b.minX, x); b.maxX = std::max(b.maxX, x);
                    b.minY = std::min(b.minY, y); b.maxY = std::max(b.maxY, y);
                }
        }
        if (!b.empty) { b.minX -= step; b.minY -= step; b.maxX += step; b.maxY += step; }
        return b;
    }

    void LensSystem::prepare(Real filmWidth, Real filmHeight) {
        m_filmW = filmWidth;
        m_filmH = filmHeight;
        m_filmRadius = Real(0.5) * std::hypot(filmWidth, filmHeight);
        const int bands = 64;
        m_pupil.assign(bands, Bounds());
#pragma omp parallel for schedule(dynamic, 1)
        for (int k = 0; k < bands; ++k)
            m_pupil[static_cast<size_t>(k)] = pupilBounds(m_filmRadius * Real(k) / Real(bands), m_filmRadius * Real(k + 1) / Real(bands));
        // The weight's normalisation: the irradiance factor at the film centre, from a 512 x 512 grid over its box
        m_centreIrradiance = Real(1);
        const Bounds& b = m_pupil[0];
        if (b.empty) return;
        const Real z = m_surfaces.back().z - m_filmZ;
        const int grid = 512;
        double sum = 0;
#pragma omp parallel for reduction(+ : sum) schedule(static)
        for (int gy = 0; gy < grid; ++gy)
            for (int gx = 0; gx < grid; ++gx) {
                const Real x = b.minX + (b.maxX - b.minX) * (Real(gx) + Real(0.5)) / Real(grid);
                const Real y = b.minY + (b.maxY - b.minY) * (Real(gy) + Real(0.5)) / Real(grid);
                Point3 o(0, 0, m_filmZ);
                Vector3 d = glm::normalize(Vector3(x, y, z));
                const Real cos4 = d.z * d.z * d.z * d.z;
                if (traceFromFilm(o, d, kLineD)) sum += double(cos4);
            }
        const Real area = (b.maxX - b.minX) * (b.maxY - b.minY);
        const Real e = Real(sum / (double(grid) * grid)) * area / (z * z);
        if (e > Real(0)) m_centreIrradiance = e;
    }

    bool LensSystem::generateRay(Real x, Real y, const Point2& u, Real nm, Point3& o, Vector3& d, Real& weight) const {
        weight = 0;
        const Real r = std::hypot(x, y);
        if (m_pupil.empty() || r >= m_filmRadius) return false;
        const Bounds& b = m_pupil[std::min(static_cast<size_t>(r / m_filmRadius * Real(m_pupil.size())), m_pupil.size() - 1)];
        if (b.empty) return false;
        // a point of the box, turned to the film point's angle (the box is for a film point on +x)
        const Real px = b.minX + (b.maxX - b.minX) * u.x, py = b.minY + (b.maxY - b.minY) * u.y;
        const Real c = r > Real(0) ? x / r : Real(1), s = r > Real(0) ? y / r : Real(0);
        const Surface& rear = m_surfaces.back();
        const Real z = rear.z - m_filmZ;
        o = Point3(x, y, m_filmZ);
        d = glm::normalize(Vector3(px * c - py * s - x, px * s + py * c - y, z));
        const Real cos4 = d.z * d.z * d.z * d.z;
        if (!traceFromFilm(o, d, nm)) return false;
        weight = cos4 * (b.maxX - b.minX) * (b.maxY - b.minY) / (z * z) / m_centreIrradiance;
        return true;
    }

    bool LensSystem::centralRay(Real x, Real y, Point3& o, Vector3& d) const {
        const Real r = std::hypot(x, y);
        if (m_pupil.empty() || r >= m_filmRadius) return false;
        // the middle of the box first, then a 7 x 7 grid over it, nearest to the middle first
        for (int ring = 0; ring <= 3; ++ring)
            for (int a = -ring; a <= ring; ++a)
                for (int b = -ring; b <= ring; ++b) {
                    if (std::max(std::abs(a), std::abs(b)) != ring) continue;
                    Real w;
                    if (generateRay(x, y, Point2(Real(0.5) + Real(a) / Real(7), Real(0.5) + Real(b) / Real(7)), kLineD, o, d, w)) return true;
                }
        return false;
    }

} // namespace rayt
