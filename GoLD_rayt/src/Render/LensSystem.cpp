// A real lens traced ray by ray (--lens); see Render/LensSystem.hpp.
#include "pch.h"
#include "Render/LensSystem.hpp"

#include "Materials/ThinFilm.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdlib>
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

        constexpr Real kAutoAperture = Real(1);   // a clear semi-diameter of 1 m or more: "auto", set by setAutoApertures

        bool ParseNumber(const std::string& token, double& value) {
            if (token == "inf" || token == "infinity" || token == "Infinity") { value = HUGE_VAL; return true; }
            char* end = nullptr;
            value = std::strtod(token.c_str(), &end);
            return end && *end == '\0' && end != token.c_str();
        }

    } // namespace

    bool LensSystem::Surface::sag(Real r, Real& z, Real& slope) const {
        // in millimetres, the units of the coefficients: z = c r^2 / (1 + s) + Sum A_i r^i, dz/dr = c r / s + Sum i A_i r^(i-1)
        const Real rm = r * Real(1e3);
        const Real c = radius == Real(0) ? Real(0) : Real(1) / (radius * Real(1e3));
        const Real arg = Real(1) - (Real(1) + conic) * c * c * rm * rm;
        if (!(arg > Real(1e-12))) return false;
        const Real s = std::sqrt(arg);
        Real zm = c * rm * rm / (Real(1) + s), dz = c * rm / s;
        Real rp = 1;   // r^(i-1), then r^i
        for (size_t i = 1; i < coef.size(); ++i) {
            dz += Real(i) * coef[i] * rp;
            rp *= rm;
            zm += coef[i] * rp;
        }
        z = zm * Real(1e-3);
        slope = dz;
        return true;
    }

    std::shared_ptr<LensSystem> LensSystem::load(const std::string& path, std::string& error) {
        std::ifstream in(path);
        if (!in) { error = "cannot open " + path; return nullptr; }
        std::vector<Surface> surfaces;
        std::vector<bool> autoDiameter;
        Real z = 0;
        double fNumberWanted = 0, imageHeightMm = 0, halfFieldDeg = 0;
        int lineNo = 0;
        for (std::string line; std::getline(in, line);) {
            ++lineNo;
            if (const size_t hash = line.find('#'); hash != std::string::npos) line.erase(hash);
            std::istringstream ss(line);
            std::vector<std::string> tok;
            for (std::string t; ss >> t;) tok.push_back(t);
            if (tok.empty()) continue;
            const std::string where = path + ", line " + std::to_string(lineNo) + ": ";
            double first = 0;
            if (!ParseNumber(tok[0], first)) {   // a keyword line
                const std::string& key = tok[0];
                if (key == "asphere") {
                    if (surfaces.empty() || surfaces.back().isStop()) { error = where + "\"asphere\" must follow the line of a lens surface"; return nullptr; }
                    if (tok.size() < 3 || tok.size() % 2 == 0) { error = where + "expected asphere followed by pairs: k K A4 a4 A6 a6 ..."; return nullptr; }
                    Surface& s = surfaces.back();
                    for (size_t i = 1; i + 1 < tok.size(); i += 2) {
                        double value = 0;
                        if (!ParseNumber(tok[i + 1], value) || !std::isfinite(value)) { error = where + "not a number: " + tok[i + 1]; return nullptr; }
                        if (tok[i] == "k") { s.conic = Real(value); continue; }
                        int power = 0;
                        if (tok[i].size() >= 2 && (tok[i][0] == 'A' || tok[i][0] == 'a')) power = std::atoi(tok[i].c_str() + 1);
                        if (power < 2 || power > 60) { error = where + "unknown coefficient " + tok[i] + " (k, or A2 ... A60)"; return nullptr; }
                        if (value != 0) s.poly.emplace_back(power, Real(value));
                    }
                } else if (key == "fnumber" || key == "image-height" || key == "half-field") {
                    double value = 0;
                    if (tok.size() != 2 || !ParseNumber(tok[1], value) || !(value > 0) || !std::isfinite(value)) { error = where + "expected " + key + " and a number above 0"; return nullptr; }
                    if (key == "half-field" && !(value < 89)) { error = where + "half-field is in degrees, below 89"; return nullptr; }
                    (key == "fnumber" ? fNumberWanted : key == "image-height" ? imageHeightMm : halfFieldDeg) = value;
                } else { error = where + "unknown keyword " + key + " (asphere, fnumber, image-height, half-field)"; return nullptr; }
                continue;
            }
            if (tok.size() < 4 || tok.size() > 5) { error = where + "expected radius, thickness, n_d, diameter [, Abbe number]"; return nullptr; }
            double v[5] = { 0, 0, 0, 0, 0 };
            const bool isAuto = tok[3] == "auto";
            for (size_t i = 0; i < tok.size(); ++i) {
                if (i == 3 && isAuto) continue;
                if (!ParseNumber(tok[i], v[i]) || (i > 0 && !std::isfinite(v[i]))) { error = where + "not a number: " + tok[i]; return nullptr; }
            }
            Surface s;
            s.flat = std::isinf(v[0]);
            s.radius = s.flat ? Real(0) : Real(v[0] * 1e-3);
            s.z = z;
            s.nd = v[2] > 1 ? Real(v[2]) : Real(1);
            s.apertureRadius = isAuto ? kAutoAperture : Real(v[3] * 0.5e-3);
            s.abbe = tok.size() == 5 ? Real(v[4]) : Real(0);
            if (!(s.apertureRadius > 0)) { error = where + "the diameter must be above 0"; return nullptr; }
            if (s.isStop() && s.nd > 1) { error = where + "a radius of 0 is the stop, in air (a flat glass surface has the radius inf)"; return nullptr; }
            surfaces.push_back(s);
            autoDiameter.push_back(isAuto);
            z -= Real(v[1] * 1e-3);
        }
        if (surfaces.empty()) { error = path + ": no surfaces"; return nullptr; }
        // A sphere smaller than its clear diameter cannot be a lens surface (an asphere's base radius can be)
        for (size_t i = 0; i < surfaces.size(); ++i) {
            const Surface& s = surfaces[i];
            if (!autoDiameter[i] && s.radius != 0 && !s.aspheric() && std::fabs(s.radius) < s.apertureRadius) {
                error = path + ", surface " + std::to_string(i + 1) + ": the radius is smaller than the semi-diameter";
                return nullptr;
            }
        }
        const bool anyAuto = std::find(autoDiameter.begin(), autoDiameter.end(), true) != autoDiameter.end();
        const auto stopIt = std::find_if(surfaces.begin(), surfaces.end(), [](const Surface& s) { return s.isStop(); });
        if (imageHeightMm > 0 && halfFieldDeg > 0) { error = path + ": give image-height or half-field, not both"; return nullptr; }
        if (anyAuto && !(imageHeightMm > 0) && !(halfFieldDeg > 0)) { error = path + ": \"auto\" clear diameters need the image height (image-height H or half-field DEG)"; return nullptr; }
        if (stopIt != surfaces.end() && stopIt->apertureRadius >= kAutoAperture && !(fNumberWanted > 0)) {
            error = path + ": an \"auto\" stop diameter needs the f-number (fnumber N)";
            return nullptr;
        }
        auto lens = std::make_shared<LensSystem>(std::move(surfaces));
        if (!(lens->focalLength() > 0)) { error = path + ": not a converging lens (no positive focal length at 587.6 nm)"; return nullptr; }
        if (fNumberWanted > 0) {   // the stop's radius for that f-number (the f-number falls as the stop opens)
            auto stop = std::find_if(lens->m_surfaces.begin(), lens->m_surfaces.end(), [](const Surface& s) { return s.isStop(); });
            if (stop == lens->m_surfaces.end()) { error = path + ": fnumber needs a stop (a surface of radius 0)"; return nullptr; }
            Real lo = Real(1e-7), hi = Real(0.5);
            stop->apertureRadius = hi;
            if (lens->fNumber() > Real(fNumberWanted)) { error = path + ": the lens cannot open to f/" + std::to_string(fNumberWanted); return nullptr; }
            for (int k = 0; k < 80; ++k) {
                stop->apertureRadius = Real(0.5) * (lo + hi);
                (lens->fNumber() > Real(fNumberWanted) ? lo : hi) = stop->apertureRadius;
            }
            stop->apertureRadius = Real(0.5) * (lo + hi);
        }
        if (halfFieldDeg > 0) {   // the image height of that field: where its chief ray meets the focal plane
            const Real h = lens->imageHeightForFieldAngle(Real(halfFieldDeg) * kPi / Real(180));
            if (!(h > 0)) { error = path + ": no chief ray of the half field " + std::to_string(halfFieldDeg) + " deg reaches the focal plane (the stop must be the first surface)"; return nullptr; }
            imageHeightMm = double(h) * 1e3;
        }
        if (imageHeightMm > 0) {
            if (anyAuto && !lens->setAutoApertures(Real(imageHeightMm * 1e-3), error)) { error = path + ": " + error; return nullptr; }
            lens->m_imageHeight = Real(imageHeightMm * 1e-3);
        }
        lens->m_wideOpen = lens->fNumber();
        return lens;
    }

    LensSystem::LensSystem(std::vector<Surface> surfaces) : m_surfaces(std::move(surfaces)) {
        for (Surface& s : m_surfaces) {
            SetCauchy(s);
            s.coef.clear();
            for (const auto& [power, a] : s.poly) {
                if (power < 1) continue;
                if (s.coef.size() <= static_cast<size_t>(power)) s.coef.resize(static_cast<size_t>(power) + 1, Real(0));
                s.coef[static_cast<size_t>(power)] += a;
            }
        }
        for (size_t i = 0; i < m_surfaces.size(); ++i)
            if (m_surfaces[i].isStop()) { m_stopIndex = static_cast<int>(i); break; }
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

    bool LensSystem::kernelParams(kernel::CameraParams<Real>& c, std::vector<kernel::LensSurface<Real>>& surfaces, std::vector<Real>& pupil,
        std::vector<Real>& bladeNormals, std::string& why) const {
        if (dispersive()) { why = "a real lens of dispersive glass"; return false; }
        if (m_pupil.empty()) { why = "a real lens that is not prepared"; return false; }
        surfaces.clear();
        for (size_t i = 0; i < m_surfaces.size(); ++i) {
            const Surface& s = m_surfaces[i];
            if (s.aspheric()) { why = "a real lens with an aspheric surface"; return false; }
            kernel::LensSurface<Real> k{};
            k.radius = s.radius;
            k.z = s.z;
            k.apertureRadius = s.apertureRadius;
            k.nBehind = indexAfter(static_cast<int>(i), kLineD);
            k.kind = s.isStop() ? kernel::kLensStop : (s.flat ? kernel::kLensFlat : kernel::kLensSphere);
            surfaces.push_back(k);
        }
        pupil.clear();
        for (const Bounds& b : m_pupil) pupil.insert(pupil.end(), { b.minX, b.minY, b.maxX, b.maxY, b.empty ? Real(1) : Real(0) });
        bladeNormals.clear();
        const bool polygon = m_blades >= 3 && !m_bladeNormals.empty();   // insideStop's test of the lens's own stop
        if (polygon)
            for (const Point2& n : m_bladeNormals) { bladeNormals.push_back(n.x); bladeNormals.push_back(n.y); }
        c.numSurfaces = static_cast<uint32_t>(m_surfaces.size());
        c.numBands = static_cast<uint32_t>(m_pupil.size());
        c.stopBlades = polygon ? static_cast<uint32_t>(m_blades) : 0u;
        c.filmZ = m_filmZ;
        c.filmW = m_filmW;
        c.filmH = m_filmH;
        c.filmRadius = m_filmRadius;
        c.centreIrradiance = m_centreIrradiance;
        c.frontPrincipalZ = m_frontPrincipal;
        c.stopScale = m_stopScale;
        c.apothemPerRadius = m_apothemPerRadius;
        return true;
    }

    void LensSystem::setStopShape(int blades, Real rotation) {
        m_blades = blades >= 3 ? blades : 0;
        m_bladeRotation = rotation;
        // the polygon's edge normals and its apothem per unit of the circle's radius, for insideStop
        m_bladeNormals.clear();
        if (m_blades < 3) return;
        const Real n = Real(m_blades);
        m_apothemPerRadius = std::sqrt(kPi / (Real(0.5) * n * std::sin(Real(2) * kPi / n))) * std::cos(kPi / n);
        for (int k = 0; k < m_blades; ++k) {
            const Real a = rotation + Real(2) * kPi * (Real(k) + Real(0.5)) / n;
            m_bladeNormals.emplace_back(std::cos(a), std::sin(a));
        }
    }

    bool LensSystem::insideStop(const Surface& s, Real x, Real y, const StopShape& stop) const {
        const Real r = s.apertureRadius * stop.scale;
        if (stop.blades < 3) return x * x + y * y <= r * r;
        if (stop.blades == m_blades && stop.rotation == m_bladeRotation && !m_bladeNormals.empty()) {   // the lens's own stop: precomputed
            const Real apothem = r * m_apothemPerRadius;
            for (const Point2& nrm : m_bladeNormals)
                if (x * nrm.x + y * nrm.y > apothem) return false;
            return true;
        }
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

    bool LensSystem::intersectAspheric(const Surface& s, const Point3& o, const Vector3& d, Real& t) const {
        // G(t) = z(t) - z_vertex + sag(rho(t)) = 0 (the surface lies at z_vertex - sag); Newton's method from the vertex plane
        if (std::fabs(d.z) < Real(1e-12)) return false;
        t = (s.z - o.z) / d.z;
        for (int it = 0; it < 50; ++it) {
            const Point3 p = o + d * t;
            const Real rho = std::sqrt(p.x * p.x + p.y * p.y);
            Real sag = 0, slope = 0;
            if (!s.sag(rho, sag, slope)) return false;
            const Real g = p.z - s.z + sag;
            const Real dg = d.z + (rho > Real(0) ? slope * (p.x * d.x + p.y * d.y) / rho : Real(0));
            if (std::fabs(dg) < Real(1e-12)) return false;
            const Real step = g / dg;
            t -= step;
            if (std::fabs(step) < Real(1e-14)) return true;
        }
        return false;
    }

    bool LensSystem::hitSurface(const Surface& s, Point3& o, const Vector3& d, bool backward, bool clip, Vector3& normal, Real& rho) const {
        const Real tMin = backward ? -HUGE_VAL : Real(1e-12);
        Real t = 0;
        if (s.aspheric() || s.flat) {
            if (s.aspheric()) { if (!intersectAspheric(s, o, d, t)) return false; }
            else {
                if (std::fabs(d.z) < Real(1e-12)) return false;
                t = (s.z - o.z) / d.z;
            }
            if (t <= tMin) return false;
            const Point3 p = o + d * t;
            rho = std::sqrt(p.x * p.x + p.y * p.y);
            if (clip && rho > s.apertureRadius) return false;
            Real sag = 0, slope = 0;
            if (!s.sag(rho, sag, slope)) return false;
            normal = rho > Real(0) ? glm::normalize(Vector3(slope * p.x / rho, slope * p.y / rho, Real(1))) : Vector3(0, 0, 1);
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
        t = ((d.z > Real(0)) == (s.radius > Real(0))) ? -b + root : -b - root;
        if (t <= tMin) return false;
        const Point3 p = o + d * t;
        rho = std::sqrt(p.x * p.x + p.y * p.y);
        if (clip && rho > s.apertureRadius) return false;
        normal = glm::normalize(p - centre);
        o = p;
        return true;
    }

    bool LensSystem::crossSurface(int i, Point3& o, Vector3& d, Real nFrom, Real nTo, const StopShape& stop, bool backward) const {
        const Surface& s = m_surfaces[static_cast<size_t>(i)];
        if (s.isStop()) {   // the aperture stop: a plane, no refraction
            if (std::fabs(d.z) < Real(1e-12)) return false;
            const Real t = (s.z - o.z) / d.z;
            if (!backward && t <= Real(0)) return false;
            const Point3 p = o + d * t;
            if (!insideStop(s, p.x, p.y, stop)) return false;
            o = p;
            return true;
        }
        Vector3 normal;
        Real rho = 0;
        if (!hitSurface(s, o, d, backward, true, normal, rho)) return false;
        if (glm::dot(normal, d) > Real(0)) normal = -normal;
        return Refract(d, normal, nFrom / nTo);
    }

    bool LensSystem::insideStopNow(Real x, Real y) const {
        if (m_stopIndex < 0) return true;
        return insideStop(m_surfaces[static_cast<size_t>(m_stopIndex)], x, y, StopShape{ m_stopScale, m_blades, m_bladeRotation });
    }

    Real LensSystem::interfaceReflectance(Real cosI, Real nIn, Real nOut, Real nm, Coating coating) {
        // A quarter wave of MgF2 at 550 nm on every glass-air surface; between two glasses (a cemented surface) bare Fresnel. With no
        // film thickness the stack's formula is Fresnel's between nIn and nOut.
        constexpr Real kMgF2 = Real(1.38), kDesignNm = Real(550);
        const bool glassAir = (nIn > Real(1)) != (nOut > Real(1));
        const Real thickness = (coating == Coating::SingleLayer && glassAir) ? kDesignNm / (Real(4) * kMgF2) : Real(0);
        return thinfilm::reflectance(cosI, nIn, kMgF2, thickness, thinfilm::Complex(nOut, Real(0)), nm);
    }

    LensSystem::ReflectanceTable LensSystem::makeReflectanceTable(int bands, const Real* bandNm, Coating coating) const {
        ReflectanceTable t;
        t.bands = bands;
        const size_t n = m_surfaces.size(), stride = ReflectanceTable::kCos + 1;
        t.r.assign(n * 2 * static_cast<size_t>(bands) * stride, Real(0));
        for (size_t k = 0; k < n; ++k)
            for (int side = 0; side < 2; ++side)
                for (int b = 0; b < bands; ++b) {
                    const Real front = indexAfter(static_cast<int>(k) - 1, bandNm[b]), back = indexAfter(static_cast<int>(k), bandNm[b]);
                    const Real nIn = side == 0 ? front : back, nOut = side == 0 ? back : front;
                    Real* row = &t.r[((k * 2 + static_cast<size_t>(side)) * static_cast<size_t>(bands) + static_cast<size_t>(b)) * stride];
                    for (size_t c = 0; c < stride; ++c)
                        row[c] = interfaceReflectance(Real(c) / Real(ReflectanceTable::kCos), nIn, nOut, bandNm[b], coating);
                }
        return t;
    }

    bool LensSystem::traceGhost(int i, int j, Point3 o, Vector3 d, Real nm, Coating coating, GhostSample& g, const ReflectanceTable* table,
        Real* bandT) const {
        const int bands = table ? table->bands : 0;
        const int n = static_cast<int>(m_surfaces.size());
        if (!(0 <= i && i < j && j < n) || m_surfaces[static_cast<size_t>(i)].isStop() || m_surfaces[static_cast<size_t>(j)].isStop()) return false;
        g = GhostSample();
        g.throughput = Real(1);
        for (int b = 0; b < bands; ++b) bandT[b] = Real(1);
        int previous = -1;
        // one crossing (or reflection) of surface k
        auto step = [&](int k, bool reflect) -> bool {
            const Surface& s = m_surfaces[static_cast<size_t>(k)];
            // a negative spacing puts the surface behind the ray's origin, along its direction of travel
            const bool backward = previous >= 0
                && (d.z < Real(0) ? s.z > m_surfaces[static_cast<size_t>(previous)].z : s.z < m_surfaces[static_cast<size_t>(previous)].z);
            previous = k;
            if (s.isStop()) {
                if (std::fabs(d.z) < Real(1e-12)) return false;
                const Real t = (s.z - o.z) / d.z;
                if (!backward && t <= Real(0)) return false;
                o = o + d * t;
                if (g.stopPasses < 3) { g.stopX[g.stopPasses] = o.x; g.stopY[g.stopPasses] = o.y; }
                ++g.stopPasses;
                return true;
            }
            Vector3 normal;
            Real rho = 0;
            if (!hitSurface(s, o, d, backward, false, normal, rho)) return false;
            g.rim = std::max(g.rim, rho / s.apertureRadius);
            if (glm::dot(normal, d) > Real(0)) normal = -normal;
            const Real cosI = -glm::dot(normal, d);
            // toward -z the ray goes from the medium in front of surface k into the one behind it; toward +z the other way
            const bool fromFront = d.z < Real(0);
            const Real nIn = fromFront ? indexAfter(k - 1, nm) : indexAfter(k, nm);
            const Real nOut = fromFront ? indexAfter(k, nm) : indexAfter(k - 1, nm);
            Real r = 0;
            if (table) {
                for (int b = 0; b < bands; ++b) {
                    const Real rb = table->at(k, fromFront, b, cosI);
                    bandT[b] *= reflect ? rb : Real(1) - rb;
                }
            }
            else r = interfaceReflectance(cosI, nIn, nOut, nm, coating);
            if (reflect) {
                d = glm::normalize(d + Real(2) * cosI * normal);
                g.throughput *= r;
                return true;
            }
            if (!Refract(d, normal, nIn / nOut)) return false;
            g.throughput *= Real(1) - r;
            return true;
        };
        for (int k = 0; k <= j; ++k) if (!step(k, k == j)) return false;
        for (int k = j - 1; k >= i; --k) if (!step(k, k == i)) return false;
        for (int k = i + 1; k < n; ++k) if (!step(k, false)) return false;
        if (!(d.z < Real(0))) return false;
        const Real t = (m_filmZ - o.z) / d.z;
        g.x = o.x + t * d.x;
        g.y = o.y + t * d.y;
        if (table && bands > 0) g.throughput = bandT[0];   // with a table, the first band's
        return g.stopPasses <= 3;
    }

    bool LensSystem::traceFromFilm(Point3& o, Vector3& d, Real nm) const {
        return traceFromFilmWith(o, d, nm, StopShape{ m_stopScale, m_blades, m_bladeRotation });
    }

    bool LensSystem::traceFromFilmWith(Point3& o, Vector3& d, Real nm, const StopShape& stop) const {
        const int n = static_cast<int>(m_surfaces.size());
        for (int i = n - 1; i >= 0; --i) {
            // a negative spacing puts surface i behind surface i + 1, which the ray has just left
            const bool backward = i + 1 < n && m_surfaces[static_cast<size_t>(i)].z < m_surfaces[static_cast<size_t>(i + 1)].z;
            if (!crossSurface(i, o, d, indexAfter(i, nm), indexAfter(i - 1, nm), stop, backward)) return false;
        }
        return true;
    }

    bool LensSystem::traceFromObject(Point3& o, Vector3& d, Real nm) const {
        const StopShape stop{ m_stopScale, m_blades, m_bladeRotation };
        for (int i = 0; i < static_cast<int>(m_surfaces.size()); ++i) {
            const bool backward = i > 0 && m_surfaces[static_cast<size_t>(i)].z > m_surfaces[static_cast<size_t>(i - 1)].z;
            if (!crossSurface(i, o, d, indexAfter(i - 1, nm), indexAfter(i, nm), stop, backward)) return false;
        }
        return true;
    }

    bool LensSystem::chiefRayHeight(Real theta, Real& y) const {
        // the chief ray of the field angle theta, coming from -x through the centre of the first-surface stop: its height on the focal plane
        if (m_surfaces.empty() || !m_surfaces.front().isStop()) return false;
        const Vector3 d0(std::sin(theta), 0, -std::cos(theta));
        Point3 o = Point3(0, 0, m_surfaces.front().z) - d0 * Real(0.01);
        Vector3 d = d0;
        if (!traceFromObject(o, d, kLineD) || d.z >= Real(0)) return false;
        y = o.x + (m_rearFocal - o.z) / d.z * d.x;
        return true;
    }

    Real LensSystem::imageHeightForFieldAngle(Real theta) const {
        Real y = 0;
        return chiefRayHeight(theta, y) && y > Real(0) ? y : Real(0);
    }

    Real LensSystem::fieldAngleForImageHeight(Real h) const {
        if (m_surfaces.empty() || !m_surfaces.front().isStop()) return Real(0);
        auto landing = [&](Real theta, Real& y) { return chiefRayHeight(theta, y); };
        Real lo = 0, hi = Real(89) * kPi / Real(180);
        for (int k = 0; k < 80; ++k) {
            const Real mid = Real(0.5) * (lo + hi);
            Real y = 0;
            (landing(mid, y) && y < h ? lo : hi) = mid;
        }
        Real y = 0;
        if (!landing(lo, y) || std::fabs(y - h) > Real(1e-6) * h) return Real(0);
        return lo;
    }

    bool LensSystem::setAutoApertures(Real imageHeight, std::string& error) {
        if (m_surfaces.empty() || !m_surfaces.front().isStop()) { error = "\"auto\" clear diameters need the stop as the first surface"; return false; }
        const Real thetaMax = fieldAngleForImageHeight(imageHeight);
        if (!(thetaMax > 0)) { error = "no chief ray reaches the image height " + std::to_string(imageHeight * 1e3) + " mm"; return false; }
        const size_t n = m_surfaces.size();
        const Real z0 = m_surfaces.front().z, rs = m_surfaces.front().apertureRadius * Real(1 - 1e-9);
        const StopShape round{ Real(1), 0, Real(0) };
        std::vector<Real> highest(n, Real(0));
        // fields from the axis to the edge; points of the stop on rings, the rim (where the heights peak) finely
        const int fields = 24, rings = 24;
        for (int f = 0; f <= fields; ++f) {
            const Real theta = thetaMax * Real(f) / Real(fields);
            const Vector3 d0(std::sin(theta), 0, -std::cos(theta));
            for (int ring = 0; ring <= rings; ++ring) {
                const int spokes = ring == 0 ? 1 : ring == rings ? 720 : 64;
                for (int k = 0; k < spokes; ++k) {
                    const Real r = rs * Real(ring) / Real(rings), a = Real(2) * kPi * Real(k) / Real(spokes);
                    Point3 o = Point3(r * std::cos(a), r * std::sin(a), z0) - d0 * Real(0.01);
                    Vector3 d = d0;
                    std::vector<Real> h(n, Real(0));
                    bool through = true;
                    for (size_t i = 0; i < n && through; ++i) {
                        const bool backward = i > 0 && m_surfaces[i].z > m_surfaces[i - 1].z;
                        through = crossSurface(static_cast<int>(i), o, d, indexAfter(static_cast<int>(i) - 1, kLineD), indexAfter(static_cast<int>(i), kLineD), round, backward);
                        h[i] = std::sqrt(o.x * o.x + o.y * o.y);
                    }
                    if (!through) continue;
                    for (size_t i = 0; i < n; ++i) highest[i] = std::max(highest[i], h[i]);
                }
            }
        }
        for (size_t i = 0; i < n; ++i) {
            Surface& s = m_surfaces[i];
            if (s.apertureRadius < kAutoAperture) continue;
            if (!(highest[i] > 0)) { error = "no ray of the bundle passes surface " + std::to_string(i + 1); return false; }
            s.apertureRadius = highest[i] * Real(1.005);   // a margin for the peaks between the samples
        }
        m_imageHeight = imageHeight;
        return true;
    }

    void LensSystem::paraxial() {
        // A ray parallel to the axis a nanometre off it, from each side: where it crosses the axis is a focal point, where its extension
        // meets the height it came in at is a principal plane. (Odd aspheric terms, A3 r^3, bend a ray in proportion to its height, not
        // its square: at a micrometre they would still move the focus of a phone lens by half a micrometre.)
        const Real h = Real(1e-9);
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
        const Real h = Real(1e-9);
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
        const bool hasStop = std::any_of(m_surfaces.begin(), m_surfaces.end(), [](const Surface& s) { return s.isStop(); });
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
