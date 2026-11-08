// Procedural triangle meshes (see Geometry/ProceduralMesh.hpp).
#include "pch.h"
#include "Geometry/ProceduralMesh.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt {

TriangleMesh MakeLatheMesh(const std::vector<Point2>& profile, int segments) {
    TriangleMesh m;
    for (size_t s = 0; s + 1 < profile.size(); ++s) {
        const Point2 a = profile[s], b = profile[s + 1];
        const Real dr = b.x - a.x, dy = b.y - a.y, len = std::hypot(dr, dy);
        if (len <= Real(0)) continue;
        const Real nr = dy / len, ny = -dr / len;
        const int base = static_cast<int>(m.p.size());
        int side = 0;
        for (const Point2& q : { a, b }) {
            for (int j = 0; j < segments; ++j) {
                const Real phi = Real(2) * constants::PI * Real(j) / Real(segments);
                m.p.push_back(Point3(q.x * std::cos(phi), q.y, q.x * std::sin(phi)));
                m.n.push_back(Vector3(nr * std::cos(phi), ny, nr * std::sin(phi)));
                // u around the axis, v along the profile: the tangent of anisotropic materials then runs around the axis (like turned metal)
                m.uv.push_back(Point2(Real(j) / Real(segments), (Real(s) + Real(side)) / Real(profile.size() - 1)));
            }
            ++side;
        }
        auto addTriangle = [&](int i0, int i1, int i2) {
            const Vector3 c = glm::cross(Vector3(m.p[i1] - m.p[i0]), Vector3(m.p[i2] - m.p[i0]));
            if (glm::length(c) < Real(1e-14)) return;                 // degenerate: two vertices on the axis
            if (glm::dot(c, m.n[i0]) < Real(0)) std::swap(i1, i2);    // counter-clockwise seen from outside
            m.indices.insert(m.indices.end(), { i0, i1, i2 });
        };
        for (int j = 0; j < segments; ++j) {
            const int j1 = (j + 1) % segments;
            addTriangle(base + j, base + j1, base + segments + j1);
            addTriangle(base + j, base + segments + j1, base + segments + j);
        }
    }
    return m;
}

TriangleMesh MakeFacetedLathe(const std::vector<Point2>& profile, int sides) {
    TriangleMesh m;
    for (size_t s = 0; s + 1 < profile.size(); ++s) {
        const Point2 a = profile[s], b = profile[s + 1];
        const Real dr = b.x - a.x, dy = b.y - a.y, len = std::hypot(dr, dy);
        if (len <= Real(0)) continue;
        const Real nr = dy / len, ny = -dr / len;   // outward normal of the profile segment
        for (int j = 0; j < sides; ++j) {
            const Real t0 = Real(2) * constants::PI * Real(j) / Real(sides), t1 = Real(2) * constants::PI * Real(j + 1) / Real(sides);
            const Real tm = Real(0.5) * (t0 + t1);
            const Vector3 ref(nr * std::cos(tm), ny, nr * std::sin(tm));
            auto P = [](const Point2& q, Real t) { return Point3(q.x * std::cos(t), q.y, q.x * std::sin(t)); };
            const Point3 p00 = P(a, t0), p01 = P(a, t1), p10 = P(b, t0), p11 = P(b, t1);
            auto addTriangle = [&](const Point3& q0, Point3 q1, Point3 q2) {
                Vector3 cr = glm::cross(Vector3(q1 - q0), Vector3(q2 - q0));
                if (glm::length(cr) < Real(1e-14)) return;                // degenerate (a vertex pair on the axis)
                if (glm::dot(cr, ref) < Real(0)) { std::swap(q1, q2); cr = -cr; }   // counter-clockwise seen from outside
                const Vector3 n = glm::normalize(cr);
                const int base = static_cast<int>(m.p.size());
                m.p.push_back(q0); m.p.push_back(q1); m.p.push_back(q2);
                m.n.push_back(n); m.n.push_back(n); m.n.push_back(n);
                m.indices.insert(m.indices.end(), { base, base + 1, base + 2 });
            };
            addTriangle(p00, p01, p11);
            addTriangle(p00, p11, p10);
        }
    }
    return m;
}

TriangleMesh MakePoolSurface(Real w, Real d, Real steep, unsigned seed, Real spacing) {
    struct Wave { Real dx, dz, k, A, phase; };
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::vector<Wave> waves;
    for (int i = 0; i < 14; ++i) {
        Wave v;
        const Real lambda = std::exp(std::log(Real(0.12)) + (std::log(Real(1.0)) - std::log(Real(0.12))) * Real(U(rng)));
        const Real ang = Real(2) * constants::PI * Real(U(rng));
        v.dx = std::cos(ang); v.dz = std::sin(ang);
        v.k = Real(2) * constants::PI / lambda;
        v.A = steep / v.k;
        v.phase = Real(2) * constants::PI * Real(U(rng));
        waves.push_back(v);
    }
    const int nx = std::max(1, static_cast<int>(std::ceil(w / spacing))), nz = std::max(1, static_cast<int>(std::ceil(d / spacing)));
    TriangleMesh m;
    for (int j = 0; j <= nz; ++j) {
        for (int i = 0; i <= nx; ++i) {
            const Real x = -w / 2 + w * Real(i) / Real(nx), z = -d / 2 + d * Real(j) / Real(nz);
            Real y = 0, nxs = 0, nzs = 0;
            for (const Wave& v : waves) {
                const Real ph = v.k * (v.dx * x + v.dz * z) + v.phase;
                y += v.A * std::cos(ph);
                nxs += v.k * v.A * v.dx * std::sin(ph);   // -dy/dx
                nzs += v.k * v.A * v.dz * std::sin(ph);   // -dy/dz
            }
            m.p.push_back(Point3(x, y, z));
            m.n.push_back(glm::normalize(Vector3(nxs, 1, nzs)));
        }
    }
    auto at = [&](int i, int j) { return j * (nx + 1) + i; };
    auto addTriangle = [&](int a, int b, int c) {
        const Vector3 cr = glm::cross(Vector3(m.p[b] - m.p[a]), Vector3(m.p[c] - m.p[a]));
        if (glm::dot(cr, m.n[a]) < Real(0)) std::swap(b, c);   // counter-clockwise seen from above
        m.indices.insert(m.indices.end(), { a, b, c });
    };
    for (int j = 0; j < nz; ++j)
        for (int i = 0; i < nx; ++i) {
            addTriangle(at(i, j), at(i + 1, j), at(i + 1, j + 1));
            addTriangle(at(i, j), at(i + 1, j + 1), at(i, j + 1));
        }
    return m;
}

TriangleMesh MakeSeaMesh(const SeaParams& P) {
    struct Wave { Real dx, dz, k, A, phase, lambda; };
    std::mt19937 rng(P.seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::normal_distribution<double> N(0.0, 1.0);
    const Real Lwind = Real(5);   // V^2 / g for V = 7 m/s
    std::vector<Wave> waves;
    for (int octave = 0; octave < 8; ++octave) {
        for (int j = 0; j < 4; ++j) {
            Wave w;
            w.lambda = Real(0.5) * std::pow(Real(2), Real(octave) + Real(U(rng)));
            w.k = Real(2) * constants::PI / w.lambda;
            const Real spread = std::clamp(Real(N(rng)) * glm::radians(Real(35)), -glm::radians(Real(80)), glm::radians(Real(80)));
            const Real ang = glm::radians(P.windDeg) + spread;
            w.dx = std::cos(ang); w.dz = std::sin(ang);
            const Real kL = w.k * Lwind;
            w.A = P.steep * std::exp(Real(-0.5) / (kL * kL)) / w.k;
            w.phase = Real(2) * constants::PI * Real(U(rng));
            waves.push_back(w);
        }
    }

    const Real Q = Real(0.5);   // Gerstner steepness: the sum of Q k A stays below 1, so crests do not loop
    const Real growth = std::pow(P.rOuter / P.rInner, Real(1) / Real(P.rings - 1));
    TriangleMesh m;
    auto addVertex = [&](Real r, Real theta) {
        const Real x0 = r * std::cos(theta), z0 = r * std::sin(theta);
        const Real spacing = std::max(r, P.rInner) * (growth - Real(1));   // radial distance to the next ring here
        Real dx = 0, dy = 0, dz = 0, nx = 0, ny = 1, nz = 0;
        for (const Wave& w : waves) {
            const Real fade = std::clamp((w.lambda / spacing - Real(2)) / Real(3), Real(0), Real(1));
            if (fade <= Real(0)) continue;
            const Real ph = w.k * (w.dx * x0 + w.dz * z0) + w.phase;
            const Real c = std::cos(ph), s = std::sin(ph), a = fade * w.A;
            dx -= Q * a * w.dx * s;  dz -= Q * a * w.dz * s;  dy += a * c;
            nx -= w.k * a * w.dx * c;  nz -= w.k * a * w.dz * c;  ny -= Q * w.k * a * s;
        }
        m.p.push_back(Point3(x0 + dx, dy, z0 + dz));
        m.n.push_back(glm::normalize(Vector3(nx, ny, nz)));
    };
    auto addTriangle = [&](int a, int b, int c) {
        const Vector3 cr = glm::cross(Vector3(m.p[b] - m.p[a]), Vector3(m.p[c] - m.p[a]));
        if (glm::length(cr) < Real(1e-14)) return;
        if (glm::dot(cr, m.n[a]) < Real(0)) std::swap(b, c);   // counter-clockwise seen from above
        m.indices.insert(m.indices.end(), { a, b, c });
    };

    addVertex(Real(0), Real(0));   // center
    for (int i = 0; i < P.rings; ++i)
        for (int j = 0; j < P.sectors; ++j)
            addVertex(P.rInner * std::pow(growth, Real(i)), Real(2) * constants::PI * Real(j) / Real(P.sectors));
    auto ring = [&](int i, int j) { return 1 + i * P.sectors + (j % P.sectors); };
    for (int j = 0; j < P.sectors; ++j) addTriangle(0, ring(0, j), ring(0, j + 1));
    for (int i = 0; i + 1 < P.rings; ++i)
        for (int j = 0; j < P.sectors; ++j) {
            addTriangle(ring(i, j), ring(i, j + 1), ring(i + 1, j + 1));
            addTriangle(ring(i, j), ring(i + 1, j + 1), ring(i + 1, j));
        }
    return m;
}

TriangleMesh MakePrismMesh(const std::vector<Point2>& xz, Real height) {
    TriangleMesh m;
    Point2 c(0, 0);
    for (const Point2& q : xz) c += q;
    c /= Real(xz.size());
    auto addTriangle = [&](const Point3& a, Point3 b, Point3 d, const Vector3& outward) {
        Vector3 cr = glm::cross(Vector3(b - a), Vector3(d - a));
        if (glm::length(cr) < Real(1e-14)) return;
        if (glm::dot(cr, outward) < Real(0)) { std::swap(b, d); cr = -cr; }
        const Vector3 n = glm::normalize(cr);
        const int base = static_cast<int>(m.p.size());
        m.p.push_back(a); m.p.push_back(b); m.p.push_back(d);
        m.n.push_back(n); m.n.push_back(n); m.n.push_back(n);
        m.indices.insert(m.indices.end(), { base, base + 1, base + 2 });
    };
    for (size_t i = 0; i < xz.size(); ++i) {
        const Point2 q0 = xz[i], q1 = xz[(i + 1) % xz.size()];
        const Point3 b0(q0.x, 0, q0.y), b1(q1.x, 0, q1.y), t0(q0.x, height, q0.y), t1(q1.x, height, q1.y);
        const Point2 mid = Real(0.5) * (q0 + q1) - c;
        const Vector3 out(mid.x, 0, mid.y);
        addTriangle(b0, b1, t1, out);
        addTriangle(b0, t1, t0, out);
        addTriangle(Point3(c.x, 0, c.y), b0, b1, Vector3(0, -1, 0));
        addTriangle(Point3(c.x, height, c.y), t0, t1, Vector3(0, 1, 0));
    }
    return m;
}

} // namespace rayt
