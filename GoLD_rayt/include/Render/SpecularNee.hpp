/**
 * @file Render/SpecularNee.hpp
 * @brief Specular manifold next event estimation (prototype): a point light seen through a closed convex flat-faced
 *        dielectric body by two refractions (enter, leave), solved deterministically instead of estimated.
 *
 * Idea (Hanika et al. 2015, "Manifold Next Event Estimation"; Zeltner et al. 2020, "Specular Manifold Sampling").
 * A caustic from a point light is a light path  L S S D  that BSDF sampling can never find (a point light cannot be hit),
 * and that a photon map only estimates by counting photons (noise, blur). Here the path is computed exactly: for a
 * receiver point x and the light y, find the two points vA, vB on the facets of the body such that Snell's law holds at both
 * (the "generalized half vector" is parallel to the facet normal), by Newton's method.
 *
 * Constraint at a refraction vertex v with neighbors a (outside, index 1) and b (inside, index n), unit directions
 * wa = (a - v)/|a - v|, wb = (b - v)/|b - v|:   (1 * wa + n * wb) . t = 0   for both tangents t of the facet.
 * That is Snell's law (the tangential parts of 1*wa and n*wb cancel). 4 unknowns (2 plane coordinates per facet), 4 equations.
 *
 * Every pair of triangles (A: the facet the receiver sees, B: the facet the light sees) is tried, with two seeds (the straight
 * line x - y projected on the facet planes, and the triangle centers); a solution counts if it lies inside both triangles and
 * the geometry is a real refraction (receiver and light outside, the chord inside), and the three segments are unoccluded.
 *
 * Contribution. The flux of the bundle that leaves the light into d(omega) arrives on the receiver surface area dA, so
 *     E = I * T * d(omega) / dA,     L = f(wo, wi) * E,
 * with T the Fresnel transmittance product of the two interfaces. d(omega)/dA is the Jacobian of the map receiver point ->
 * direction leaving the light, computed by finite differences of two re-solved nearby paths (warm started Newton). (For a free
 * path it reduces to cos(theta) / d^2, the usual point-light irradiance.) The 1/eta^2 radiance scaling of refraction cancels
 * between entering and leaving.
 *
 * Dispersion (spectral mode, Abbe number): the path depends on the wavelength, so one of the 4 path wavelengths is chosen at
 * random, the chain is solved for its IOR and the value is written to that component times 4 (an unbiased estimate of the
 * 4-component sum).
 *
 * Scope of this prototype: point lights, ONE closed convex body made of flat triangles (outward winding), exactly two
 * refractions, no total internal reflection chains (R, TIR), no absorption. Other lights (area, environment) would need the light
 * sampled and the same solve per sample (not done). Do not combine with the caustic photon map in the same render
 * (the caustic would be counted twice).
 */

#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <memory>
#include <vector>

#include "Core/Types.hpp"
#include "Core/Math.hpp"
#include "Core/Ray.hpp"
#include "Core/Sampling.hpp"
#include "Core/Constants.hpp"
#include "Core/Fresnel.hpp"
#include "Core/Interaction.hpp"
#include "Core/SpectrumUtils.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Geometry/Frame.hpp"
#include "Geometry/TriangleMesh.hpp"
#include "Scene/Scene.hpp"
#include "Lights/PointLight.hpp"
#include "Materials/Dielectric.hpp"
#include "Materials/Material.hpp"

namespace rayt {

    class SpecularNee {
    public:
        SpecularNee(std::shared_ptr<const TriangleMesh> mesh, std::shared_ptr<const Dielectric> glass, const Scene& scene)
            : m_mesh(std::move(mesh)), m_glass(std::move(glass))
        {
            // Coplanar triangles with the same outward normal form one facet (a flat face): the solver works per facet, which keeps the
            // number of (receiver facet, light facet) pairs small.
            for (size_t i = 0; i < m_mesh->numTriangles(); ++i) {
                const Point3 q0 = m_mesh->p[static_cast<size_t>(m_mesh->indices[3 * i])];
                const Point3 q1 = m_mesh->p[static_cast<size_t>(m_mesh->indices[3 * i + 1])];
                const Point3 q2 = m_mesh->p[static_cast<size_t>(m_mesh->indices[3 * i + 2])];
                const Vector3 c = glm::cross(Vector3(q1 - q0), Vector3(q2 - q0));
                if (glm::length(c) < Real(1e-14)) continue;   // degenerate
                const Vector3 n = glm::normalize(c);          // outward (counter-clockwise winding)
                Facet* host = nullptr;
                for (Facet& F : m_facets)
                    if (glm::dot(n, F.n) > Real(1) - Real(1e-9) && std::abs(glm::dot(n, Vector3(q0 - F.p0))) < Real(1e-9)) { host = &F; break; }
                if (!host) {
                    Facet f;
                    f.p0 = q0;
                    f.n = n;
                    f.t1 = glm::normalize(Vector3(q1 - q0));
                    f.t2 = glm::cross(f.n, f.t1);
                    m_facets.push_back(f);
                    host = &m_facets.back();
                }
                auto coords = [&](const Point3& P) { return Point2(glm::dot(Vector3(P - host->p0), host->t1), glm::dot(Vector3(P - host->p0), host->t2)); };
                host->tris.push_back({ coords(q0), coords(q1), coords(q2) });
            }
            for (Facet& F : m_facets) {
                Point2 c(0, 0);
                for (const auto& t : F.tris) c += (t[0] + t[1] + t[2]) / Real(3);
                F.center = c / Real(F.tris.size());
                for (const auto& t : F.tris) for (const Point2& q : t) F.radius = std::max(F.radius, glm::length(Point2(q - F.center)));
            }
            for (const auto& l : scene.lights())
                if (const auto* pl = dynamic_cast<const PointLight*>(l.get())) m_lights.push_back(pl);
        }

        size_t numFacets() const { return m_facets.size(); }
        size_t numLights() const { return m_lights.size(); }
        long long queries() const { return m_queries.load(); }        ///< receiver points evaluated
        long long solutions() const { return m_solutions.load(); }    ///< solved paths added

        /// Where the candidate (receiver facet, light facet) pairs drop out, for diagnostics.
        struct Counters {
            std::atomic<long long> pairs{ 0 };        ///< pairs tried by the solver
            std::atomic<long long> converged{ 0 };    ///< Newton reached the tolerance
            std::atomic<long long> insideFacets{ 0 }; ///< ... and the vertices lie in their facets
            std::atomic<long long> refraction{ 0 };   ///< ... and it is a real refraction (orientation)
            std::atomic<long long> unblocked{ 0 };    ///< ... and both outside segments are free
            std::atomic<long long> focusOk{ 0 };      ///< ... and the Jacobian (focusing) is finite and positive
        };
        const Counters& counters() const { return m_c; }

        /**
         * @brief Radiance toward `wo` at the non-specular point `rec`, from every point light through the body (light-like Spectrum).
         */
        Spectrum sample(const Scene& scene, const SurfaceInteraction& rec, const Vector3& wo, const SampledWavelengths& lambda) const {
            Spectrum total(0.0);
            if (m_lights.empty() || m_facets.empty()) return total;
            ++m_queries;

            int k = -1;   // the wavelength component carried by a dispersive solve
#if defined(RAYT_SPECTRAL_MODE)
            if (m_glass->isDispersive()) k = std::min(static_cast<int>(sampling::Random() * Real(kSpectralSamples)), kSpectralSamples - 1);
#endif
            const Real n = m_glass->iorAtNm(lambda[k >= 0 ? k : 0]);
            const Point3 x = rec.p;
            const frame::Frame fr(rec.gn);
            const Vector3 tx = fr.localToWorld(Vector3(1, 0, 0)), ty = fr.localToWorld(Vector3(0, 1, 0));

            for (const PointLight* light : m_lights) {
                const Point3 y = light->position();
                const Spectrum I = light->intensity(lambda);
                for (size_t ia = 0; ia < m_facets.size(); ++ia) {
                    const Facet& A = m_facets[ia];
                    if (glm::dot(A.n, Vector3(x - A.p0)) <= Real(1e-9)) continue;           // the receiver must see the facet from outside
                    for (size_t ib = 0; ib < m_facets.size(); ++ib) {
                        if (ib == ia) continue;
                        const Facet& B = m_facets[ib];
                        if (glm::dot(B.n, Vector3(y - B.p0)) <= Real(1e-9)) continue;        // so must the light
                        Point2 a, b;
                        ++m_c.pairs;
                        if (!solvePair(x, y, A, B, n, a, b)) continue;
                        const Point3 vA = A.p0 + A.t1 * a.x + A.t2 * a.y;
                        const Point3 vB = B.p0 + B.t1 * b.x + B.t2 * b.y;

                        // visibility of the two outside segments
                        if (blocked(scene, rec, x, vA) || blockedFrom(scene, vB, y)) continue;
                        ++m_c.unblocked;

                        const Real dOmegaDA = focusing(x, tx, ty, y, A, B, n, a, b, vB);
                        if (!(dOmegaDA > Real(0)) || !std::isfinite(dOmegaDA)) continue;
                        ++m_c.focusOk;

                        const Vector3 wx = glm::normalize(Vector3(x - vA));     // at vA: toward the receiver
                        const Vector3 wi = -wx;                                 // at the receiver: the direction the light arrives from (toward vA)
                        const Vector3 wy = glm::normalize(Vector3(y - vB));     // at vB: toward the light
                        const Real T = (Real(1) - fresnel::fresnelDielectric(std::abs(glm::dot(wx, A.n)), Real(1), n))
                                     * (Real(1) - fresnel::fresnelDielectric(std::abs(glm::dot(wy, B.n)), Real(1), n));

                        const Spectrum f = rec.matPtr->eval(rec, wo, wi, lambda);
                        if (isBlack(f)) continue;
                        ++m_solutions;
                        if (k >= 0) {
#if defined(RAYT_SPECTRAL_MODE)
                            Spectrum c(0.0);
                            c[k] = Real(kSpectralSamples) * f[k] * I[k] * T * dOmegaDA;   // one wavelength, times 4: unbiased for the sum
                            total += c;
#endif
                        }
                        else total += f * (I * (T * dOmegaDA));
                    }
                }
            }
            return total;
        }

    private:
        struct Facet {
            Point3 p0;
            Vector3 n, t1, t2;                       ///< outward normal and an orthonormal tangent frame (p0 is the origin of the plane coordinates)
            std::vector<std::array<Point2, 3>> tris; ///< the facet's triangles in (t1, t2) coordinates
            Point2 center;                           ///< mean of the triangle centers
            Real radius = 0;                         ///< largest distance from the center to a corner
        };

        std::shared_ptr<const TriangleMesh> m_mesh;
        std::shared_ptr<const Dielectric> m_glass;
        std::vector<Facet> m_facets;
        std::vector<const PointLight*> m_lights;
        mutable std::atomic<long long> m_queries{ 0 };
        mutable std::atomic<long long> m_solutions{ 0 };
        mutable Counters m_c;

        using Vec4 = glm::vec<4, Real, glm::defaultp>;
        using Mat4 = glm::mat<4, 4, Real, glm::defaultp>;

        /// The four Snell conditions for the chain x - vA - vB - y with the unknown plane coordinates a (on A) and b (on B).
        static void residual(const Point3& x, const Point3& y, const Facet& A, const Facet& B, Real n,
            const Point2& a, const Point2& b, Real R[4])
        {
            const Point3 vA = A.p0 + A.t1 * a.x + A.t2 * a.y;
            const Point3 vB = B.p0 + B.t1 * b.x + B.t2 * b.y;
            const Vector3 wx = glm::normalize(Vector3(x - vA));
            const Vector3 wB = glm::normalize(Vector3(vB - vA));   // inside the body, seen from vA
            const Vector3 wy = glm::normalize(Vector3(y - vB));
            const Vector3 hA = wx + wB * n;                         // 1 * wx + n * wB
            const Vector3 hB = wB * (-n) + wy;                      // n * (back toward vA) + 1 * wy
            R[0] = glm::dot(hA, A.t1);
            R[1] = glm::dot(hA, A.t2);
            R[2] = glm::dot(hB, B.t1);
            R[3] = glm::dot(hB, B.t2);
        }

        /// Newton's method from (a, b) (numerical Jacobian). True when the residual is below 1e-13. A start near a solution converges in a few
        /// steps, so a run that has not got below 1e-4 after 6 steps (or leaves the neighborhood of the body) is given up.
        static bool newton(const Point3& x, const Point3& y, const Facet& A, const Facet& B, Real n, Point2& a, Point2& b) {
            for (int it = 0; it < 12; ++it) {
                Real R[4];
                residual(x, y, A, B, n, a, b, R);
                const Real err = std::max(std::max(std::abs(R[0]), std::abs(R[1])), std::max(std::abs(R[2]), std::abs(R[3])));
                if (!std::isfinite(err)) return false;
                if (err < Real(1e-13)) return true;
                if (it >= 6 && err > Real(1e-4)) return false;
                if (glm::length(Point2(a - A.center)) > Real(5) * A.radius + Real(0.1) || glm::length(Point2(b - B.center)) > Real(5) * B.radius + Real(0.1)) return false;

                const Real u0[4] = { a.x, a.y, b.x, b.y };
                Mat4 J(0.0);
                const Real du = Real(1e-7);
                for (int j = 0; j < 4; ++j) {
                    Real u[4] = { u0[0], u0[1], u0[2], u0[3] };
                    u[j] += du;
                    Real Rp[4];
                    residual(x, y, A, B, n, Point2(u[0], u[1]), Point2(u[2], u[3]), Rp);
                    for (int i = 0; i < 4; ++i) J[j][i] = (Rp[i] - R[i]) / du;   // glm is column-major: J[column][row]
                }
                const Vec4 step = glm::inverse(J) * Vec4(R[0], R[1], R[2], R[3]);
                if (!std::isfinite(step.x + step.y + step.z + step.w)) return false;
                const Real lim = Real(0.02);   // never move more than 2 cm in one step
                a.x -= std::clamp(step.x, -lim, lim); a.y -= std::clamp(step.y, -lim, lim);
                b.x -= std::clamp(step.z, -lim, lim); b.y -= std::clamp(step.w, -lim, lim);
            }
            return false;
        }

        /// Is the point (plane coordinates) inside one of the facet's triangles?
        static bool inside(const Facet& F, const Point2& p) {
            const Real eps = Real(1e-9);
            auto edge = [&](const Point2& u, const Point2& v) { return (v.x - u.x) * (p.y - u.y) - (v.y - u.y) * (p.x - u.x); };
            for (const auto& t : F.tris) {
                const Real e0 = edge(t[0], t[1]), e1 = edge(t[1], t[2]), e2 = edge(t[2], t[0]);
                if ((e0 >= -eps && e1 >= -eps && e2 >= -eps) || (e0 <= eps && e1 <= eps && e2 <= eps)) return true;
            }
            return false;
        }

        /// Refracted direction of the unit direction d (travelling) at a surface with unit normal nrm on the side d comes from, relative IOR eta.
        static bool refractDir(const Vector3& d, const Vector3& nrm, Real eta, Vector3& out) {
            const Real cosi = -glm::dot(d, nrm);
            const Real k = Real(1) - eta * eta * (Real(1) - cosi * cosi);
            if (k < Real(0)) return false;
            out = d * eta + nrm * (eta * cosi - std::sqrt(k));
            return true;
        }

        /// Solves the chain through facet A (receiver side) and facet B (light side); the solution must lie in both facets.
        bool solvePair(const Point3& x, const Point3& y, const Facet& A, const Facet& B, Real n, Point2& a, Point2& b) const {
            // seed 1: aim at the middle of facet A from the receiver, refract there (Snell), and see where that ray meets the plane of B;
            // seed 2: the facet centers
            Point2 seed[2][2];
            bool haveSnell = false;
            {
                const Point3 vA0 = A.p0 + A.t1 * A.center.x + A.t2 * A.center.y;
                Vector3 dIn;
                if (refractDir(glm::normalize(Vector3(vA0 - x)), A.n, Real(1) / n, dIn)) {
                    const Real den = glm::dot(B.n, dIn);
                    if (den > Real(1e-9)) {   // the ray must be leaving through B
                        const Real s = glm::dot(B.n, Vector3(B.p0 - vA0)) / den;
                        if (s > Real(0)) {
                            const Point3 P = vA0 + dIn * s;
                            seed[0][0] = A.center;
                            seed[0][1] = Point2(glm::dot(Vector3(P - B.p0), B.t1), glm::dot(Vector3(P - B.p0), B.t2));
                            haveSnell = true;
                        }
                    }
                }
            }
            seed[1][0] = A.center;
            seed[1][1] = B.center;

            for (int s = haveSnell ? 0 : 1; s < 2; ++s) {
                a = seed[s][0];
                b = seed[s][1];
                if (!newton(x, y, A, B, n, a, b)) continue;
                ++m_c.converged;
                if (!inside(A, a) || !inside(B, b)) continue;
                ++m_c.insideFacets;
                const Point3 vA = A.p0 + A.t1 * a.x + A.t2 * a.y;
                const Point3 vB = B.p0 + B.t1 * b.x + B.t2 * b.y;
                const Vector3 chord = vB - vA;
                if (glm::length(chord) < Real(1e-6)) continue;
                const Vector3 wx = glm::normalize(Vector3(x - vA)), wB = glm::normalize(chord), wy = glm::normalize(Vector3(y - vB));
                // a real refraction: both ends outside, the chord inside the body
                if (glm::dot(wx, A.n) <= Real(1e-9) || glm::dot(wB, A.n) >= Real(-1e-9)) continue;
                if (glm::dot(wy, B.n) <= Real(1e-9) || glm::dot(wB, B.n) <= Real(1e-9)) continue;
                ++m_c.refraction;
                return true;
            }
            return false;
        }

        /// d(omega)/dA: how a patch of the receiver surface maps to directions leaving the light (finite differences of two re-solved paths).
        static Real focusing(const Point3& x, const Vector3& tx, const Vector3& ty, const Point3& y, const Facet& A, const Facet& B,
            Real n, const Point2& a, const Point2& b, const Point3& vB0)
        {
            const Real d = Real(1e-6);
            Point2 a1 = a, b1 = b, a2 = a, b2 = b;
            if (!newton(x + tx * d, y, A, B, n, a1, b1)) return Real(0);
            if (!newton(x + ty * d, y, A, B, n, a2, b2)) return Real(0);
            const Point3 vB1 = B.p0 + B.t1 * b1.x + B.t2 * b1.y;
            const Point3 vB2 = B.p0 + B.t1 * b2.x + B.t2 * b2.y;
            const Vector3 w0 = glm::normalize(Vector3(vB0 - y)), w1 = glm::normalize(Vector3(vB1 - y)), w2 = glm::normalize(Vector3(vB2 - y));
            return glm::length(glm::cross(w1 - w0, w2 - w0)) / (d * d);
        }

        /// Is the segment receiver -> vA blocked (the ray starts at the receiver offset along its normal, and stops just before vA)?
        static bool blocked(const Scene& scene, const SurfaceInteraction& rec, const Point3& from, const Point3& to) {
            const Vector3 toward = glm::normalize(Vector3(to - from));
            const Point3 o = rec.OffsetRayOrigin(toward);
            const Vector3 d = to - o;
            const Real dist = glm::length(d);
            if (dist <= Real(0)) return false;
            Ray r(o, d / dist, Real(0), rec.medium);
            r.tMax = dist * (Real(1) - Real(1e-4));
            return scene.occluded(r);
        }

        static bool blockedFrom(const Scene& scene, const Point3& from, const Point3& to) {
            const Vector3 d = to - from;
            const Real dist = glm::length(d);
            if (dist <= Real(0)) return false;
            Ray r(from, d / dist, Real(1e-5), nullptr);
            r.tMax = dist * (Real(1) - Real(1e-4));
            return scene.occluded(r);
        }
    };

} // namespace rayt
