/**
 * @file Render/SpecularNee.hpp
 * @brief Specular manifold next event estimation (prototype): a point light seen through a flat-faced dielectric body by two
 *        refractions (enter, leave), solved deterministically instead of estimated.
 *
 * Idea (Hanika et al. 2015, "Manifold Next Event Estimation"; Zeltner et al. 2020, "Specular Manifold Sampling").
 * A caustic from a point light is a light path  L S S D  that BSDF sampling can never find (a point light cannot be hit),
 * and that a photon map only estimates by counting photons (noise, blur). Here the path is computed exactly: for a
 * receiver point x and the light y, find the two points vA, vB on the facets of a body such that Snell's law holds at both
 * (the "generalized half vector" is parallel to the facet normal), by Newton's method.
 *
 * Constraint at a refraction vertex v with neighbors a (outside, index 1) and b (inside, index n), unit directions
 * wa = (a - v)/|a - v|, wb = (b - v)/|b - v|:   (1 * wa + n * wb) . t = 0   for both tangents t of the facet.
 * That is Snell's law (the tangential parts of 1*wa and n*wb cancel). 4 unknowns (2 plane coordinates per facet), 4 equations.
 *
 * Every pair of facets (A: the facet the receiver sees, B: the facet the light sees) of a body is tried, with two seeds (the
 * receiver's ray refracted at the middle of A, and the facet centers); a solution counts if it lies inside both facets and the
 * geometry is a real refraction (receiver and light outside, the chord inside), the chord A - B runs inside the glass without meeting
 * any surface (which makes hollow vessels work: a chord across the cavity is rejected, the chord across a wall is accepted), and the
 * two outside segments are unoccluded.
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
 * Scope: point lights; any number of bodies, each made of flat triangles (outward winding), convex or hollow; exactly two
 * refractions through the SAME body (a wall, a prism, a block). Not handled: light that passes through two bodies or through a body and
 * a water filling (four or more refractions), reflection and total-internal-reflection chains, curved surfaces, area / environment
 * lights (they would need the light sampled and the same solve per sample). Combine it with the caustic photon map for those chains:
 * give the map the materials() of the bodies (CausticPhotonParams::skipTransmissionPairs) so that it does not store the T T chains
 * of these bodies a second time.
 */

#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <ostream>
#include <string>
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
#include "Render/SpecularBody.hpp"

namespace rayt {

    class SpecularNee {
    public:
        /// A flat-faced body (the mesh, with outward counter-clockwise triangles) and its glass (Render/SpecularBody.hpp).
        using Body = SpecularBody;

        SpecularNee(const std::vector<Body>& bodies, const Scene& scene) {
            for (const Body& b : bodies) {
                BodyData bd;
                bd.glass = b.glass;
                buildFacets(*b.mesh, bd.facets);
                if (b.glass->isDispersive()) m_anyDispersive = true;
                m_bodies.push_back(std::move(bd));
            }
            for (const auto& l : scene.lights())
                if (const auto* pl = dynamic_cast<const PointLight*>(l.get())) m_lights.push_back(pl);
        }

        /// One body (a prism).
        SpecularNee(std::shared_ptr<const TriangleMesh> mesh, std::shared_ptr<const Dielectric> glass, const Scene& scene)
            : SpecularNee(std::vector<Body>{ Body{ std::move(mesh), std::move(glass) } }, scene) {}

        size_t numFacets() const { size_t n = 0; for (const BodyData& b : m_bodies) n += b.facets.size(); return n; }
        size_t numLights() const { return m_lights.size(); }
        long long queries() const { return m_queries.load(); }        ///< receiver points evaluated
        long long solutions() const { return m_solutions.load(); }    ///< solved paths added

        /// The glass materials of the bodies (for CausticPhotonParams::skipTransmissionPairs).
        std::vector<const Material*> materials() const {
            std::vector<const Material*> m;
            for (const BodyData& b : m_bodies) m.push_back(b.glass.get());
            return m;
        }

        /// Where the candidate (receiver facet, light facet) pairs drop out, for diagnostics.
        struct Counters {
            std::atomic<long long> pairs{ 0 };        ///< pairs tried by the solver
            std::atomic<long long> converged{ 0 };    ///< Newton reached the tolerance
            std::atomic<long long> insideFacets{ 0 }; ///< ... and the vertices lie in their facets
            std::atomic<long long> refraction{ 0 };   ///< ... and it is a real refraction (orientation)
            std::atomic<long long> chordFree{ 0 };    ///< ... and the chord between the two vertices stays inside the glass
            std::atomic<long long> unblocked{ 0 };    ///< ... and both outside segments are free
            std::atomic<long long> focusOk{ 0 };      ///< ... and the Jacobian (focusing) is finite and positive
            std::atomic<long long> landTried[7]{};    ///< pairs by where the Snell seed lands on B (in B radii: <1, <2, <4, <8, <16, more, no seed)
            std::atomic<long long> landOk[7]{};       ///< ... of which gave a real refraction
        };
        const Counters& counters() const { return m_c; }

        /**
         * @brief Radiance toward `wo` at the non-specular point `rec`, from every point light through the bodies (light-like Spectrum).
         */
        Spectrum sample(const Scene& scene, const SurfaceInteraction& rec, const Vector3& wo, const SampledWavelengths& lambda) const {
            Spectrum total(0.0);
            if (m_lights.empty() || m_bodies.empty()) return total;
            ++m_queries;

            int k = -1;   // the wavelength component carried by a dispersive solve
#if defined(RAYT_SPECTRAL_MODE)
            if (m_anyDispersive) k = std::min(static_cast<int>(sampling::Random() * Real(kSpectralSamples)), kSpectralSamples - 1);
#endif
            const Point3 x = rec.p;
            const frame::Frame fr(rec.gn);
            const Vector3 tx = fr.localToWorld(Vector3(1, 0, 0)), ty = fr.localToWorld(Vector3(0, 1, 0));

            for (const PointLight* light : m_lights) {
                const Point3 y = light->position();
                const Spectrum I = light->intensity(lambda);
                for (const BodyData& body : m_bodies) {
                    const Real n = body.glass->iorAtNm(lambda[k >= 0 ? k : 0]);
                    const std::vector<Facet>& facets = body.facets;
                    for (size_t ia = 0; ia < facets.size(); ++ia) {
                        const Facet& A = facets[ia];
                        if (glm::dot(A.n, Vector3(x - A.p0)) <= Real(1e-9)) continue;           // the receiver must see the facet from outside
                        for (size_t ib = 0; ib < facets.size(); ++ib) {
                            if (ib == ia) continue;
                            const Facet& B = facets[ib];
                            if (glm::dot(B.n, Vector3(y - B.p0)) <= Real(1e-9)) continue;        // so must the light
                            Point2 a, b;
                            ++m_c.pairs;
                            if (!solvePair(x, y, A, B, n, a, b)) continue;
                            const Point3 vA = A.p0 + A.t1 * a.x + A.t2 * a.y;
                            const Point3 vB = B.p0 + B.t1 * b.x + B.t2 * b.y;

                            // the chord inside the glass must not meet a surface (a chord across a cavity does), and the two outside
                            // segments must be free
                            if (chordBlocked(scene, vA, vB)) continue;
                            ++m_c.chordFree;
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
            }
            return total;
        }

        /// Diagnostics: for the receiver x, every (receiver facet, light facet) pair that has a solution found either by the normal seeds or by a
        /// dense grid of Newton starting points (so a solution the normal seeds miss shows up as exhaustive > 0, std = 0).
        void diagnose(const Scene& scene, const SurfaceInteraction& rec, std::ostream& os) const {
            const Point3 x = rec.p;
            for (const PointLight* light : m_lights) {
                const Point3 y = light->position();
                for (const BodyData& body : m_bodies) {
                    const Real n = body.glass->iorAtNm(Real(550));
                    const auto& facets = body.facets;
                    for (size_t ia = 0; ia < facets.size(); ++ia) {
                        const Facet& A = facets[ia];
                        if (glm::dot(A.n, Vector3(x - A.p0)) <= Real(1e-9)) continue;
                        for (size_t ib = 0; ib < facets.size(); ++ib) {
                            if (ib == ia) continue;
                            const Facet& B = facets[ib];
                            if (glm::dot(B.n, Vector3(y - B.p0)) <= Real(1e-9)) continue;
                            Point2 a, b;
                            const bool stdFound = solvePair(x, y, A, B, n, a, b);
                            int found = 0;
                            Point2 fa(0, 0), fb(0, 0);
                            for (int i = -2; i <= 2; ++i) for (int j = -2; j <= 2; ++j) for (int k = -2; k <= 2; ++k) for (int l = -2; l <= 2; ++l) {
                                Point2 a1 = A.center + Point2(i, j) * (A.radius * Real(0.5));
                                Point2 b1 = B.center + Point2(k, l) * (B.radius * Real(0.5));
                                if (!newton(x, y, A, B, n, a1, b1) || !inside(A, a1) || !inside(B, b1)) continue;
                                const Point3 vA = A.p0 + A.t1 * a1.x + A.t2 * a1.y, vB = B.p0 + B.t1 * b1.x + B.t2 * b1.y;
                                const Vector3 ch = vB - vA;
                                if (glm::length(ch) < Real(1e-6)) continue;
                                const Vector3 wx = glm::normalize(Vector3(x - vA)), wB = glm::normalize(ch), wy = glm::normalize(Vector3(y - vB));
                                if (glm::dot(wx, A.n) <= Real(1e-9) || glm::dot(wB, A.n) >= Real(-1e-9)) continue;
                                if (glm::dot(wy, B.n) <= Real(1e-9) || glm::dot(wB, B.n) <= Real(1e-9)) continue;
                                ++found; fa = a1; fb = b1;
                            }
                            if (!stdFound && found == 0) continue;
                            const Point3 vA = A.p0 + A.t1 * (stdFound ? a.x : fa.x) + A.t2 * (stdFound ? a.y : fa.y);
                            const Point3 vB = B.p0 + B.t1 * (stdFound ? b.x : fb.x) + B.t2 * (stdFound ? b.y : fb.y);
                            const std::string chord = std::string("  blocked: chord ") + (chordBlocked(scene, vA, vB) ? "YES" : "no") + ", receiver-A "
                                + (blocked(scene, rec, x, vA) ? "YES" : "no") + ", B-light " + (blockedFrom(scene, vB, y) ? "YES" : "no");
                            os << "    pair A" << ia << " (n " << A.n.x << " " << A.n.y << " " << A.n.z << ") B" << ib << " (n " << B.n.x << " " << B.n.y << " " << B.n.z
                               << "): normal seeds " << (stdFound ? "FOUND" : "no") << ", grid of starts found " << found << " (vA y " << vA.y << ", vB y " << vB.y << ")" << chord << "\n";
                        }
                    }
                }
            }
        }

    private:
        struct Facet {
            Point3 p0;
            Vector3 n, t1, t2;                       ///< outward normal and an orthonormal tangent frame (p0 is the origin of the plane coordinates)
            std::vector<std::array<Point2, 3>> tris; ///< the facet's triangles in (t1, t2) coordinates
            Point2 center;                           ///< mean of the triangle centers
            Real radius = 0;                         ///< largest distance from the center to a corner
        };

        struct BodyData {
            std::shared_ptr<const Dielectric> glass;
            std::vector<Facet> facets;
        };

        std::vector<BodyData> m_bodies;
        std::vector<const PointLight*> m_lights;
        bool m_anyDispersive = false;
        mutable std::atomic<long long> m_queries{ 0 };
        mutable std::atomic<long long> m_solutions{ 0 };
        mutable Counters m_c;

        using Vec4 = glm::vec<4, Real, glm::defaultp>;
        using Mat4 = glm::mat<4, 4, Real, glm::defaultp>;

        /// Coplanar triangles with the same outward normal form one facet (a flat face): the solver works per facet, which keeps the number of
        /// (receiver facet, light facet) pairs small.
        static void buildFacets(const TriangleMesh& mesh, std::vector<Facet>& facets) {
            for (size_t i = 0; i < mesh.numTriangles(); ++i) {
                const Point3 q0 = mesh.p[static_cast<size_t>(mesh.indices[3 * i])];
                const Point3 q1 = mesh.p[static_cast<size_t>(mesh.indices[3 * i + 1])];
                const Point3 q2 = mesh.p[static_cast<size_t>(mesh.indices[3 * i + 2])];
                const Vector3 c = glm::cross(Vector3(q1 - q0), Vector3(q2 - q0));
                if (glm::length(c) < Real(1e-14)) continue;   // degenerate
                const Vector3 n = glm::normalize(c);          // outward (counter-clockwise winding)
                Facet* host = nullptr;
                for (Facet& F : facets)
                    if (glm::dot(n, F.n) > Real(1) - Real(1e-9) && std::abs(glm::dot(n, Vector3(q0 - F.p0))) < Real(1e-9)) { host = &F; break; }
                if (!host) {
                    Facet f;
                    f.p0 = q0;
                    f.n = n;
                    f.t1 = glm::normalize(Vector3(q1 - q0));
                    f.t2 = glm::cross(f.n, f.t1);
                    facets.push_back(f);
                    host = &facets.back();
                }
                auto coords = [&](const Point3& P) { return Point2(glm::dot(Vector3(P - host->p0), host->t1), glm::dot(Vector3(P - host->p0), host->t2)); };
                host->tris.push_back({ coords(q0), coords(q1), coords(q2) });
            }
            for (Facet& F : facets) {
                Point2 c(0, 0);
                for (const auto& t : F.tris) c += (t[0] + t[1] + t[2]) / Real(3);
                F.center = c / Real(F.tris.size());
                for (const auto& t : F.tris) for (const Point2& q : t) F.radius = std::max(F.radius, glm::length(Point2(q - F.center)));
            }
        }

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

        /// The 2-D formulation: the ray from the receiver through the point a of facet A is refracted exactly (Snell) at A, followed to the plane
        /// of B, and R is what is left of Snell's law at that point of B (the tangential parts of n * (direction in the glass) - (direction to the
        /// light)). R = 0 is the solution. b receives the point of B in plane coordinates. False if the ray does not reach the plane of B.
        static bool traceFromA(const Point3& x, const Point3& y, const Facet& A, const Facet& B, Real n, const Point2& a, Real R[2], Point2& b) {
            const Point3 vA = A.p0 + A.t1 * a.x + A.t2 * a.y;
            Vector3 dIn;
            if (!refractDir(glm::normalize(Vector3(vA - x)), A.n, Real(1) / n, dIn)) return false;
            const Real den = glm::dot(B.n, dIn);
            if (den <= Real(1e-9)) return false;
            const Real s = glm::dot(B.n, Vector3(B.p0 - vA)) / den;
            if (!(s > Real(0))) return false;
            const Point3 vB = vA + dIn * s;
            const Vector3 h = glm::normalize(Vector3(y - vB)) - dIn * n;
            R[0] = glm::dot(h, B.t1);
            R[1] = glm::dot(h, B.t2);
            b = Point2(glm::dot(Vector3(vB - B.p0), B.t1), glm::dot(Vector3(vB - B.p0), B.t2));
            return std::isfinite(R[0] + R[1]);
        }

        /// Newton's method on traceFromA (numerical 2x2 Jacobian). True when the residual is below 1e-12; a and b then hold the solution.
        static bool newton2(const Point3& x, const Point3& y, const Facet& A, const Facet& B, Real n, Point2& a, Point2& b) {
            for (int it = 0; it < 10; ++it) {
                Real R[2];
                if (!traceFromA(x, y, A, B, n, a, R, b)) return false;
                const Real err = std::max(std::abs(R[0]), std::abs(R[1]));
                if (err < Real(1e-12)) return true;
                if (it >= 5 && err > Real(1e-4)) return false;
                if (glm::length(Point2(a - A.center)) > Real(3) * A.radius + Real(0.02)) return false;
                const Real du = Real(1e-7);
                Real Ra[2], Rb[2]; Point2 tmp;
                if (!traceFromA(x, y, A, B, n, Point2(a.x + du, a.y), Ra, tmp)) return false;
                if (!traceFromA(x, y, A, B, n, Point2(a.x, a.y + du), Rb, tmp)) return false;
                const Real j00 = (Ra[0] - R[0]) / du, j10 = (Ra[1] - R[1]) / du;   // d R / d a.x
                const Real j01 = (Rb[0] - R[0]) / du, j11 = (Rb[1] - R[1]) / du;   // d R / d a.y
                const Real det = j00 * j11 - j01 * j10;
                if (!(std::abs(det) > Real(1e-18))) return false;
                const Real sx = (j11 * R[0] - j01 * R[1]) / det, sy = (-j10 * R[0] + j00 * R[1]) / det;
                const Real lim = Real(0.02);
                a.x -= std::clamp(sx, -lim, lim);
                a.y -= std::clamp(sy, -lim, lim);
            }
            return false;
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

            int bin = 6;
            if (haveSnell) {
                const Real rl = glm::length(Point2(seed[0][1] - B.center)) / std::max(B.radius, Real(1e-9));
                bin = rl < 1 ? 0 : rl < 2 ? 1 : rl < 4 ? 2 : rl < 8 ? 3 : rl < 16 ? 4 : 5;
            }
            ++m_c.landTried[bin];

            // Measured on the faceted tumbler: 94 % of the pairs have no Snell seed (the ray refracted at the middle of A does not travel toward
            // B at all) and none of them ever produced a solution, so they are not worth a Newton run.
            static const bool noPrune = [] {   // diagnostics: the old behavior
                char* v = nullptr; size_t len = 0;
                const bool set = (_dupenv_s(&v, &len, "RAYT_MNEE_NO_PRUNE") == 0 && v != nullptr);
                std::free(v);
                return set;
            }();
            if (!haveSnell && !noPrune) return false;

            // The 2-D formulation (unknown: the point on A) from the middle of A and, if that fails, from four points around it. The 4-D Newton
            // iteration on (a, b) was found to miss thin, nearly parallel walls (an outer and an inner face of a tumbler), where the basin is tiny.
            static const Real offs[5][2] = { { 0, 0 }, { 0.5, 0 }, { -0.5, 0 }, { 0, 0.5 }, { 0, -0.5 } };
            for (int s = 0; s < 5; ++s) {
                a = A.center + Point2(offs[s][0], offs[s][1]) * A.radius;
                if (!newton2(x, y, A, B, n, a, b)) continue;
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
                ++m_c.landOk[bin];
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

        /// Does the chord vA -> vB (inside the glass) meet a surface before vB? (Then it crosses air, e.g. the cavity of a vessel.)
        static bool chordBlocked(const Scene& scene, const Point3& vA, const Point3& vB) {
            const Vector3 d = vB - vA;
            const Real dist = glm::length(d);
            if (dist <= Real(0)) return false;
            Ray r(vA, d / dist, Real(1e-6), nullptr);
            r.tMax = dist * (Real(1) - Real(1e-3));
            return scene.occluded(r);
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
