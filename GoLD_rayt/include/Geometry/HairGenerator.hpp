/**
 * @file Geometry/HairGenerator.hpp
 * @brief Procedural strands for HairCurves: a hanging lock of hair, and hair that falls over a head.
 *
 * Not a grooming system: just enough structure for the fibers to look like hair, not like noise.
 *   - A strand is grown segment by segment: its direction turns toward gravity ("droop", radians per metre
 *     of arc length) and sways with a smooth wave whose phase depends on the root position, so
 *     neighbouring strands move together (frizz adds a little independent wobble).
 *   - Strands form clumps (the locks you see in real hair): every strand follows a guide path of its clump,
 *     more and more toward the tip, so the clumps separate by gaps while the fibers inside one stay together.
 *   - Lock: the clumps share one wavy centerline and sway around it each in their own way; the whole bunch fans out
 *     toward the tip, the ends are ragged and a few stray fibers leave the bunch (flyaways).
 *   - Head hair: roots on the cap of a sphere (a hairline cut out of the forehead), initial direction away from the
 *     scalp, then gravity; a strand that would enter the head sphere is pushed back onto its surface.
 * All lengths are in metres. A real hair is 40-60 micrometres in radius; the strand radius is a parameter because a
 * rendering often uses fibers a few times thicker (fewer strands, less noise).
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Geometry/HairCurves.hpp"

namespace rayt::hairgen {

    struct LockParams {
        int strands = 3000;
        Real length = Real(0.30);          ///< m
        int segments = 24;
        Real rootRadius = Real(0.00008);   ///< fiber radius at the root (m)
        Real tipRadius = Real(0.00005);
        Real width = Real(0.012);          ///< radius of the bunch at the root (m)
        Real tipWidth = Real(1.4);         ///< width of the bunch at the tip relative to the root (1 = a cylinder)
        int clumps = 9;                    ///< number of sub-locks
        Real waveAmp = Real(0.018);        ///< sideways amplitude of the shared centerline (m)
        Real waveLength = Real(0.16);      ///< m
        Real twist = Real(0.5);            ///< turns of the clumps around the centerline over the whole length
        Real frizz = Real(0.0006);         ///< amplitude of the independent wobble of each strand (m)
        Real flyaway = Real(0.03);         ///< fraction of the strands that wander off the bunch
        Point3 root = Point3(0.0);         ///< the top of the lock (it hangs down the -y axis)
        unsigned seed = 1;
    };

    inline void MakeLock(HairCurves& hair, const LockParams& P) {
        std::mt19937 rng(P.seed);
        std::uniform_real_distribution<Real> U(0, 1);
        const Real TWO_PI = Real(2) * constants::PI;
        const int n = std::max(P.segments, 2);
        const int nc = std::max(P.clumps, 1);
        std::vector<Point3> pts(static_cast<size_t>(n) + 1);

        // Shared centerline c(s): hangs down -y, sways along x and z (different phases) so it is a gentle 3-D wave.
        auto center = [&](Real s) {
            const Real k = TWO_PI / P.waveLength;
            return P.root + Vector3(P.waveAmp * std::sin(k * s), -s, P.waveAmp * Real(0.6) * std::sin(k * s * Real(0.7) + Real(1.0)));
        };

        // Each clump: where it is rooted, how wide it is, and how it sways away from the centerline.
        struct Clump { Real angle, dist, swayAmpX, swayAmpZ, swayFreq, phaseX, phaseZ, curl; };
        std::vector<Clump> clumps(static_cast<size_t>(nc));
        for (Clump& c : clumps) {
            c.angle = TWO_PI * U(rng);
            c.dist = P.width * std::sqrt(U(rng)) * Real(0.8);
            c.swayAmpX = P.width * Real(0.9) * (Real(0.3) + U(rng));
            c.swayAmpZ = P.width * Real(0.9) * (Real(0.3) + U(rng));
            c.swayFreq = TWO_PI / (P.waveLength * (Real(0.5) + Real(0.9) * U(rng)));
            c.phaseX = TWO_PI * U(rng); c.phaseZ = TWO_PI * U(rng);
            c.curl = (U(rng) - Real(0.5)) * Real(2);
        }
        const Real clumpRadius = P.width * Real(1.1) / std::sqrt(Real(nc));

        for (int i = 0; i < P.strands; ++i) {
            const Clump& c = clumps[static_cast<size_t>(i) % static_cast<size_t>(nc)];
            const bool stray = U(rng) < P.flyaway;
            // offset of the strand inside its clump (kept for the whole length, shrinking a little toward the tip)
            const Real rr = clumpRadius * std::sqrt(U(rng)), th = TWO_PI * U(rng);
            const Real f1 = Real(6) + Real(8) * U(rng), f2 = Real(14) + Real(20) * U(rng);   // wobble frequencies (rad / m)
            const Real p1 = TWO_PI * U(rng), p2 = TWO_PI * U(rng), p3 = TWO_PI * U(rng), p4 = TWO_PI * U(rng);
            const Real lenScale = stray ? Real(0.75) + Real(0.35) * U(rng) : Real(0.90) + Real(0.10) * U(rng);   // ragged ends
            const Real strayAngle = TWO_PI * U(rng), strayDist = P.width * (Real(1.5) + Real(4) * U(rng));
            for (int k = 0; k <= n; ++k) {
                const Real u = Real(k) / Real(n);
                const Real s = u * P.length * lenScale;
                const Point3 cl = center(s);
                // frame around the centerline: the tangent is nearly -y, so x and z span the cross-section
                const Vector3 T = glm::normalize(center(s + Real(1e-3)) - center(std::max(s - Real(1e-3), Real(0))));
                Vector3 ex = glm::normalize(Vector3(1, 0, 0) - T * T.x);
                const Vector3 ez = glm::cross(T, ex);
                const Real fan = Real(1) + (P.tipWidth - Real(1)) * u;
                // the clump center: rooted at (angle, dist), fanning out, turning around the axis, swaying on its own
                const Real ang = c.angle + P.twist * TWO_PI * u * (Real(0.5) + c.curl);
                Real cx = c.dist * fan * std::cos(ang) + c.swayAmpX * u * std::sin(c.swayFreq * s + c.phaseX);
                Real cz = c.dist * fan * std::sin(ang) + c.swayAmpZ * u * std::sin(c.swayFreq * s * Real(0.8) + c.phaseZ);
                // the strand inside the clump
                const Real shrink = Real(1) - Real(0.3) * u;
                cx += rr * shrink * std::cos(th + P.twist * TWO_PI * u); cz += rr * shrink * std::sin(th + P.twist * TWO_PI * u);
                // independent wobble, growing from the root (the root is fixed)
                const Real wob = P.frizz * std::sqrt(u) * (stray ? Real(4) : Real(1));
                cx += wob * (std::sin(f1 * s + p1) + Real(0.5) * std::sin(f2 * s + p2));
                cz += wob * (std::sin(f1 * s + p3) + Real(0.5) * std::sin(f2 * s + p4));
                if (stray) { cx += strayDist * u * u * std::cos(strayAngle); cz += strayDist * u * u * std::sin(strayAngle); }
                pts[static_cast<size_t>(k)] = cl + ex * cx + ez * cz;
            }
            hair.addStrand(pts, P.rootRadius, P.tipRadius);
        }
    }

    struct HeadHairParams {
        Point3 headCenter = Point3(0.0);
        Real headRadius = Real(0.09);
        int strands = 40000;
        Real length = Real(0.28);
        int segments = 26;
        Real rootRadius = Real(0.00008);
        Real tipRadius = Real(0.00005);
        Real droop = Real(9.0);            ///< how fast a strand turns toward gravity (radians per metre of length)
        Real outward = Real(0.55);         ///< weight of the scalp normal in the initial direction (the rest follows the head down)
        Real waveAmp = Real(5.0);          ///< sway of the direction (radians per metre)
        Real waveFreq = Real(28.0);        ///< radians per metre of arc length
        Real frizz = Real(1.2);            ///< independent sway (radians per metre)
        int clumps = 500;                  ///< guide paths; every strand follows the nearest one
        Real clumpStrength = Real(0.65);   ///< how closely the strands of a clump stay together at the tip (0 = not at all, 1 = on the guide)
        Real sweep = Real(1.0);            ///< how strongly hair near the face is combed back and to the sides (0 = straight down)
        Real maxPolarDeg = Real(105);      ///< the cap of hair reaches this angle from the crown (up axis)
        Real hairlineHeight = Real(0.55);  ///< forehead cut: no roots where (dir . front) > 0.5 and dir.y below this
        unsigned seed = 7;
    };

    inline void MakeHeadHair(HairCurves& hair, const HeadHairParams& P) {
        std::mt19937 rng(P.seed);
        std::uniform_real_distribution<Real> U(0, 1);
        const Real TWO_PI = Real(2) * constants::PI;
        const int n = std::max(P.segments, 2);
        const Real ds = P.length / Real(n);
        const Real cosMax = std::cos(P.maxPolarDeg * constants::PI / Real(180));
        const Vector3 down(0, -1, 0);

        // A root direction on the cap, or false if the random direction falls outside it / on the forehead.
        auto randomRoot = [&](Vector3& dir) {
            const Real z = Real(1) - Real(2) * U(rng), ph = TWO_PI * U(rng);
            const Real rxy = std::sqrt(std::max(Real(0), Real(1) - z * z));
            dir = Vector3(rxy * std::cos(ph), z, rxy * std::sin(ph));   // y is up; the face looks along +z
            if (dir.y < cosMax) return false;
            if (dir.z > Real(0.5) && dir.y < P.hairlineHeight + Real(0.25) * (dir.z - Real(0.5))) return false;
            return true;
        };

        // Grows one strand from a root direction; the phases make its sway.
        auto grow = [&](const Vector3& dir, Real phase1, Real phase2, Real fphase1, Real fphase2, Real fr, Real frizz, int nSeg, std::vector<Point3>& pts) {
            Point3 p = P.headCenter + dir * (P.headRadius + Real(0.0002));
            // combing direction: down, and near the face swept back and to the sides (a part at x = 0), so the face stays free
            const Real front = std::max(dir.z, Real(0));
            const Vector3 comb = down + Vector3(0, 0, -1) * (P.sweep * Real(1.2) * front) + Vector3(dir.x >= 0 ? Real(1) : Real(-1), 0, 0) * (P.sweep * Real(1.0) * front * front);
            Vector3 tang = comb - dir * glm::dot(comb, dir);   // along the scalp, away from the crown
            const Real tl = glm::length(tang);
            tang = tl > Real(1e-6) ? tang / tl : Vector3(dir.z, 0, -dir.x);
            Vector3 d = glm::normalize(dir * P.outward + tang * (Real(1) - P.outward) + Vector3(0, -0.05, 0));
            pts.assign(static_cast<size_t>(nSeg) + 1, p);
            for (int k = 1; k <= nSeg; ++k) {
                const Real s = Real(k) * ds;
                d += down * (P.droop * ds);
                Vector3 side = glm::cross(d, Vector3(0, 1, 0));
                const Real sl = glm::length(side);
                side = sl > Real(1e-6) ? side / sl : Vector3(1, 0, 0);
                const Vector3 up2 = glm::cross(side, d);
                d += side * (ds * (P.waveAmp * std::sin(P.waveFreq * s + phase1) + frizz * std::sin(fr * s + fphase1)));
                d += up2 * (ds * (P.waveAmp * Real(0.5) * std::sin(P.waveFreq * s * Real(0.8) + phase2) + frizz * std::sin(fr * s * Real(1.3) + fphase2)));
                d = glm::normalize(d);
                p += d * ds;
                const Vector3 rel = p - P.headCenter;   // keep the strand outside the head and slide it along the surface
                const Real dist = glm::length(rel);
                const Real minDist = P.headRadius + Real(0.0006);
                if (dist < minDist) {
                    p = P.headCenter + rel * (minDist / std::max(dist, Real(1e-9)));
                    const Vector3 nrm = glm::normalize(p - P.headCenter);
                    d = glm::normalize(d - nrm * std::min(glm::dot(d, nrm), Real(0)));
                }
                pts[static_cast<size_t>(k)] = p;
            }
        };

        // The guides (smooth paths, no frizz): the phases depend on the root direction so that nearby guides move together.
        struct Guide { Vector3 dir; std::vector<Point3> pts; };
        std::vector<Guide> guides;
        const int nGuides = std::max(P.clumps, 1);
        while (static_cast<int>(guides.size()) < nGuides) {
            Guide g;
            if (!randomRoot(g.dir)) continue;
            grow(g.dir, Real(9.0) * g.dir.x + Real(7.0) * g.dir.z, Real(8.0) * g.dir.z - Real(6.0) * g.dir.x, 0, 0, 1, 0, n, g.pts);
            guides.push_back(std::move(g));
        }

        std::vector<Point3> pts;
        int made = 0;
        while (made < P.strands) {
            Vector3 dir;
            if (!randomRoot(dir)) continue;
            const Real phase1 = Real(9.0) * dir.x + Real(7.0) * dir.z + Real(0.6) * (U(rng) - Real(0.5));   // neighbours share the phase
            const Real phase2 = Real(8.0) * dir.z - Real(6.0) * dir.x + Real(0.6) * (U(rng) - Real(0.5));
            const Real fphase1 = TWO_PI * U(rng), fphase2 = TWO_PI * U(rng);
            const Real fr = Real(20) + Real(25) * U(rng);
            const Real lenScale = Real(0.8) + Real(0.2) * U(rng);
            const int nSeg = std::max(2, int(std::lround(n * lenScale)));
            grow(dir, phase1, phase2, fphase1, fphase2, fr, P.frizz, nSeg, pts);

            // pull the strand toward the guide of the nearest root (more and more toward the tip)
            size_t best = 0;
            Real bestDot = -2;
            for (size_t gi = 0; gi < guides.size(); ++gi) {
                const Real dd = glm::dot(guides[gi].dir, dir);
                if (dd > bestDot) { bestDot = dd; best = gi; }
            }
            const Guide& g = guides[best];
            for (int k = 1; k <= nSeg; ++k) {
                const Real u = Real(k) / Real(n);
                const Real t = std::clamp((u - Real(0.1)) / Real(0.7), Real(0), Real(1));
                const Real w = P.clumpStrength * t * t * (Real(3) - Real(2) * t);   // smoothstep
                Point3 q = pts[static_cast<size_t>(k)] + (g.pts[static_cast<size_t>(k)] - pts[static_cast<size_t>(k)]) * w;
                const Vector3 rel = q - P.headCenter;
                const Real dist = glm::length(rel), minDist = P.headRadius + Real(0.0006);
                if (dist < minDist) q = P.headCenter + rel * (minDist / std::max(dist, Real(1e-9)));
                pts[static_cast<size_t>(k)] = q;
            }
            hair.addStrand(pts, P.rootRadius, P.tipRadius);
            ++made;
        }
    }

} // namespace rayt::hairgen
