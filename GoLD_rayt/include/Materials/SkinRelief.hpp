/**
 * @file Materials/SkinRelief.hpp
 * @brief The micro-relief of skin (furrows and pores) as a procedural height field, used to bend the shading normal.
 *
 * Real skin is not smooth: a network of fine furrows (the skin lines) divides it into small polygons about half a
 * millimetre wide, and pores sit in some of them. The pattern is cellular, so it is built from a 3-D Worley (cellular) noise
 * in world space (no UV map is needed, the pattern is the same on every object that uses the material):
 *   - the distance to the nearest and the second nearest feature point, F1 and F2, of the jittered grid of cells;
 *   - furrows: a groove where F2 - F1 is small (the cell boundaries): depth * exp(-((F2 - F1) / width)^2); the depth varies from
 *     groove to groove (a random factor per pair of neighbouring cells, the same from both sides, so the groove stays continuous);
 *   - two scales: the fine network and a coarser one (cells `coarseScale` times larger, `coarseWeight` times as deep): real skin has
 *     a few deep primary lines and many shallow secondary ones;
 *   - `stretch` elongates the cells along the x axis (skin lines run in a preferred direction);
 *   - pores: in a fraction of the cells a round dimple around the feature point: depth * exp(-(F1 / radius)^2).
 * The height h(p) (metres, negative = a dent) changes the normal, n' = normalize(n - grad_t h), the way a bump map does: the
 * geometry is not displaced, so the silhouette stays smooth and at grazing angles there are the usual artefacts of bump mapping
 * (shadowing inside the dents is not computed). The gradient is taken by central differences along two tangents.
 *
 * Sizes are in metres. Typical values: cell 0.5-0.8 mm, furrows 30-50 um deep, pores 60-100 um deep and 70 um in radius.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "Core/Types.hpp"
#include "Core/Math.hpp"

namespace rayt {

    class SkinRelief {
    public:
        struct Params {
            Real cell = Real(0.0006);          ///< size of a furrow cell (m)
            Real furrowDepth = Real(0.00004);  ///< depth of the furrows (m)
            Real furrowWidth = Real(0.12);     ///< width of the furrows in units of the cell
            Real stretch = Real(1.0);          ///< cells are this much longer along x (1 = round)
            Real coarseScale = Real(3.0);      ///< the coarse network has cells this many times larger ...
            Real coarseWeight = Real(0.8);     ///< ... and furrows this many times as deep as the fine one (0 = fine network only)
            Real poreDepth = Real(0.00008);    ///< depth of the pores (m)
            Real poreRadius = Real(0.00007);   ///< radius of the pores (m)
            Real poreFraction = Real(0.35);    ///< fraction of the cells that hold a pore
            unsigned seed = 1;
        };

        explicit SkinRelief(const Params& p) : m_p(p) {}

        const Params& params() const { return m_p; }

        /// Height of the relief at p (m): 0 on the smooth skin, negative in a furrow or a pore.
        Real height(const Point3& p) const {
            Point3 q = p / m_p.cell;
            q.x /= m_p.stretch;
            Real f1, f2;
            unsigned id, id2;
            worley(q, m_p.seed, f1, f2, id, id2);
            const Real e = (f2 - f1) / m_p.furrowWidth;
            Real h = -m_p.furrowDepth * (Real(0.3) + Real(0.7) * Real(0.5) * (unit(id, 0u) + unit(id2, 0u))) * std::exp(-e * e);
            if (m_p.coarseWeight > Real(0)) {
                Real g1, g2;
                unsigned jd, jd2;
                worley(q / m_p.coarseScale, m_p.seed + 7919u, g1, g2, jd, jd2);
                const Real ec = (g2 - g1) / (Real(0.7) * m_p.furrowWidth);
                h -= m_p.furrowDepth * m_p.coarseWeight * (Real(0.3) + Real(0.7) * Real(0.5) * (unit(jd, 0u) + unit(jd2, 0u))) * std::exp(-ec * ec);
            }
            if (unit(id, 0u) < m_p.poreFraction) {
                const Real r = f1 * m_p.cell / m_p.poreRadius;
                h -= m_p.poreDepth * std::exp(-r * r);
            }
            return h;
        }

        /// The bent normal at p for the smooth normal n. `wo` (towards the viewer) keeps it from turning away from the viewer.
        Vector3 perturb(const Point3& p, const Vector3& n, const Vector3& wo) const {
            const Vector3 axis = std::abs(n.y) < Real(0.9) ? Vector3(0, 1, 0) : Vector3(1, 0, 0);
            const Vector3 t = glm::normalize(glm::cross(n, axis));
            const Vector3 b = glm::cross(n, t);
            const Real eps = Real(0.04) * std::min(m_p.poreRadius, m_p.furrowWidth * m_p.cell);   // a few micrometres
            const Real dt = (height(p + t * eps) - height(p - t * eps)) / (Real(2) * eps);
            const Real db = (height(p + b * eps) - height(p - b * eps)) / (Real(2) * eps);
            Vector3 np = glm::normalize(n - t * dt - b * db);
            const Real c = glm::dot(np, wo);
            if (c < Real(0.05)) np = glm::normalize(np + wo * (Real(0.05) - c));   // never turn the normal away from the viewer
            return np;
        }

    private:
        static unsigned hash(int x, int y, int z, unsigned seed, unsigned salt) {
            std::uint32_t h = static_cast<std::uint32_t>(x) * 0x8da6b343u ^ static_cast<std::uint32_t>(y) * 0xd8163841u
                ^ static_cast<std::uint32_t>(z) * 0xcb1ab31fu ^ seed * 0x9e3779b9u ^ salt * 0x85ebca6bu;
            h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16;
            return h;
        }
        static Real unit(unsigned h, unsigned) { return Real(h) / Real(4294967296.0); }

        /// F1, F2 and the ids (hashes) of the nearest and the second nearest cell, for the point q in cell units.
        static void worley(const Point3& q, unsigned seed, Real& f1, Real& f2, unsigned& id, unsigned& id2) {
            const int ix = static_cast<int>(std::floor(q.x)), iy = static_cast<int>(std::floor(q.y)), iz = static_cast<int>(std::floor(q.z));
            f1 = f2 = Real(1e9);
            id = id2 = 0;
            for (int dz = -1; dz <= 1; ++dz)
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dx = -1; dx <= 1; ++dx) {
                        const int cx = ix + dx, cy = iy + dy, cz = iz + dz;
                        const Point3 f(Real(cx) + unit(hash(cx, cy, cz, seed, 0u), 0u),
                            Real(cy) + unit(hash(cx, cy, cz, seed, 1u), 0u),
                            Real(cz) + unit(hash(cx, cy, cz, seed, 2u), 0u));
                        const Vector3 d = q - f;
                        const Real dist = std::sqrt(glm::dot(d, d));
                        if (dist < f1) { f2 = f1; id2 = id; f1 = dist; id = hash(cx, cy, cz, seed, 3u); }
                        else if (dist < f2) { f2 = dist; id2 = hash(cx, cy, cz, seed, 3u); }
                    }
        }

        Params m_p;
    };

} // namespace rayt
