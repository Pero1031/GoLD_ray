/**
 * @file Materials/SkinPigment.hpp
 * @brief The unevenness of skin color as a procedural field: freckles, uneven pigment, and redness.
 *
 * Real skin is never one color. The field gives an RGB factor F(p) in (0, 1] (p in world space, no UV map needed) by which the
 * reflectance of the skin is multiplied:
 *   - freckles: small brown spots, each with its own size and darkness, irregular at the edge (domain-warped), and clustered:
 *     the chance of a spot in a cell follows a low-frequency noise (sun-exposed areas have many, others none). Spots are put
 *     on the jittered points of a 3-D grid, so no two are alike and none is a grid;
 *   - uneven pigment: broad, soft, faint patches (a fractal noise);
 *   - redness: a very low-frequency flush that lowers green and blue (blood absorbs there).
 * Pores darken the skin too (a shadowed, slightly red dot): the caller multiplies in SkinRelief::poreMask.
 *
 * The darkening of a spot is taken from the colors of a tan spot against fair skin (about 0.58 / 0.48 / 0.42 of the skin around).
 * All numbers are plausible values chosen by eye, not measurements. Lengths are in metres.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "Core/Types.hpp"
#include "Materials/DensityField.hpp"

namespace rayt {

    class SkinPigment {
    public:
        struct Params {
            Real freckleCell = Real(0.0025);     ///< grid size of the freckle spots (m)
            Real freckleFraction = Real(0.5);    ///< chance of a spot in a cell where the clustering noise is at its maximum
            Real freckleRadius = Real(0.0008);   ///< mean radius of a spot (m); each is 0.5-1.5 times this
            Real freckleStrength = Real(0.8);    ///< darkening of the darkest spots (0..1)
            Real clusterFreq = Real(40);         ///< features per metre of the clustering noise (about 2.5 cm patches)
            Real blotchFreq = Real(150);         ///< features per metre of the uneven pigment
            Real blotchStrength = Real(0.15);    ///< its darkening at the most (0..1)
            Real rednessFreq = Real(35);         ///< features per metre of the flush
            Real rednessStrength = Real(0.5);    ///< 0 = none, 1 = a strong flush
            unsigned seed = 1;
        };

        explicit SkinPigment(const Params& p)
            : m_p(p)
            , m_cluster(p.clusterFreq, 2, Real(0.35), Real(0.65), p.seed + 11u)
            , m_blotch(p.blotchFreq, 2, Real(0.3), Real(0.75), p.seed + 23u)
            , m_flush(p.rednessFreq, 1, Real(0.25), Real(0.75), p.seed + 37u)
            , m_warpX(Real(2500), 1, Real(0), Real(1), p.seed + 41u)
            , m_warpY(Real(2500), 1, Real(0), Real(1), p.seed + 43u)
            , m_warpZ(Real(2500), 1, Real(0), Real(1), p.seed + 47u) {}

        const Params& params() const { return m_p; }

        /// The multiplier of the reflectance at p, per RGB channel, in (0, 1].
        Vector3 factor(const Point3& p) const {
            const Real fr = freckles(p);
            const Real bl = m_p.blotchStrength > Real(0) ? m_blotch(p) * m_p.blotchStrength : Real(0);
            const Real rd = m_p.rednessStrength > Real(0) ? m_flush(p) * m_p.rednessStrength : Real(0);
            // a tan spot: (0.58, 0.48, 0.42) of the skin at full strength; blotches the same hue, fainter; flush takes green and blue
            const Vector3 f = (Vector3(1.0) - fr * Vector3(0.42, 0.52, 0.58))
                * (Vector3(1.0) - bl * Vector3(0.20, 0.28, 0.33))
                * (Vector3(1.0) - rd * Vector3(0.0, 0.16, 0.10));
            return glm::clamp(f, Vector3(0.05), Vector3(1.0));
        }

    private:
        static std::uint32_t hash(int x, int y, int z, std::uint32_t seed, std::uint32_t salt) {
            std::uint32_t h = static_cast<std::uint32_t>(x) * 0x8da6b343u ^ static_cast<std::uint32_t>(y) * 0xd8163841u
                ^ static_cast<std::uint32_t>(z) * 0xcb1ab31fu ^ seed * 0x9e3779b9u ^ salt * 0x85ebca6bu;
            h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16;
            return h;
        }
        static Real unit(std::uint32_t h) { return Real(h) / Real(4294967296.0); }

        /// Darkening of the freckles at p in 0..1: the strongest spot there (spots are soft Gaussians, combined as 1 - prod(1 - s)).
        Real freckles(const Point3& p) const {
            if (!(m_p.freckleStrength > Real(0)) || !(m_p.freckleFraction > Real(0))) return Real(0);
            // irregular edges: the position is moved by a little noise (up to ~0.3 of a spot radius)
            const Real w = Real(0.6) * m_p.freckleRadius;
            const Point3 q = p + Vector3(m_warpX(p) - Real(0.5), m_warpY(p) - Real(0.5), m_warpZ(p) - Real(0.5)) * w;
            const Point3 g = q / m_p.freckleCell;
            const int ix = static_cast<int>(std::floor(g.x)), iy = static_cast<int>(std::floor(g.y)), iz = static_cast<int>(std::floor(g.z));
            Real keep = Real(1);
            for (int dz = -1; dz <= 1; ++dz)
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dx = -1; dx <= 1; ++dx) {
                        const int cx = ix + dx, cy = iy + dy, cz = iz + dz;
                        const Point3 c((Real(cx) + unit(hash(cx, cy, cz, m_p.seed, 0u))) * m_p.freckleCell,
                            (Real(cy) + unit(hash(cx, cy, cz, m_p.seed, 1u))) * m_p.freckleCell,
                            (Real(cz) + unit(hash(cx, cy, cz, m_p.seed, 2u))) * m_p.freckleCell);
                        const Real u = unit(hash(cx, cy, cz, m_p.seed, 3u));
                        if (u >= m_p.freckleFraction || u >= m_p.freckleFraction * m_cluster(c)) continue;   // clustered (the cheap test first)
                        const Real radius = m_p.freckleRadius * (Real(0.5) + unit(hash(cx, cy, cz, m_p.seed, 4u)));
                        const Real strength = m_p.freckleStrength * (Real(0.4) + Real(0.6) * unit(hash(cx, cy, cz, m_p.seed, 5u)));
                        const Vector3 d = q - c;
                        const Real r2 = glm::dot(d, d) / (radius * radius);
                        if (r2 > Real(9)) continue;
                        keep *= Real(1) - strength * std::exp(-std::pow(r2, Real(1.5)));   // soft, but with an edge near the radius
                    }
            return Real(1) - keep;
        }

        Params m_p;
        NoiseDensity m_cluster, m_blotch, m_flush, m_warpX, m_warpY, m_warpZ;
    };

} // namespace rayt
