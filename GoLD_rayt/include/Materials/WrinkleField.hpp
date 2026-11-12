/**
 * @file Materials/WrinkleField.hpp
 * @brief Skin micro wrinkles and pores as a graph simulation (after Weiss et al. 2023, "Graph-Based Synthesis for Skin Micro Wrinkles").
 *
 * On real skin the finest wrinkles (the micro-relief) connect the pores; they begin and end at pores, form rectangles along the direction
 * of the stress (the "flow" of the skin) with shallower diagonal ones in between, and their density and depth change from place to place.
 * The method here (a CPU version of the paper's, for a flat tile) has three steps, repeated for several levels from coarse to fine:
 *
 *   1. Graph. Pores are nodes, Poisson-disk sampled with the average distance `spacing` (a finer level keeps the nodes of the coarser
 *      one and adds new ones between them). Every node is joined to its k nearest neighbours; the edges are candidate wrinkles.
 *   2. Simulation, inspired by how pedestrians wear trails into grass. Repeatedly a random node is chosen and one of its edges is
 *      chosen with a probability
 *          P(h) ~ max(eps, 1 - a_dist (d_h / d_min - 1) + p_cont + p_pref - a_cross p_cross)
 *      where d_h is the length of the edge (d_min: the shortest edge at the node), p_pref prefers edges along the primary orientation
 *      (a Gaussian of the angle plus the base strength a_s), p_cont rewards a continuation of a strong edge on the other side of the node
 *      and punishes edges that run alongside one (a_cont, a_sim), and p_cross punishes edges that cross strong ones. The chosen edge is
 *      carved: W += dt a_deposit (1 - W). After each round all weights decay: W *= 1 - dt a_decay. Strong edges are those that win again
 *      and again.
 *   3. Drawing. Each edge becomes a groove with a width proportional to its weight, a sideways wobble (so the line is not straight)
 *      and a profile with a valley and a small bulge beside it; each node becomes a pore (deeper where wrinkles meet). Fragments are
 *      blended by a smooth maximum (the paper's MellowMax) so that grooves join without creases.
 *
 * The result is a height map of one tile (a torus: the pattern repeats without a seam), `tileMm` wide. It is looked up in world space by
 * a blend of the three axis-aligned projections (triplanar), so any shape can use it without UV coordinates. Heights are in metres
 * (negative = a valley). The same pass gives a pore mask for the skin color.
 *
 * Differences from the paper: no UV map and no CUDA (a flat tile, the simulation runs on the CPU); the distance term is relative to the
 * shortest edge of the node; the weights are normalised by the largest weight of the level before drawing; one pass of drawing for all
 * levels. The formula of the selection probability was read from a text extraction of the paper in which some symbols were lost; the form
 * above is the one used here.
 *
 * The construction (graph, simulation, drawing) lives in src/Materials/WrinkleField.cpp; only the per-hit lookup is in this header.
 */

#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"

namespace rayt {

    class WrinkleField {
    public:
        struct Level {
            Real spacingMm = Real(1.5);     ///< average distance of the pores (a_d)
            Real angleDeg = Real(0);        ///< primary orientation of the wrinkles (a_theta), in the plane of the tile
            Real base = Real(0.2);          ///< a_s: 0 = only the primary direction, large = every direction alike
            Real distance = Real(0.7);      ///< a_dist: how much shorter edges are preferred
            Real continuation = Real(0.3);  ///< a_cont: reward for continuing a strong edge
            Real similarity = Real(0.3);    ///< a_sim: penalty for running alongside a strong edge
            Real crossing = Real(0.5);      ///< a_cross: penalty for crossing a strong edge
            Real deposit = Real(0.25);      ///< a_deposit: how fast a chosen edge is carved
            Real width = Real(0.15);        ///< a_wrinkle-width: groove width = this * spacing * weight
            Real poreWidth = Real(0.1);     ///< a_pore-width: pore radius = this * spacing * pore strength
            Real poreMin = Real(0.1);       ///< a_dmin: the weakest pore strength
            Real poreBlend = Real(0.5);     ///< a_blend: 0 = a pore as strong as its strongest edge, 1 = as the sum of its edges
            Real poreSkew = Real(0);        ///< a_skew: pores are elongated along the primary orientation by 1 + 2 * this
            Real wobble = Real(0.15);       ///< sideways wobble of a wrinkle, in units of the spacing
            Real wobbleFreq = Real(2);      ///< wobbles per spacing
            Real depthUm = Real(30);        ///< depth of the deepest wrinkle of this level (um)
            bool bulge = true;              ///< a small ridge beside each groove (the finest level may do without)
            bool darkPores = true;          ///< the pores of this level darken the skin color (the pore mask); the finest level has too many
        };

        struct Params {
            Real tileMm = Real(24);         ///< width of the tile (mm); the pattern repeats with this period
            int resolution = 2048;          ///< texels per side
            unsigned seed = 1;
            int steps = 300;                ///< rounds of the simulation per level
            Real depthScale = Real(1);      ///< multiplies all depths
            Real densityVariation = Real(0.5);  ///< 0 = pores evenly spread; larger: regions with fewer and more pores (a periodic noise)
            Real carveThreshold = Real(0.3);    ///< edges whose weight is below this fraction of the strongest one of their level are not drawn
            Real preference = Real(2);      ///< weight of the preference for the primary orientation in the selection of an edge
            Real selectivity = Real(2.5);   ///< the selection probability is the score to this power (more: the strong edges win more clearly)
            std::vector<Level> levels;      ///< empty: the default skin (skinDefault())
        };

        /// Three levels for a region with a weak preferred direction (a cheek): primary wrinkles with a slight direction along x, secondary ones across,
        /// and the finest, almost isotropic network. The finest level is wide enough to show at a pixel of 0.04 mm.
        static std::vector<Level> skinDefault();

        /// Builds the tile (prints two lines of progress; a few seconds at 2048^2).
        explicit WrinkleField(const Params& p);

        const Params& params() const { return m_p; }
        int resolution() const { return m_p.resolution; }
        const std::vector<float>& heightMap() const { return m_height; }   ///< metres, row-major, y * resolution + x
        const std::vector<float>& poreMap() const { return m_pore; }       ///< 0..1

        /// Height (m) at the world point p on a surface with the normal n (negative in a valley).
        Real height(const Point3& p, const Vector3& n) const { return triplanar(m_height, p, n); }

        /// 0..1: how much of a pore is at p.
        Real poreMask(const Point3& p, const Vector3& n) const { return triplanar(m_pore, p, n); }

    private:
        struct V2 { Real x = 0, y = 0; };
        struct HalfEdge { int edge = 0, to = 0; Real theta = 0, dist = 0; };
        struct Graph {
            std::vector<V2> nodes;
            std::vector<std::array<int, 2>> edges;
            std::vector<V2> edgeVec;                       // from edges[e][0] to edges[e][1]
            std::vector<std::vector<HalfEdge>> adj;        // per node
            std::vector<std::vector<int>> crossing;        // per edge: the edges that cross it
        };
        struct Canvas {
            std::vector<float> accV, accB, pore;
            std::vector<std::uint8_t> count;
        };

        // ------------------------------------------------------------------------------------------------------------------------------
        // lookup (per hit: stays in the header)
        // ------------------------------------------------------------------------------------------------------------------------------
        Real texel(const std::vector<float>& map, Real xmm, Real ymm) const {
            const int R = m_p.resolution;
            const Real u = xmm / m_p.tileMm * Real(R) - Real(0.5), v = ymm / m_p.tileMm * Real(R) - Real(0.5);
            const Real fu = std::floor(u), fv = std::floor(v);
            const Real a = u - fu, b = v - fv;
            auto wrap = [R](Real f) { int i = static_cast<int>(std::fmod(f, Real(R))); return i < 0 ? i + R : i; };
            const int x0 = wrap(fu), y0 = wrap(fv), x1 = (x0 + 1) % R, y1 = (y0 + 1) % R;
            const Real h00 = map[static_cast<size_t>(y0) * R + x0], h10 = map[static_cast<size_t>(y0) * R + x1];
            const Real h01 = map[static_cast<size_t>(y1) * R + x0], h11 = map[static_cast<size_t>(y1) * R + x1];
            return (h00 * (1 - a) + h10 * a) * (1 - b) + (h01 * (1 - a) + h11 * a) * b;
        }

        Real triplanar(const std::vector<float>& map, const Point3& p, const Vector3& n) const {
            Real wx = std::abs(n.x), wy = std::abs(n.y), wz = std::abs(n.z);
            wx *= wx; wx *= wx; wy *= wy; wy *= wy; wz *= wz; wz *= wz;   // |n|^4: a soft choice of the dominant axis
            const Real sum = wx + wy + wz;
            if (!(sum > Real(0))) return Real(0);
            const Real mx = p.x * Real(1000), my = p.y * Real(1000), mz = p.z * Real(1000);   // metres to mm
            Real h = 0;
            if (wx > 0) h += wx * texel(map, my, mz);
            if (wy > 0) h += wy * texel(map, mx, mz);
            if (wz > 0) h += wz * texel(map, mx, my);
            return h / sum;
        }

        // ------------------------------------------------------------------------------------------------------------------------------
        // construction (src/Materials/WrinkleField.cpp)
        // ------------------------------------------------------------------------------------------------------------------------------
        Real wrapDelta(Real d) const;
        V2 delta(const V2& a, const V2& b) const;
        static Real angDiff2Pi(Real a, Real b);
        static Real angDiffPi(Real a, Real b);
        Real density(const V2& q) const;
        void addNodes(std::vector<V2>& nodes, Real r, std::mt19937& rng) const;
        Graph makeGraph(const std::vector<V2>& nodes, int k, Real spacing) const;
        std::vector<Real> simulate(const Graph& g, const Level& L, std::mt19937& rng) const;
        static void profile(Real x, bool bulge, Real& valley, Real& ridge);
        static Real hash1(std::uint32_t a, std::uint32_t b);
        static Real noise1(Real x, std::uint32_t seed);
        void flush(Canvas& cv, int x0, int y0, int w, int h, const std::vector<float>& lv, const std::vector<float>& lb, Real scale) const;
        void drawEdge(Canvas& cv, const V2& a, const V2& d, Real strength, const Level& L, Real depthNorm, std::uint32_t seed) const;
        void drawPore(Canvas& cv, const V2& c, Real strength, const Level& L, Real depthNorm) const;
        void build();

        Params m_p;
        std::vector<float> m_height, m_pore;
        std::vector<std::array<Real, 3>> m_waves;
        bool m_quiet = false;
    };

} // namespace rayt
