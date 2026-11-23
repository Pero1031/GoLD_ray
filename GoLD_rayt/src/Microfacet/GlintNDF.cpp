/**
 * @file Microfacet/GlintNDF.cpp
 * @brief The glinty GGX distribution of Kemppinen et al. 2025 (see the header and docs/GLINTS.md).
 */

#include "pch.h"

#include "Microfacet/GlintNDF.hpp"

#include <algorithm>
#include <cmath>

#include "Core/Constants.hpp"
#include "Core/Math.hpp"

namespace rayt {

    namespace {

        /// splitmix64: a well-mixed 64-bit value from a 64-bit state.
        inline std::uint64_t Mix(std::uint64_t z) {
            z += 0x9E3779B97F4A7C15ull;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
            return z ^ (z >> 31);
        }

        /// The random numbers of one facet: seeded by its cell (two spatial and two angular indices), its level and the arrangement.
        struct FacetRandom {
            std::uint64_t s;
            FacetRandom(std::int64_t ix, std::int64_t iy, std::int64_t ax, std::int64_t ay, int level, std::uint32_t seed, std::uint32_t stream) {
                s = Mix(static_cast<std::uint64_t>(ix) ^ (static_cast<std::uint64_t>(stream) << 48));
                s = Mix(s ^ static_cast<std::uint64_t>(iy));
                s = Mix(s ^ static_cast<std::uint64_t>(ax));
                s = Mix(s ^ static_cast<std::uint64_t>(ay));
                s = Mix(s ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(level)) << 32 | seed));
            }
            Real next() {
                s = Mix(s);
                return static_cast<Real>(s >> 11) * Real(1.0 / 9007199254740992.0);   // 53 bits: [0, 1)
            }
        };

        /// The mass of a 1D Gaussian (mean x, standard deviation sigma) between a and b.
        inline Real GaussMass(Real a, Real b, Real x, Real sigma) {
            const Real k = Real(1) / (sigma * Real(1.4142135623730951));
            return Real(0.5) * (std::erf((b - x) * k) - std::erf((a - x) * k));
        }

    }   // namespace

    /// One evaluation's constants: the query in grid units (level 0 cells of side 1) and in the disc.
    struct GlintNDF::Query {
        Real gx, gy;            ///< the texture point
        Real sxx, sxy, syy;     ///< the filter's covariance
        Real ixx, ixy, iyy;     ///< its inverse
        Real norm;              ///< 1 / (2 pi sqrt(det))
        Real hx, hy;            ///< half extents of the box that holds kFilterReach standard deviations
        Real sigmaMax;          ///< the standard deviation along the major axis
        Point2 xa;              ///< T(wh)
        Real sigmaA;            ///< the lobe's standard deviation in the disc
    };

    GlintNDF::GlintNDF(const GlintParams& params) : m_p(params) {
        m_p.alphaX = std::max(m_p.alphaX, Real(1e-4));
        m_p.alphaY = std::max(m_p.alphaY, Real(1e-4));
        m_p.microRoughness = std::max(m_p.microRoughness, Real(1e-6));
        m_p.density = std::max(m_p.density, Real(1e-12));
        m_p.angularCells = std::clamp(m_p.angularCells, 1, 8);
        m_p.rouletteSoftness = std::clamp(m_p.rouletteSoftness, Real(0), Real(0.5));
    }

    Real GlintNDF::smoothD(const Vector3& wh) const {
        if (!(wh.z > Real(0))) return 0;
        // D = det M / (pi |M w|^4) with M = diag(1/ax, 1/ay, 1): the GGX distribution
        const Real x = wh.x / m_p.alphaX, y = wh.y / m_p.alphaY;
        const Real n2 = x * x + y * y + wh.z * wh.z;
        return Real(1) / (constants::PI * m_p.alphaX * m_p.alphaY * n2 * n2);
    }

    Point2 GlintNDF::toDisc(const Vector3& wh) const {
        Vector3 s(wh.x / m_p.alphaX, wh.y / m_p.alphaY, wh.z);
        s /= glm::length(s);
        const Real k = Real(0.5) / std::sqrt(Real(1) + s.z);
        return Point2(Real(0.5) + s.x * k, Real(0.5) + s.y * k);
    }

    Vector3 GlintNDF::fromDisc(const Point2& u) const {
        const Real px = Real(2) * (u.x - Real(0.5)), py = Real(2) * (u.y - Real(0.5));
        const Real r2 = std::min(px * px + py * py, Real(1));
        const Real k = std::sqrt(Real(2) - r2);
        const Vector3 m(m_p.alphaX * px * k, m_p.alphaY * py * k, Real(1) - r2);
        const Real len = glm::length(m);
        return len > Real(0) ? m / len : Vector3(0, 0, 1);
    }

    Real GlintNDF::jacobian(const Vector3& wh) const {
        const Real x = wh.x / m_p.alphaX, y = wh.y / m_p.alphaY;
        const Real n = std::sqrt(x * x + y * y + wh.z * wh.z);
        return Real(1) / (Real(8) * m_p.alphaX * m_p.alphaY * n * n * n);
    }

    Real GlintNDF::level(const GlintFootprint& fp) const {
        const Real tr = (fp.xx + fp.yy) * m_p.density, det = (fp.xx * fp.yy - fp.xy * fp.xy) * m_p.density * m_p.density;
        const Real lmax = Real(0.5) * (tr + std::sqrt(std::max(tr * tr - Real(4) * det, Real(0))));
        if (!(lmax > Real(0))) return -std::numeric_limits<Real>::infinity();
        return std::log2(kCellsPerSigma * std::sqrt(lmax));
    }

    Real GlintNDF::levelSum(const Query& q, int level, Real weight, bool exhaustive, std::uint32_t stream) const {
        const Real cell = std::ldexp(Real(1), level);                  // side of a spatial cell
        const std::int64_t resA = std::int64_t(1) << level;            // angular cells per side
        const Real invResA = Real(1) / static_cast<Real>(resA);
        // the angular block summed facet by facet: the cells nearest to T(wh), within the unit square
        std::int64_t a0x = 0, a1x = resA - 1, a0y = 0, a1y = resA - 1;
        if (!exhaustive) {
            const Real back = Real(0.5) * static_cast<Real>(m_p.angularCells - 1);
            a0x = static_cast<std::int64_t>(std::floor(q.xa.x * static_cast<Real>(resA) - back));
            a0y = static_cast<std::int64_t>(std::floor(q.xa.y * static_cast<Real>(resA) - back));
            a1x = std::min(a0x + m_p.angularCells - 1, resA - 1);
            a1y = std::min(a0y + m_p.angularCells - 1, resA - 1);
            a0x = std::max<std::int64_t>(a0x, 0);
            a0y = std::max<std::int64_t>(a0y, 0);
        }
        const Real reach = exhaustive ? Real(1.5) : Real(1);
        const std::int64_t i0x = static_cast<std::int64_t>(std::floor((q.gx - reach * q.hx) / cell));
        const std::int64_t i1x = static_cast<std::int64_t>(std::floor((q.gx + reach * q.hx) / cell));
        const std::int64_t i0y = static_cast<std::int64_t>(std::floor((q.gy - reach * q.hy) / cell));
        const std::int64_t i1y = static_cast<std::int64_t>(std::floor((q.gy + reach * q.hy) / cell));
        const Real invVarA = Real(1) / (q.sigmaA * q.sigmaA);
        const Real normA = invVarA / (Real(2) * constants::PI);
        const Real eps = m_p.rouletteSoftness;
        Real sum = 0;
        for (std::int64_t iy = i0y; iy <= i1y; ++iy)
            for (std::int64_t ix = i0x; ix <= i1x; ++ix)
                for (std::int64_t ay = a0y; ay <= a1y; ++ay)
                    for (std::int64_t ax = a0x; ax <= a1x; ++ax) {
                        FacetRandom rnd(ix, iy, ax, ay, level, m_p.seed, stream);
                        // the roulette between the two levels: on average the level's weight
                        Real keep = 1;
                        const Real r = rnd.next();
                        if (!exhaustive) {
                            if (eps > Real(0)) {
                                const Real lo = std::max(Real(0), r - eps), hi = std::min(Real(1), r + eps);
                                const Real t = std::clamp((weight - lo) / (hi - lo), Real(0), Real(1));
                                keep = t * t * (Real(3) - Real(2) * t);
                            }
                            else keep = r < weight ? Real(1) : Real(0);
                            if (keep == Real(0)) continue;
                        }
                        // the facet's normal: a point of its angular cell, a facet only inside the disc
                        const Real ux = (static_cast<Real>(ax) + rnd.next()) * invResA, uy = (static_cast<Real>(ay) + rnd.next()) * invResA;
                        if (math::sqr(ux - Real(0.5)) + math::sqr(uy - Real(0.5)) >= Real(0.25)) continue;
                        const Real dax = q.xa.x - ux, day = q.xa.y - uy;
                        const Real ea = Real(0.5) * (dax * dax + day * day) * invVarA;
                        if (ea > Real(40)) continue;
                        // its place: a point of its spatial cell
                        const Real dx = q.gx - (static_cast<Real>(ix) + rnd.next()) * cell, dy = q.gy - (static_cast<Real>(iy) + rnd.next()) * cell;
                        const Real es = Real(0.5) * (q.ixx * dx * dx + Real(2) * q.ixy * dx * dy + q.iyy * dy * dy);
                        if (es > Real(40)) continue;
                        sum += keep * normA * q.norm * std::exp(-ea - es);
                    }
        if (exhaustive) return sum;
        // the other facets by their expectation: the lobe's mass over the unit square minus the block
        const Real bx0 = static_cast<Real>(a0x) * invResA, bx1 = static_cast<Real>(a1x + 1) * invResA;
        const Real by0 = static_cast<Real>(a0y) * invResA, by1 = static_cast<Real>(a1y + 1) * invResA;
        const Real whole = GaussMass(0, 1, q.xa.x, q.sigmaA) * GaussMass(0, 1, q.xa.y, q.sigmaA);
        const Real block = a1x >= a0x && a1y >= a0y ? GaussMass(bx0, bx1, q.xa.x, q.sigmaA) * GaussMass(by0, by1, q.xa.y, q.sigmaA) : Real(0);
        return sum + weight * std::max(whole - block, Real(0));
    }

    namespace {
        /// Fills the query's spatial part; false for a footprint without area.
        template <class Q>
        bool SpatialQuery(Q& q, const Point2& uv, const GlintFootprint& fp, Real density) {
            const Real res = std::sqrt(density);
            q.gx = uv.x * res;
            q.gy = uv.y * res;
            q.sxx = fp.xx * density;
            q.sxy = fp.xy * density;
            q.syy = fp.yy * density;
            const Real tr = q.sxx + q.syy;
            if (!(tr > Real(0))) return false;
            // a footprint seen edge-on has no width: give it a little, so that the filter stays a function
            Real det = q.sxx * q.syy - q.sxy * q.sxy;
            if (det < Real(1e-8) * tr * tr) {
                const Real add = Real(1e-4) * tr;
                q.sxx += add;
                q.syy += add;
                det = q.sxx * q.syy - q.sxy * q.sxy;
            }
            q.ixx = q.syy / det;
            q.iyy = q.sxx / det;
            q.ixy = -q.sxy / det;
            q.norm = Real(1) / (Real(2) * constants::PI * std::sqrt(det));
            q.hx = GlintNDF::kFilterReach * std::sqrt(q.sxx);
            q.hy = GlintNDF::kFilterReach * std::sqrt(q.syy);
            const Real t2 = q.sxx + q.syy;
            q.sigmaMax = std::sqrt(Real(0.5) * (t2 + std::sqrt(std::max(t2 * t2 - Real(4) * det, Real(0)))));
            return true;
        }
    }   // namespace

    Real GlintNDF::ratio(const Vector3& wh, const Point2& uv, const GlintFootprint& fp, Real weight, std::uint32_t stream) const {
        if (!(wh.z > Real(0)) || !(weight > Real(0))) return 0;
        Query q;
        if (!SpatialQuery(q, uv, fp, m_p.density)) return weight;
        const Real lambda = std::log2(kCellsPerSigma * q.sigmaMax);
        if (!(lambda < Real(kMaxLevel))) return weight;
        q.xa = toDisc(wh);
        q.sigmaA = m_p.microRoughness * std::sqrt(jacobian(wh));
        if (lambda <= Real(0)) return levelSum(q, 0, weight, false, stream);
        const int l0 = static_cast<int>(std::floor(lambda));
        const Real w1 = lambda - static_cast<Real>(l0);
        Real sum = levelSum(q, l0, weight * (Real(1) - w1), false, stream);
        if (w1 > Real(0)) sum += levelSum(q, l0 + 1, weight * w1, false, stream);
        return sum;
    }

    // ---------------------------------------------------------------------------------------------------------------------------
    // Importance sampling (docs/GLINTS.md, section 3). The paper's procedure has no density that can be written down exactly, so
    // this is another one with the same aim, to put samples inside the facets' lobes:
    //   1. a level, with the probabilities of the blend between the two levels;
    //   2. an angular cell of that level, uniformly over the unit square (probability 4^-level);
    //   3. a facet of that cell among the spatial cells under the pixel, in proportion to its weight c_i (the pixel filter at the
    //      facet, times the roulette that keeps it);
    //   4. a point of the disc from the facet's proposal: its Gaussian lobe cut to the square window of half width
    //      h = min(4 sigma, half an angular cell) around the facet (uniform in the window where it is narrower than sigma).
    // The window keeps the density exact and cheap: a point can only come from facets of the 2 x 2 angular cells nearest to it,
    //      p(x) = Sum_levels w_l 4^-l Sum_cells (1 / C_cell) Sum_i c_i q_i(x),   C_cell = Sum_i c_i,
    // the same cells the evaluation visits. An empty cell or a point outside the disc gives no sample: the density integrates to
    // less than 1, which an estimator that divides by it does not mind. The far facets' expectation is not sampled here (the
    // smooth distribution's sampling covers it).
    // ---------------------------------------------------------------------------------------------------------------------------

    namespace {
        struct Facet {
            Real c = 0;        ///< the facet's weight under the pixel (the filter's Gaussian without its constant, times the roulette)
            Real ux = 0, uy = 0;
        };

        /// The proposal of a facet at (ux, uy): its standard deviation and its window.
        struct Proposal {
            Real sigma, h;
            bool uniform;
            Real density(Real dx, Real dy) const {
                if (std::abs(dx) > h || std::abs(dy) > h) return 0;
                if (uniform) return Real(1) / (Real(4) * h * h);
                const Real m = std::erf(h / (sigma * Real(1.4142135623730951)));
                return std::exp(Real(-0.5) * (dx * dx + dy * dy) / (sigma * sigma)) / (Real(2) * constants::PI * sigma * sigma * m * m);
            }
        };
    }   // namespace

    Real GlintNDF::jacobianAtDisc(const Point2& u) const {
        const Real px = Real(2) * (u.x - Real(0.5)), py = Real(2) * (u.y - Real(0.5));
        const Real r2 = std::min(px * px + py * py, Real(1));
        const Real k2 = Real(2) - r2;
        // |M^-1 s|^3 / (8 ax ay) for the hemisphere normal s = (px k, py k, 1 - r2)
        const Real n2 = (math::sqr(m_p.alphaX * px) + math::sqr(m_p.alphaY * py)) * k2 + math::sqr(Real(1) - r2);
        return n2 * std::sqrt(n2) / (Real(8) * m_p.alphaX * m_p.alphaY);
    }

    namespace {
        /// The facet of one cell (spatial ix, iy; angular ax, ay), as the evaluation generates it; false where there is none under the
        /// pixel (dropped by the roulette, outside the disc, beyond the filter).
        template <class Q>
        bool MakeFacet(const Q& q, std::int64_t ix, std::int64_t iy, std::int64_t ax, std::int64_t ay, int level, Real weight, Real eps,
            std::uint32_t seed, std::uint32_t stream, Facet& f) {
            const Real cell = std::ldexp(Real(1), level), invResA = std::ldexp(Real(1), -level);
            FacetRandom rnd(ix, iy, ax, ay, level, seed, stream);
            const Real r = rnd.next();
            Real keep;
            if (eps > Real(0)) {
                const Real lo = std::max(Real(0), r - eps), hi = std::min(Real(1), r + eps);
                const Real t = std::clamp((weight - lo) / (hi - lo), Real(0), Real(1));
                keep = t * t * (Real(3) - Real(2) * t);
            }
            else keep = r < weight ? Real(1) : Real(0);
            if (keep == Real(0)) return false;
            f.ux = (static_cast<Real>(ax) + rnd.next()) * invResA;
            f.uy = (static_cast<Real>(ay) + rnd.next()) * invResA;
            if (math::sqr(f.ux - Real(0.5)) + math::sqr(f.uy - Real(0.5)) >= Real(0.25)) return false;
            const Real dx = q.gx - (static_cast<Real>(ix) + rnd.next()) * cell, dy = q.gy - (static_cast<Real>(iy) + rnd.next()) * cell;
            const Real es = Real(0.5) * (q.ixx * dx * dx + Real(2) * q.ixy * dx * dy + q.iyy * dy * dy);
            if (es > Real(40)) return false;
            f.c = keep * std::exp(-es);
            return true;
        }
    }   // namespace

    Real GlintNDF::levelPdf(const Query& q, int level, Real weight, std::uint32_t stream, const Point2& x) const {
        const Real cell = std::ldexp(Real(1), level);
        const std::int64_t resA = std::int64_t(1) << level;
        const Real invResA = Real(1) / static_cast<Real>(resA);
        const std::int64_t a0x = static_cast<std::int64_t>(std::floor(x.x * static_cast<Real>(resA) - Real(0.5)));
        const std::int64_t a0y = static_cast<std::int64_t>(std::floor(x.y * static_cast<Real>(resA) - Real(0.5)));
        const std::int64_t i0x = static_cast<std::int64_t>(std::floor((q.gx - q.hx) / cell)), i1x = static_cast<std::int64_t>(std::floor((q.gx + q.hx) / cell));
        const std::int64_t i0y = static_cast<std::int64_t>(std::floor((q.gy - q.hy) / cell)), i1y = static_cast<std::int64_t>(std::floor((q.gy + q.hy) / cell));
        Real p = 0;
        for (std::int64_t ay = std::max<std::int64_t>(a0y, 0); ay <= std::min(a0y + 1, resA - 1); ++ay)
            for (std::int64_t ax = std::max<std::int64_t>(a0x, 0); ax <= std::min(a0x + 1, resA - 1); ++ax) {
                Real total = 0, here = 0;
                for (std::int64_t iy = i0y; iy <= i1y; ++iy)
                    for (std::int64_t ix = i0x; ix <= i1x; ++ix) {
                        Facet f;
                        if (!MakeFacet(q, ix, iy, ax, ay, level, weight, m_p.rouletteSoftness, m_p.seed, stream, f)) continue;
                        total += f.c;
                        const Real sigma = m_p.microRoughness * std::sqrt(jacobianAtDisc(Point2(f.ux, f.uy)));
                        const Real h = std::min(Real(4) * sigma, Real(0.5) * invResA);
                        here += f.c * Proposal{ sigma, h, h < sigma }.density(x.x - f.ux, x.y - f.uy);
                    }
                if (total > Real(0)) p += here / total;
            }
        return p * invResA * invResA;
    }

    bool GlintNDF::levelSample(const Query& q, int level, Real weight, std::uint32_t stream, const std::function<Real()>& rnd, Point2& x) const {
        const Real cell = std::ldexp(Real(1), level);
        const std::int64_t resA = std::int64_t(1) << level;
        const Real invResA = Real(1) / static_cast<Real>(resA);
        const std::int64_t ax = std::min(static_cast<std::int64_t>(rnd() * static_cast<Real>(resA)), resA - 1);
        const std::int64_t ay = std::min(static_cast<std::int64_t>(rnd() * static_cast<Real>(resA)), resA - 1);
        const std::int64_t i0x = static_cast<std::int64_t>(std::floor((q.gx - q.hx) / cell)), i1x = static_cast<std::int64_t>(std::floor((q.gx + q.hx) / cell));
        const std::int64_t i0y = static_cast<std::int64_t>(std::floor((q.gy - q.hy) / cell)), i1y = static_cast<std::int64_t>(std::floor((q.gy + q.hy) / cell));
        // the cell's facets under the pixel: one of them in proportion to its weight (a reservoir of one)
        Real total = 0;
        Facet chosen;
        const Real pick = rnd();
        // two passes keep the choice a plain inversion of the cumulative weights
        for (int pass = 0; pass < 2; ++pass) {
            Real run = 0;
            bool done = false;
            for (std::int64_t iy = i0y; iy <= i1y && !done; ++iy)
                for (std::int64_t ix = i0x; ix <= i1x && !done; ++ix) {
                    Facet f;
                    if (!MakeFacet(q, ix, iy, ax, ay, level, weight, m_p.rouletteSoftness, m_p.seed, stream, f)) continue;
                    if (pass == 0) total += f.c;
                    else {
                        run += f.c;
                        chosen = f;
                        if (run >= pick * total) done = true;
                    }
                }
            if (pass == 0 && !(total > Real(0))) return false;
        }
        const Real sigma = m_p.microRoughness * std::sqrt(jacobianAtDisc(Point2(chosen.ux, chosen.uy)));
        const Real h = std::min(Real(4) * sigma, Real(0.5) * invResA);
        Real dx, dy;
        if (h < sigma) {
            dx = (Real(2) * rnd() - Real(1)) * h;
            dy = (Real(2) * rnd() - Real(1)) * h;
        }
        else {
            // the Gaussian, drawn again while it falls outside the window (at least 46 % fall inside)
            int tries = 0;
            do {
                const Real r = sigma * std::sqrt(Real(-2) * std::log(std::max(Real(1) - rnd(), Real(1e-300))));
                const Real phi = Real(2) * constants::PI * rnd();
                dx = r * std::cos(phi);
                dy = r * std::sin(phi);
            } while ((std::abs(dx) > h || std::abs(dy) > h) && ++tries < 200);
            if (tries >= 200) return false;
        }
        x = Point2(chosen.ux + dx, chosen.uy + dy);
        return math::sqr(x.x - Real(0.5)) + math::sqr(x.y - Real(0.5)) < Real(0.25);
    }

    bool GlintNDF::samplable(const GlintFootprint& fp) const {
        const Real l = level(fp);
        return std::isfinite(l) && l < kSampleLevel;
    }

    bool GlintNDF::sample(const Point2& uv, const GlintFootprint& fp, Real weight, std::uint32_t stream, const std::function<Real()>& rnd, Vector3& wh) const {
        Query q;
        if (!(weight > Real(0)) || !SpatialQuery(q, uv, fp, m_p.density)) return false;
        const Real lambda = std::log2(kCellsPerSigma * q.sigmaMax);
        if (!(lambda < Real(kMaxLevel))) return false;
        Point2 x;
        if (lambda <= Real(0)) {
            if (!levelSample(q, 0, weight, stream, rnd, x)) return false;
        }
        else {
            const int l0 = static_cast<int>(std::floor(lambda));
            const Real w1 = lambda - static_cast<Real>(l0);
            const bool upper = rnd() < w1;
            if (!levelSample(q, upper ? l0 + 1 : l0, weight * (upper ? w1 : Real(1) - w1), stream, rnd, x)) return false;
        }
        wh = fromDisc(x);
        return wh.z > Real(0);
    }

    Real GlintNDF::pdf(const Vector3& wh, const Point2& uv, const GlintFootprint& fp, Real weight, std::uint32_t stream) const {
        if (!(wh.z > Real(0)) || !(weight > Real(0))) return 0;
        Query q;
        if (!SpatialQuery(q, uv, fp, m_p.density)) return 0;
        const Real lambda = std::log2(kCellsPerSigma * q.sigmaMax);
        if (!(lambda < Real(kMaxLevel))) return 0;
        const Point2 x = toDisc(wh);
        Real p;
        if (lambda <= Real(0)) p = levelPdf(q, 0, weight, stream, x);
        else {
            const int l0 = static_cast<int>(std::floor(lambda));
            const Real w1 = lambda - static_cast<Real>(l0);
            p = (Real(1) - w1) * levelPdf(q, l0, weight * (Real(1) - w1), stream, x);
            if (w1 > Real(0)) p += w1 * levelPdf(q, l0 + 1, weight * w1, stream, x);
        }
        return p * jacobian(wh);   // from area in the disc to solid angle of the half vector
    }

    Real GlintNDF::ratioExhaustive(const Vector3& wh, const Point2& uv, const GlintFootprint& fp, int level) const {
        if (!(wh.z > Real(0))) return 0;
        Query q;
        if (!SpatialQuery(q, uv, fp, m_p.density)) return 1;
        q.xa = toDisc(wh);
        q.sigmaA = m_p.microRoughness * std::sqrt(jacobian(wh));
        return levelSum(q, std::clamp(level, 0, 10), 1, true, 0);
    }

} // namespace rayt
