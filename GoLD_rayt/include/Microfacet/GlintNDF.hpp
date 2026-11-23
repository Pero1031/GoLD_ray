/**
 * @file Microfacet/GlintNDF.hpp
 * @brief A glinty normal distribution: a GGX distribution made of a finite number of facets (docs/GLINTS.md).
 *
 * After Kemppinen, Paulin, Thonat, Thiery, Lehtinen and Boubekeur 2025, "Evaluating and Sampling Glinty NDFs in Constant Time"
 * (ACM TOG 44(6)), written from the paper's text.
 *
 * A facet is a point in four dimensions: two for its place on the surface (texture space) and two for its normal. The normals live in
 * a disc U inscribed in the unit square, where the map T(w) = L(M w / |M w|) puts the micro-normal w: M = diag(1/ax, 1/ay, 1) turns
 * the GGX distribution into the normals of a hemisphere, L(x, y, z) = (x, y) / sqrt(1 + z) is Lambert's equal-area projection. One
 * facet sits at a random place in every cell of a four-dimensional grid; the cell's index seeds the random numbers, so nothing is
 * stored. A facet answers the half vector with a Gaussian lobe in U and the pixel with the pixel's filter (a Gaussian in texture
 * space), and the distribution is the smooth D(w) times the sum of the facets:
 *
 *     D*(x, w) = D(w) x Sum_i  N_a(T(w) - mu_i)  f(x_i - x)        (in grid units, one facet per unit of area x area)
 *
 * so that it averages to D(w) over the surface. The grid adapts to the pixel: level l has spatial cells of side 2^l and 2^l x 2^l
 * angular cells, the same number of facets at every level, and the level is chosen from the footprint so that it covers a few
 * spatial cells. Within them, the facets of the angular cells nearest to T(w) are summed one by one and the others are replaced by
 * their expectation (the lobe's integral over the unit square minus that block: error functions). Between two levels each facet is
 * kept or dropped by a roulette with the level's weight (not blended), so the glints keep their brightness.
 *
 * Where this differs from the paper is listed in docs/GLINTS.md (the normalization fixed by the average, the lobe's width from the
 * map's true Jacobian, a spatial neighbourhood that covers four standard deviations of the filter, an exact roulette by default).
 * The model does not conserve energy exactly (the paper says so): lobes near the disc's rim lose or gain a little.
 */

#pragma once

#include <cstdint>

#include "Core/Types.hpp"

namespace rayt {

    /// The pixel filter on the surface: the covariance of a Gaussian in texture space (uv units squared).
    struct GlintFootprint {
        Real xx = 0, xy = 0, yy = 0;
    };

    struct GlintParams {
        Real alphaX = Real(0.3);          ///< GGX roughness of the whole distribution
        Real alphaY = Real(0.3);
        Real microRoughness = Real(0.01); ///< the width of one facet's lobe: its standard deviation in radians of the half vector
        Real density = Real(1e6);         ///< facets per unit area of texture space
        int angularCells = 2;             ///< angular cells per side summed facet by facet (the paper's k = 4 is 2 x 2)
        Real rouletteSoftness = 0;        ///< the paper's epsilon (0.1 there, for animation); 0 is the exact roulette
        std::uint32_t seed = 0;           ///< another arrangement of the facets
    };

    class GlintNDF {
    public:
        explicit GlintNDF(const GlintParams& params);

        const GlintParams& params() const { return m_p; }

        /// The smooth GGX distribution the facets add up to.
        Real smoothD(const Vector3& wh) const;

        /// The glinty distribution at the texture point uv for the half vector wh (local frame, z the normal), seen through a pixel
        /// whose filter has the covariance `footprint` on the surface.
        Real D(const Vector3& wh, const Point2& uv, const GlintFootprint& footprint) const { return smoothD(wh) * ratio(wh, uv, footprint); }

        /// D* / D: 1 on average over the surface. With `weight` below 1 every facet is kept by a roulette of that probability (and the
        /// expectation scaled): `weight` on average, for a surface that mixes several arrangements (`stream` tells them apart).
        Real ratio(const Vector3& wh, const Point2& uv, const GlintFootprint& footprint, Real weight = Real(1), std::uint32_t stream = 0) const;

        /// For the checks: every facet of one level summed one by one (no roulette, no expectation); D* / D of that level alone.
        Real ratioExhaustive(const Vector3& wh, const Point2& uv, const GlintFootprint& footprint, int level) const;

        /// The level (real-valued) the footprint asks for; negative when the pixel is smaller than a cell of level 0.
        Real level(const GlintFootprint& footprint) const;

        // --- the map T between micro-normals and the disc of radius 1/2 centred in the unit square -------------------------------

        Point2 toDisc(const Vector3& wh) const;
        Vector3 fromDisc(const Point2& u) const;
        /// Area in the disc per solid angle of the micro-normal: det M / (8 |M w|^3).
        Real jacobian(const Vector3& wh) const;

        static constexpr int kMaxLevel = 20;        ///< beyond it a pixel holds so many facets that the distribution is the smooth one
        static constexpr Real kCellsPerSigma = 4;   ///< a cell of the chosen level is 2 to 4 standard deviations of the filter wide
        static constexpr Real kFilterReach = 4;     ///< the spatial cells within this many standard deviations are visited

    private:
        struct Query;
        Real levelSum(const Query& q, int level, Real weight, bool exhaustive, std::uint32_t stream) const;

        GlintParams m_p;
    };

} // namespace rayt
