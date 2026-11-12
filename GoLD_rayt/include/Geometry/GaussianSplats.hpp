/**
 * @file Geometry/GaussianSplats.hpp
 * @brief A cloud of 3D Gaussians (3D Gaussian Splatting) traced by rays, as emissive and absorbing particles.
 *
 * Each particle is an anisotropic Gaussian with an opacity and a view-dependent color (spherical harmonics up to degree 3), as
 * the training of 3D Gaussian Splatting (Kerbl et al. 2023) writes it to point_cloud.ply. Along a ray a particle responds with
 *     alpha = opacity * exp(-1/2 |x(t*) - mu|^2_Sigma)
 * at the point t* of the ray that is closest to the center in the metric of the covariance (the ray-Gaussian response of
 * 3D Gaussian Ray Tracing, Moenne-Loccoz et al. 2024). The particles a ray meets are sorted by t* and composited front to back:
 *     L = sum_i T_i alpha_i c_i,   T_i = prod_{j<i} (1 - alpha_j),
 * and the ray goes on with the transmittance T = prod (1 - alpha_i). There is no scattering and no reflectance: c_i is the radiance
 * the capture saw, so the cloud is a backdrop and a light that paths see when they fly through it (also mirrored in a metal), and
 * it shadows what is behind it by T. A BVH over the particles' boxes (where alpha can still reach alphaMin) finds a ray's particles.
 *
 * Colors: the stored coefficients give the pixel values of the training images (sRGB-encoded as a rule); `linearize` decodes them to
 * linear radiance and `intensity` scales them. The spherical harmonics are evaluated in the cloud's own frame (the direction is
 * turned back by the scene rotation).
 *
 * The implementation (loading, the BVH, the ray response and the compositing) is in src/Geometry/GaussianSplats.cpp.
 */

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "Core/Types.hpp"
#include "Core/AABB.hpp"
#include "Core/Ray.hpp"

namespace rayt {

    /// One particle as read from the file: in the cloud's own coordinates, with the raw (activated) parameters.
    struct SplatSource {
        Point3 mean = Point3(0.0);
        Vector3 scale = Vector3(1.0);       ///< standard deviations along the particle's axes (exp of the stored value)
        glm::dvec4 rot = glm::dvec4(1, 0, 0, 0);   ///< unit quaternion (w, x, y, z): the particle's axes in the cloud's frame
        Real opacity = 1;                   ///< 0..1 (sigmoid of the stored value)
        std::vector<float> sh;              ///< 3 * (degree + 1)^2 coefficients: f_dc (3), then f_rest channel by channel
    };

    struct SplatParams {
        Real alphaMin = Real(1) / Real(255);   ///< a particle's response below this is ignored (and bounds its box)
        Real intensity = 1;                   ///< multiplies the radiance
        bool linearize = true;                ///< decode the sRGB-encoded colors to linear radiance
        int degree = 3;                       ///< spherical harmonics used (at most what the file has)
        Real scale = 1;                       ///< the scene transform: scale, then rotation, then translation
        glm::dmat3 rotation = glm::dmat3(1.0);
        Vector3 translate = Vector3(0.0);
        bool hasKeep = false, hasRemove = false;   ///< particles kept only inside `keep` / dropped inside `remove` (world space, by center)
        AABB keep, remove;
    };

    class GaussianSplats {
    public:
        /// Builds the cloud (transform, culling, BVH) from the particles. `degreeInFile` is what the coefficients contain.
        GaussianSplats(const std::vector<SplatSource>& sources, int degreeInFile, const SplatParams& params);

        size_t count() const { return m_particles.size(); }
        int degree() const { return m_degree; }
        const AABB& bounds() const { return m_bounds; }
        size_t numNodes() const { return m_nodes.size(); }

        /**
         * @brief What a ray sees of the cloud between r.tMin and tMax: the composited radiance (linear RGB) and the transmittance left.
         * @param tMax  distance of the surface the ray hits (infinite: nothing), the particles beyond it are not seen.
         */
        void gather(const Ray& r, Real tMax, Vector3& radiance, Real& transmittance) const;

        /// The transmittance alone (shadow rays): the order of the particles does not matter.
        Real transmittance(const Ray& r, Real tMax) const;

        /// The particles with alpha >= alphaMin along the ray (unsorted), by brute force: for the checks.
        void collectBruteForce(const Ray& r, Real tMax, std::vector<std::pair<Real, Real>>& tAlpha) const;

        /// The view-dependent color (linear RGB radiance) of a particle for a world direction.
        Vector3 colorForDirection(size_t index, const Vector3& dirWorld) const;

        /**
         * @brief Reads point_cloud.ply of 3D Gaussian Splatting: x y z, opacity, scale_0..2, rot_0..3, f_dc_0..2, f_rest_0..N (N = 8, 23 or 44).
         * The stored values are activated here (sigmoid, exp, normalization). `degreeOut` is the spherical harmonics degree found.
         */
        static bool loadPly(const std::string& filename, std::vector<SplatSource>& out, int& degreeOut, std::string* outError = nullptr);

        /// Rotation matrix of a unit quaternion (w, x, y, z); its columns are the rotated axes.
        static glm::dmat3 quatToMatrix(const glm::dvec4& q);

        /// The 3D Gaussian Splatting color for spherical harmonics coefficients (3 per band term) and a unit direction: max(0, SH(dir) + 0.5).
        static Vector3 evalSh(const float* sh, int degree, const Vector3& d);

    private:
        static constexpr Real kStopT = Real(1e-3);   // a path with this little transmittance left stops compositing (as 3DGRT does)
        static constexpr int kLeafSize = 4;

        struct Particle {
            Point3 mean;
            float rotT[3][3];           // R^T (rows = the particle's axes): x_local = R^T (x - mean)
            glm::vec<3, float> invScale;
            float opacity;
        };
        struct Hit { Real t; Real alpha; int index; };
        struct Node {
            AABB bounds;
            int left = -1, right = -1;   // inner node
            int start = 0, count = 0;    // leaf: particles m_order[start .. start + count)
        };

        static bool inside(const AABB& b, const Point3& p);
        bool respond(const Particle& p, const Ray& r, Real tMax, Real& tOut, Real& alphaOut) const;
        Vector3 color(int index, const Vector3& dirLocal) const;
        void collect(const Ray& r, Real tMax, std::vector<Hit>& hits) const;
        static bool slab(const AABB& b, const Ray& r, const Vector3& invDir, Real tMax);
        void buildBvh(const std::vector<AABB>& boxes);
        int buildNode(const std::vector<AABB>& boxes, int start, int end);

        SplatParams m_params;
        int m_degree = 0, m_coeffs = 1;
        glm::dmat3 m_rotT = glm::dmat3(1.0);   // world direction -> the cloud's frame (for the harmonics)
        std::vector<Particle> m_particles;
        std::vector<float> m_sh;                // 3 * m_coeffs per particle
        std::vector<Node> m_nodes;
        std::vector<int> m_order;
        AABB m_bounds;
    };

} // namespace rayt
