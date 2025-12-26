#pragma once

#include <algorithm>
#include <cmath>
#include <vector>
#include <memory>

#include "Core/Types.hpp"
#include "Core/Constants.hpp"
#include "Core/Image.hpp"
#include "Core/Distribution2D.hpp"

namespace rayt {

    /**
     * @class EnvMap
     * @brief A class representing an Image-Based Lighting (IBL) environment map.
     *
     * This class handles an equirectangular environment map (latitude-longitude format).
     * It provides functionality to sample radiance from a given direction using bilinear interpolation.
     */
    class EnvMap {
    public:
        /**
         * @brief Constructs an environment map from an image.
         * @param image The source image (equirectangular projection).
         */
        explicit EnvMap(Image image) noexcept
            : m_img(std::move(image)) {
            if (m_img.isValid()) {
                buildDistribution();
            }
        }

        /**
         * @brief Evaluates the environment radiance from a specific direction.
         *
         * Converts the direction vector to UV coordinates and samples the texture
         * using bilinear interpolation.
         *
         * @param dir The normalized direction vector pointing to the environment.
         * @return The radiance color (Vector3) at that direction.
         * Returns black if the image is invalid.
         */
        Vector3 eval(const Vector3& dir) const {
            if (!m_img.isValid()) {
                return Vector3(0.0);
            }

            Real u, v;

            // Ensure direction is normalized before calculating UV
            dirToUV(glm::normalize(dir), u, v);
            return sampleBilinear(static_cast<float>(u), static_cast<float>(v));
        }

        /**
         * @brief Importance sample the environment map (for NEE).
         * @param u 2D uniform random in [0,1)^2
         * @param wi [out] sampled direction (world)
         * @param pdfW [out] pdf w.r.t. solid angle (sr^-1)
         * @return Le radiance from the sampled direction
         */
        Vector3 sample(const Point2& u, Vector3& wi, Real& pdfW) const {
            pdfW = Real(0);
            if (!m_img.isValid() || !m_dist) return Vector3(0.0);

            // 1) Sample UV from 2D distribution (pdf in uv-domain)
            Point2 uvImg;
            float pdfUV_f = 0.0f;
            m_dist->sampleContinuous(u, uvImg, pdfUV_f); // uv in [0,1), pdfUV is density on [0,1]^2
            Real pdfUV = Real(pdfUV_f);

            if (pdfUV <= Real(0)) return Vector3(0.0);

            Real uSph = uvImg.x;
            Real vSph = Real(1) - uvImg.y;

            // 2) Convert UV -> direction
            wi = uvToDir(uSph, vSph);

            // 3) Evaluate radiance from that direction
            Vector3 Le = sampleBilinear((float)uSph, (float)vSph);

            // 4) Convert pdf_uv -> pdf_omega
            // For lat-long map: dω = 2π^2 sinθ du dv
            Real sinTheta = sinThetaFromV(vSph);
            if (sinTheta <= Real(0)) {
                pdfW = Real(0);
                return Vector3(0.0);
            }

            pdfW = pdfUV / (Real(2) * constants::PI * constants::PI * sinTheta);
            return Le;
        }

        /**
         * @brief Pdf of sampling direction wi by EnvMap::sample (w.r.t solid angle).
         */
        Real pdf(const Vector3& wi) const {
            if (!m_img.isValid() || !m_dist) return Real(0);

            //------------------------------------------------------------------------
            // EnvMap が期待するv球面とDistribution2D が扱う v（画像)が上下反転
            // していたので直した今後改善が必要
            //------------------------------------------------------------------------

            // Direction -> UV
            Real uSph, vSph;
            dirToUV(glm::normalize(wi), uSph, vSph);

            Real vImg = Real(1) - vSph;

            // pdf in uv-domain
            float pdfUV_f = m_dist->pdf(Point2((float)uSph, (float)vImg));
            Real pdfUV = Real(pdfUV_f);
            if (pdfUV <= Real(0)) return Real(0);

            // Convert to solid angle pdf
            Real sinTheta = sinThetaFromV(vSph);
            if (sinTheta <= Real(0)) return Real(0);

            return pdfUV / (Real(2) * constants::PI * constants::PI * sinTheta);
        }

    private:
        Image m_img;

        // 2D importance distribution (built from luminance * sinθ)
        std::unique_ptr<Distribution2D> m_dist;

        static Real luminance(const Vector3& rgb) {
            // Rec.709 / sRGB luminance weights (common choice)
            return Real(0.2126) * rgb.x + Real(0.7152) * rgb.y + Real(0.0722) * rgb.z;
        }

        // v in [0,1], where v = 1 - theta/pi  (per your dirToUV)
        static Real sinThetaFromV(Real v) {
            // theta = pi * (1 - v)
            Real theta = constants::PI * (Real(1) - std::clamp(v, Real(0), Real(1)));
            return std::sin(theta);
        }

        // Inverse of dirToUV consistent with your mapping:
        // u = (phi + pi) / (2pi), v = 1 - theta/pi
        static Vector3 uvToDir(Real u, Real v) {
            Real phi = (std::clamp(u, Real(0), Real(1)) * Real(2) * constants::PI) - constants::PI;
            Real theta = constants::PI * (Real(1) - std::clamp(v, Real(0), Real(1)));

            Real sinTheta = std::sin(theta);
            Real cosTheta = std::cos(theta);
            Real cosPhi = std::cos(phi);
            Real sinPhi = std::sin(phi);

            // This matches dirToUV's atan2(z,x) and y-up:
            // x = sinθ cosφ, y = cosθ, z = sinθ sinφ
            return Vector3(sinTheta * cosPhi, cosTheta, sinTheta * sinPhi);
        }

        /**
         * @brief Converts a 3D direction vector to spherical UV coordinates.
         *
         * Assumes a Y-up coordinate system.
         * Maps (0, 1, 0) to V=1.0 and (0, -1, 0) to V=0.0.
         * Maps atan2(z, x) to U.
         *
         * @param d The normalized direction vector.
         * @param u [out] The resulting U coordinate [0, 1].
         * @param v [out] The resulting V coordinate [0, 1].
         */
        static void dirToUV(const Vector3& d, Real& u, Real& v) {
            const Real theta = std::acos(std::clamp(d.y, Real(-1), Real(1)));
            const Real phi = std::atan2(d.z, d.x);

            u = (phi + constants::PI) / (Real(2) * constants::PI);
            v = Real(1) - (theta / constants::PI);
        }

        /**
         * @brief Fetches a pixel color with boundary handling.
         *
         * Wraps horizontally (longitude) and clamps vertically (latitude).
         *
         * @param x The x-coordinate (column).
         * @param y The y-coordinate (row).
         * @return The pixel color.
         */
        Vector3 texel(int x, int y) const {
            const int w = m_img.width();
            const int h = m_img.height();

            x = (x % w + w) % w;              // wrap horizontally
            y = std::clamp(y, 0, h - 1);      // clamp vertically

            return m_img.at(x, y);
        }

        /**
         * @brief Samples the texture using bilinear interpolation.
         *
         * @param u The U coordinate [0, 1].
         * @param v The V coordinate [0, 1].
         * @return The interpolated color.
         */
        Vector3 sampleBilinear(float u, float v) const {
            u = u - std::floor(u);            // wrap [0,1)
            v = std::clamp(v, 0.0f, 1.0f);

            // Map UV to pixel coordinates (center aligned)
            // Note: (1.0 - v) flips V because image origin (0,0) is top-left,
            // but spherical V=1 usually corresponds to the top (North Pole).
            const float x = u * m_img.width() - 0.5f;
            const float y = (1.0f - v) * m_img.height() - 0.5f;

            const int x0 = static_cast<int>(std::floor(x));
            const int y0 = static_cast<int>(std::floor(y));
            const int x1 = x0 + 1;
            const int y1 = y0 + 1;

            const float tx = x - std::floor(x);
            const float ty = y - std::floor(y);

            const Vector3 c00 = texel(x0, y0);
            const Vector3 c10 = texel(x1, y0);
            const Vector3 c01 = texel(x0, y1);
            const Vector3 c11 = texel(x1, y1);

            // Bilinear interpolation
            return glm::mix(
                glm::mix(c00, c10, tx),
                glm::mix(c01, c11, tx),
                ty
            );
        }

        void buildDistribution() {
            const int w = m_img.width();
            const int h = m_img.height();
            if (w <= 0 || h <= 0) return;

            // Build weights = luminance * sin(theta)
            std::vector<float> weights;
            weights.resize(size_t(w) * size_t(h));

            for (int y = 0; y < h; ++y) {
                // v coordinate at pixel center (match sampleBilinear convention)
                // Using pixel centers helps stability.
                Real v = (Real(y) + Real(0.5)) / Real(h);
                // Note: sampleBilinear uses (1 - v) for y mapping; our v here is the "spherical v"
                // consistent with dirToUV: v=1 at north pole.
                // Image row y=0 is top, so spherical v at row y is:
                Real vSph = Real(1) - v;

                Real sinTheta = sinThetaFromV(vSph);
                float sinT = static_cast<float>(sinTheta);

                for (int x = 0; x < w; ++x) {
                    Vector3 rgb = m_img.at(x, y);
                    Real lum = luminance(rgb);
                    // Clamp negative luminance (just in case)
                    float wt = static_cast<float>(std::max(lum, Real(0))) * sinT;
                    weights[size_t(y) * size_t(w) + size_t(x)] = wt;
                }
            }

            m_dist = std::make_unique<Distribution2D>(weights.data(), w, h);
        }

    };

} // namespace rayt