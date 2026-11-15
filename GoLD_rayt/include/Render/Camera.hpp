/**
 * @file Renderer/Camera.hpp
 * @brief Physically-based camera model with Depth of Field (DoF).
 * * Implements a thin-lens camera system. Supports adjustable Field of View (FOV),
 * focal distance, and aperture size to simulate realistic optical effects
 * like defocus blur.
 */

#pragma once

#include <algorithm>
#include <cmath>

#include "Core/Core.hpp"
#include "Core/Ray.hpp"
#include "Core/Math.hpp"
#include "Core/Sampling.hpp"

namespace rayt {

    struct CameraSample {
        // Film sample in normalized [0,1]^2 (pixel + jitter / NDC)
        Point2 pFilm{ 0.5, 0.5 };

        // Lens sample in [0,1)^2
        Point2 pLens{ 0.5, 0.5 };

        // Time sample in [0,1)
        Real time{ 0.0 };

        // Wavelength of the path (nm). Only the chromatic aberrations (Camera::setAberrations) depend on it;
        // the reference wavelength is the one the lens is focused and scaled for.
        Real lambdaNm{ 550.0 };
    };

    /// The wavelength at which the lens has no chromatic aberration (nm).
    constexpr Real kReferenceWavelengthNm = Real(550);

    /**
     * @brief A physically-based camera.
     * * The camera defines the transformation from film coordinates to world-space rays.
     * It uses an orthonormal basis (u, v, w) to orient the view and supports
     * lens sampling for depth-of-field effects.
     */
    class Camera {
    public:
        /**
         * @brief Constructs    a camera with full control over orientation and optics.
         * @param lookFrom      The position of the camera in world space.
         * @param lookAt        The point the camera is looking at.
         * @param vUp           The "up" vector for the world (used to stabilize the camera).
         * @param vfov          Vertical Field of View in degrees.
         * @param aspect        Aspect ratio (Width / Height).
         * @param aperture      Lens diameter (0.0 results in a sharp pinhole camera).
         * @param focusDist     The distance from the lens to the plane of perfect focus.
         * @param shutterOpen   Shutter open time
         * @param shutterClose  Shutter close time
         * @param medium        Medium that the camera is inside (optional)
         */
        Camera(const Point3 lookFrom,
            const Point3 lookAt,
            const Vector3 vUp,
            Real vfovDeg,
            Real aspect,
            Real aperture,
            Real focusDist,
            Real shutterOpen = 0.0,
            Real shutterClose = 1.0,
            const Medium* medium = nullptr)
            : m_origin(lookFrom)
            , m_lensRadius(aperture* Real(0.5))
            , m_focusDist(std::max(focusDist, Real(0)))
            , m_shutterOpen(shutterOpen)
            , m_shutterClose(std::max(shutterClose, shutterOpen))
            , m_medium(medium) {

            m_lensRadius = aperture / 2.0;

            // Convert vertical FOV from degrees to radians and compute half-height
            // Real theta = rayt::math::toRadians(vfovDeg);
            // Real h = std::tan(theta / 2.0);

            // Viewport dimensions in the focal plane
            // Real viewportHeight = 2.0 * h;
            // Real viewportWidth = aspect * viewportHeight;

            // Construct the Camera Coordinate System (Right-Handed)
            // m_w: Points directly away from the target (View direction = -m_w)
            m_w = glm::normalize(lookFrom - lookAt);

            // m_u: "Right" vector (perpendicular to world-up and view direction)
            m_u = glm::normalize(glm::cross(vUp, m_w));

            // m_v: "Up" vector within the camera's local frame
            m_v = glm::cross(m_w, m_u);

            // --- Perspective viewport on the focal plane ---
            const Real theta = rayt::math::toRadians(vfovDeg);
            const Real halfHeight = std::tan(theta * Real(0.5));
            const Real viewportHeight = Real(2) * halfHeight;
            const Real viewportWidth = aspect * viewportHeight;

            // Calculate the view vectors scaled by focus distance to reach the focal plane
            m_horizontal = focusDist * viewportWidth * m_u;
            m_vertical = focusDist * viewportHeight * m_v;

            // The film rectangle at unit depth, for the radial distortion (its radius is measured in half diagonals)
            m_viewportW = viewportWidth;
            m_viewportH = viewportHeight;
            m_invHalfDiag2 = Real(4) / (viewportWidth * viewportWidth + viewportHeight * viewportHeight);

            // Compute the lower-left corner of the projected image plane in world space
            m_lowerLeftCorner = m_origin
                - m_horizontal * Real(0.5)
                - m_vertical * Real(0.5)
                - focusDist * m_w;
        }

        /**
         * @brief Generate a camera ray from a sample (PBRT-style).
         * * Uses thin-lens model if lensRadius > 0.
         */
        Ray generateRay(const CameraSample& cs) const {
            // --- time (motion blur) ---
            const Real time = math::lerp(m_shutterOpen, m_shutterClose, math::saturate(cs.time));

            // --- film sample (NDC) ---
            const Real s = math::saturate(cs.pFilm.x);
            const Real t = math::saturate(cs.pFilm.y);

            // --- where on the film this ray comes from: radial distortion and lateral chromatic aberration ---
            // The pixel at (s, t) looks along the direction of the ideal (pinhole) film point (s', t') = centre + (s, t - centre) * g,
            // g = (1 + k1 r^2 + k2 r^4) * (1 + lateral * (lambda - 550 nm) / 100 nm), r = distance from the centre in half diagonals.
            Real sd = s, td = t;
            if (m_distorted || m_lateralCA != Real(0)) {
                const Real dx = s - Real(0.5), dy = t - Real(0.5);
                const Real xw = dx * m_viewportW, yh = dy * m_viewportH;
                const Real r2 = (xw * xw + yh * yh) * m_invHalfDiag2;
                Real g = Real(1) + m_k1 * r2 + m_k2 * r2 * r2;
                g *= Real(1) + m_lateralCA * (cs.lambdaNm - kReferenceWavelengthNm) * Real(0.01);
                sd = Real(0.5) + dx * g;
                td = Real(0.5) + dy * g;
            }

            // --- point on focal plane ---
            Point3 target = m_lowerLeftCorner + sd * m_horizontal + td * m_vertical;
            // Axial chromatic aberration: the plane of focus of this wavelength lies nearer or farther (the chief ray is not bent)
            if (m_axialCA != Real(0) && m_lensRadius > Real(0))
                target = m_origin + (target - m_origin) * (focusDistanceFor(cs.lambdaNm) / m_focusDist);

            // --- lens sample (DoF) ---
            Vector3 offset(0);
            if (m_lensRadius > Real(0)) {
                // A round aperture, or a regular polygon (blades >= 3): out-of-focus highlights take its shape.
                const Vector3 rd = (m_blades >= 3) ? sampleBlades(cs.pLens) : m_lensRadius * sampling::UniformSampleDisk(cs.pLens);
                offset = m_u * rd.x + m_v * rd.y;
            }

            const Point3 rayOrigin = m_origin + offset;
            // The BSDFs (reflect / refract / GGX) assume unit-length directions, and the camera ray
            // becomes `wo` at the first hit, so it must be normalized (its raw length is ~ focusDist).
            const Vector3 rayDir = glm::normalize(target - rayOrigin);

            Ray r(rayOrigin, rayDir, constants::RAY_EPSILON, m_medium);
            r.time = time;
            return r;
        }

        /**
         * @brief Generates a ray for a specific pixel coordinate.
         * * @param s      Normalized horizontal coordinate on the film [0, 1].
         * @param t      Normalized vertical coordinate on the film [0, 1].
         * @param uLens  2D random sample for lens/aperture sampling (for DoF).
         * @return Ray   A world-space ray originating from the lens toward the focal plane.
         */
        Ray getRay(Real s, Real t, const Point2& uLens, Real timeSample = 0.0, Real lambdaNm = kReferenceWavelengthNm) const {

            CameraSample cs;
            cs.pFilm = Point2(s, t);
            cs.pLens = uLens;
            cs.time = timeSample;
            cs.lambdaNm = lambdaNm;
            return generateRay(cs);

            /*Vector3 offset(0.0);

            // Defocus Blur (Depth of Field) calculation
            if (m_lensRadius > 0.0) {
                // Map the 2D random sample to a point on a circular lens/aperture
                Vector3 rd = m_lensRadius * sampling::UniformSampleDisk(uLens);

                // Offset the ray origin within the disk of the lens
                offset = m_u * rd.x + m_v * rd.y;
            }

            // Target point on the plane of perfect focus
            Point3 target = m_lowerLeftCorner + s * m_horizontal + t * m_vertical;

            // The ray originates from the sampled point on the lens
            Point3 rayOrigin = m_origin + offset;

            // The direction points from the lens offset to the target on the focal plane
            Vector3 rayDirection = glm::normalize(target - rayOrigin);

            // Pass 'medium' (PBRT style)
            return Ray(rayOrigin, rayDirection);*/
        }

        /**
         * @brief Projects a world point onto the film of the pinhole view (the lens is ignored): the inverse of generateRay().
         *
         * @param s, t     film coordinates in [0, 1) (as in generateRay)
         * @param cosCam   cosine of the angle between the direction to the point and the view axis
         * @param dist     distance from the camera position to the point
         * @return false if the point is behind the camera or outside the film.
         */
        bool project(const Point3& p, Real& s, Real& t, Real& cosCam, Real& dist) const {
            const Vector3 d = p - m_origin;
            const Real z = glm::dot(d, -m_w);   // depth along the view axis
            if (z <= Real(0)) return false;
            dist = glm::length(d);
            cosCam = z / dist;
            const Point3 hit = m_origin + d * (m_focusDist / z);   // the point on the focal plane
            const Vector3 q = hit - m_lowerLeftCorner;
            s = glm::dot(q, m_horizontal) / glm::dot(m_horizontal, m_horizontal);
            t = glm::dot(q, m_vertical) / glm::dot(m_vertical, m_vertical);
            return s >= Real(0) && s < Real(1) && t >= Real(0) && t < Real(1);
        }

        /**
         * @brief Gives the aperture the shape of a regular polygon with `blades` corners (3 or more; fewer: a circle again).
         *
         * The polygon has the area of the round aperture of the same f-number, so the size of the blur stays the same and only its
         * shape changes. `rotation` (radians) turns it about the view axis.
         */
        void setApertureShape(int blades, Real rotation = Real(0)) {
            m_blades = (blades >= 3) ? blades : 0;
            m_bladeRotation = rotation;
            if (m_blades >= 3) {
                const Real n = Real(m_blades);
                // area of a regular n-gon of circumradius R: (n / 2) R^2 sin(2 pi / n); equal to pi r^2
                m_bladeRadius = m_lensRadius * std::sqrt(constants::PI / (Real(0.5) * n * std::sin(constants::TWO_PI / n)));
            }
        }

        /**
         * @brief The lens aberrations (all 0: an ideal lens).
         *
         * @param k1, k2      radial distortion, r' = r (1 + k1 r^2 + k2 r^4) with r in half diagonals of the film: the pixel at radius r sees
         *                    what the ideal picture shows at r'. k1 > 0 is barrel distortion (straight lines bow outwards), k1 < 0 pincushion.
         * @param lateralCA   lateral chromatic aberration: the picture of wavelength lambda is scaled by 1 + lateralCA (lambda - 550 nm) / 100 nm
         *                    about the centre (0.004: red and blue differ by 0.8 % in size, 8 pixels at the corner of a 2000 px wide picture).
         *                    > 0: the red picture is the smaller one.
         * @param axialCA     axial chromatic aberration: the focal length of a thin lens varies as 1 + axialCA (lambda - 550 nm) / 100 nm
         *                    (0.01 is an uncorrected singlet of crown glass, 0.0005 an achromat), so every wavelength has its own plane of focus.
         *                    Needs the focal length (metres) and an aperture; the film stays where the lens is focused for 550 nm.
         *
         * The chromatic ones need the wavelength of the path (CameraSample::lambdaNm) and so the spectral build.
         */
        void setAberrations(Real k1, Real k2, Real lateralCA, Real axialCA, Real focalLength) {
            m_k1 = k1;
            m_k2 = k2;
            m_distorted = (k1 != Real(0) || k2 != Real(0));
            m_lateralCA = lateralCA;
            m_axialCA = axialCA;
            m_focalLength = focalLength;
            if (m_axialCA != Real(0)) {
                // the film distance that focuses the reference wavelength at m_focusDist (thin lens equation)
                if (focalLength > Real(0) && m_focusDist > focalLength) m_filmDist = focalLength * m_focusDist / (m_focusDist - focalLength);
                else m_axialCA = Real(0);   // focused at or inside the focal length: no meaningful film distance
            }
        }

        /// True when the rays depend on the wavelength (lateral or axial chromatic aberration is on).
        bool chromatic() const { return m_lateralCA != Real(0) || (m_axialCA != Real(0) && m_lensRadius > Real(0)); }

        /// True when the distortion maps the film onto itself with a fold (r' stops growing with r somewhere on the film).
        bool distortionFolds() const {
            for (int i = 0; i <= 64; ++i) {
                const Real r2 = Real(i) / Real(64) * Real(i) / Real(64);
                if (Real(1) + Real(3) * m_k1 * r2 + Real(5) * m_k2 * r2 * r2 <= Real(0)) return true;
            }
            return false;
        }

        /// Area of the film rectangle on the focal plane (world units), and the focus distance: the geometry of the pixel footprint.
        Real focalPlaneArea() const { return glm::length(m_horizontal) * glm::length(m_vertical); }
        Real focusDistance() const { return m_focusDist; }
        Real lensRadius() const { return m_lensRadius; }

        // Accessors
        const Point3& origin() const { return m_origin; }
        const Medium* medium() const { return m_medium; }

    private:
        /// The distance of the plane of focus for light of this wavelength: the film is fixed, the focal length changes with the wavelength.
        Real focusDistanceFor(Real lambdaNm) const {
            const Real f = m_focalLength * (Real(1) + m_axialCA * (lambdaNm - kReferenceWavelengthNm) * Real(0.01));
            if (!(f > Real(0))) return m_focusDist;
            const Real gap = m_filmDist - f;
            return (gap > f * Real(1e-6)) ? f * m_filmDist / gap : Real(1e6);   // film at or inside the focal length: focused at infinity
        }

        /// A uniform point of the polygonal aperture: pick one of the n triangles (centre, corner k, corner k+1), then a point in it.
        Vector3 sampleBlades(const Point2& u) const {
            const Real scaled = u.x * Real(m_blades);
            const int k = std::min(static_cast<int>(scaled), m_blades - 1);
            const Real a = scaled - Real(k);                    // position along the outer edge, uniform again
            const Real r = std::sqrt(u.y);                      // distance factor towards the edge
            const Real step = constants::TWO_PI / Real(m_blades);
            const Real a0 = m_bladeRotation + step * Real(k);
            const Real a1 = a0 + step;
            const Real x = m_bladeRadius * r * ((Real(1) - a) * std::cos(a0) + a * std::cos(a1));
            const Real y = m_bladeRadius * r * ((Real(1) - a) * std::sin(a0) + a * std::sin(a1));
            return Vector3(x, y, Real(0));
        }

        Point3 m_origin;
        Point3 m_lowerLeftCorner;
        Vector3 m_horizontal;
        Vector3 m_vertical;

        // Internal Orthonormal Basis
        Vector3 m_u, m_v, m_w;

        // Thin lens
        Real m_lensRadius = 0.0;
        Real m_focusDist = 1.0;
        int m_blades = 0;           // 0: a round aperture
        Real m_bladeRotation = 0.0;
        Real m_bladeRadius = 0.0;   // circumradius of the polygon

        // Aberrations (setAberrations)
        Real m_viewportW = 1.0, m_viewportH = 1.0;   // the film at unit depth
        Real m_invHalfDiag2 = 1.0;                   // 1 / (half diagonal)^2 of it
        Real m_k1 = 0.0, m_k2 = 0.0;
        bool m_distorted = false;
        Real m_lateralCA = 0.0;
        Real m_axialCA = 0.0;
        Real m_focalLength = 0.0;                    // metres
        Real m_filmDist = 0.0;                       // lens to film, metres

        // Motion blur (shutter interval)
        Real m_shutterOpen = 0.0;
        Real m_shutterClose = 1.0;

        // Medium around camera
        const Medium* m_medium = nullptr;
    };

} // namespace rayt