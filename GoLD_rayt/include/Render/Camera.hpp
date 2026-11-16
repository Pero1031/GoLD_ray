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
            Real shutterClose = 0.0,
            const Medium* medium = nullptr)
            : m_origin(lookFrom)
            , m_lensRadius(aperture* Real(0.5))
            , m_focusDist(std::max(focusDist, Real(0)))
            , m_vUp(vUp)
            , m_shutterOpen(shutterOpen)
            , m_shutterClose(std::max(shutterClose, shutterOpen))
            , m_medium(medium) {

            m_lensRadius = aperture / 2.0;

            // --- Perspective viewport on the focal plane ---
            const Real theta = rayt::math::toRadians(vfovDeg);
            const Real halfHeight = std::tan(theta * Real(0.5));
            const Real viewportHeight = Real(2) * halfHeight;
            const Real viewportWidth = aspect * viewportHeight;

            // The film rectangle at unit depth, for the radial distortion (its radius is measured in half diagonals)
            m_viewportW = viewportWidth;
            m_viewportH = viewportHeight;
            m_invHalfDiag2 = Real(4) / (viewportWidth * viewportWidth + viewportHeight * viewportHeight);

            // The basis and the image plane at the focus distance
            orient(lookFrom, glm::normalize(lookFrom - lookAt));
        }

        /**
         * @brief The exposure: the shutter is open from time 0 to `seconds` (ray times are in seconds).
         *
         * Every camera ray gets a time in that interval; the moving objects (Scene/MovingObject.hpp) and the moving camera
         * (setMotion) are where they are at that time. 0: no exposure time, nothing blurs and no extra random number is drawn.
         */
        void setShutter(Real seconds) {
            m_shutterOpen = Real(0);
            m_shutterClose = std::max(seconds, Real(0));
        }
        bool hasShutter() const { return m_shutterClose > m_shutterOpen; }
        Real shutterSeconds() const { return m_shutterClose - m_shutterOpen; }

        /**
         * @brief The camera moves during the exposure.
         *
         * @param velocity   metres per second (the view direction stays unless a pan is set)
         * @param yawRate    degrees per second: the view turns to the right (negative: to the left) about the up vector
         * @param pitchRate  degrees per second: the view turns up (negative: down)
         *
         * The pose at time 0 is the one the camera was built with. The focus distance does not change.
         * project() (used by light tracing) keeps the pose at time 0.
         */
        void setMotion(const Vector3& velocity, Real yawRate, Real pitchRate) {
            m_velocity = velocity;
            m_yawRate = yawRate * constants::DEG_TO_RAD;
            m_pitchRate = pitchRate * constants::DEG_TO_RAD;
            m_moving = (glm::dot(velocity, velocity) > Real(0)) || yawRate != Real(0) || pitchRate != Real(0);
        }
        bool moving() const { return m_moving; }

        /**
         * @brief Generate a camera ray from a sample (PBRT-style).
         * * Uses thin-lens model if lensRadius > 0.
         */
        Ray generateRay(const CameraSample& cs, Real* weight = nullptr) const {
            // --- time (motion blur) ---
            const Real time = math::lerp(m_shutterOpen, m_shutterClose, math::saturate(cs.time));

            // A moving camera: the same camera at its pose of that moment (a copy: a few dozen numbers, next to a whole path).
            if (m_moving) {
                Camera at(*this);
                at.m_moving = false;
                at.placeAt(*this, time);
                return at.generateRay(cs, weight);   // the same time again
            }

            // --- the point on the focal plane that the film sample looks at (distortion, chromatic aberrations) ---
            const Point3 target = focalTarget(cs.pFilm.x, cs.pFilm.y, cs.lambdaNm);

            // --- lens sample (DoF) ---
            Vector3 rd(0);
            Vector3 offset(0);
            if (m_lensRadius > Real(0)) {
                // A round aperture, or a regular polygon (blades >= 3): out-of-focus highlights take its shape.
                rd = (m_blades >= 3) ? sampleBlades(cs.pLens) : m_lensRadius * sampling::UniformSampleDisk(cs.pLens);
                offset = m_u * rd.x + m_v * rd.y;
            }

            const Point3 rayOrigin = m_origin + offset;
            // The BSDFs (reflect / refract / GGX) assume unit-length directions, and the camera ray
            // becomes `wo` at the first hit, so it must be normalized (its raw length is ~ focusDist).
            const Vector3 rayDir = glm::normalize(target - rayOrigin);

            // --- vignetting: how much of the light of this ray the camera lets through (0 = blocked) ---
            if (weight) {
                Real w = Real(1);
                if (m_naturalVignetting > Real(0)) {
                    // natural vignetting: the film is lit in proportion to cos^4 of the angle off the axis (relative to the centre)
                    const Real cosTheta = glm::dot(glm::normalize(target - m_origin), -m_w);
                    w *= std::pow(std::max(cosTheta, Real(0)), Real(4) * m_naturalVignetting);
                }
                if (m_barrelRadius > Real(0) && m_lensRadius > Real(0)) {
                    // optical vignetting: the opening of the lens barrel, `m_barrelDist` in front of the stop, cuts the bundle of an oblique
                    // ray; what is left is the overlap of the two circles, and the out-of-focus highlights become cat eyes
                    const Real dz = glm::dot(rayDir, -m_w);
                    if (dz <= Real(0)) w = Real(0);
                    else {
                        const Real k = m_barrelDist / dz;
                        const Real lu = rd.x + glm::dot(rayDir, m_u) * k;
                        const Real lv = rd.y + glm::dot(rayDir, m_v) * k;
                        if (lu * lu + lv * lv > m_barrelRadius * m_barrelRadius) w = Real(0);
                    }
                }
                *weight = w;
            }

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
        Ray getRay(Real s, Real t, const Point2& uLens, Real timeSample = 0.0, Real lambdaNm = kReferenceWavelengthNm, Real* weight = nullptr) const {

            CameraSample cs;
            cs.pFilm = Point2(s, t);
            cs.pLens = uLens;
            cs.time = timeSample;
            cs.lambdaNm = lambdaNm;
            return generateRay(cs, weight);

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
         * @brief The ray through the centre of the lens toward film position (s, t), at time 0 and the reference wavelength.
         *
         * The pinhole view of the camera (the distortion included, no defocus, no motion): what a pixel sees as geometry, e.g. its depth.
         */
        Ray centerRay(Real s, Real t) const {
            const Point3 target = focalTarget(s, t, kReferenceWavelengthNm);
            Ray r(m_origin, glm::normalize(target - m_origin), constants::RAY_EPSILON, m_medium);
            r.time = m_shutterOpen;
            return r;
        }

        /// The direction the camera looks in (unit length; at time 0).
        Vector3 viewAxis() const { return -m_w; }

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

        /**
         * @brief Vignetting (all 0: none).
         *
         * @param barrelRadius  radius (metres) of an opening of the lens barrel `barrelDist` metres in front of the aperture stop: it cuts the light
         *                      of oblique rays, which darkens the corners and turns out-of-focus highlights there into cat eyes. It clips nothing
         *                      once the stop is smaller than it (the lens stopped down). 0: no such opening.
         * @param natural       natural vignetting, the cos^4 law of the angle off the axis: 0 = none, 1 = the full law (a power in between scales the exponent).
         *
         * Reported through the `weight` of generateRay() / getRay(): 0 when the ray is blocked, else the factor of the natural vignetting.
         */
        void setVignetting(Real barrelRadius, Real barrelDist, Real natural) {
            m_barrelRadius = std::max(barrelRadius, Real(0));
            m_barrelDist = std::max(barrelDist, Real(0));
            m_naturalVignetting = std::max(natural, Real(0));
        }

        /// The share of the (round) aperture the barrel leaves open for the film position (s, t) in [0,1]^2, from the chief ray (1 = all of it).
        Real apertureFraction(Real s, Real t) const {
            if (!(m_barrelRadius > Real(0)) || !(m_lensRadius > Real(0))) return Real(1);
            const Point3 target = focalTarget(s, t, kReferenceWavelengthNm);
            const Vector3 chief = glm::normalize(target - m_origin);
            const Real dz = glm::dot(chief, -m_w);
            if (dz <= Real(0)) return Real(0);
            const Real k = m_barrelDist / dz;
            const Real d = std::hypot(glm::dot(chief, m_u) * k, glm::dot(chief, m_v) * k);   // shift of the barrel opening against the stop
            const Real r = m_lensRadius, R = m_barrelRadius;
            if (d >= r + R) return Real(0);
            if (d <= std::fabs(R - r)) return (R >= r) ? Real(1) : (R * R) / (r * r);
            const Real area = r * r * std::acos((d * d + r * r - R * R) / (Real(2) * d * r))
                + R * R * std::acos((d * d + R * R - r * r) / (Real(2) * d * R))
                - Real(0.5) * std::sqrt((-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R));
            return area / (constants::PI * r * r);
        }

        /// True when some ray can be blocked or dimmed (vignetting is on): callers then pass a weight to getRay() and apply it.
        bool vignetted() const { return m_naturalVignetting > Real(0) || (m_barrelRadius > Real(0) && m_lensRadius > Real(0)); }

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
        /// The image's horizontal axis (to the right) in world space: the reference direction of the Stokes vectors (Core/Polarization.hpp).
        const Vector3& right() const { return m_u; }
        const Medium* medium() const { return m_medium; }

    private:
        /// Places the camera at `from`, looking along -w: the orthonormal basis and the image plane at the focus distance.
        void orient(const Point3& from, const Vector3& w) {
            m_origin = from;
            // m_w points directly away from the target (the view direction is -m_w)
            m_w = w;
            // m_u: "right" (perpendicular to the world up and the view direction), m_v: "up" in the camera's frame
            m_u = glm::normalize(glm::cross(m_vUp, m_w));
            m_v = glm::cross(m_w, m_u);
            // the view vectors scaled by the focus distance reach the focal plane
            m_horizontal = m_focusDist * m_viewportW * m_u;
            m_vertical = m_focusDist * m_viewportH * m_v;
            m_lowerLeftCorner = m_origin - m_horizontal * Real(0.5) - m_vertical * Real(0.5) - m_focusDist * m_w;
        }

        /// Puts this camera where `start` (the pose at time 0) is at `time` seconds: moved by its velocity, then turned by the pan rates.
        void placeAt(const Camera& start, Real time) {
            Vector3 d = -start.m_w;   // the view direction
            if (start.m_yawRate != Real(0)) d = rotateAbout(d, glm::normalize(start.m_vUp), -start.m_yawRate * time);   // to the right is clockwise seen from above
            if (start.m_pitchRate != Real(0)) d = rotateAbout(d, start.m_u, start.m_pitchRate * time);
            orient(start.m_origin + start.m_velocity * time, -glm::normalize(d));
        }

        /// Rodrigues' rotation of x about the unit axis k by `angle` (radians).
        static Vector3 rotateAbout(const Vector3& x, const Vector3& k, Real angle) {
            const Real c = std::cos(angle), s = std::sin(angle);
            return x * c + glm::cross(k, x) * s + k * (glm::dot(k, x) * (Real(1) - c));
        }

        /// The point of the focal plane that film position (s, t) looks at: radial distortion, lateral and axial chromatic aberration applied.
        Point3 focalTarget(Real sFilm, Real tFilm, Real lambdaNm) const {
            const Real s = math::saturate(sFilm);
            const Real t = math::saturate(tFilm);

            // The pixel at (s, t) looks along the direction of the ideal (pinhole) film point (s', t') = centre + (s, t - centre) * g,
            // g = (1 + k1 r^2 + k2 r^4) * (1 + lateral * (lambda - 550 nm) / 100 nm), r = distance from the centre in half diagonals.
            Real sd = s, td = t;
            if (m_distorted || m_lateralCA != Real(0)) {
                const Real dx = s - Real(0.5), dy = t - Real(0.5);
                const Real xw = dx * m_viewportW, yh = dy * m_viewportH;
                const Real r2 = (xw * xw + yh * yh) * m_invHalfDiag2;
                Real g = Real(1) + m_k1 * r2 + m_k2 * r2 * r2;
                g *= Real(1) + m_lateralCA * (lambdaNm - kReferenceWavelengthNm) * Real(0.01);
                sd = Real(0.5) + dx * g;
                td = Real(0.5) + dy * g;
            }

            Point3 target = m_lowerLeftCorner + sd * m_horizontal + td * m_vertical;
            // Axial chromatic aberration: the plane of focus of this wavelength lies nearer or farther (the chief ray is not bent)
            if (m_axialCA != Real(0) && m_lensRadius > Real(0))
                target = m_origin + (target - m_origin) * (focusDistanceFor(lambdaNm) / m_focusDist);
            return target;
        }

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

        // Vignetting (setVignetting)
        Real m_barrelRadius = 0.0;
        Real m_barrelDist = 0.0;
        Real m_naturalVignetting = 0.0;

        // The up vector the basis is built from (also for the moving camera)
        Vector3 m_vUp = Vector3(0, 1, 0);

        // Motion blur: the shutter interval in seconds, and the motion of the camera in it (setShutter, setMotion)
        Real m_shutterOpen = 0.0;
        Real m_shutterClose = 0.0;
        bool m_moving = false;
        Vector3 m_velocity = Vector3(0.0);
        Real m_yawRate = 0.0, m_pitchRate = 0.0;   // radians per second

        // Medium around camera
        const Medium* m_medium = nullptr;
    };

} // namespace rayt