#pragma once

/**
 * @file Camera.hpp
 * @brief Physically-based camera model with Depth of Field (DoF).
 * * Implements a thin-lens camera system. Supports adjustable Field of View (FOV),
 * focal distance, and aperture size to simulate realistic optical effects
 * like defocus blur.
 */

#include "Core/Core.hpp"
#include "Core/Ray.hpp"
#include "Core/Math.hpp"
#include "Core/Sampling.hpp"

namespace rayt {

    /**
     * @brief A physically-based camera.
     * * The camera defines the transformation from film coordinates to world-space rays.
     * It uses an orthonormal basis (u, v, w) to orient the view and supports
     * lens sampling for depth-of-field effects.
     */
    class Camera {
    public:
        /**
         * @brief Constructs a camera with full control over orientation and optics.
         * @param lookFrom   The position of the camera in world space.
         * @param lookAt     The point the camera is looking at.
         * @param vUp        The "up" vector for the world (used to stabilize the camera).
         * @param vfov       Vertical Field of View in degrees.
         * @param aspect     Aspect ratio (Width / Height).
         * @param aperture   Lens diameter (0.0 results in a sharp pinhole camera).
         * @param focusDist  The distance from the lens to the plane of perfect focus.
         */
        Camera(Point3 lookFrom,
            Point3 lookAt,
            Vector3 vUp,
            Real vfov,
            Real aspect,
            Real aperture,
            Real focusDist) {

            m_lensRadius = aperture / 2.0;

            // Convert vertical FOV from degrees to radians and compute half-height
            Real theta = rayt::math::toRadians(vfov);
            Real h = std::tan(theta / 2.0);

            // Viewport dimensions in the focal plane
            Real viewportHeight = 2.0 * h;
            Real viewportWidth = aspect * viewportHeight;

            // Construct the Camera Coordinate System (Right-Handed)
            // m_w: Points directly away from the target (View direction = -m_w)
            m_w = glm::normalize(lookFrom - lookAt);

            // m_u: "Right" vector (perpendicular to world-up and view direction)
            m_u = glm::normalize(glm::cross(vUp, m_w));

            // m_v: "Up" vector within the camera's local frame
            m_v = glm::cross(m_w, m_u);

            m_origin = lookFrom;

            // Calculate the view vectors scaled by focus distance to reach the focal plane
            m_horizontal = focusDist * viewportWidth * m_u;
            m_vertical = focusDist * viewportHeight * m_v;

            // Compute the lower-left corner of the projected image plane in world space
            m_lowerLeftCorner = m_origin
                - m_horizontal / 2.0
                - m_vertical / 2.0
                - focusDist * m_w;
        }

        /**
         * @brief Generates a ray for a specific pixel coordinate.
         * * @param s      Normalized horizontal coordinate on the film [0, 1].
         * @param t      Normalized vertical coordinate on the film [0, 1].
         * @param uLens  2D random sample for lens/aperture sampling (for DoF).
         * @return Ray   A world-space ray originating from the lens toward the focal plane.
         */
        Ray getRay(Real s, Real t, const Point2& uLens) const {

            Vector3 offset(0.0);

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
            return Ray(rayOrigin, rayDirection);
        }

    private:
        Point3 m_origin;
        Point3 m_lowerLeftCorner;
        Vector3 m_horizontal;
        Vector3 m_vertical;

        // Internal Orthonormal Basis
        Vector3 m_u, m_v, m_w;

        Real m_lensRadius;
    };

} // namespace rayt