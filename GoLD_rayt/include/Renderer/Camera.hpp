#pragma once

#include "Core/Core.hpp"
#include "Core/Ray.hpp"
#include "Core/Utils.hpp"

namespace rayt {

    // A physically-based camera model.
    // Supports:
    // 1. Adjustable Field of View (FOV)
    // 2. Depth of Field (Defocus Blur) via aperture size
    // 3. Focus Distance setting
    class Camera {
    public:
        // vfov: Vertical Field of View in degrees
        // aspect: Width / Height
        // aperture: Lens diameter (0.0 = Pinhole camera, sharp everywhere)
        // focusDist: Distance to the plane of perfect focus
        Camera(Point3 lookFrom,
            Point3 lookAt,
            Vector3 vUp,
            Real vfov,
            Real aspect,
            Real aperture,
            Real focusDist) {

            m_lensRadius = aperture / 2.0;

            // Convert FOV to radians
            Real theta = Utils::toRadians(vfov);
            Real h = std::tan(theta / 2.0);

            Real viewportHeight = 2.0 * h;
            Real viewportWidth = aspect * viewportHeight;

            // Build an Orthonormal Basis (u, v, w) for the camera
            // w points opposite to the view direction (Right-Handed System)
            m_w = glm::normalize(lookFrom - lookAt);
            m_u = glm::normalize(glm::cross(vUp, m_w));
            m_v = glm::cross(m_w, m_u);

            m_origin = lookFrom;

            // Calculate the "horizontal" and "vertical" vectors of the projection plane.
            // These are scaled by focusDist to project rays to the correct focal plane.
            m_horizontal = focusDist * viewportWidth * m_u;
            m_vertical = focusDist * viewportHeight * m_v;

            // Calculate the lower-left corner of the image plane
            m_lowerLeftCorner = m_origin
                - m_horizontal / 2.0
                - m_vertical / 2.0
                - focusDist * m_w;
        }

        // Generates a ray for a given film coordinate (s, t).
        // s, t: Normalized coordinates [0, 1] on the film.
        Ray getRay(Real s, Real t) const {
            // For Depth of Field:
            // Sample a random point on the lens disk if aperture > 0.
            Vector3 rd = m_lensRadius * Utils::randomInUnitDisk();
            Vector3 offset = m_u * rd.x + m_v * rd.y;

            // Determine the target point on the focus plane
            Point3 target = m_lowerLeftCorner + s * m_horizontal + t * m_vertical;

            // The ray starts from the lens offset (not pure center) and goes to the target.
            return Ray(m_origin + offset, glm::normalize(target - (m_origin + offset)));
        }

    private:
        Point3 m_origin;
        Point3 m_lowerLeftCorner;
        Vector3 m_horizontal;
        Vector3 m_vertical;

        // Camera Basis Vectors
        Vector3 m_u, m_v, m_w;

        Real m_lensRadius;
    };

} // namespace rayt