/**
 * @file include/Render/LensSystem.hpp
 * @brief A real lens: the spherical surfaces and the aperture stop of a published prescription, traced ray by ray (--lens). docs/CAMERA.md
 *
 * A camera ray starts on the film, crosses every surface from the back to the front (Snell's law at each, a dispersive glass bending each
 * wavelength its own way) and is lost when it meets a rim, the stop or total internal reflection. Distortion, field curvature, vignetting,
 * cat's-eye bokeh, focus breathing and chromatic aberration are not parameters here: they follow from the prescription.
 *
 * Lens frame: the optical axis is z, pointing at the object; the front vertex is at z = 0 and the film behind the lens at negative z.
 * Lengths in metres (the files are in millimetres). Kolb, Mitchell and Hanrahan 1995; pbrt's RealisticCamera.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Core/Types.hpp"

namespace rayt {

    class LensSystem {
    public:
        struct Surface {
            Real radius = 0;            ///< > 0: centre of curvature on the film side; 0: the aperture stop (a plane)
            Real z = 0;                 ///< vertex on the axis
            Real apertureRadius = 0;    ///< clear semi-diameter
            Real nd = 1;                ///< index of the medium behind the surface (toward the film) at 587.6 nm; 1 = air
            Real abbe = 0;              ///< Abbe number V_d of that medium (0: not dispersive)
            Real cauchyA = 0, cauchyB = 0;
        };

        /**
         * @brief Reads a prescription: one surface per line from the front (object side) to the back, in millimetres:
         *   radius  thickness  n_d  clear-diameter  [Abbe number]
         * `thickness` is the distance to the next surface (the last one's is not used); n_d = 0 or 1 is air; a radius of 0 is the aperture
         * stop. '#' starts a comment. (The column layout of pbrt's lens files, plus the optional Abbe number.)
         */
        static std::shared_ptr<LensSystem> load(const std::string& path, std::string& error);

        /// Surfaces given directly (front to back, metres; z of the first = 0); the checks build lenses this way.
        explicit LensSystem(std::vector<Surface> surfaces);

        // ---- paraxial properties at the d line (587.6 nm), from rays a micrometre off the axis ----
        Real focalLength() const { return m_focal; }
        Real frontPrincipalZ() const { return m_frontPrincipal; }
        Real rearPrincipalZ() const { return m_rearPrincipal; }
        Real frontFocalZ() const { return m_frontFocal; }
        Real rearFocalZ() const { return m_rearFocal; }
        Real focalLengthAt(Real nm) const;   ///< the paraxial focal length at another wavelength (spectral build; else the d line's)

        /// Focal length over the diameter of the entrance pupil (the widest bundle parallel to the axis that passes), for the current stop.
        Real fNumber() const;
        Real wideOpenFNumber() const { return m_wideOpen; }
        /// Closes the stop to f-number N (no wider than the lens allows); returns the f-number set.
        Real setFNumber(Real n);
        /// A stop of `blades` straight blades (a regular polygon of the circle's area, turned by `rotation` radians); fewer than 3: round.
        void setStopShape(int blades, Real rotation) { m_blades = blades >= 3 ? blades : 0; m_bladeRotation = rotation; }

        /// Moves the film so that a point `distance` in front of the front principal plane is in focus (Gauss's lens equation with the
        /// principal planes). False, and focus at infinity, when the distance is not beyond the focal length.
        bool focus(Real distance);
        Real filmZ() const { return m_filmZ; }

        /// The film (metres) and the exit pupil table: call after the stop and the focus are set, before rendering.
        void prepare(Real filmWidth, Real filmHeight);
        Real filmWidth() const { return m_filmW; }
        Real filmHeight() const { return m_filmH; }

        /**
         * @brief A camera ray from the film point (x, y) (lens frame) through the point u of the exit pupil's bounds, at `nm`.
         * @param weight  the film's irradiance factor of the sample, cos^4 theta x the bounds' area / (film to rear vertex)^2, divided by
         *                its value at the film centre: 1 on average at the centre, less off the axis (natural and optical vignetting).
         * @return false when the ray is blocked (a rim, the stop, total internal reflection).
         */
        bool generateRay(Real x, Real y, const Point2& u, Real nm, Point3& o, Vector3& d, Real& weight) const;
        /// A ray from (x, y) through the middle of the exit pupil (nearly the chief ray), at the d line: the geometry of a pixel (depth, mask).
        bool centralRay(Real x, Real y, Point3& o, Vector3& d) const;

        /// True when the glass is dispersive and the build is spectral: every wavelength takes its own path.
        bool dispersive() const;
        const std::vector<Surface>& surfaces() const { return m_surfaces; }

        /// Traces a ray in place, from the film toward the object / from the object toward the film. False if it is blocked.
        bool traceFromFilm(Point3& o, Vector3& d, Real nm) const;
        bool traceFromObject(Point3& o, Vector3& d, Real nm) const;

    private:
        struct Bounds { Real minX = 0, minY = 0, maxX = 0, maxY = 0; bool empty = true; };
        struct StopShape { Real scale = 1; int blades = 0; Real rotation = 0; };   ///< the stop a trace tests against

        Real indexAfter(int i, Real nm) const;   ///< the medium behind surface i (toward the film); 1 before the front surface (i = -1)
        bool crossSurface(int i, Point3& o, Vector3& d, Real nFrom, Real nTo, const StopShape& stop) const;
        bool insideStop(const Surface& s, Real x, Real y, const StopShape& stop) const;
        bool traceFromFilmWith(Point3& o, Vector3& d, Real nm, const StopShape& stop) const;
        void paraxial();
        Real entrancePupilRadius() const;
        Bounds pupilBounds(Real r0, Real r1) const;

        std::vector<Surface> m_surfaces;
        Real m_focal = 0, m_frontPrincipal = 0, m_rearPrincipal = 0, m_frontFocal = 0, m_rearFocal = 0;
        Real m_wideOpen = 0;
        Real m_stopScale = 1;
        int m_blades = 0;
        Real m_bladeRotation = 0;
        Real m_filmZ = 0;
        Real m_filmW = 0.036, m_filmH = 0.024;
        std::vector<Bounds> m_pupil;   ///< per band of film radius, for a film point on the +x axis
        Real m_filmRadius = 0;          ///< the half diagonal the table covers
        Real m_centreIrradiance = 1;    ///< the weight's normalisation
    };

} // namespace rayt
