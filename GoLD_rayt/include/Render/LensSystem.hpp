/**
 * @file include/Render/LensSystem.hpp
 * @brief A real lens: the spherical, flat and aspheric surfaces and the aperture stop of a published prescription, traced ray by ray
 *        (--lens). docs/CAMERA.md
 *
 * A camera ray starts on the film, crosses every surface from the back to the front (Snell's law at each, a dispersive glass bending each
 * wavelength its own way) and is lost when it meets a rim, the stop or total internal reflection. Distortion, field curvature, vignetting,
 * cat's-eye bokeh, focus breathing and chromatic aberration are not parameters here: they follow from the prescription.
 *
 * An aspheric surface has the sag of the lens patents, z(r) = c r^2 / (1 + sqrt(1 - (1 + k) c^2 r^2)) + Sum_i A_i r^i (c = 1 / radius,
 * k the conic constant, any powers i >= 2, odd ones too), measured from the vertex toward the film. A ray meets it where Newton's method
 * on the ray parameter, started on the vertex plane, converges; the normal comes from the derivative of the sag.
 *
 * Lens frame: the optical axis is z, pointing at the object; the first surface's vertex is at z = 0 and the film behind the lens at
 * negative z. A negative spacing (a stop behind the first lens's vertex, as phone lenses have it) puts the next surface in front of the
 * previous one; the trace then follows the ray's line backward to it, as sequential lens-design programs do.
 * Lengths in metres (the files are in millimetres). Kolb, Mitchell and Hanrahan 1995; pbrt's RealisticCamera.
 */

#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Core/Types.hpp"
#include "Kernel/KernelScene.hpp"

namespace rayt {

    class LensSystem {
    public:
        struct Surface {
            Real radius = 0;            ///< > 0: centre of curvature on the film side; 0: the aperture stop (a plane), or a flat surface
            bool flat = false;          ///< radius 0 is a flat refracting surface (a filter, a cover glass), not the stop
            Real z = 0;                 ///< vertex on the axis
            Real apertureRadius = 0;    ///< clear semi-diameter
            Real nd = 1;                ///< index of the medium behind the surface (toward the film) at 587.6 nm; 1 = air
            Real abbe = 0;              ///< Abbe number V_d of that medium (0: not dispersive)
            Real cauchyA = 0, cauchyB = 0;
            Real conic = 0;             ///< conic constant k of an aspheric surface (0: the base is a sphere)
            std::vector<std::pair<int, Real>> poly;   ///< aspheric terms (power i, A_i) of the sag, in millimetres as the patents give them
            std::vector<Real> coef;     ///< the same terms by power (index i), filled by the LensSystem constructor

            bool isStop() const { return radius == Real(0) && !flat; }
            bool aspheric() const { return conic != Real(0) || !poly.empty(); }
            /// The sag z(r) (metres, positive toward the film) and its slope dz/dr at the distance r (metres) from the axis; false where
            /// the conic is undefined (1 - (1 + k) c^2 r^2 < 0).
            bool sag(Real r, Real& z, Real& slope) const;
        };

        /**
         * @brief Reads a prescription: one surface per line from the front (object side) to the back, in millimetres:
         *   radius  thickness  n_d  clear-diameter  [Abbe number]
         * `thickness` is the distance to the next surface (the last one's is not used, and it may be negative); n_d = 0 or 1 is air; a
         * radius of 0 is the aperture stop, `inf` a flat refracting surface. '#' starts a comment. (The column layout of pbrt's lens files,
         * plus the optional Abbe number.) Further lines, for the lenses of phone-camera patents:
         *   asphere k K A4 a4 A6 a6 ...   the surface of the line above is aspheric: conic constant K and the terms A_i (any i >= 2, in the
         *                                 patents' units of millimetres; terms left out are 0)
         *   fnumber N                     the stop's diameter is set so that the lens has this f-number (its diameter column may be `auto`)
         *   image-height H                half the image diagonal (mm) the lens is designed for: the film becomes the picture's rectangle of
         *                                 that half diagonal (instead of 36 x 24 mm), and every `auto` clear diameter is set to pass the full
         *                                 bundle of the stop to every point of that image (the stop must then be the first surface)
         *   half-field DEG                instead of image-height: the image height is where the chief ray of this field angle meets
         *                                 the focal plane (for patents that give the field of view but not the image height)
         */
        static std::shared_ptr<LensSystem> load(const std::string& path, std::string& error);

        /// Surfaces given directly (front to back, metres; z of the first = 0); the checks build lenses this way.
        explicit LensSystem(std::vector<Surface> surfaces);

        // ---- paraxial properties at the d line (587.6 nm), from rays a nanometre off the axis ----
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
        void setStopShape(int blades, Real rotation);

        /// Moves the film so that a point `distance` in front of the front principal plane is in focus (Gauss's lens equation with the
        /// principal planes). False, and focus at infinity, when the distance is not beyond the focal length.
        bool focus(Real distance);
        Real filmZ() const { return m_filmZ; }
        /// After prepare(): the projected solid angle (sr) of the rays that reach the film centre through the lens, at the d line: the film's
        /// irradiance there is the scene's radiance times this (pi / (4 N_w^2) for a lens that meets the sine condition). Physical exposure.
        Real centreIrradiance() const { return m_centreIrradiance; }

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

        /// This lens as the kernel renderer models it (Kernel/Camera.hpp): the surfaces (with the index behind each), the exit pupil's
        /// bounds per band (minX, minY, maxX, maxY, empty), the stop's edge normals, and the film, focus and stop in `c`. Call after
        /// prepare(). False, with the reason, for an aspheric surface or dispersive glass.
        bool kernelParams(kernel::CameraParams<Real>& c, std::vector<kernel::LensSurface<Real>>& surfaces, std::vector<Real>& pupil,
            std::vector<Real>& bladeNormals, std::string& why) const;
        const std::vector<Surface>& surfaces() const { return m_surfaces; }
        /// Half the image diagonal the lens is designed for (metres; `image-height` of the file), or 0: the 35 mm film is used.
        Real imageHeight() const { return m_imageHeight; }

        /**
         * @brief Sets every clear semi-diameter at or above 1 m (the loader's "auto") to the largest height (plus 0.5 %) at which a ray passes that
         *        surface on its way from the stop, which must be the first surface, to an image point within `imageHeight` of the axis
         *        (the lens focused at infinity, the d line): the stop's full bundle reaches the whole image, no surface vignettes it.
         *        Also sets imageHeight(). False (and the reason) when the stop is not the first surface or no ray gets through.
         */
        bool setAutoApertures(Real imageHeight, std::string& error);
        /// The angle (radians) between the axis and the object-side ray through the centre of the first-surface stop that lands at the
        /// height h of the focal plane (the chief ray of that image point, the d line); 0 when it cannot be found.
        Real fieldAngleForImageHeight(Real h) const;
        /// The inverse: where the chief ray of the field angle theta (radians) meets the focal plane (metres); 0 when it does not.
        Real imageHeightForFieldAngle(Real theta) const;

        /// Traces a ray in place, from the film toward the object / from the object toward the film. False if it is blocked.
        bool traceFromFilm(Point3& o, Vector3& d, Real nm) const;
        bool traceFromObject(Point3& o, Vector3& d, Real nm) const;

        // ---- ghosts (Render/LensFlare): light reflected twice inside the lens ----

        /// The anti-reflection coating of the glass-air surfaces: none (bare glass, about 4 % per surface) or a single layer of MgF2
        /// (n = 1.38) a quarter wave thick at 550 nm, the classic coating (about 1.3 % at 550 nm, more toward blue and red).
        enum class Coating { Uncoated, SingleLayer };

        /// One ray of a ghost, traced without the clear apertures: where it meets the film and what would clip it, so that the flare
        /// renderer can test the stop and the rims between its rays.
        struct GhostSample {
            Real x = 0, y = 0;          ///< where it meets the film plane (lens frame, metres)
            Real throughput = 0;        ///< the two reflectances times the transmittances (1 - R) of every other crossing
            int stopPasses = 0;         ///< how often it crossed the stop's plane (1 to 3), and where
            Real stopX[3] = { 0, 0, 0 }, stopY[3] = { 0, 0, 0 };
            Real rim = 0;               ///< the largest height / clear radius at the lens surfaces it crossed (above 1: clipped)
        };

        /// The reflectances of every surface's interface (from the front and from behind), per wavelength band, over the cosine of
        /// incidence: the flare's rays meet the surfaces millions of times, and a coated surface's reflectance costs a thin-film formula.
        struct ReflectanceTable {
            static constexpr int kCos = 512;   ///< intervals of the cosine (linear interpolation)
            int bands = 0;
            std::vector<Real> r;               ///< [surface][side: 0 from the front, 1 from behind][band][kCos + 1]
            Real at(int surface, bool fromFront, int band, Real cosI) const {
                const Real x = std::clamp(cosI, Real(0), Real(1)) * Real(kCos);
                const int i = std::min(static_cast<int>(x), kCos - 1);
                const Real f = x - Real(i);
                const size_t base = ((static_cast<size_t>(surface) * 2 + (fromFront ? 0 : 1)) * static_cast<size_t>(bands) + static_cast<size_t>(band))
                    * static_cast<size_t>(kCos + 1) + static_cast<size_t>(i);
                return r[base] + (r[base + 1] - r[base]) * f;
            }
        };
        /// The table for the wavelengths bandNm[0 .. bands) (each band's own indices of refraction in the spectral build).
        ReflectanceTable makeReflectanceTable(int bands, const Real* bandNm, Coating coating) const;

        /**
         * @brief The ghost of the surfaces i < j: light from the object side crosses the surfaces to j, is reflected there back toward
         *        the front, crosses back to i, is reflected there toward the film, and crosses the rest to the film (the film as focused).
         *        The reflectances and transmittances are the unpolarized ones of each interface (Fresnel, or the coated stack of a glass-air
         *        surface, Materials/ThinFilm.hpp). The path is traced at `nm`; its throughput (the two reflectances times 1 - R of every
         *        other crossing) at the wavelength `nm` goes to g.throughput, or, with a table, at each of its bands to bandT (the RGB build
         *        traces one path for all bands: its glass bends every wavelength alike while the coating still colours the reflections).
         *        False if the ray misses a surface or is totally reflected where it should cross.
         */
        bool traceGhost(int i, int j, Point3 o, Vector3 d, Real nm, Coating coating, GhostSample& g, const ReflectanceTable* table = nullptr,
            Real* bandT = nullptr) const;
        /// Whether a point of the stop's plane (lens frame) passes the stop as it is set (f-number, blades).
        bool insideStopNow(Real x, Real y) const;
        /// The unpolarized reflectance of an interface for light from the medium of index nIn into nOut at the cosine cosI: Fresnel, or
        /// through the coating when one side is air and the other glass.
        static Real interfaceReflectance(Real cosI, Real nIn, Real nOut, Real nm, Coating coating);

    private:
        struct Bounds { Real minX = 0, minY = 0, maxX = 0, maxY = 0; bool empty = true; };
        struct StopShape { Real scale = 1; int blades = 0; Real rotation = 0; };   ///< the stop a trace tests against

        Real indexAfter(int i, Real nm) const;   ///< the medium behind surface i (toward the film); 1 before the front surface (i = -1)
        /// Moves the ray to surface i and refracts it. `backward`: the surface lies behind the ray's origin (a negative spacing), so the
        /// intersection may be at a negative ray parameter.
        bool crossSurface(int i, Point3& o, Vector3& d, Real nFrom, Real nTo, const StopShape& stop, bool backward) const;
        /// Moves the ray to the (non-stop) surface s: the point, its distance from the axis and the normal there (facing either way).
        /// `clip`: false where the point lies outside the clear aperture.
        bool hitSurface(const Surface& s, Point3& o, const Vector3& d, bool backward, bool clip, Vector3& normal, Real& rho) const;
        /// The ray parameter where the ray meets the aspheric (or flat) surface i, by Newton's method from its vertex plane; false if it
        /// does not converge or leaves the conic's domain.
        bool intersectAspheric(const Surface& s, const Point3& o, const Vector3& d, Real& t) const;
        bool chiefRayHeight(Real theta, Real& y) const;
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
        Real m_imageHeight = 0;         ///< half the design's image diagonal (0: not given)
        int m_stopIndex = -1;           ///< the stop surface (-1: none)
        std::vector<Point2> m_bladeNormals;   ///< the polygonal stop's edge normals (setStopShape)
        Real m_apothemPerRadius = 1;          ///< its apothem per unit of the round stop's radius (the same area)
    };

} // namespace rayt
