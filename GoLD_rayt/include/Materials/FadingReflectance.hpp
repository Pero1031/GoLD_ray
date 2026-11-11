/**
 * @file Materials/FadingReflectance.hpp
 * @brief A measured reflectance that changes with the hours of exposure, with the exposure varying over the surface.
 *
 * The reflectance at a point is R(lambda, hours(p)) of a measured series (IO/FadingSeries.hpp): the linear blend of the two
 * measured rounds around that time. The time comes from an exposure fraction f(p) in [0, 1]:
 *     hours(p) = unexposed + f(p) * (exposed - unexposed)
 * and f is the product of whichever of these factors are set (each is 1 where it is not set):
 *   - a map     : a gray image over the surface's uv (white = fully exposed)
 *   - a mat     : a border of the uv square that was covered (a mat board, a frame) and so did not fade
 *   - a gradient: exposure 1 at one point, 0 at another, linear in between along their axis (light from one side; works on any shape)
 *   - a noise   : uneven fading, 1 - amount * noise(p) (works on any shape)
 *
 * Where the exposure comes from is the scene author's input: the measurements give the reflectance against time only.
 */

#pragma once

#include <algorithm>
#include <memory>
#include <vector>

#include "Core/Types.hpp"
#include "Core/Interaction.hpp"
#include "Core/Texture.hpp"
#include "Color/SpectralConversion.hpp"
#include "IO/FadingSeries.hpp"
#include "Materials/DensityField.hpp"

namespace rayt {

    /// The exposure fraction over a surface (see the file comment).
    struct ExposureField {
        std::shared_ptr<const Texture2D> map;   ///< gray image; its sRGB value is the fraction

        bool hasMat = false;                    ///< covered border, in uv: left, right, bottom, top widths (the uv square has side 1)
        Real matLeft = 0, matRight = 0, matBottom = 0, matTop = 0;
        Real matSoft = Real(0.002);             ///< width of the soft edge (uv)

        bool hasGradient = false;
        Point3 bright = Point3(0.0);            ///< exposure 1 here ...
        Point3 dark = Point3(0, 1, 0);          ///< ... and 0 here
        Real gradientExponent = 1;

        std::shared_ptr<const DensityField> noise;   ///< uneven fading: factor 1 - noiseAmount * noise(p)
        Real noiseAmount = 0;

        bool any() const { return map || hasMat || hasGradient || noise; }

        static Real smooth(Real x) { x = std::clamp(x, Real(0), Real(1)); return x * x * (Real(3) - Real(2) * x); }

        /// 1 inside the uncovered rectangle, 0 in the covered border, a soft step across each edge.
        Real matAt(Real u, Real v) const {
            const Real s = std::max(matSoft, Real(1e-6));
            auto edge = [&](Real x, Real inset) { return smooth((x - inset) / s + Real(0.5)); };
            return edge(u, matLeft) * edge(Real(1) - u, matRight) * edge(v, matBottom) * edge(Real(1) - v, matTop);
        }

        Real at(const SurfaceInteraction& rec) const {
            Real f = 1;
            if (map) {
                const Real c = map->rgb(rec.uv.x, rec.uv.y).x;   // back to the 8-bit value: what the author painted
                f *= std::clamp(c <= Real(0.0031308) ? c * Real(12.92) : Real(1.055) * std::pow(c, Real(1) / Real(2.4)) - Real(0.055), Real(0), Real(1));
            }
            if (hasMat) f *= matAt(rec.uv.x, rec.uv.y);
            if (hasGradient) {
                const Vector3 axis = dark - bright;
                const Real len2 = std::max(glm::dot(axis, axis), Real(1e-12));
                const Real t = std::clamp(Real(1) - glm::dot(rec.p - bright, axis) / len2, Real(0), Real(1));
                f *= gradientExponent == Real(1) ? t : std::pow(t, gradientExponent);
            }
            if (noise) f *= std::clamp(Real(1) - noiseAmount * (*noise)(rec.p), Real(0), Real(1));
            return f;
        }
    };

    /// R(lambda, hours(p)) of a measured series over a surface.
    class FadingReflectance {
    public:
        FadingReflectance(io::FadingSeries series, ExposureField exposure, Real exposedHours, Real unexposedHours)
            : m_series(std::move(series)), m_exposure(std::move(exposure)), m_exposed(exposedHours), m_unexposed(unexposedHours) {
            for (const auto& c : m_series.curves) m_rgb.push_back(color::reflectanceCurveToRgb(c));
        }

        Real hoursAt(const SurfaceInteraction& rec) const { return m_unexposed + m_exposure.at(rec) * (m_exposed - m_unexposed); }

        /// Reflectance at one wavelength (nm) after `hours`.
        Real reflectance(Real hours, Real lambda_nm) const {
            size_t lo, hi; Real t;
            m_series.bracket(hours, lo, hi, t);
            const Real a = std::max(m_series.curves[lo].evaluate(lambda_nm, color::Extrapolate::Clamp), Real(0));
            if (lo == hi) return a;
            const Real b = std::max(m_series.curves[hi].evaluate(lambda_nm, color::Extrapolate::Clamp), Real(0));
            return a + (b - a) * t;
        }

        /// The color under D65 after `hours` (what the RGB build uses; linear in the reflectance, so it is the blend of the rounds' colors).
        Vector3 rgb(Real hours) const {
            size_t lo, hi; Real t;
            m_series.bracket(hours, lo, hi, t);
            return lo == hi ? m_rgb[lo] : m_rgb[lo] * (Real(1) - t) + m_rgb[hi] * t;
        }

        const io::FadingSeries& series() const { return m_series; }
        const ExposureField& exposure() const { return m_exposure; }

    private:
        io::FadingSeries m_series;
        ExposureField m_exposure;
        Real m_exposed, m_unexposed;
        std::vector<Vector3> m_rgb;   ///< D65 color of every round
    };

} // namespace rayt
