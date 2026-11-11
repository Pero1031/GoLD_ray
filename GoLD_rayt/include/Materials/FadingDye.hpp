/**
 * @file Materials/FadingDye.hpp
 * @brief Fading of a dyed paper from its color alone: a one-parameter Kubelka-Munk model, with the exposure varying over the surface.
 *
 * K/S = (1 - R)^2 / (2 R) of a dyed paper is the K/S of the bare paper plus the dye's share D = max(0, K/S_color - K/S_paper).
 * Fading lowers the dye's concentration:   K/S(t) = K/S_paper + c(t) D.    A color brighter than the bare paper has no dye share and
 * does not change (so white stays white). c(t) is a table (hours, c) with c(0) = 1, from tools/data/fading_dye_model.py.
 *
 * The spectral build turns the color (RGB or a texel of an image) into the engine's smooth reflectance spectrum and fades it per
 * wavelength. The RGB build fades each channel as if it were a reflectance (an approximation).
 *
 * Only dyed paper was tested against measurements (the dye series of the lab's fading experiment: the color predicted from the RGB of the new
 * paper alone, left out of the fit of c(t), is off by about 5 CIEDE2000); pigments and dyed origami do not follow this model.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "Core/Types.hpp"
#include "Core/Interaction.hpp"
#include "Core/Texture.hpp"
#include "Materials/FadingReflectance.hpp"

namespace rayt {

    class DyeFading {
    public:
        struct Point { Real hours, c; };

        explicit DyeFading(std::vector<Point> table, Real paperReflectance = Real(0.9))
            : m_table(std::move(table)), m_ksPaper(ks(paperReflectance)), m_paper(paperReflectance) {}

        /// Reads "hours,c" rows (a header and comments are skipped); the hours must increase.
        static bool load(const std::string& file, Real paperReflectance, std::shared_ptr<const DyeFading>& out, std::string* err = nullptr) {
            std::ifstream in(file);
            if (!in) { if (err) *err = "cannot open " + file; return false; }
            std::vector<Point> table;
            std::string line;
            while (std::getline(in, line)) {
                if (line.empty() || line[0] == '#') continue;
                const size_t comma = line.find(',');
                if (comma == std::string::npos) continue;
                try {
                    size_t i0 = 0, i1 = 0;
                    const double h = std::stod(line.substr(0, comma), &i0);
                    const double c = std::stod(line.substr(comma + 1), &i1);
                    if (i0 == 0 || i1 == 0) continue;
                    table.push_back({ Real(h), Real(c) });
                }
                catch (...) { continue; }   // the header
            }
            if (table.empty()) { if (err) *err = "no \"hours,c\" rows in " + file; return false; }
            for (size_t i = 1; i < table.size(); ++i)
                if (!(table[i].hours > table[i - 1].hours)) { if (err) *err = "the hours must increase from row to row in " + file; return false; }
            out = std::make_shared<const DyeFading>(std::move(table), paperReflectance);
            return true;
        }

        Real firstHours() const { return m_table.front().hours; }
        Real lastHours() const { return m_table.back().hours; }
        Real paper() const { return m_paper; }

        /// The concentration after `hours` of exposure: linear between the rows, held at the ends (nothing is extrapolated).
        Real concentration(Real hours) const {
            if (hours <= m_table.front().hours) return m_table.front().c;
            if (hours >= m_table.back().hours) return m_table.back().c;
            size_t hi = 1;
            while (m_table[hi].hours < hours) ++hi;
            const Point& a = m_table[hi - 1];
            const Point& b = m_table[hi];
            return a.c + (b.c - a.c) * (hours - a.hours) / (b.hours - a.hours);
        }

        /// The reflectance (one wavelength or one channel) of a color r0 after `hours`.
        Real fade(Real r0, Real hours) const {
            const Real c = concentration(hours);
            if (c >= Real(1)) return r0;
            const Real d = ks(std::clamp(r0, Real(0.02), Real(0.98))) - m_ksPaper;
            if (d <= Real(0)) return r0;   // brighter than the bare paper: no dye
            return reflectance(m_ksPaper + c * d);
        }

        static Real ks(Real r) { r = std::clamp(r, Real(0.02), Real(0.98)); return (Real(1) - r) * (Real(1) - r) / (Real(2) * r); }
        static Real reflectance(Real k) { return Real(1) + k - std::sqrt(k * k + Real(2) * k); }

    private:
        std::vector<Point> m_table;
        Real m_ksPaper, m_paper;
    };

    /// A color (constant or an image) that fades with an exposure that varies over the surface.
    class DyeFadingSurface {
    public:
        DyeFadingSurface(std::shared_ptr<const DyeFading> model, Vector3 rgb, std::shared_ptr<const Texture2D> image, ExposureField exposure,
                         Real exposedHours, Real unexposedHours, Real fastness)
            : m_model(std::move(model)), m_rgb(rgb), m_image(std::move(image)), m_exposure(std::move(exposure)),
              m_exposed(exposedHours), m_unexposed(unexposedHours), m_fastness(fastness) {}

        Vector3 baseRgb(const SurfaceInteraction& rec) const { return m_image ? m_image->rgb(rec.uv.x, rec.uv.y) : m_rgb; }

        /// The time the model sees at this point: the hours of exposure times the fastness.
        Real hoursAt(const SurfaceInteraction& rec) const { return (m_unexposed + m_exposure.at(rec) * (m_exposed - m_unexposed)) * m_fastness; }

        Real fade(Real r0, Real modelHours) const { return m_model->fade(r0, modelHours); }

        /// The faded color as RGB (the RGB build and the denoiser's albedo): every channel faded as a reflectance.
        Vector3 rgbAt(const SurfaceInteraction& rec) const {
            const Vector3 c = baseRgb(rec);
            const Real h = hoursAt(rec);
            return Vector3(m_model->fade(c.x, h), m_model->fade(c.y, h), m_model->fade(c.z, h));
        }

        const DyeFading& model() const { return *m_model; }

    private:
        std::shared_ptr<const DyeFading> m_model;
        Vector3 m_rgb;
        std::shared_ptr<const Texture2D> m_image;
        ExposureField m_exposure;
        Real m_exposed, m_unexposed, m_fastness;
    };

} // namespace rayt
