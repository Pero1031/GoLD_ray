/**
 * @file Color/EmissionSpectrum.hpp
 * @brief Spectral power distribution of a light source (point / area light emission).
 *
 * Three ways to describe a light:
 *  - EmissionSpectrum(RgbColor)          : the old behavior, RGB x D65 (white = D65).
 *  - EmissionSpectrum::Blackbody(K)      : Planck radiator (candle 1850 K, tungsten A 2856 K, ...).
 *  - EmissionSpectrum::FromCurve(curve)  : any measured / generated SPD (see IO/SpectralCurveLoader.hpp).
 *  - EmissionSpectrum::FromCurves(...)   : a weighted sum of SPDs (several LEDs lit together).
 *
 * Brightness convention: curve-based lights are scaled so that their luminance Y equals the
 * `luminance` argument. Y = 1 is exactly as bright as RgbColor(1) (= D65, RGB (1,1,1)), so
 * switching a light between D65 and another SPD changes its color but not its brightness.
 *
 * RGB mode: the SPD is reduced once to the linear sRGB color the film sensor would record for it
 * (the same sensor that SpectrumToRGB() uses), so RGB and spectral mode agree.
 *
 * A non-D65 light has a color cast in the image (as with a real camera at fixed white balance).
 * setWhiteBalance(light) removes it at the film (CAT16 adaptation to D65, see Color/ChromaticAdaptation.hpp).
 */

#pragma once

#include <algorithm>
#include <vector>

#include "Core/Types.hpp"
#include "Core/SampledWavelengths.hpp"
#include "Color/SpectralCurve.hpp"
#include "Color/Illuminants.hpp"
#include "Color/SpectralConversion.hpp"

namespace rayt {

    class EmissionSpectrum {
    public:
        EmissionSpectrum() = default;

        /// RGB x D65 emission (implicit on purpose: existing RgbColor call sites keep working).
        EmissionSpectrum(const RgbColor& rgb) : m_rgb(rgb.rgb()) {}

        /**
         * @brief Emission with an arbitrary SPD (relative values; zero outside the curve's range).
         * @param luminance Y of the resulting light (1 = as bright as RgbColor(1)).
         */
        static EmissionSpectrum FromCurve(const color::SpectralCurve& curve, Real luminance = Real(1)) {
            return FromCurves({ &curve }, { Real(1) }, luminance);
        }

        /**
         * @brief Emission of several sources lit together: the sum of the curves, each times its weight (a mix of LEDs).
         *        Only the proportions of the weights matter (the result is scaled to `luminance`); the weights are >= 0.
         */
        static EmissionSpectrum FromCurves(const std::vector<const color::SpectralCurve*>& curves, const std::vector<Real>& weights,
                                           Real luminance = Real(1)) {
            EmissionSpectrum e;
            std::vector<Real> lut(static_cast<size_t>(kLutLastNm - kLutFirstNm + 1), Real(0));
            for (size_t k = 0; k < curves.size() && k < weights.size(); ++k)
                for (int l = kLutFirstNm; l <= kLutLastNm; ++l)
                    lut[static_cast<size_t>(l - kLutFirstNm)] += weights[k] * std::max(curves[k]->evaluate(Real(l), color::Extrapolate::Zero), Real(0));

            // Sensor RGB of the unscaled SPD (1 nm trapezoid over 380..780, as the sensor's white point).
            Vector3 rgb(0.0);
            for (int l = 380; l <= 780; ++l) {
                const Real w = (l == 380 || l == 780) ? Real(0.5) : Real(1);
                rgb += color::rgbSensorWeights(Real(l)) * (lut[static_cast<size_t>(l - kLutFirstNm)] * w);
            }
            const Real Y = Real(0.2126) * rgb.x + Real(0.7152) * rgb.y + Real(0.0722) * rgb.z;
            if (Y <= Real(0)) return e;   // no visible power: black light

            const Real scale = luminance / Y;
            for (Real& v : lut) v *= scale;
            e.m_rgb = rgb * scale;
            e.m_lut = std::move(lut);
            return e;
        }

        /// Planck radiator of the given temperature (K).
        static EmissionSpectrum Blackbody(Real kelvin, Real luminance = Real(1)) {
            color::SpectralCurve c;
            for (int l = kLutFirstNm; l <= kLutLastNm; l += 5) c.add(Real(l), color::evaluateBlackbody(Real(l), kelvin));
            c.finalize();
            return FromCurve(c, luminance);
        }

        /// CIE standard illuminant A (tungsten, blackbody 2856 K).
        static EmissionSpectrum IlluminantA(Real luminance = Real(1)) { return Blackbody(Real(2856), luminance); }

        /// Emission at the wavelengths of the current path (RGB mode: the RGB color).
        Spectrum operator()(const SampledWavelengths& lambda) const {
#if defined(RAYT_SPECTRAL_MODE)
            if (m_lut.empty()) return RgbColor(m_rgb).illuminant(lambda);
            Spectrum s;
            for (int i = 0; i < kSpectralSamples; ++i) s[i] = lutAt(lambda[i]);
            return s;
#else
            (void)lambda;
            return m_rgb;
#endif
        }

        /// Linear sRGB color of this light as recorded by the film sensor.
        const Vector3& rgb() const { return m_rgb; }

        EmissionSpectrum operator*(Real s) const {
            EmissionSpectrum e = *this;
            e.m_rgb *= s;
            for (Real& v : e.m_lut) v *= s;
            return e;
        }

    private:
        static constexpr int kLutFirstNm = 360;
        static constexpr int kLutLastNm = 830;

        Real lutAt(Real lambda_nm) const {
            const Real x = std::clamp(lambda_nm, Real(kLutFirstNm), Real(kLutLastNm)) - Real(kLutFirstNm);
            const int i = std::min(static_cast<int>(x), static_cast<int>(m_lut.size()) - 2);
            const Real f = x - Real(i);
            return m_lut[static_cast<size_t>(i)] * (Real(1) - f) + m_lut[static_cast<size_t>(i) + 1] * f;
        }

        Vector3 m_rgb = Vector3(1.0);
        std::vector<Real> m_lut;   ///< 1 nm table (empty = RGB x D65)
    };

    namespace color {

        /**
         * @brief Makes `light` render neutral: the film adapts from the light's white to D65 with CAT16 (luminance preserved).
         *
         * The light's white is its XYZ (from the sensor RGB the film would record); the destination is D65 at the same Y.
         * Unlike a division of the RGB channels this also works for a light whose sRGB blue is negative (the candle).
         * Returns false (and leaves the film unchanged) when the light has no visible power.
         */
        inline bool setWhiteBalance(const EmissionSpectrum& light, Real degree = Real(1)) {
            const Matrix3 toXYZ = filmToXYZ();
            const Vector3 src = toXYZ * light.rgb();
            const Vector3 d65 = toXYZ * Vector3(1.0);   // the film's (1,1,1) is D65
            if (!(src.y > Real(0)) || !cat16WhiteValid(src)) return false;
            const Vector3 dst = d65 * (src.y / d65.y);
            WhiteBalance& wb = whiteBalance();
            wb.matrix = glm::inverse(toXYZ) * cat16Matrix(src, dst, degree) * toXYZ;
            wb.active = true;
            return true;
        }

        inline void resetWhiteBalance() { whiteBalance() = WhiteBalance(); }

    } // namespace color

} // namespace rayt
