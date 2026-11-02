/**
 * @file Color/ComplexIorTable.hpp
 * @brief Tabulated complex index of refraction (eta, k) of a conductor as a function of wavelength.
 *
 * This is the canonical spectral description of a metal: every other input (RefractiveIndex.info
 * CSV, 3-column tables, separate n / k files, normal-incidence reflectance, ...) is converted into
 * one of these tables, and RoughConductor / the RGB reduction only ever see this class.
 *
 * Internally eta(l) and k(l) are two SpectralCurve objects (Color/SpectralCurve.hpp).
 *
 * Notes:
 * - Wavelengths are in nanometers (nm).
 * - Out-of-range wavelengths are clamped to the nearest endpoint (no extrapolation), which is
 *   safer for optical constants than linear extrapolation.
 * - Fill with add() (any order) and call finalize() once; file parsing lives in the IO loaders.
 */

#pragma once

#include <initializer_list>
#include <utility>
#include <vector>

#include "Core/Precision.hpp"
#include "Color/SpectralCurve.hpp"

namespace rayt::color {

    class ComplexIorTable {
    public:
        /// One measured / fitted sample: eta = real part, k = extinction coefficient.
        struct Sample {
            Real lambda_nm = Real(0);
            Real eta = Real(1);
            Real k = Real(0);
        };

        /// Appends a sample. Call finalize() after the last one.
        void add(Real lambda_nm, Real eta, Real k) {
            m_eta.add(lambda_nm, eta);
            m_k.add(lambda_nm, k);
        }

        /// Replaces the content with the given samples (already finalized).
        void assign(std::initializer_list<Sample> samples) {
            m_eta = SpectralCurve();
            m_k = SpectralCurve();
            for (const Sample& s : samples) add(s.lambda_nm, s.eta, s.k);
            finalize();
        }

        /// Sorts by wavelength and merges duplicate wavelengths.
        void finalize() {
            m_eta.finalize();
            m_k.finalize();
        }

        /// Builds a table on the union of both curves' wavelengths (each curve clamped at its ends).
        static ComplexIorTable fromCurves(const SpectralCurve& eta, const SpectralCurve& k) {
            ComplexIorTable t;
            for (const auto& s : eta.samples()) t.add(s.lambda_nm, eta.evaluate(s.lambda_nm, Extrapolate::Clamp), k.evaluate(s.lambda_nm, Extrapolate::Clamp));
            for (const auto& s : k.samples()) t.add(s.lambda_nm, eta.evaluate(s.lambda_nm, Extrapolate::Clamp), k.evaluate(s.lambda_nm, Extrapolate::Clamp));
            t.finalize();   // identical wavelengths from both loops are merged
            return t;
        }

        bool empty() const { return m_eta.empty(); }
        size_t size() const { return m_eta.size(); }

        const SpectralCurve& etaCurve() const { return m_eta; }
        const SpectralCurve& kCurve() const { return m_k; }

        /// Linear interpolation, clamped at both ends. Returns (eta, k); an empty table is (1, 0).
        std::pair<Real, Real> evaluateNK(Real lambda_nm) const {
            if (m_eta.empty()) return { Real(1), Real(0) };
            return { m_eta.evaluate(lambda_nm, Extrapolate::Clamp), m_k.evaluate(lambda_nm, Extrapolate::Clamp) };
        }

    private:
        SpectralCurve m_eta;
        SpectralCurve m_k;
    };

} // namespace rayt::color
