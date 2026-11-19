// A diffuse surface that fluoresces (Materials/Fluorescent.hpp, docs/FLUORESCENCE.md).
#include "pch.h"

#include "Materials/Fluorescent.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

#include "Color/Illuminants.hpp"
#include "Color/SpectralConversion.hpp"
#include "Color/SpectralCurve.hpp"
#include "Core/Sampling.hpp"
#include "Geometry/Frame.hpp"

namespace rayt {

    int BispectralMatrix::diagonalColumn(int r) const {
        const Real x = (em(r) - exFirst) / exStep;
        const int c = static_cast<int>(std::lround(x));
        return (c >= 0 && c < nEx && std::abs(x - Real(c)) < Real(1e-3)) ? c : -1;
    }

    namespace {
        bool fail(std::string* err, const std::string& msg) { if (err) *err = msg; return false; }

        /// True when every step of v is `step` (within 1e-3 nm).
        bool evenlySpaced(const std::vector<Real>& v, Real& step) {
            if (v.size() < 2) return false;
            step = v[1] - v[0];
            for (size_t i = 1; i < v.size(); ++i) if (std::abs(v[i] - v[i - 1] - step) > Real(1e-3)) return false;
            return step > Real(0);
        }
    }

    bool LoadBfcMatrix(const std::string& path, BispectralMatrix& out, BispectralCleaning* cleaning, std::string* err) {
        std::ifstream in(path);
        if (!in) return fail(err, "cannot open " + path);
        out = BispectralMatrix();
        std::string line;
        int lineNo = 0;
        std::vector<Real> ex, em;
        std::vector<std::vector<Real>> rows;
        bool header = false;
        while (std::getline(in, line)) {
            ++lineNo;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (lineNo == 3 && !line.empty() && line[0] == ';') { out.title = line.substr(1); continue; }   // ";BCRbk1wr. BFC": file name
            if (lineNo == 4 && !line.empty() && line[0] == ';') { out.title = line.substr(1); continue; }   // ";BCRA black Tile ...": the sample
            if (line.rfind("r:c:", 0) == 0) {
                std::istringstream s(line.substr(4));
                Real v;
                while (s >> v) ex.push_back(v);
                header = true;
                continue;
            }
            if (!header) continue;
            std::istringstream s(line);
            Real w;
            if (!(s >> w)) continue;
            std::vector<Real> row;
            Real v;
            while (s >> v) row.push_back(v);
            if (row.empty()) continue;
            if (row.size() != ex.size()) return fail(err, path + ": line " + std::to_string(lineNo) + " has " + std::to_string(row.size()) +
                                                         " values for " + std::to_string(ex.size()) + " excitation wavelengths");
            em.push_back(w);
            rows.push_back(std::move(row));
        }
        if (!header) return fail(err, path + ": no \"r:c:\" row of excitation wavelengths (not a BFC-450 matrix file)");
        if (!evenlySpaced(ex, out.exStep) || !evenlySpaced(em, out.emStep)) return fail(err, path + ": the wavelengths are not evenly spaced");
        out.exFirst = ex.front();
        out.emFirst = em.front();
        out.nEx = static_cast<int>(ex.size());
        out.nEm = static_cast<int>(em.size());
        out.d.assign(static_cast<size_t>(out.nEm) * static_cast<size_t>(out.nEx), Real(0));
        BispectralCleaning c;
        auto antiStokes = [&](int r, int k) { return out.ex(k) >= out.em(r) + Real(1.5) * out.exStep; };   // excitation 20 nm or more above
        double sq = 0;
        int nq = 0;
        for (int r = 0; r < out.nEm; ++r)
            for (int k = 0; k < out.nEx; ++k)
                if (antiStokes(r, k)) { const double v = rows[static_cast<size_t>(r)][static_cast<size_t>(k)]; sq += v * v; ++nq; }
        c.noise = nq > 0 ? std::sqrt(sq / nq) : 0.0;
        const Real floor = Real(3 * c.noise);
        for (int r = 0; r < out.nEm; ++r) {
            const int dc = out.diagonalColumn(r);
            for (int k = 0; k < out.nEx; ++k) {
                Real v = rows[static_cast<size_t>(r)][static_cast<size_t>(k)];
                if (antiStokes(r, k)) {
                    if (v > 0) { ++c.antiStokes; c.antiStokesSum += v; }
                    v = 0;
                }
                if (v < 0) { ++c.negative; c.negativeSum += v; v = 0; }
                if (k != dc && v > 0 && v < floor) { ++c.belowNoise; c.belowNoiseSum += v; v = 0; }
                out.d[static_cast<size_t>(r) * static_cast<size_t>(out.nEx) + static_cast<size_t>(k)] = v;
            }
        }
        if (cleaning) *cleaning = c;
        return true;
    }

    // ------------------------------------------------------------------------------------------------ the material

    Fluorescent::Fluorescent(std::shared_ptr<const BispectralMatrix> m, bool fluorescence) : m_m(std::move(m)), m_fluo(fluorescence) {
        const BispectralMatrix& M = *m_m;
        m_R.assign(static_cast<size_t>(M.nEm), Real(0));
        m_F.assign(static_cast<size_t>(M.nEm), Real(0));
        for (int r = 0; r < M.nEm; ++r) {
            const int dc = M.diagonalColumn(r);
            for (int c = 0; c < M.nEx; ++c) {
                if (c == dc) m_R[static_cast<size_t>(r)] = M.at(r, c);
                else m_F[static_cast<size_t>(r)] += M.at(r, c);
            }
        }

        // The color under D65 with its ultraviolet: the total radiance factor beta(lo) = R(lo) + Sum_c D(lo, c) Ebar(c) / E(lo), as a
        // reflectance curve on 380..780 nm, reduced to RGB as any measured reflectance is.
        color::SpectralCurve beta;
        for (int l = 380; l <= 780; ++l) {
            const Real lo = Real(l);
            Real b = reflectance(lo);
            const Real e = color::evaluateD65(lo);
            int r;
            Real t;
            if (m_fluo && e > Real(0) && rowOf(lo, r, t)) {
                for (int c = 0; c < M.nEx; ++c) {
                    const Real w = columnWeight(r, t, c);
                    if (w <= Real(0)) continue;
                    Real ebar = 0;   // D65 averaged over the band of column c (1 nm steps)
                    for (int k = 0; k < 10; ++k) ebar += color::evaluateD65Full(M.ex(c) - Real(0.5) * M.exStep + (Real(k) + Real(0.5)) * M.exStep / Real(10));
                    b += w * (ebar / Real(10)) / e;
                }
            }
            beta.add(lo, b);
        }
        beta.finalize();
        m_rgbD65 = color::reflectanceCurveToRgb(beta);
    }

    bool Fluorescent::rowOf(Real lambdaO, int& r, Real& t) const {
        const BispectralMatrix& M = *m_m;
        const Real x = (lambdaO - M.emFirst) / M.emStep;
        if (x < Real(0) || x > Real(M.nEm - 1)) return false;
        r = std::min(static_cast<int>(x), M.nEm - 2);
        t = x - Real(r);
        return true;
    }

    Real Fluorescent::columnWeight(int r, Real t, int c) const {
        const BispectralMatrix& M = *m_m;
        const Real a = (c == M.diagonalColumn(r)) ? Real(0) : M.at(r, c);
        const Real b = (c == M.diagonalColumn(r + 1)) ? Real(0) : M.at(r + 1, c);
        return (Real(1) - t) * a + t * b;
    }

    Real Fluorescent::reflectance(Real lambdaO) const {
        const BispectralMatrix& M = *m_m;
        const Real x = std::clamp((lambdaO - M.emFirst) / M.emStep, Real(0), Real(M.nEm - 1));
        const int r = std::min(static_cast<int>(x), M.nEm - 2);
        const Real t = x - Real(r);
        return std::max(Real(0), (Real(1) - t) * m_R[static_cast<size_t>(r)] + t * m_R[static_cast<size_t>(r) + 1]);
    }

    Real Fluorescent::fluorescence(Real lambdaO) const {
        int r;
        Real t;
        if (!m_fluo || !rowOf(lambdaO, r, t)) return Real(0);
        return (Real(1) - t) * m_F[static_cast<size_t>(r)] + t * m_F[static_cast<size_t>(r) + 1];
    }

    Real Fluorescent::sampleExcitation(Real lambdaO, Real u) const {
        const BispectralMatrix& M = *m_m;
        int r;
        Real t;
        if (!rowOf(lambdaO, r, t)) return lambdaO;
        const Real total = fluorescence(lambdaO);
        if (!(total > Real(0))) return lambdaO;
        Real target = u * total;
        int last = -1;
        for (int c = 0; c < M.nEx; ++c) {
            const Real w = columnWeight(r, t, c);
            if (w <= Real(0)) continue;
            last = c;
            if (target < w) {   // inside this band: the remainder of u places the wavelength uniformly in it
                const Real v = std::clamp(target / w, Real(0), Real(1) - Real(1e-6));
                return M.ex(c) + (v - Real(0.5)) * M.exStep;
            }
            target -= w;
        }
        return last >= 0 ? M.ex(last) : lambdaO;   // rounding: the last band
    }

    Real Fluorescent::excitationDensity(Real lambdaO, Real lambdaI) const {
        const BispectralMatrix& M = *m_m;
        int r;
        Real t;
        if (!m_fluo || !rowOf(lambdaO, r, t)) return Real(0);
        const int c = static_cast<int>(std::floor((lambdaI - M.exFirst) / M.exStep + Real(0.5)));
        if (c < 0 || c >= M.nEx) return Real(0);
        return columnWeight(r, t, c) / M.exStep;
    }

    bool Fluorescent::shiftsWavelength() const {
#if defined(RAYT_SPECTRAL_MODE)
        return m_fluo;
#else
        return false;
#endif
    }

    // ------------------------------------------------------------------------------------------------ sampling with the light

    // Share of the excitation samples drawn from the matrix row alone when a light spectrum is known; the rest follow row x light.
    // The weight row / mixture is then at most 1 / kRowShare = 4 wherever the given spectrum is wrong (a light it leaves out). With
    // the RIT samples under a black light, one sample's relative deviation at the vertex is 3 to 6 with the row alone, 1.3 to 1.5 with
    // a share of 1/2 and 0.9 to 1.2 with 1/4.
    constexpr Real kRowShare = Real(0.25);

    void Fluorescent::setExcitationLight(const std::function<Real(Real)>& spectrum) {
        const BispectralMatrix& M = *m_m;
        m_cellsPerBand = std::max(1, static_cast<int>(std::lround(M.exStep)));   // cells of about 1 nm
        m_cellWidth = M.exStep / Real(m_cellsPerBand);
        m_cellOrigin = M.exFirst - Real(0.5) * M.exStep;
        m_cell.assign(static_cast<size_t>(M.nEx) * static_cast<size_t>(m_cellsPerBand), Real(0));
        m_bandLight.assign(static_cast<size_t>(M.nEx), Real(0));
        Real total = 0;
        for (size_t j = 0; j < m_cell.size(); ++j) {
            const Real v = std::max(Real(0), spectrum(m_cellOrigin + (Real(j) + Real(0.5)) * m_cellWidth));
            m_cell[j] = std::isfinite(v) ? v : Real(0);
            m_bandLight[j / static_cast<size_t>(m_cellsPerBand)] += m_cell[j];
            total += m_cell[j];
        }
        m_hasLight = total > Real(0);
    }

    Real Fluorescent::lightNorm(int r, Real t) const {
        Real s = 0;
        for (int c = 0; c < m_m->nEx; ++c) s += columnWeight(r, t, c) * m_bandLight[static_cast<size_t>(c)];
        return s;
    }

    Real Fluorescent::sampleExcitationWithLight(Real lambdaO, Real u) const {
        const BispectralMatrix& M = *m_m;
        int r;
        Real t;
        if (!m_hasLight || !rowOf(lambdaO, r, t)) return lambdaO;
        const Real norm = lightNorm(r, t);
        if (!(norm > Real(0))) return lambdaO;
        Real target = u * norm;
        for (int c = 0; c < M.nEx; ++c) {
            const Real w = columnWeight(r, t, c) * m_bandLight[static_cast<size_t>(c)];
            if (w <= Real(0)) continue;
            if (target < w || c == M.nEx - 1) {   // this band; then a cell in proportion to the light, uniform inside it
                Real left = std::clamp(target / w, Real(0), Real(1) - Real(1e-6)) * m_bandLight[static_cast<size_t>(c)];
                const size_t j0 = static_cast<size_t>(c) * static_cast<size_t>(m_cellsPerBand);
                for (int k = 0; k < m_cellsPerBand; ++k) {
                    const Real s = m_cell[j0 + static_cast<size_t>(k)];
                    if (s <= Real(0)) continue;
                    if (left < s || k == m_cellsPerBand - 1) {
                        const Real v = std::clamp(left / s, Real(0), Real(1) - Real(1e-6));
                        return m_cellOrigin + (Real(j0 + static_cast<size_t>(k)) + v) * m_cellWidth;
                    }
                    left -= s;
                }
            }
            target -= w;
        }
        return lambdaO;
    }

    Real Fluorescent::lightDensity(Real lambdaO, Real lambdaI) const {
        const BispectralMatrix& M = *m_m;
        int r;
        Real t;
        if (!m_hasLight || !m_fluo || !rowOf(lambdaO, r, t)) return Real(0);
        const int j = static_cast<int>(std::floor((lambdaI - m_cellOrigin) / m_cellWidth));
        if (j < 0 || j >= static_cast<int>(m_cell.size())) return Real(0);
        const Real norm = lightNorm(r, t);
        if (!(norm > Real(0))) return Real(0);
        const int c = j / m_cellsPerBand;
        (void)M;
        return columnWeight(r, t, c) * m_cell[static_cast<size_t>(j)] / (m_cellWidth * norm);
    }

    Real Fluorescent::mixtureDensity(Real lambdaO, Real lambdaI) const {
        const Real F = fluorescence(lambdaO);
        if (!(F > Real(0))) return Real(0);
        const Real pRow = excitationDensity(lambdaO, lambdaI) / F;
        int r;
        Real t;
        const bool light = m_hasLight && rowOf(lambdaO, r, t) && lightNorm(r, t) > Real(0);
        return light ? kRowShare * pRow + (Real(1) - kRowShare) * lightDensity(lambdaO, lambdaI) : pRow;
    }

    void Fluorescent::reemit(SurfaceInteraction& rec, SampledWavelengths& lambda) const {
#if defined(RAYT_SPECTRAL_MODE)
        static_assert(kSpectralSamples <= 4, "SurfaceInteraction::reemitAlbedo holds four wavelengths");
        // Each wavelength the path still carries (one after a dispersive surface) chooses on its own: the direction a diffuse
        // surface sends light in does not depend on the wavelength, so all of them go on along the same path.
        const int n = lambda.secondaryTerminated() ? 1 : kSpectralSamples;
#else
        const int n = 1;
#endif
        for (int i = 0; i < n; ++i) {
            const Real lo = lambda.lambda[static_cast<size_t>(i)];
            const Real R = reflectance(lo), F = fluorescence(lo);
            rec.reemitAlbedo[static_cast<size_t>(i)] = R + F;
            if (!(F > Real(0)) || sampling::Random() * (R + F) < R) continue;   // reflection: the wavelength stays
            int r;
            Real t;
            const bool light = m_hasLight && rowOf(lo, r, t) && lightNorm(r, t) > Real(0);
            const Real li = (light && sampling::Random() >= kRowShare) ? sampleExcitationWithLight(lo, sampling::Random())
                                                                         : sampleExcitation(lo, sampling::Random());
            const Real pRow = excitationDensity(lo, li) / F;
            const Real q = light ? kRowShare * pRow + (Real(1) - kRowShare) * lightDensity(lo, li) : pRow;
            rec.reemitAlbedo[static_cast<size_t>(i)] = q > Real(0) ? (R + F) * pRow / q : Real(0);
            lambda.lambda[static_cast<size_t>(i)] = li;
        }
        for (int i = n; i < 4; ++i) rec.reemitAlbedo[static_cast<size_t>(i)] = rec.reemitAlbedo[0];
    }

    Spectrum Fluorescent::eval(const SurfaceInteraction& rec, const Vector3& /*wo*/, const Vector3& wi, const SampledWavelengths& lambda,
                               TransportMode /*mode*/) const {
        if (glm::dot(rec.gn, wi) <= Real(0)) return Spectrum(0.0);
#if defined(RAYT_SPECTRAL_MODE)
        if (rec.reemitAlbedo[0] >= Real(0)) {
            Spectrum s;
            for (int i = 0; i < kSpectralSamples; ++i) s[i] = rec.reemitAlbedo[static_cast<size_t>(i)] / constants::PI;
            return s;
        }
#endif
#if defined(RAYT_SPECTRAL_MODE)
        // Not prepared by reemit() (a renderer other than the path tracer): the reflectance alone, no fluorescence.
        Spectrum s;
        for (int i = 0; i < kSpectralSamples; ++i) s[i] = reflectance(lambda[i]) / constants::PI;
        return s;
#else
        return RgbColor(m_rgbD65).reflectance(lambda) * (Real(1) / constants::PI);
#endif
    }

    std::optional<BSDFSample> Fluorescent::sample(const SurfaceInteraction& rec, const Vector3& wo, const Point2& u,
                                                  const SampledWavelengths& lambda, TransportMode mode) const {
        BSDFSample s;
        const frame::Frame onb(rec.n);
        const Vector3 local = sampling::CosineSampleHemisphere(u);
        s.wi = onb.localToWorld(local);
        if (glm::dot(rec.gn, s.wi) <= Real(0)) return std::nullopt;
        s.pdf = local.z / constants::PI;
        s.f = eval(rec, wo, s.wi, lambda, mode);
        s.flags = BxDFFlags::Diffuse | BxDFFlags::Reflection;
        return s;
    }

    Real Fluorescent::pdf(const SurfaceInteraction& rec, const Vector3& /*wo*/, const Vector3& wi, const SampledWavelengths& /*lambda*/) const {
        if (glm::dot(rec.gn, wi) <= Real(0)) return Real(0);
        const Real c = glm::dot(rec.n, wi);
        return c > Real(0) ? c / constants::PI : Real(0);
    }

    Vector3 Fluorescent::aovAlbedo(const SurfaceInteraction& /*rec*/) const {
        return glm::clamp(m_rgbD65, Vector3(0.0), Vector3(1.0));
    }

} // namespace rayt
