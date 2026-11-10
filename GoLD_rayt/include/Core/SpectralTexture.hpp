#pragma once
/**
 * @file Core/SpectralTexture.hpp
 * @brief A reflectance texture whose texels are spectra: R(u, v, lambda) (header-only).
 *
 * A measured spectral image (e.g. a 2D spectroradiometer shot of a cloth divided by a white plate) is
 * stored compactly as a mean spectrum plus K basis spectra and K coefficients per texel:
 *     R(u, v, lambda) = mean(lambda) + sum_k c_k(u, v) * basis_k(lambda)
 * The coefficients are bilinearly filtered with repeat wrapping (as Texture2D), the spectrum is
 * linearly interpolated between its wavelength samples and held constant outside their range, and the
 * result is clamped to [0, 1]. For the RGB build the D65 color of every texel is computed once at load.
 *
 * File format ".sptex" (little-endian), written by tools/data/hsm_texture.py:
 *   char[4] "SPTX", uint32 version = 1, uint32 width, uint32 height, uint32 K, uint32 nLambda,
 *   float32 lambda0_nm, float32 step_nm,
 *   float32 mean[nLambda], float32 basis[K][nLambda], float32 coef[height][width][K]  (row 0 = top).
 */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "Core/Types.hpp"
#include "Color/Illuminants.hpp"
#include "Color/SpectralConversion.hpp"

namespace rayt {

    class SpectralTexture {
    public:
        /// Loads a .sptex file. Returns nullptr and sets *err if the file cannot be read.
        static std::shared_ptr<SpectralTexture> load(const std::string& file, std::string* err = nullptr) {
            std::ifstream in(file, std::ios::binary);
            if (!in) { if (err) *err = "cannot open"; return nullptr; }
            char magic[4];
            uint32_t version = 0, w = 0, h = 0, k = 0, nl = 0;
            float l0 = 0, step = 0;
            in.read(magic, 4);
            in.read(reinterpret_cast<char*>(&version), 4);
            in.read(reinterpret_cast<char*>(&w), 4);
            in.read(reinterpret_cast<char*>(&h), 4);
            in.read(reinterpret_cast<char*>(&k), 4);
            in.read(reinterpret_cast<char*>(&nl), 4);
            in.read(reinterpret_cast<char*>(&l0), 4);
            in.read(reinterpret_cast<char*>(&step), 4);
            if (!in || std::memcmp(magic, "SPTX", 4) != 0 || version != 1) { if (err) *err = "not a .sptex file (version 1)"; return nullptr; }
            if (w == 0 || h == 0 || k == 0 || k > 64 || nl < 2 || nl > 4096 || step <= 0) { if (err) *err = "bad header"; return nullptr; }
            auto t = std::make_shared<SpectralTexture>();
            t->m_w = int(w); t->m_h = int(h); t->m_k = int(k); t->m_nl = int(nl); t->m_l0 = l0; t->m_step = step;
            t->m_mean.resize(nl);
            t->m_basis.resize(size_t(k) * nl);
            t->m_coef.resize(size_t(w) * h * k);
            in.read(reinterpret_cast<char*>(t->m_mean.data()), std::streamsize(t->m_mean.size() * 4));
            in.read(reinterpret_cast<char*>(t->m_basis.data()), std::streamsize(t->m_basis.size() * 4));
            in.read(reinterpret_cast<char*>(t->m_coef.data()), std::streamsize(t->m_coef.size() * 4));
            if (!in) { if (err) *err = "truncated file"; return nullptr; }
            t->computeRgb();
            return t;
        }

        int width() const { return m_w; }
        int height() const { return m_h; }
        int basisCount() const { return m_k; }
        Real firstNm() const { return Real(m_l0); }
        Real lastNm() const { return Real(m_l0 + m_step * (m_nl - 1)); }

        /// Reflectance at texture coordinates (u, v) and one wavelength (nm). Bilinear in (u, v), linear in lambda.
        Real reflectance(Real u, Real v, Real lambda_nm) const {
            Texel tx[4]; Real tw[4];
            footprint(u, v, tx, tw);
            Real x = (lambda_nm - Real(m_l0)) / Real(m_step);
            x = std::clamp(x, Real(0), Real(m_nl - 1));
            const int i0 = std::min(int(x), m_nl - 2);
            const Real f = x - Real(i0);
            Real r = 0;
            for (int t = 0; t < 4; ++t) {
                const float* c = &m_coef[(size_t(tx[t].y) * m_w + tx[t].x) * m_k];
                Real a = m_mean[i0], b = m_mean[i0 + 1];
                for (int kk = 0; kk < m_k; ++kk) {
                    a += c[kk] * m_basis[size_t(kk) * m_nl + i0];
                    b += c[kk] * m_basis[size_t(kk) * m_nl + i0 + 1];
                }
                r += tw[t] * (a + (b - a) * f);
            }
            return std::clamp(r, Real(0), Real(1));
        }

        /// Linear sRGB of the texel spectrum under D65 (bilinear). What the RGB build uses.
        Vector3 rgb(Real u, Real v) const {
            Texel tx[4]; Real tw[4];
            footprint(u, v, tx, tw);
            Vector3 c(0.0);
            for (int t = 0; t < 4; ++t) {
                const float* p = &m_rgb[(size_t(tx[t].y) * m_w + tx[t].x) * 3];
                c += tw[t] * Vector3(p[0], p[1], p[2]);
            }
            return c;
        }

    private:
        struct Texel { int x, y; };

        void footprint(Real u, Real v, Texel tx[4], Real tw[4]) const {
            const Real x = (u - std::floor(u)) * Real(m_w) - Real(0.5);
            const Real y = (v - std::floor(v)) * Real(m_h) - Real(0.5);
            const int x0 = int(std::floor(x)), y0 = int(std::floor(y));
            const Real fx = x - Real(x0), fy = y - Real(y0);
            int t = 0;
            for (int dy = 0; dy < 2; ++dy)
                for (int dx = 0; dx < 2; ++dx, ++t) {
                    tx[t] = { ((x0 + dx) % m_w + m_w) % m_w, ((y0 + dy) % m_h + m_h) % m_h };
                    tw[t] = (dx ? fx : Real(1) - fx) * (dy ? fy : Real(1) - fy);
                }
        }

        /// The spectrum of one texel at a wavelength (held outside the sampled range), no filtering.
        Real texelReflectance(int x, int y, Real lambda_nm) const {
            Real s = (lambda_nm - Real(m_l0)) / Real(m_step);
            s = std::clamp(s, Real(0), Real(m_nl - 1));
            const int i0 = std::min(int(s), m_nl - 2);
            const Real f = s - Real(i0);
            const float* c = &m_coef[(size_t(y) * m_w + x) * m_k];
            Real a = m_mean[i0], b = m_mean[i0 + 1];
            for (int kk = 0; kk < m_k; ++kk) {
                a += c[kk] * m_basis[size_t(kk) * m_nl + i0];
                b += c[kk] * m_basis[size_t(kk) * m_nl + i0 + 1];
            }
            return std::clamp(a + (b - a) * f, Real(0), Real(1));
        }

        /// D65 color of every texel: Integral( R(l) D65(l) s_c(l) dl ), like color::reflectanceCurveToRgb.
        void computeRgb() {
            m_rgb.assign(size_t(m_w) * m_h * 3, 0.0f);
            for (int y = 0; y < m_h; ++y)
                for (int x = 0; x < m_w; ++x) {
                    Vector3 rgb(0.0);
                    for (int l = 380; l <= 780; ++l) {
                        const Real trapezoid = (l == 380 || l == 780) ? Real(0.5) : Real(1);
                        rgb += color::rgbSensorWeights(Real(l)) * (texelReflectance(x, y, Real(l)) * color::evaluateD65(Real(l)) * trapezoid);
                    }
                    float* p = &m_rgb[(size_t(y) * m_w + x) * 3];
                    p[0] = float(rgb.x); p[1] = float(rgb.y); p[2] = float(rgb.z);
                }
        }

        int m_w = 0, m_h = 0, m_k = 0, m_nl = 0;
        float m_l0 = 380.0f, m_step = 1.0f;
        std::vector<float> m_mean, m_basis, m_coef, m_rgb;
    };

} // namespace rayt
