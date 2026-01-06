#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <glm/glm.hpp>

#include "Core/Types.hpp"   // Real, Float, etc.

namespace rayt {

    enum class SpectrumMode {
        RGB,
        Spectral
    };

    // ---- 設定（ここだけ触れば切り替わる） ----
    constexpr SpectrumMode kSpectrumMode = SpectrumMode::RGB;
    // constexpr SpectrumMode kSpectrumMode = SpectrumMode::Spectral;

    constexpr int kSpectralSamples = 30; // Spectral時のみ有効（例: 30 or 36）

    // -------------------------
    // 1) RGB 実装（軽い）
    // -------------------------
    class RGBSpectrum {
    public:
        glm::vec<3, Real, glm::defaultp> v;

        RGBSpectrum(Real x = 0) : v(x, x, x) {}
        RGBSpectrum(Real r, Real g, Real b) : v(r, g, b) {}

        bool isBlack() const { return v.x == 0 && v.y == 0 && v.z == 0; }

        RGBSpectrum& operator+=(const RGBSpectrum& o) { v += o.v; return *this; }
        RGBSpectrum& operator*=(const RGBSpectrum& o) { v *= o.v; return *this; }
        RGBSpectrum& operator*=(Real s) { v *= s; return *this; }

        friend RGBSpectrum operator+(RGBSpectrum a, const RGBSpectrum& b) { return a += b; }
        friend RGBSpectrum operator*(RGBSpectrum a, const RGBSpectrum& b) { return a *= b; }
        friend RGBSpectrum operator*(RGBSpectrum a, Real s) { return a *= s; }
        friend RGBSpectrum operator*(Real s, RGBSpectrum a) { return a *= s; }

        glm::vec3 toRGB() const { return glm::vec3((float)v.x, (float)v.y, (float)v.z); }
    };

    // -------------------------
    // 2) Spectral 実装（正確）
    // -------------------------
    template<int N>
    class SampledSpectrum {
    public:
        std::array<Real, N> c{};

        SampledSpectrum(Real x = 0) { c.fill(x); }

        bool isBlack() const {
            for (Real x : c) if (x != 0) return false;
            return true;
        }

        SampledSpectrum& operator+=(const SampledSpectrum& o) {
            for (int i = 0; i < N; ++i) c[i] += o.c[i];
            return *this;
        }
        SampledSpectrum& operator*=(const SampledSpectrum& o) {
            for (int i = 0; i < N; ++i) c[i] *= o.c[i];
            return *this;
        }
        SampledSpectrum& operator*=(Real s) {
            for (int i = 0; i < N; ++i) c[i] *= s;
            return *this;
        }

        friend SampledSpectrum operator+(SampledSpectrum a, const SampledSpectrum& b) { return a += b; }
        friend SampledSpectrum operator*(SampledSpectrum a, const SampledSpectrum& b) { return a *= b; }
        friend SampledSpectrum operator*(SampledSpectrum a, Real s) { return a *= s; }
        friend SampledSpectrum operator*(Real s, SampledSpectrum a) { return a *= s; }

        // TODO: 波長サンプル＋CIE等色関数でXYZ→RGB変換
        glm::vec3 toRGB() const {
            // 暫定: 平均をグレーにする（まず動かす用）
            Real m = 0;
            for (Real x : c) m += x;
            m /= (Real)N;
            return glm::vec3((float)m);
        }
    };

    // -------------------------
    // 3) 外部公開する Spectrum 型
    // -------------------------
    using Spectrum =
        std::conditional_t<kSpectrumMode == SpectrumMode::RGB,
        RGBSpectrum,
        SampledSpectrum<kSpectralSamples>>;

} // namespace rayt