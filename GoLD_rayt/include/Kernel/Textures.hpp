#pragma once
/**
 * @file Kernel/Textures.hpp
 * @brief Image textures of the kernel renderer (docs/GPU.md).
 *
 * Core/Texture.hpp (Texture2D::rgb: 8-bit sRGB RGBA, decoded by a 256-entry table) and Core/SpectralTexture.hpp (SpectralTexture::rgb:
 * the linear RGB of each texel's spectrum under D65, what the RGB build uses) step by step: (0, 0) is the top-left corner, the
 * coordinates repeat, texel centres sit at half-integer positions and the four nearest texels are blended bilinearly. With T = double
 * the values are the engine's bits; with T = float they are what the GPU computes.
 */

#include "Kernel/KernelMath.hpp"
#include "Kernel/KernelScene.hpp"

namespace rayt::kernel {

    /// A scene's texels in the precision of the evaluation (KScene::texels8, texelsF, srgbLut).
    template <class T> struct TextureArrays {
        const uint8_t* texels8;
        const T* texelsF;
        const T* srgbLut;
    };

    /// The linear RGB of texture `tex` at (u, v).
    template <class T> RAYT_HD Vec3<T> textureRgb(const KTexture& tex, const TextureArrays<T>& a, T u, T v) {
        const int w = int(tex.width), h = int(tex.height);
        if (w <= 0 || h <= 0) return Vec3<T>(0);
        const T x = (u - kfloor(u)) * T(w) - T(0.5);   // repeat, texel centres at +0.5
        const T y = (v - kfloor(v)) * T(h) - T(0.5);
        const int x0 = int(kfloor(x)), y0 = int(kfloor(y));
        const T fx = x - T(x0), fy = y - T(y0);
        if (tex.kind == kTexSrgb8) {
            T r = T(0), g = T(0), b = T(0);
            for (int dy = 0; dy < 2; ++dy)
                for (int dx = 0; dx < 2; ++dx) {
                    const T wt = (dx ? fx : T(1) - fx) * (dy ? fy : T(1) - fy);
                    const int px = ((x0 + dx) % w + w) % w;
                    const int py = ((y0 + dy) % h + h) % h;
                    const uint8_t* t = a.texels8 + tex.offset + (size_t(py) * size_t(w) + size_t(px)) * 4;
                    r += wt * a.srgbLut[t[0]];
                    g += wt * a.srgbLut[t[1]];
                    b += wt * a.srgbLut[t[2]];
                }
            return Vec3<T>(r, g, b);
        }
        Vec3<T> c(0);
        for (int dy = 0; dy < 2; ++dy)
            for (int dx = 0; dx < 2; ++dx) {
                const T wt = (dx ? fx : T(1) - fx) * (dy ? fy : T(1) - fy);
                const int px = ((x0 + dx) % w + w) % w;
                const int py = ((y0 + dy) % h + h) % h;
                const T* p = a.texelsF + tex.offset + (size_t(py) * size_t(w) + size_t(px)) * 3;
                c += wt * Vec3<T>(p[0], p[1], p[2]);
            }
        return c;
    }

} // namespace rayt::kernel
