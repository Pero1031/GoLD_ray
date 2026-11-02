/**
 * @file Core/RenderMode.hpp
 * @brief Compile-time switch: RGB rendering or spectral rendering.
 *
 * This is THE place that decides how light is represented, in the same spirit as
 * "Real = double / float" in Core/Precision.hpp (which includes this header, so every file
 * that knows Real also knows the mode).
 *
 *   RAYT_ENABLE_SPECTRAL 0 : RGB rendering (default).
 *       Every light path carries (R, G, B); Spectrum is a plain RGB triplet. Fast.
 *   RAYT_ENABLE_SPECTRAL 1 : Spectral rendering.
 *       Every light path carries 4 sampled wavelengths (SampledSpectrumT<4>). Glass dispersion
 *       (Dielectric Abbe number) and per-wavelength metal Fresnel become available. Slower
 *       to converge on caustics, but the film and the saved images are still sRGB.
 *
 * HOW TO SWITCH: change the number below, then rebuild everything
 * (Build > Rebuild Solution). Changing the mode changes the Spectrum type in every file, so a
 * normal incremental build is not enough.
 *
 * The mode can also be forced from the compiler command line with /DRAYT_SPECTRAL_MODE
 * (experiments, scripts). The code itself only tests  #if defined(RAYT_SPECTRAL_MODE).
 *
 * @note This header is intentionally dependency-free (see Precision.hpp).
 */

#pragma once

// 0 = RGB rendering, 1 = spectral rendering   <-- the one line to edit
#define RAYT_ENABLE_SPECTRAL 0

#if RAYT_ENABLE_SPECTRAL && !defined(RAYT_SPECTRAL_MODE)
#define RAYT_SPECTRAL_MODE 1
#endif

namespace rayt {

    /// Human-readable name of the mode this binary was compiled for (print it at startup).
#if defined(RAYT_SPECTRAL_MODE)
    inline constexpr const char* kRenderModeName = "SPECTRAL (4 wavelengths per path)";
#else
    inline constexpr const char* kRenderModeName = "RGB";
#endif

} // namespace rayt
