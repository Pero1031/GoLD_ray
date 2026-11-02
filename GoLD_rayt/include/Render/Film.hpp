/**
 * @file include/Render/Film.hpp
 * @brief Film class for storing per-pixel accumulated radiance (sum buffer).
 *
 * Film represents the camera sensor buffer used by the renderer.
 * It stores high-dynamic-range (HDR) radiance values per pixel as an
 * accumulation (sum) over samples/passes. Any normalization (1/spp),
 * tone mapping, and display encoding (gamma/sRGB) are performed at output time
 * (e.g., UI display or file saving).
 */

#pragma once

#include <string>
#include <vector>
#include <memory>
#include <algorithm>

#include "Core/Types.hpp" 

namespace rayt {

    /**
    * @brief Accumulation buffer for rendered radiance.
    *
    * The internal pixel buffer stores the sum of radiance contributions:
    *   sum(x,y) = Σ_k L_k(x,y)
    *
    * To obtain the unbiased estimate (average), scale by 1/spp at output:
    *   avg = sum / spp
    *
    * This design is convenient for progressive rendering where the application
    * controls the current sample count.
    */
    class Film {
    public:
        /**
         * @brief Initializes the film with a given resolution.
         * @param width  Horizontal resolution in pixels.
         * @param height Vertical resolution in pixels.
         */
        Film(int width, int height);

        /**
         * @brief Sets the spectral radiance for a specific pixel coordinate.
         * * In advanced integrators, this may be extended to an 'addSample' method
         * to facilitate multi-sample accumulation and reconstruction filtering.
         * @param x        Pixel x-coordinate.
         * @param y        Pixel y-coordinate.
         * @param radiance The physical intensity (radiance) to store.
         */
        void setPixel(int x, int y, const Vector3& radiance);

        /// Adds radiance to the accumulation buffer at (x,y) (progressive rendering).
        void addPixel(int x, int y, const Vector3& radiance);

        /// Returns pointer to the accumulation (sum) buffer (linear sRGB, HDR).
        const Vector3* getData() const { return m_pixels.data(); }

        // clear buffer
        void clear();

        // ---- Auxiliary buffers (AOVs) for denoising -------------------------------------------
        // Lighting-free color (albedo) and surface normal of what each pixel sees, each averaged over
        // the pixel's samples. Off by default; the integrator fills them when they are enabled.

        void enableAOVs();
        bool hasAOVs() const { return !m_albedo.empty(); }

        /// Overwrites the AOVs of a pixel (albedo in [0,1], normal in [-1,1]; (0,0,0) = none).
        void setAOV(int x, int y, const Vector3& albedo, const Vector3& normal);

        const Vector3* getAlbedoData() const { return m_albedo.data(); }
        const Vector3* getNormalData() const { return m_normal.data(); }

        /// .hdr: linear floats; .png/.bmp/.jpg: albedo sRGB-encoded.
        void saveAlbedo(const std::string& filename) const;

        /// .hdr: floats in [-1,1]; .png/.bmp/.jpg: remapped to [0,1] (n * 0.5 + 0.5), no encoding.
        void saveNormal(const std::string& filename) const;

        /**
         * @brief Returns the width of the film in pixels.
         */
        int width() const { return m_width; }

        /**
         * @brief Returns the height of the film in pixels.
         */
        int height() const { return m_height; }

        /**
        * @brief Saves the film to a file.
        *
        * @note If the film stores accumulated sums (recommended), use the overload with
        *       an explicit scale (e.g., scale = 1/spp) for correct exposure.
        */
        void save(const std::string& filename) const;

        /**
        * @brief Saves the film to a file with an explicit scale applied.
        *
        * Typical usage:
        *   - Progressive rendering: scale = 1 / currentSpp
        *   - HDR analysis: scale = 1 (raw sum) or 1/spp depending on intent
         */
        void save(const std::string& filename, float scale) const;

    private:
        int m_width;
        int m_height;

        /**
         * @brief Internal buffer for raw pixel data (linear sRGB, floating point HDR).
         * * In spectral mode the integrator converts each path's sampled spectrum to RGB
         * (SpectrumToRGB) before it reaches the film, so the film is always RGB.
         */
        std::vector<Vector3> m_pixels;

        std::vector<Vector3> m_albedo;   ///< empty unless enableAOVs() was called
        std::vector<Vector3> m_normal;

        /// Writes an AOV buffer; `encodeSrgb` for color, `remapSigned` for normals (see saveAlbedo / saveNormal).
        void saveBuffer(const std::string& filename, const std::vector<Vector3>& data,
            bool encodeSrgb, bool remapSigned) const;
    };

} // namespace rayt