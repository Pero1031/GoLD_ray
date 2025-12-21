#pragma once

#include <string>
#include <vector>
#include <memory>
#include <algorithm>

#include "Core/Core.hpp"  // Include core definitions (Spectrum, Real, etc.)

namespace rayt {

    // The Film class models the sensing device in a simulated camera.
    // It stores the spectral radiance values accumulated for each pixel.
    class Film {
    public:
        // Constructor: Initializes the film with a specific resolution.
        Film(int width, int height);

        // Sets the spectral radiance for a specific pixel.
        // In a full PBRT style, this might be 'addSample' to accumulate multiple rays.
        void setPixel(int x, int y, const Spectrum& radiance);

        // Returns the width of the film in pixels.
        int width() const { return m_width; }

        // Returns the height of the film in pixels.
        int height() const { return m_height; }

        // Saves the current film data to an image file.
        // Automatically detects the format based on the file extension:
        // - ".hdr": Saves as Radiance HDR (high precision, linear data).
        // - ".png" / ".jpg": Applies tone mapping and gamma correction, then saves as 8-bit LDR.
        void save(const std::string& filename) const;

    private:
        int m_width;
        int m_height;

        // Stores raw radiance values.
        // Using 'Spectrum' (float/double) preserves the physical intensity of light,
        // allowing for proper post-processing or HDR output.
        std::vector<Spectrum> m_pixels;
    };

} // namespace rayt