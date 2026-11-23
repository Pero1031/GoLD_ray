/**
 * @file Render/LensFlare.hpp
 * @brief Ghosts of a real lens (--ghosts; docs/CAMERA.md): light that the lens surfaces reflect twice on its way to the film.
 *
 * Every glass surface reflects a little of the light that crosses it (Fresnel: about 4 % for bare glass, about 1 % through a single
 * anti-reflection coating). Light reflected at surface j back toward the front and again at an earlier surface i heads for the film once
 * more, out of focus and displaced: a ghost, one for every pair i < j, each an image of the stop (round, or the blades' polygon) or a
 * crescent of it, strung along the line through the picture's centre and the light. Light reflected four times or more is left out (it
 * is weaker by another factor of R^2).
 *
 * The ghosts are traced from the lights, not from the camera: a camera path rarely takes the two reflections, and never reaches a point
 * light. For every light (points on an area light), every pair and a set of wavelengths, a grid of parallel rays covers the lens's
 * entrance; each is traced through the ghost's path without the clear apertures (LensSystem::traceGhost), and the grid's cells become
 * triangles on the film. A triangle spreads the light of its cell over its area (the irradiance follows from the area); within it, the
 * points where its rays crossed the stop and the rims are interpolated and tested per sub-pixel, so the ghosts keep the stop's sharp
 * edges without noise (Hullin, Eisemann, Seidel and Lee 2011, "Physically-Based Real-Time Lens Flare Rendering", traced exactly here
 * instead of their polynomial fits). The result is a picture in the film's units (scene radiance at the centre of the picture, as the
 * camera's rays give it), added to the rendered picture.
 *
 * Approximations: a light is a parallel beam over the lens (its direction from the lens's centre; the lens is millimetres across);
 * the light reaches the lens or not as a whole (one shadow ray); the reflectances are those of unpolarized light at each interface
 * (the polarization picked up at the first reflection is not carried to the second). An environment map casts the ghosts of its sun
 * only (EnvMap::brightSpot: the texels within 2 deg of the brightest above 1 % of it, as one beam), not of the rest of its sky.
 */

#pragma once

#include <utility>
#include <vector>

#include "Core/Types.hpp"
#include "Render/LensSystem.hpp"

namespace rayt {

    class Scene;
    class Camera;

    struct LensFlareSettings {
        LensSystem::Coating coating = LensSystem::Coating::SingleLayer;
        int grid = 192;          ///< rays per side of the entrance grid, per ghost, light point and wavelength
        int wavelengths = 8;     ///< wavelength bins over the visible (both builds: the coating's colour; the spectral build: dispersion too)
        int lightSamples = 4;    ///< points per area light (point and spot lights: 1)
        std::pair<int, int> only = { -1, -1 };   ///< one ghost only (its two surfaces), or all
    };

    struct LensFlareStats {
        int lights = 0;          ///< light points that reach the lens
        int ghosts = 0;          ///< surface pairs
        double rays = 0;         ///< rays traced
        double seconds = 0;
        Real peak = 0;           ///< the brightest pixel of the flare (film units)
    };

    /**
     * @brief The ghosts of the scene's lights through the camera's real lens (prepared: stop, focus, film), as a picture of
     *        width x height in linear sRGB (film units), rows from the top like Film. All zero without a real lens.
     */
    std::vector<Vector3> RenderLensFlare(const Scene& scene, const Camera& camera, int width, int height, const LensFlareSettings& settings,
        LensFlareStats* stats = nullptr);

    /**
     * @brief The same ghosts as point samples: about `rays` rays per light point and ghost (one in each of sqrt(rays)^2 strata of the
     *        entrance, jittered), each landing in one pixel with its exact clipping. Noisy, but independent of the triangles: a check of
     *        RenderLensFlare (--scene camera).
     */
    std::vector<Vector3> RenderLensFlareBySampling(const Scene& scene, const Camera& camera, int width, int height, const LensFlareSettings& settings,
        int rays);

} // namespace rayt
