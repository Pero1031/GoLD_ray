/**
 * @file include/Render/MaskDepth.hpp
 * @brief An object mask and a depth map from camera rays alone (--mask, --depth-map).
 *
 * A pass of its own after the render, which follows only the first surface each camera ray meets (no lighting). It works with every
 * way the picture is made (path tracing, photon passes, --lvc) and draws its own random numbers, so the picture of a --seed run stays
 * the same with or without it.
 *
 * - Mask: the share of the pixel that surfaces of the chosen materials cover, from n x n stratified samples per pixel through the
 *   camera's lens and shutter (the edges are anti-aliased, defocused and motion-blurred as in the picture). A translucent layer
 *   (opacity < 1) stops a sample with the probability of its opacity, as in the picture; vignetting is not counted (a dark corner is
 *   still covered).
 * - Depth: the distance along the view axis from the camera to the first surface seen through the pixel centre by the pinhole view
 *   (the distortion included, the shutter at its opening); 0 where the ray meets no surface. A translucent layer counts from opacity 1/2.
 *
 * Invisible boundaries (index-matched surfaces of a fog) are crossed; glass, water and mirrors are surfaces. Gaussian splats are not.
 */

#pragma once

#include <vector>

namespace rayt {

    class Scene;
    class Camera;
    class Film;
    class Material;

    /// Fills the film's mask (if Film::enableMask was called) and depth (Film::enableDepth). `samplesPerSide`: n of the n x n mask samples.
    void RenderMaskAndDepth(const Scene& scene, const Camera& camera, Film& film, const std::vector<const Material*>& maskMaterials, int samplesPerSide);

} // namespace rayt
