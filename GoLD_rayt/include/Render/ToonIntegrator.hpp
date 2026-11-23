/**
 * @file Render/ToonIntegrator.hpp
 * @brief An experimental non-photorealistic mode (--toon): cel shading with outlines, from the same scene, camera and film
 *        as the path tracer.
 *
 * Every pixel is split into s x s sub-samples (s from the sample count, at most 4). Each sub-sample follows its camera ray to
 * the first surface (reflected by perfect mirrors; smooth glass is seen straight through, tinted) and shades it:
 * - the light's share, sum(strength * max(N.L, 0) * visible) / sum(strength), with N.L from a point light or a 6 x 6 grid on an
 *   area light and one shadow ray to the grid point that lights the surface most (a hard shadow); a scene without lights of its
 *   own is lit by a key light from the brightest part of the environment (else from above the camera);
 * - a plain material: its base color (Material::aovAlbedo, textures included) in two tones, lit and shadow, split at a threshold
 *   of the share;
 * - a material with a ToonLook (an MMD model, --scene pmx): the model's own recipe, computed as MMD does in display (sRGB)
 *   values: (diffuse x light + ambient) x texture, the sphere map (multiplied or added, looked up with the normal as the camera
 *   sees it: a Lit-Sphere / MatCap), times the bottom color of the toon ramp blended to white by saturate(3 x share) (MMD's
 *   self-shadow formula), plus the Blinn-Phong specular. Overlay layers (alpha blended over what lies behind, e.g. a hair
 *   highlight drawn by an additive sphere map), face normals taken from an ellipsoid round the head (so that the nose and the
 *   eye sockets cast no blotchy shade) and eyebrows seen through the fringe are options of the look;
 * - a rim of light on the lit sub-samples within about 2.5 pixels of a silhouette (screen space).
 * Outlines: a plain material draws them where neighbouring sub-samples see different materials, a depth step (the neighbour's
 * point off the tangent plane by more than a few footprints) or a crease; a ToonLook draws them at depth steps only (the
 * silhouettes, as MMD's inverted hull does), in its edge color and width, and not at all with its edge flag off. The sub-samples
 * then average to the pixel, so the lines are anti-aliased.
 *
 * The result is a display-referred picture in [0, 1] (save it with --tonemap clip; the CLI does this itself).
 */

#pragma once

#include <memory>
#include <unordered_map>

#include "Core/Texture.hpp"
#include "Render/Integrator.hpp"

namespace rayt {

    /// How a material looks in the toon mode when it has its own recipe (an MMD material). Colors are display (sRGB) values, as in the file.
    struct ToonLook {
        Vector3 diffuse = Vector3(1.0);
        Vector3 ambient = Vector3(0.0);
        Vector3 specular = Vector3(0.0);
        Real shininess = 0;
        Real alpha = 1;
        std::shared_ptr<const Texture2D> texture;            ///< the color texture (null: none)
        std::shared_ptr<const Texture2D> sphere;             ///< the sphere map (null: none)
        int sphereMode = 0;                                  ///< 1 multiply, 2 add (0 or 3: not used)
        std::shared_ptr<const Texture2D> toon;               ///< the toon ramp (its bottom row is the shadow color; null: toonShadow)
        Vector3 toonShadow = Vector3(0.8);                   ///< the shadow color without a ramp
        bool edge = true;                                    ///< draws an outline at its silhouettes
        Vector3 edgeColor = Vector3(0.0);
        Real edgeSize = 1;                                   ///< times ToonSettings::lineWidth
        bool overlay = false;                                ///< alpha blended over what lies behind it; casts no shadow and draws no outline
        bool faceNormals = false;                            ///< normals from the ellipsoid (faceCenter, faceRadii) instead of the mesh
        Vector3 faceCenter = Vector3(0.0);
        Vector3 faceRadii = Vector3(1.0);
        bool hair = false;                                   ///< hair: an eyebrow just behind it shows through
        bool showThroughHair = false;                        ///< an eyebrow
    };
    using ToonLooks = std::unordered_map<const Material*, ToonLook>;

    struct ToonSettings {
        Real threshold = Real(0.2);                          ///< plain materials: lit where the light's share exceeds this
        Real softness = Real(0.005);                         ///< half width of the ramp between the two tones (the sub-samples anti-alias the edge)
        Vector3 shadowTone = Vector3(0.50, 0.53, 0.70);      ///< plain materials: the base color is multiplied by this in the shadow (a cool shade)
        Real rim = Real(0.35);                               ///< strength of the rim light (0: none)
        Real lineWidth = Real(1.0);                          ///< outline width in pixels (0: no outlines)
        Real creaseDeg = Real(50);                           ///< plain materials: normals further apart than this draw a line
        Real lightColor = Real(154.0 / 255.0);               ///< ToonLook: the light color of MMD (its default 154 / 255)
        Real hairSeeThrough = Real(0.04);                    ///< ToonLook: an eyebrow at most this far (scene units) behind the hair shows through
        std::shared_ptr<const ToonLooks> looks;              ///< the materials with their own recipe (null: none)
    };

    class ToonIntegrator : public Integrator {
    public:
        ToonIntegrator(std::shared_ptr<Camera> camera, int spp, const ToonSettings& settings = {})
            : m_camera(std::move(camera)), m_spp(spp), m_set(settings) {}

        void render(const Scene& scene, Film& film) override;

    private:
        std::shared_ptr<Camera> m_camera;
        int m_spp;
        ToonSettings m_set;
    };

} // namespace rayt
