/**
 * @file Preview/MmdPreview.hpp
 * @brief A real-time preview of an MMD (PMX) model in the window application: OpenGL rasterization with the toon recipe.
 *
 * The same picture as the toon render (--toon, Render/ToonIntegrator.hpp) as far as rasterization allows, at interactive
 * rates: MMD's self-shadow formula in display (sRGB) values with the model's own toon ramp, sphere maps (looked up with the
 * normal as the camera sees it), Blinn-Phong specular, face normals from an ellipsoid round the head, and the '+' overlay
 * layers alpha blended on top; shadows from a shadow map of one directional key light, outlines by the inverted hull (the
 * edge-flagged materials drawn again with their faces turned inside out and pushed out along the normals by a width in
 * pixels, as MMD does). The geometry is the one the ray tracer builds (cli::PmxToMesh: metres, +Y up, feet on y = 0), so a
 * camera of the preview can be handed to the ray tracer.
 *
 * Not drawn: the eyebrows seen through the fringe and the rim light of the toon render.
 */

#pragma once

#include <glad/glad.h>

#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Core/Types.hpp"
#include "Render/ToonIntegrator.hpp"

namespace rayt::preview {

    /// An orbit camera: it looks at `target` from `distance` away, turned by yaw and pitch (degrees).
    struct OrbitCamera {
        Vector3 target = Vector3(0, 1, 0);
        Real yaw = 0;          ///< 0 = in front of the model (on +z), positive to the right
        Real pitch = 5;        ///< positive = from above
        Real distance = 4;
        Real vfov = 35;        ///< vertical field of view, degrees

        Vector3 eye() const;
        Vector3 forward() const;   ///< unit, from the eye to the target
    };

    struct PreviewOptions {
        Real lightAzimuth = 30;       ///< degrees, 0 = the light comes from the front (+z), positive = from the model's left
        Real lightElevation = 50;     ///< degrees above the horizon
        bool shadows = true;
        bool edges = true;
        bool floor = true;
        Real edgeScale = 1;           ///< times the model's edge sizes
        Vector3 background = Vector3(0.62, 0.66, 0.74);
    };

    class MmdPreview {
    public:
        MmdPreview() = default;
        ~MmdPreview();
        MmdPreview(const MmdPreview&) = delete;
        MmdPreview& operator=(const MmdPreview&) = delete;

        /**
         * @brief Loads a model, optionally posed by a pose file (bones and the expression entries of a .vpd).
         * Needs a current OpenGL context. Replaces the model loaded before.
         * @return false (with a message in `error`) if the model could not be read
         */
        bool load(const std::filesystem::path& pmx, const std::filesystem::path& pose, std::string* error);

        bool loaded() const { return m_indexCount > 0; }

        /// Draws the model into the preview's own framebuffer of this size (rebuilt when the size changes).
        void render(int width, int height, const OrbitCamera& camera, const PreviewOptions& options);

        /// The color texture of the last render (RGBA8, display values), for ImGui::Image.
        GLuint texture() const { return m_colorTexture; }

        /// A camera that frames the whole model (as the ray tracer's PMX scene does).
        OrbitCamera framingCamera() const;

        Real modelHeight() const { return m_height; }
        const std::string& modelName() const { return m_name; }
        /// Milliseconds of the last load: reading, posing, converting, uploading.
        double loadMilliseconds() const { return m_loadMs; }

    private:
        struct DrawRange {
            int first = 0, count = 0;     ///< in indices
            ToonLook look;
            GLuint texture = 0, sphere = 0, toon = 0;
            bool doubleSided = false;
            bool hidden = false;          ///< alpha 0: not drawn
        };

        void release();
        bool ensurePrograms(std::string* error);
        void ensureTargets(int width, int height);
        GLuint uploadTexture(const std::shared_ptr<const Texture2D>& tex, bool clampToEdge);

        // model
        std::string m_name;
        Real m_height = 1;
        Vector3 m_lo = Vector3(0.0), m_hi = Vector3(1.0);
        GLuint m_vao = 0, m_vbo = 0, m_ebo = 0;
        int m_indexCount = 0;
        bool m_ccwFront = true;           ///< the winding the model's front faces have after the conversion
        std::vector<DrawRange> m_ranges;
        std::map<const Texture2D*, GLuint> m_glTextures;
        std::vector<std::shared_ptr<const Texture2D>> m_keepTextures;
        double m_loadMs = 0;

        // floor
        GLuint m_floorVao = 0, m_floorVbo = 0;

        // programs
        GLuint m_mainProgram = 0, m_edgeProgram = 0, m_shadowProgram = 0;

        // targets
        int m_w = 0, m_h = 0;
        GLuint m_msaaFbo = 0, m_msaaColor = 0, m_msaaDepth = 0;
        GLuint m_resolveFbo = 0, m_colorTexture = 0;
        GLuint m_shadowFbo = 0, m_shadowTexture = 0;
        static constexpr int kShadowSize = 2048;
    };

} // namespace rayt::preview
