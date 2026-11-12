/**
 * @file Render/SpecularBody.hpp
 * @brief A flat-faced glass body that the specular manifold NEE (Render/SpecularNee.hpp) can solve: its mesh and its glass.
 *
 * A scene lists these bodies (Cli/Scene.hpp, SceneSetup::manifoldBodies) without needing the prototype's header.
 */

#pragma once

#include <memory>

namespace rayt {

    struct TriangleMesh;
    class Dielectric;

    struct SpecularBody {
        std::shared_ptr<const TriangleMesh> mesh;
        std::shared_ptr<const Dielectric> glass;
    };

} // namespace rayt
