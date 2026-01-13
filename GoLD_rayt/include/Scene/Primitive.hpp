/**
* @file Scene/Primitive.hpp
*/

#pragma once

#include <memory>

#include "Geometry/Shape.hpp"
#include "Materials/Material.hpp"
#include "Lights/Light.hpp"

namespace rayt {

    /**
     * @brief Binds geometry (shape) with shading data (material) and optional area light.
     *
     * Shape is responsible only for intersection and geometric data.
     * Primitive assigns material / area light to the SurfaceInteraction.
     */
    class Primitive final : public Hittable {
    public:
        Primitive(std::shared_ptr<Shape> shape,
            std::shared_ptr<Material> material,
            std::shared_ptr<Light> areaLight = nullptr)
            : m_shape(std::move(shape))
            , m_material(std::move(material))
            , m_areaLight(std::move(areaLight)) {}

        bool hit(const Ray& r, SurfaceInteraction& rec) const override {
            if (!m_shape->hit(r, rec)) return false;

            rec.matPtr = m_material.get();

            // Optional: only set if it is actually an area light.
            // (If you want stricter, check m_areaLight->type() == LightType::Area)
            rec.areaLight = m_areaLight.get();

            return true;
        }

        AABB bounds() const override {
            return m_shape->bounds();
        }

        const Material* material() const { return m_material.get(); }
        const Light* areaLight() const { return m_areaLight.get(); }

    private:
        std::shared_ptr<Shape> m_shape;
        std::shared_ptr<Material> m_material;
        std::shared_ptr<Light> m_areaLight; // optional
    };

} // namespace rayt