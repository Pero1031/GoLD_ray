#pragma once

#include "Material.hpp"

namespace rayt {

    class DiffuseLight : public Material {
    public:
        // color: 光の色と強さ (例: (10, 10, 10) なら非常に明るい白)
        DiffuseLight(const Spectrum& color) : m_emit(color) {}

        // 散乱しない (光を反射するのではなく、光源そのものなので false)
        virtual bool scatter(const Ray& r_in, const SurfaceInteraction& rec,
            Spectrum& attenuation, Ray& scattered) const override {
            return false;
        }

        // 発光する！
        virtual Spectrum emitted(const SurfaceInteraction& rec) const override {
            return m_emit;
        }

    private:
        Spectrum m_emit;
    };

}