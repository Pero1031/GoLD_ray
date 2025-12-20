#pragma once

#include "Core.hpp"
#include "Material.hpp"
#include "Utils.hpp"

namespace rayt {

    class Lambertian : public Material {
    public:
        // albedo: (R, G, B)
        Lambertian(const Spectrum& albedo) : m_albedo(albedo) {}

        virtual bool scatter(const Ray& r_in, const SurfaceInteraction& rec,
            Spectrum& attenuation, Ray& scattered) const override {
            // Lambertian reflection direction: normal + random vector
            Vector3 scatter_direction = rec.n + Utils::randomUnitVector();

            // ランダムベクトルがたまたま法線と真逆でゼロベクトルになるのを防ぐ
            if (glm::length(scatter_direction) < 1e-8)
                scatter_direction = rec.n;

            scattered = Ray(rec.p, glm::normalize(scatter_direction));
            attenuation = m_albedo;
            return true;
        }

    private:
        Spectrum m_albedo;
    };

}