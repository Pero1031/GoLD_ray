#pragma once
#include "IO/ImageLoader.hpp"
#include <glm/glm.hpp>
#include <algorithm>
#include <cmath>

namespace rayt {

    class EnvMap {
    public:
        explicit EnvMap(io::Image image)
            : m_img(std::move(image)) {}

        // 環境からの放射輝度
        glm::vec3 eval(const glm::vec3& dir) const {
            float u, v;
            dirToUV(glm::normalize(dir), u, v);
            return sampleBilinear(u, v);
        }

    private:
        io::Image m_img;

        static void dirToUV(const glm::vec3& d, float& u, float& v) {
            constexpr float PI = 3.14159265358979323846f;
            float theta = std::acos(std::clamp(d.y, -1.0f, 1.0f));
            float phi = std::atan2(d.z, d.x);
            u = (phi + PI) / (2.0f * PI);
            v = 1.0f - (theta / PI);
        }

        glm::vec3 texel(int x, int y) const {
            x = (x % m_img.width + m_img.width) % m_img.width; // wrap
            y = std::clamp(y, 0, m_img.height - 1);
            return m_img.at(x, y);
        }

        glm::vec3 sampleBilinear(float u, float v) const {
            if (!m_img.isValid()) return glm::vec3(0.0f);

            u = u - std::floor(u);            // wrap
            v = std::clamp(v, 0.0f, 1.0f);

            float x = u * m_img.width - 0.5f;
            float y = (1.0f - v) * m_img.height - 0.5f;

            int x0 = static_cast<int>(std::floor(x));
            int y0 = static_cast<int>(std::floor(y));
            int x1 = x0 + 1;
            int y1 = y0 + 1;

            float tx = x - std::floor(x);
            float ty = y - std::floor(y);

            glm::vec3 c00 = texel(x0, y0);
            glm::vec3 c10 = texel(x1, y0);
            glm::vec3 c01 = texel(x0, y1);
            glm::vec3 c11 = texel(x1, y1);

            return glm::mix(
                glm::mix(c00, c10, tx),
                glm::mix(c01, c11, tx),
                ty
            );
        }
    };

} // namespace rayt
