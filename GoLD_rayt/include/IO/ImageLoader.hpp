#pragma once

#include <string>
#include <vector>

#include <glm/vec3.hpp>

namespace rayt::io {

    // 線形RGB画像（HDR/LDR共通の内部表現）
    struct Image {
        int width = 0;
        int height = 0;
        std::vector<glm::vec3> pixels; // row-major: y*width + x

        bool isValid() const {
            return width > 0 && height > 0 && !pixels.empty();
        }

        const glm::vec3& at(int x, int y) const {
            return pixels[y * width + x];
        }
    };

    // 拡張子で自動判別（現時点では .hdr のみ対応）
    Image loadImage(const std::string& filename);

    // 明示的に HDR を読む
    Image loadHDR(const std::string& filename);

} //namespace rayt::io