#include "IO/ImageLoader.hpp"

#include <stdexcept>
#include <algorithm>
#include <cctype>

#include "stb_image.h"

namespace rayt::io {

    static std::string toLower(std::string s) {
        for (char& c : s) {
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        }
        return s;
    }

    static std::string getExt(const std::string& path) {
        auto p = path.find_last_of('.');
        if (p == std::string::npos) return "";
        return toLower(path.substr(p + 1));
    }

    Image loadHDR(const std::string& filename) {
        int w = 0, h = 0, n = 0;

        // float RGB として読み込む
        float* data = stbi_loadf(filename.c_str(), &w, &h, &n, 3);
        if (!data) {
            throw std::runtime_error(
                "Failed to load HDR: " + filename + " (" + stbi_failure_reason() + ")"
            );
        }

        Image img;
        img.width = w;
        img.height = h;
        img.pixels.resize(static_cast<size_t>(w) * h);

        for (int i = 0; i < w * h; ++i) {
            img.pixels[i] = glm::vec3(
                data[3 * i + 0],
                data[3 * i + 1],
                data[3 * i + 2]
            );
        }

        stbi_image_free(data);
        return img;
    }

    Image loadImage(const std::string& filename) {
        const std::string ext = getExt(filename);

        if (ext == "hdr") {
            return loadHDR(filename);
        }

        throw std::runtime_error("Unsupported image format: " + ext);
    }

} // namespace rayt::io
