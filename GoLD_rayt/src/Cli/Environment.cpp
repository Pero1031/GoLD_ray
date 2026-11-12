// The environment (map and light) a scene wants: the --env map, a uniform white, a procedural sky, or nothing.
#include "pch.h"
#include "Cli/Scenes.hpp"
#include "IO/ImageLoader.hpp"
#include "Lights/EnvLight.hpp"

namespace rayt::cli {

Environment MakeEnvironment(const CliOptions& opt, EnvKind kind) {
    Environment env;

    if (kind == EnvKind::Uniform) {
        // Uniform white environment: radiance (1,1,1) in every direction.
        env.map = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        env.light = std::make_shared<EnvLight>(env.map, 1.0);
        std::cout << "[EnvMap] Uniform white environment (furnace test).\n";
    }
    else if (kind == EnvKind::ProceduralSky) {
        // Procedural sky (equirectangular, row 0 = zenith): horizon-to-zenith blue gradient, a warm glow around a low sun and the
        // sun disc itself (angular radius 1.2 degrees, radiance 4000: irradiance about 5.5, versus about 2 from the sky). The disc
        // is larger and dimmer than the real sun on purpose: a smaller, brighter one makes the caustics on the seabed (found by
        // random bounces only) far too spiky. Below the horizon the sky continues as the horizon color, seen in tilted waves.
        const Real sunAz = (opt.camYaw < 1e8 ? glm::radians(opt.camYaw) : Real(0)) + glm::radians(Real(15));
        const Real sunEl = glm::radians(opt.seaSunElev);
        env.peakDir = Vector3(std::cos(sunEl) * std::sin(sunAz), std::sin(sunEl), std::cos(sunEl) * std::cos(sunAz));
        const int W = 2048, H = 1024;
        std::vector<Vector3> px(static_cast<size_t>(W) * static_cast<size_t>(H));
        const Vector3 horizon(0.85, 0.88, 0.95), zenith(0.14, 0.32, 0.75);
        for (int y = 0; y < H; ++y) {
            const Real theta = constants::PI * (Real(y) + Real(0.5)) / Real(H);
            for (int x = 0; x < W; ++x) {
                const Real phi = Real(2) * constants::PI * (Real(x) + Real(0.5)) / Real(W) - constants::PI;
                const Vector3 d(std::sin(theta) * std::cos(phi), std::cos(theta), std::sin(theta) * std::sin(phi));
                Vector3 c = (d.y >= 0) ? glm::mix(horizon, zenith, std::pow(std::min(d.y, Real(1)), Real(0.4))) : horizon * Real(0.8);
                const Real angle = std::acos(std::clamp(glm::dot(d, env.peakDir), Real(-1), Real(1)));
                c += Vector3(1.0, 0.62, 0.30) * Real(1.2) * std::exp(-angle * angle / Real(0.125));
                if (angle < glm::radians(Real(1.2))) c = Vector3(4000.0, 3700.0, 3200.0);
                px[static_cast<size_t>(y) * static_cast<size_t>(W) + static_cast<size_t>(x)] = c;
            }
        }
        env.map = std::make_shared<EnvMap>(Image(W, H, std::move(px)));
        env.light = std::make_shared<EnvLight>(env.map, 1.0);
        std::cout << "[Sea] procedural sky: sun at elevation " << opt.seaSunElev << " deg, azimuth " << glm::degrees(sunAz) << " deg\n";
    }
    else if (kind == EnvKind::None) {
        std::cout << "[EnvMap] none (spd scene: black background, one sphere light).\n";
    }
    else try {
        std::cout << "[EnvMap] Loading: " << opt.env << std::endl;
        auto envImg = io::loadHDR(opt.env);
        if (kind == EnvKind::HdrWithPeak) {
            // Brightest spot of the map (the "sun"): where the caustic comes from. Same mapping as EnvMap::uvToDir (row 0 = north pole).
            const int w = envImg.width(), h = envImg.height();
            auto dirOf = [&](int x, int y) {
                const Real theta = constants::PI * (Real(y) + Real(0.5)) / Real(h);
                const Real phi = Real(2) * constants::PI * (Real(x) + Real(0.5)) / Real(w) - constants::PI;
                return Vector3(std::sin(theta) * std::cos(phi), std::cos(theta), std::sin(theta) * std::sin(phi));
            };
            auto lum = [&](int x, int y) { const Vector3& c = envImg.at(x, y); return 0.2126 * c.x + 0.7152 * c.y + 0.0722 * c.z; };
            int bx = 0, by = 0; double best = -1;
            for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x) if (lum(x, y) > best) { best = lum(x, y); bx = x; by = y; }
            env.peakDir = dirOf(bx, by);
            double total = 0, inCone = 0;   // power (radiance x solid angle) of the whole map and within 6 degrees of the peak
            for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x) {
                const double p = lum(x, y) * std::sin(constants::PI * (Real(y) + Real(0.5)) / Real(h));
                total += p;
                if (glm::dot(dirOf(x, y), env.peakDir) > std::cos(glm::radians(6.0))) inCone += p;
            }
            std::cout << "[Water] brightest spot of " << opt.env << ": radiance " << best << ", elevation "
                << glm::degrees(std::asin(env.peakDir.y)) << " deg, azimuth " << glm::degrees(std::atan2(env.peakDir.x, env.peakDir.z))
                << " deg (about +y, 0 = +z); " << 100.0 * inCone / total << " % of the map's power lies within 6 deg of it\n";
        }
        env.map = std::make_shared<EnvMap>(std::move(envImg));
        std::cout << "[EnvMap] Success.\n";

        // the light of the map (scale 1.0)
        env.light = std::make_shared<EnvLight>(env.map, 1.0);
        std::cout << "[EnvLight] Created from environment map.\n";
    }
    catch (const std::exception& e) {
        std::cerr << "[EnvMap] Failed: " << e.what() << "\n";
        std::cerr << "[EnvMap] Fallback to black background.\n";
        env.map = nullptr;
        env.light = nullptr;
    }
    return env;
}

} // namespace rayt::cli
