// Film output checks (--scene film): the exposure and the tone map of the written pictures (--exposure, --tonemap), the mask and the
// depth map (--mask, --depth-map).
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Color/ColorSpace.hpp"
#include "Color/ColorTransform.hpp"
#include "Geometry/Quad.hpp"
#include "Geometry/Sphere.hpp"
#include "Materials/Dielectric.hpp"
#include "Materials/Lambertian.hpp"
#include "Color/Illuminants.hpp"
#include "Color/SpectralConversion.hpp"
#include "Core/Image.hpp"
#include "IO/EnvMap.hpp"
#include "Lights/EnvLight.hpp"
#include "Render/Camera.hpp"
#include "Render/Film.hpp"
#include "Render/Integrator.hpp"
#include "Render/MaskDepth.hpp"
#include "Render/Sensor.hpp"
#include "IO/DngWriter.hpp"
#include "Scene/Aggregate.hpp"
#include "Scene/Scene.hpp"

#include "stb_image.h"

#include <filesystem>
#include <fstream>
#include <iterator>

namespace rayt::cli {

int RunFilmTests() {
    bool allPassed = true;
    auto report = [&](const std::string& name, bool ok, const std::string& detail = "") {
        std::cout << "[Film]   " << name << (detail.empty() ? "" : "   (" + detail + ")") << "   " << (ok ? "OK" : "FAIL") << "\n";
        allPassed = allPassed && ok;
    };
    auto near = [](Real a, Real b, Real tol) { return std::fabs(a - b) <= tol; };

    // (a) The default is the old picture: Reinhard per channel after 2^EV, then the sRGB curve.
    {
        bool ok = true;
        for (const Real ev : { Real(0), Real(1), Real(-2) })
            for (const Real x : { Real(0), Real(0.25), Real(1), Real(4), Real(100) }) {
                const Real e = x * std::pow(Real(2), ev);
                const Real expected = color::linearToSRGB(e / (Real(1) + e));
                ok = ok && near(color::toDisplaySRGB_Reinhard(Vector3(x), ev)[0], expected, Real(1e-12));
            }
        report("a  reinhard: x 2^EV / (1 + x 2^EV), then the sRGB curve", ok);
    }

    // (b) Clip: linear up to white, white above it; negative and non-finite values are black.
    {
        const auto d = [](Real x, Real ev) { return color::toDisplaySRGB_Clip(Vector3(x, 0, 0), ev)[0]; };
        const bool ok = near(d(0.25, 1), color::linearToSRGB(Real(0.5)), Real(1e-12))
            && near(d(0.6, -1), color::linearToSRGB(Real(0.3)), Real(1e-12))
            // white: 1.055 - 0.055 of the sRGB curve is 1 - 1e-16 in doubles (still code 255 / 65535)
            && near(d(1.0, 0), 1, 1e-12) && near(d(3.0, 0), 1, 1e-12) && near(d(0.6, 1), 1, 1e-12)
            && d(-0.2, 0) == Real(0) && d(std::numeric_limits<Real>::quiet_NaN(), 0) == Real(0)
            && d(std::numeric_limits<Real>::infinity(), 0) == Real(0);
        report("b  clip: 0.25 at +1 EV is 0.5, white and above is 1, negative / NaN / inf is 0", ok);
    }

    // A film of four pixels, each channel a power of two of the others, and the files it writes.
    Film film(4, 1);
    const Vector3 px[4] = { Vector3(0.125, 0.0625, 0.03125), Vector3(0.5, 0.25, 0.125), Vector3(1.0, 0.5, 0.25), Vector3(4.0, 2.0, 1.0) };
    for (int x = 0; x < 4; ++x) film.setPixel(x, 0, px[x]);
    const std::filesystem::path dir = std::filesystem::temp_directory_path();
    const std::string png8 = (dir / "rayt_film_check_8.png").string();
    const std::string png16 = (dir / "rayt_film_check_16.png").string();
    const std::string pfm = (dir / "rayt_film_check.pfm").string();

    // (c) The 8-bit and 16-bit PNG hold the clipped mapping at +1 EV (each value within one code of the formula).
    {
        film.setDisplay(Real(1), Film::ToneMap::Clip);
        film.setPng16(false);
        film.save(png8);
        film.setPng16(true);
        film.save(png16);
        film.setPng16(false);
        bool ok = true;
        int w = 0, h = 0, n = 0;
        unsigned char* b8 = stbi_load(png8.c_str(), &w, &h, &n, 3);
        ok = ok && b8 && w == 4 && h == 1;
        unsigned short* b16 = stbi_load_16(png16.c_str(), &w, &h, &n, 3);
        ok = ok && b16 && w == 4 && h == 1;
        for (int i = 0; ok && i < 4; ++i)
            for (int c = 0; c < 3; ++c) {
                const Real v = color::linearToSRGB(std::min(px[i][c] * Real(2), Real(1)));
                ok = ok && near(Real(b8[i * 3 + c]), v * Real(255), Real(1)) && near(Real(b16[i * 3 + c]), v * Real(65535), Real(1));
            }
        if (b8) stbi_image_free(b8);
        if (b16) stbi_image_free(b16);
        report("c  .png (8 and 16 bit) at +1 EV, clip: the formula within one code", ok);
    }

    // (d) The .pfm ignores the display settings: it holds the radiance as rendered.
    {
        film.setDisplay(Real(3), Film::ToneMap::Clip);
        film.save(pfm);
        std::ifstream f(pfm, std::ios::binary);
        std::string magic;
        int w = 0, h = 0;
        double scale = 0;
        f >> magic >> w >> h >> scale;
        f.get();   // the single whitespace after the header
        std::vector<float> data(12);
        f.read(reinterpret_cast<char*>(data.data()), static_cast<std::streamsize>(data.size() * sizeof(float)));
        bool ok = f.good() && magic == "PF" && w == 4 && h == 1;
        for (int i = 0; ok && i < 4; ++i)
            for (int c = 0; c < 3; ++c) ok = ok && near(Real(data[static_cast<size_t>(i) * 3 + c]), px[i][c], Real(1e-6));
        report("d  .pfm at +3 EV, clip: the radiance unchanged", ok);
    }

    // (e) The clip report: at +1 EV the brightest channels are 0.25, 1, 2 and 8 times white, so two of four pixels clip and the brightest
    //     is 3 stops above white; at -2 EV nothing clips and the brightest sits at white exactly (white itself does not count as clipped).
    {
        film.setDisplay(Real(1), Film::ToneMap::Clip);
        const Film::WhiteClip a = film.whiteClip();
        film.setDisplay(Real(-2), Film::ToneMap::Clip);
        const Film::WhiteClip b = film.whiteClip();
        const Film::WhiteClip c = film.whiteClip(0.5f);   // the 1/spp scale counts like the exposure
        const bool ok = near(a.share, 0.5, 1e-12) && near(std::log2(a.brightest), 3, 1e-9)
            && b.share == Real(0) && near(b.brightest, 1, 1e-9) && near(std::log2(c.brightest), -1, 1e-9);
        std::cout << "[Film]      +1 EV: " << a.share * 100 << " % above white, brightest " << std::log2(a.brightest) << " stops; -2 EV: "
            << b.share * 100 << " %, " << std::log2(b.brightest) << " stops\n";
        report("e  clip report: 2 of 4 pixels and +3 stops at +1 EV, none and 0 stops at -2 EV", ok);
    }

    std::error_code ec;
    for (const std::string& p : { png8, png16, pfm }) std::filesystem::remove(p, ec);

    // ---- The mask and the depth (--mask, --depth-map; Render/MaskDepth.hpp) ----
    // A pinhole camera at the origin looks down -z (65 x 65 pixels, so the centre pixel looks exactly along the axis, 40 degrees).
    // A sphere of radius 1 at 5 m (material "ball") in front of a wall at 10 m (material "wall") that fills the picture.
    const int size = 65;
    const Real vfov = 40;
    const auto ballMat = std::make_shared<Lambertian>(RgbColor(Vector3(0.5)));
    const auto wallMat = std::make_shared<Lambertian>(RgbColor(Vector3(0.5)));
    auto makeScene = [&](const std::shared_ptr<Material>& pane, Real paneZ) {
        auto world = std::make_shared<Aggregate>();
        world->add(std::make_shared<Sphere>(Point3(0, 0, -5), 1.0, ballMat));
        world->add(std::make_shared<Quad>(Point3(-20, -20, -10), Vector3(40, 0, 0), Vector3(0, 40, 0), wallMat));
        if (pane) world->add(std::make_shared<Quad>(Point3(-5, -5, paneZ), Vector3(10, 0, 0), Vector3(0, 10, 0), pane));
        auto scene = std::make_shared<Scene>(world);
        scene->finalize("uniform");
        return scene;
    };
    const Camera camera(Point3(0.0), Point3(0, 0, -1), Vector3(0, 1, 0), vfov, 1.0, 0.0, 1.0);
    auto run = [&](const Scene& scene, const std::vector<const Material*>& mats, Film& f) {
        f.enableMask();
        f.enableDepth();
        RenderMaskAndDepth(scene, camera, f, mats, 16);
    };
    auto maskSum = [](const Film& f) { double s = 0; for (const float c : f.mask()) s += c; return s; };
    const int centre = (size / 2) * size + size / 2;   // index of the centre pixel (row-major, top row first)

    // (f) Depth: 4 m at the centre (the sphere's front), 10 m at every corner (the wall: the distance along the axis, not along the ray).
    {
        const auto scene = makeScene(nullptr, 0);
        Film f(size, size);
        run(*scene, { ballMat.get() }, f);
        const auto& d = f.depth();
        const bool ok = near(d[centre], 4, 1e-5) && near(d[0], 10, 1e-5) && near(d[size - 1], 10, 1e-5)
            && near(d[static_cast<size_t>(size) * size - 1], 10, 1e-5) && near(d[static_cast<size_t>(size - 1) * size], 10, 1e-5);
        std::cout << "[Film]      depth: centre " << d[centre] << ", corners " << d[0] << " " << d[size - 1] << "\n";
        report("f  depth: the sphere's front at 4, the wall at 10 in every corner (along the view axis)", ok);

        // (g) Mask of the ball: its silhouette is a circle of radius tan(asin(1/5)) / tan(20 deg) half-heights; the coverage adds up to
        //     its area in pixels, the centre is covered, the wall is not.
        const Real rPx = (Real(1) / std::sqrt(Real(24))) / std::tan(vfov * constants::DEG_TO_RAD / 2) * Real(size) / 2;
        const Real area = constants::PI * rPx * rPx;
        const double sum = maskSum(f);
        std::cout << "[Film]      mask of the ball: " << sum << " pixels, the silhouette " << area << "\n";
        report("g  mask: the coverage adds up to the silhouette's area (0.5 %), 1 at the centre, 0 on the wall",
            near(Real(sum), area, area * Real(0.005)) && f.mask()[centre] == 1.0f && f.mask()[0] == 0.0f);
    }

    // (h) Two materials: the ball and the wall together cover the whole picture.
    {
        const auto scene = makeScene(nullptr, 0);
        Film f(size, size);
        run(*scene, { ballMat.get(), wallMat.get() }, f);
        report("h  mask of two materials: the union covers every pixel", maskSum(f) == double(size) * size);
    }

    // (i) A translucent layer in front (opacity 0.25 at 2 m): it covers a quarter of every pixel on average, and the depth looks through it;
    //     at opacity 0.75 the depth stops at it.
    {
        auto layer = std::make_shared<Lambertian>(RgbColor(Vector3(0.5)));
        layer->setTranslucent(0.25);
        Film f(size, size);
        run(*makeScene(layer, -2), { layer.get() }, f);
        const double mean = maskSum(f) / (double(size) * size);
        layer->setTranslucent(0.75);
        Film g(size, size);
        run(*makeScene(layer, -2), { layer.get() }, g);
        std::cout << "[Film]      translucent layer: mean coverage " << mean << " (0.25), depth at the centre " << f.depth()[centre] << " (4), then "
            << g.depth()[centre] << " (2)\n";
        report("i  translucent layer: coverage 0.25 on average, the depth through it below opacity 1/2 and at it above",
            near(Real(mean), 0.25, 0.005) && near(f.depth()[centre], 4, 1e-5) && near(g.depth()[centre], 2, 1e-5));
    }

    // (j) An index-matched boundary (n = 1) is crossed; a glass pane (n = 1.5) is a surface.
    {
        Film f(size, size), g(size, size);
        run(*makeScene(std::make_shared<Dielectric>(1.0, 0.0), -3), { ballMat.get() }, f);
        run(*makeScene(std::make_shared<Dielectric>(1.5, 0.0), -3), { ballMat.get() }, g);
        std::cout << "[Film]      depth at the centre: index-matched pane " << f.depth()[centre] << " (4), glass " << g.depth()[centre] << " (3)\n";
        report("j  index-matched boundary crossed, glass a surface (depth and mask)",
            near(f.depth()[centre], 4, 1e-5) && f.mask()[centre] == 1.0f && near(g.depth()[centre], 3, 1e-5) && g.mask()[centre] == 0.0f);
    }

    // (k) The files: the mask as a one-channel 8- and 16-bit PNG, the depth as a one-channel PFM ("Pf"), both read back.
    {
        const auto scene = makeScene(nullptr, 0);
        Film f(size, size);
        run(*scene, { ballMat.get() }, f);
        const std::string m8 = (dir / "rayt_mask_check_8.png").string(), m16 = (dir / "rayt_mask_check_16.png").string();
        const std::string dp = (dir / "rayt_depth_check.pfm").string();
        f.saveMask(m8);
        f.setPng16(true);
        f.saveMask(m16);
        f.saveDepth(dp);
        int w = 0, h = 0, n8 = 0, n16 = 0;
        unsigned char* b8 = stbi_load(m8.c_str(), &w, &h, &n8, 0);
        unsigned short* b16 = stbi_load_16(m16.c_str(), &w, &h, &n16, 0);
        bool ok = b8 && b16 && n8 == 1 && n16 == 1 && w == size && h == size;
        for (size_t i = 0; ok && i < f.mask().size(); ++i)
            ok = near(Real(b8[i]), Real(f.mask()[i]) * 255, 0.5) && near(Real(b16[i]), Real(f.mask()[i]) * 65535, 0.5);
        if (b8) stbi_image_free(b8);
        if (b16) stbi_image_free(b16);
        std::ifstream in(dp, std::ios::binary);
        std::string magic;
        int pw = 0, ph = 0;
        double scale = 0;
        in >> magic >> pw >> ph >> scale;
        in.get();
        std::vector<float> rows(static_cast<size_t>(size) * size);
        in.read(reinterpret_cast<char*>(rows.data()), static_cast<std::streamsize>(rows.size() * sizeof(float)));
        ok = ok && in.good() && magic == "Pf" && pw == size && ph == size;
        for (int y = 0; ok && y < size; ++y)   // PFM rows bottom to top
            for (int x = 0; x < size; ++x)
                ok = ok && rows[static_cast<size_t>(size - 1 - y) * size + x] == f.depth()[static_cast<size_t>(y) * size + x];
        report("k  files: mask .png (8 / 16 bit, one channel), depth .pfm (one channel), read back", ok);
        for (const std::string& p : { m8, m16, dp }) std::filesystem::remove(p, ec);
    }

    // ---- The spectral image (--spectral-image; the spectral build) ----
#if defined(RAYT_SPECTRAL_MODE)
    {
        // A uniform sky of RGB 0.5 (as light: 0.5 x D65), alone or around a red ball, seen by a pinhole; 32 x 32 pixels, 256 spp.
        // 2048 spp: a 10 nm band then gets about 1e5 wavelength samples (noise near 0.4 %, the band's worst of 40 near 1 %), and a fixed
        // seed makes the result the same at every run (at 256 spp the worst band wandered between 1 and 3 % from run to run).
        const int n = 32, spp = 2048;
        sampling::SetRunSeed(20261117);
        auto camera = std::make_shared<Camera>(Point3(0, 0, 3), Point3(0, 0, 0), Vector3(0, 1, 0), 30.0, 1.0, 0.0, 3.0);
        auto render = [&](bool ball, Real first, Real step, int bands) {
            auto world = std::make_shared<Aggregate>();
            if (ball) world->add(std::make_shared<Sphere>(Point3(0.0), 0.6, std::make_shared<Lambertian>(RgbColor(Vector3(0.7, 0.15, 0.1)))));
            auto scene = std::make_shared<Scene>(world);
            scene->setEnvLight(std::make_shared<EnvLight>(std::make_shared<EnvMap>(Image(16, 8, std::vector<Vector3>(16 * 8, Vector3(0.5)))), 1.0));
            scene->finalize("uniform");
            auto film = std::make_unique<Film>(n, n);
            film->enableSpectral(first, step, bands);
            PathIntegrator(camera, 8, spp).render(*scene, *film);
            return film;
        };
        auto bandMean = [&](const Film& f, int k) {   // the band's value averaged over the picture
            double s = 0;
            for (int i = 0; i < n * n; ++i) s += f.spectral()[static_cast<size_t>(k) * n * n + i];
            return s / (n * n);
        };

        // (l) The sky seen directly: every band is 0.5 x the mean of D65 over it (the picture's average, 2 %).
        {
            const auto f = render(false, 380, 10, 40);
            double worst = 0;
            for (int k = 0; k < 40; ++k) {
                double d65 = 0;
                for (int s = 0; s < 100; ++s) d65 += color::evaluateD65(Real(380 + 10 * k + (s + 0.5) * 0.1)) / 100;
                worst = std::max(worst, std::fabs(bandMean(*f, k) / (0.5 * d65) - 1));
            }
            std::cout << "[Film]      sky: largest deviation of a 10 nm band from 0.5 x D65 " << worst * 100 << " %\n";
            report("l  spectral image of a uniform sky: each band is 0.5 x the band's mean of D65 (2 %)", worst < 0.02);
        }

        // (m) A red ball under the sky, in 5 nm bands: the spectral image weighted by the sensor curves over each band gives the RGB
        //     picture rendered from the same paths (the picture's average, 0.5 %).
        {
            const auto f = render(true, 380, 5, 80);
            Vector3 fromBands(0.0), rgb(0.0);
            for (int k = 0; k < 80; ++k) {
                Vector3 w(0.0);   // the sensor curves integrated over the band
                for (int s = 0; s < 50; ++s) w += color::rgbSensorWeights(Real(380 + 5 * k + (s + 0.5) * 0.1)) * Real(0.1);
                fromBands += w * Real(bandMean(*f, k));
            }
            for (int i = 0; i < n * n; ++i) rgb += f->getData()[i];
            rgb /= Real(n * n);
            const Real worst = std::max({ std::fabs(fromBands.x / rgb.x - 1), std::fabs(fromBands.y / rgb.y - 1), std::fabs(fromBands.z / rgb.z - 1) });
            std::cout << "[Film]      red ball: RGB from the bands " << fromBands.x << " " << fromBands.y << " " << fromBands.z << ", the RGB picture "
                << rgb.x << " " << rgb.y << " " << rgb.z << "\n";
            report("m  spectral image of a red ball, converted with the sensor curves, gives the RGB picture (0.5 %)", worst < 0.005);

            // (n) The ENVI files: the header states the size and the band centres, the data reads back as written.
            const std::string stem = (dir / "rayt_spectral_check").string();
            f->saveSpectral(stem);
            std::ifstream hdr(stem + "_spectral.hdr");
            const std::string text((std::istreambuf_iterator<char>(hdr)), std::istreambuf_iterator<char>());
            std::ifstream img(stem + "_spectral.img", std::ios::binary);
            std::vector<float> data(f->spectral().size());
            img.read(reinterpret_cast<char*>(data.data()), static_cast<std::streamsize>(data.size() * sizeof(float)));
            const bool ok = text.rfind("ENVI", 0) == 0 && text.find("samples = 32\n") != std::string::npos && text.find("bands = 80\n") != std::string::npos
                && text.find("interleave = bsq") != std::string::npos && text.find("wavelength = {382.5, 387.5,") != std::string::npos
                && img.good() && data == f->spectral();
            report("n  spectral image files (ENVI header and float data) read back", ok);
            hdr.close();
            img.close();
            for (const char* ext : { "_spectral.hdr", "_spectral.img" }) std::filesystem::remove(stem + ext, ec);
        }
    }

    // ---- The camera response (--camera-response; Color/CameraResponse.hpp) ----
    {
        // Integrals over 380..780 nm with the film's rule (1 nm, trapezoid)
        auto integrate = [](auto&& weight, auto&& spectrum) {
            Vector3 sum(0.0);
            for (int l = 380; l <= 780; ++l) sum += weight(Real(l)) * (spectrum(Real(l)) * ((l == 380 || l == 780) ? Real(0.5) : Real(1)));
            return sum;
        };
        auto cie = [](Real l) { return color::rgbSensorWeights(l); };
        auto camera = [](Real l) { return color::cameraFilmWeight(color::cameraResponse(), l); };
        std::string err;

        // (o) A "camera" whose sensitivities are the CIE film's own weights: its white balance and colour matrix are the identity.
        {
            std::vector<Real> nm;
            std::vector<Vector3> rgb;
            for (int l = 360; l <= 830; ++l) { nm.push_back(Real(l)); rgb.push_back(color::rgbSensorWeights(Real(l))); }
            const bool set = color::setCameraResponse(nm, rgb, "CIE", err);
            const color::CameraResponse& cam = color::cameraResponse();
            double worstWeight = 0, worstMatrix = 0;
            for (int l = 360; l <= 830; ++l) worstWeight = std::max(worstWeight, double(glm::length(camera(Real(l)) - cie(Real(l)))));
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) worstMatrix = std::max(worstMatrix, std::fabs(double(cam.matrix[j][i]) - (i == j ? 1.0 : 0.0)));
            report("o  camera response of the CIE film's own curves: identity white balance and matrix, the same weights",
                set && worstWeight < 1e-9 && worstMatrix < 1e-9 && cam.fitRms < 1e-9);
            color::resetCameraResponse();
        }

        // (p) A camera of three Gaussian sensitivities (R 600, G 535, B 455 nm): D65 becomes (1, 1, 1), a flat grey stays grey, and the raw
        //     file's value of D65 is the integral of the curves (G = 1).
        auto gauss = [](Real l, Real mu, Real sigma) { return std::exp(-0.5 * (l - mu) * (l - mu) / (sigma * sigma)); };
        std::vector<Real> nm;
        std::vector<Vector3> rgb;
        for (int l = 380; l <= 780; l += 5) {
            nm.push_back(Real(l));
            rgb.push_back(Vector3(gauss(l, 600, 35), 0.9 * gauss(l, 535, 35), 0.7 * gauss(l, 455, 25)));
        }
        const bool set = color::setCameraResponse(nm, rgb, "gaussian", err);
        {
            const color::CameraResponse& cam = color::cameraResponse();
            const Vector3 white = integrate(camera, [](Real l) { return color::evaluateD65(l); });
            const Vector3 grey = integrate(camera, [](Real l) { return Real(0.3) * color::evaluateD65(l); });
            const Vector3 raw = color::filmToCameraRaw(Vector3(1.0));
            const Vector3 rawDirect = integrate([&](Real l) { return cam.curves[static_cast<size_t>(int(l) - color::CameraResponse::kFirstNm)]; },
                [](Real l) { return color::evaluateD65(l); });
            const Vector3 colour(0.7, 0.2, 0.05);
            const Vector3 back = color::cameraRawToFilm(color::filmToCameraRaw(colour));   // raw and back (the sensor's way)
            report("p  camera response of Gaussian curves: D65 -> (1, 1, 1), grey stays grey, the raw of D65 as integrated, raw -> film undoes film -> raw",
                set && glm::length(white - Vector3(1.0)) < 1e-9 && glm::length(grey - Vector3(0.3)) < 1e-9
                && glm::length(raw - rawDirect / rawDirect.y) < 1e-9 && glm::length(back - colour) < 1e-9,
                "fit RMS " + std::to_string(cam.fitRms) + ", raw of D65 " + std::to_string(raw.x) + " 1 " + std::to_string(raw.z));
        }

        // (q) Metamerism: a light and the same light plus a "metameric black" (a wave with no component along the three CIE curves) give
        //     the same CIE film RGB but not the same camera RGB.
        {
            std::vector<Vector3> basis;   // the CIE curves, orthonormalized for the film's integral
            auto dot = [&](auto&& f, auto&& g) {
                double s = 0;
                for (int l = 380; l <= 780; ++l) s += double(f(Real(l))) * double(g(Real(l))) * ((l == 380 || l == 780) ? 0.5 : 1.0);
                return s;
            };
            std::vector<std::vector<double>> e;
            for (int c = 0; c < 3; ++c) {
                std::vector<double> v(401);
                for (int l = 380; l <= 780; ++l) v[l - 380] = double(cie(Real(l))[c]);
                for (const auto& u : e) {
                    const double p = dot([&](Real l) { return v[int(l) - 380]; }, [&](Real l) { return u[int(l) - 380]; });
                    for (int k = 0; k < 401; ++k) v[k] -= p * u[k];
                }
                const double n = std::sqrt(dot([&](Real l) { return v[int(l) - 380]; }, [&](Real l) { return v[int(l) - 380]; }));
                for (double& x : v) x /= n;
                e.push_back(v);
            }
            std::vector<double> black(401);
            for (int l = 380; l <= 780; ++l) black[l - 380] = std::sin(2 * constants::PI * (l - 380) / 60.0);
            for (const auto& u : e) {
                const double p = dot([&](Real l) { return black[int(l) - 380]; }, [&](Real l) { return u[int(l) - 380]; });
                for (int k = 0; k < 401; ++k) black[k] -= p * u[k];
            }
            auto light1 = [](Real l) { return color::evaluateD65(l) * color::rgbToReflectance(Vector3(0.5, 0.4, 0.3), l); };
            auto light2 = [&](Real l) { return light1(l) + Real(25) * Real(black[int(l) - 380]); };   // stays positive
            const Vector3 cie1 = integrate(cie, light1), cie2 = integrate(cie, light2);
            const Vector3 cam1 = integrate(camera, light1), cam2 = integrate(camera, light2);
            const double cieDiff = glm::length(cie1 - cie2) / glm::length(cie1), camDiff = glm::length(cam1 - cam2) / glm::length(cam1);
            report("q  metamers: equal for the CIE film, different through the camera",
                cieDiff < 1e-9 && camDiff > 1e-3, "relative differences " + std::to_string(cieDiff) + " and " + std::to_string(camDiff));
        }
        color::resetCameraResponse();
    }
#else
    std::cout << "[Film]   l-q  spectral image and camera response: need the spectral build (skipped)\n";
#endif

    // ---- The sensor (--sensor; Render/Sensor.hpp), the default model: 50000 e- full well at ISO 100, read noise 1.5 + 3 e-, 14 bits ----
    {
        const SensorModel model;
        // N draws of one raw value; the electrons read back from the DN, (DN - black) / K: their mean and variance
        auto draw = [&](Real raw, const SensorSetting& s, int n, std::uint64_t seed, double& mean, double& var, int& atTop) {
            std::mt19937_64 rng(seed);
            double sum = 0, sum2 = 0;
            atTop = 0;
            for (int i = 0; i < n; ++i) {
                const int dn = SenseDn(raw, s, rng);
                if (dn == s.maxDn) ++atTop;
                const double e = (double(dn) - s.black) / double(s.dnPerElectron);
                sum += e;
                sum2 += e * e;
            }
            mean = sum / n;
            var = sum2 / n - mean * mean;
        };
        const int n = 100000;

        // (r) Photon transfer at the base ISO: the mean electrons are the light's (within 4 standard errors) and their variance is the mean
        //     plus the read noise squared plus the rounding of the ADC (1/12 DN^2), within 3 % (the variance of a variance of 1e5 draws:
        //     0.45 %). The slope of the variance against the mean, both in DN, is the conversion gain K (the photon transfer curve).
        {
            const SensorSetting s = MakeSensorSetting(model, 0);
            bool ok = std::fabs(s.whiteElectrons - 50000) < 1e-9;
            double worst = 0, dnMean[2] = { 0, 0 }, dnVar[2] = { 0, 0 };
            int k = 0;
            for (const Real raw : { Real(0), Real(0.0005), Real(0.005), Real(0.05), Real(0.5) }) {
                double mean, var;
                int top;
                draw(raw, s, n, 11 + static_cast<std::uint64_t>(k), mean, var, top);
                const double want = double(raw * s.whiteElectrons);
                const double wantVar = want + double(s.readElectrons * s.readElectrons) + 1.0 / (12.0 * double(s.dnPerElectron * s.dnPerElectron));
                ok = ok && std::fabs(mean - want) < 4 * std::sqrt(wantVar / n) && top == 0;
                worst = std::max(worst, std::fabs(var / wantVar - 1));
                if (raw >= Real(0.05)) { dnMean[k - 3] = mean * double(s.dnPerElectron); dnVar[k - 3] = var * double(s.dnPerElectron * s.dnPerElectron); }
                ++k;
            }
            const double slope = (dnVar[1] - dnVar[0]) / (dnMean[1] - dnMean[0]);
            ok = ok && worst < 0.03 && std::fabs(slope / double(s.dnPerElectron) - 1) < 0.03;
            report("r  photon transfer: electrons as the light, variance = mean + read^2 (3 %), its slope in DN the gain K", ok,
                "variance off by up to " + std::to_string(100 * worst) + " %, slope " + std::to_string(slope) + " DN/e- against "
                + std::to_string(double(s.dnPerElectron)));
        }

        // (s) The ISO amplifies, the light stays: 20 electrons' worth of light at ISO 100 (raw 0.0004) and at ISO 1600 (raw 0.0064, the
        //     picture 16 times brighter) give 20 electrons both times. The signal-to-noise ratio is 20 / sqrt(20 + read^2 + rounding): 3.53
        //     at ISO 100 (read noise 3.35 e-) and 4.24 at ISO 1600 (1.51 e-), so raising the ISO helps a dark picture; within 3 %.
        {
            const SensorSetting s100 = MakeSensorSetting(model, 100), s1600 = MakeSensorSetting(model, 1600);
            double m100, v100, m1600, v1600;
            int top;
            draw(Real(0.0004), s100, n, 21, m100, v100, top);
            draw(Real(0.0064), s1600, n, 22, m1600, v1600, top);
            auto snrOf = [](const SensorSetting& s) {
                return 20 / std::sqrt(20 + double(s.readElectrons * s.readElectrons) + 1.0 / (12.0 * double(s.dnPerElectron * s.dnPerElectron)));
            };
            const double want100 = snrOf(s100), want1600 = snrOf(s1600);
            const double snr100 = m100 / std::sqrt(v100), snr1600 = m1600 / std::sqrt(v1600);
            report("s  ISO 100 and 1600, the same light: 20 electrons both times, SNR 20 / sqrt(20 + read^2) (3 %)",
                std::fabs(m100 - 20) < 0.1 && std::fabs(m1600 - 20) < 0.1 && std::fabs(snr100 / want100 - 1) < 0.03 && std::fabs(snr1600 / want1600 - 1) < 0.03,
                "electrons " + std::to_string(m100) + " / " + std::to_string(m1600) + ", SNR " + std::to_string(snr100) + " / " + std::to_string(snr1600)
                + " against " + std::to_string(want100) + " / " + std::to_string(want1600));
        }

        // (t) The ends of the range. Above white at the base ISO every pixel is at the top DN. At ISO 50 the white would need 100000 e-:
        //     the well is full at 0.5 of it (the mean of raw 0.6 reads 0.5, 0.1 %), and 0.4 is still below it. No light reads 0 on average
        //     (the black level keeps the noise below it), and negative light counts as none.
        {
            const SensorSetting s100 = MakeSensorSetting(model, 100), s50 = MakeSensorSetting(model, 50);
            double m, v, m6, m4, m0, v0, mNeg;
            int top, top6, top0;
            draw(Real(1.2), s100, 1000, 31, m, v, top);
            draw(Real(0.6), s50, 20000, 32, m6, v, top6);
            draw(Real(0.4), s50, 20000, 33, m4, v, top0);
            draw(Real(0), s100, n, 34, m0, v0, top0);
            draw(Real(-1), s100, 1000, 35, mNeg, v, top0);
            const double r6 = m6 / double(s50.whiteElectrons), r4 = m4 / double(s50.whiteElectrons);
            const bool ok = top == 1000 && top6 == 20000 && std::fabs(r6 - 0.5) < 5e-4 && std::fabs(r4 - 0.4) < 5e-4
                && std::fabs(m0) < 4 * std::sqrt(v0 / n) + 0.02 && std::fabs(mNeg) < 0.5;
            report("t  ends of the range: above white all at the top; ISO 50 full at 0.5 of white; no light reads 0, negative light none", ok,
                "raw 0.6 / 0.4 at ISO 50 read " + std::to_string(r6) + " / " + std::to_string(r4) + ", no light " + std::to_string(m0) + " e-");
        }

        // (u) A film through a sensor without the mosaic (ApplySensor): the same seed draws the same picture and another seed another; a
        //     grey of 2500 electrons keeps its mean radiance (0.5 %) and gets the shot noise's spread (5 %); the white above the clip is
        //     the white.
        {
            const int side = 128;
            auto makeFilm = [&](Real value) { Film f(side, side); for (int y = 0; y < side; ++y) for (int x = 0; x < side; ++x) f.setPixel(x, y, Vector3(value)); return f; };
            SensorModel plain = model;
            plain.mosaic = false;
            const SensorSetting s = MakeSensorSetting(plain, 0);
            const Real gain = Real(4);   // the film holds 0.0125: the sensor sees 0.05 (2500 e-)
            Film a = makeFilm(Real(0.0125)), b = makeFilm(Real(0.0125)), c = makeFilm(Real(0.0125)), top = makeFilm(Real(1));
            ApplySensor(a, plain, s, gain, 7);
            ApplySensor(b, plain, s, gain, 7);
            ApplySensor(c, plain, s, gain, 8);
            const Real clipped = ApplySensor(top, plain, s, gain, 7);
            bool same = true, differ = false;
            double sum = 0, sum2 = 0;
            for (int i = 0; i < side * side; ++i) {
                same = same && a.getData()[i] == b.getData()[i];
                differ = differ || a.getData()[i] != c.getData()[i];
                const double g = double(a.getData()[i].y);
                sum += g;
                sum2 += g * g;
            }
            const double mean = sum / (side * side), sd = std::sqrt(sum2 / (side * side) - mean * mean);
            const double wantSd = 0.0125 * std::sqrt(2500.0 + double(s.readElectrons * s.readElectrons)) / 2500.0;
            const double topValue = double(top.getData()[0].y) * gain;
            report("u  a film through the sensor: same seed same picture, the mean kept (0.5 %), the spread of the shot noise (5 %), the white clipped",
                same && differ && std::fabs(mean / 0.0125 - 1) < 0.005 && std::fabs(sd / wantSd - 1) < 0.05 && clipped == Real(1)
                && std::fabs(topValue - double(Real(s.maxDn - s.black) / (s.dnPerElectron * s.whiteElectrons))) < 1e-6,
                "mean off by " + std::to_string(100 * (mean / 0.0125 - 1)) + " %, spread " + std::to_string(sd / wantSd) + " of the expected, white reads "
                + std::to_string(topValue));
        }

        // (v) Demosaicing (bilinear and Malvar-He-Cutler) gives back a uniform colour everywhere, the edges included, and a linear ramp
        //     inside (2 pixels from the edges): both interpolate a plane exactly. At a slanted grey edge (64 x 64) they give the false
        //     colour that colour-demosaicing 0.2.6 (Python, its own implementation of both) gives: RMS of each channel's distance from
        //     the pixel's grey 0.013633 and 0.006779 (Malvar halves it), and at pixel (30, 20) (0.8, 0.725, 0.65) and (0.8, 0.81125, 0.779375).
        {
            const int w = 64, h = 64;
            std::vector<Real> flat(w * h), ramp(w * h), edge(w * h);
            const Vector3 colour(0.3, 0.5, 0.2);
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; ++x) {
                    const int i = y * w + x;
                    flat[i] = colour[BayerChannel(x, y)];
                    ramp[i] = Real(0.1) + Real(0.01) * x + Real(0.005) * y;
                    edge[i] = Real(0.2) + Real(0.6) * std::clamp((Real(x) - 32 - Real(0.25) * (Real(y) - 32)) * Real(0.7) + Real(0.5), Real(0), Real(1));
                }
            bool ok = true;
            double chroma[2] = { 0, 0 };
            Vector3 probe[2];
            for (int k = 0; k < 2; ++k) {
                const DemosaicMethod method = k == 0 ? DemosaicMethod::Bilinear : DemosaicMethod::Malvar;
                const std::vector<Vector3> f = DemosaicBayer(flat, w, h, method), r = DemosaicBayer(ramp, w, h, method), e = DemosaicBayer(edge, w, h, method);
                double sq = 0;
                int count = 0;
                for (int y = 0; y < h; ++y)
                    for (int x = 0; x < w; ++x) {
                        const int i = y * w + x;
                        ok = ok && glm::length(f[i] - colour) < 1e-6;
                        if (x >= 2 && x < w - 2 && y >= 2 && y < h - 2) ok = ok && glm::length(r[i] - Vector3(ramp[i])) < 1e-6;
                        if (x >= 4 && x < w - 4 && y >= 4 && y < h - 4) {
                            const Vector3 d = e[i] - Vector3((e[i].x + e[i].y + e[i].z) / 3);
                            sq += double(glm::dot(d, d));
                            count += 3;
                        }
                    }
                chroma[k] = std::sqrt(sq / count);
                probe[k] = e[20 * w + 30];
            }
            ok = ok && std::fabs(chroma[0] - 0.013633) < 2e-6 && std::fabs(chroma[1] - 0.006779) < 2e-6
                && glm::length(probe[0] - Vector3(0.8, 0.725, 0.65)) < 1e-6 && glm::length(probe[1] - Vector3(0.8, 0.81125, 0.779375)) < 1e-6;
            report("v  demosaicing: a uniform colour and a ramp come back exactly; a grey edge's false colour as colour-demosaicing gives it", ok,
                "RMS false colour " + std::to_string(chroma[0]) + " (bilinear), " + std::to_string(chroma[1]) + " (Malvar)");
        }

        // (w) A film through the sensor with the mosaic: the recorded DN holds R at (0, 0), G at (1, 0) and (0, 1), B at (1, 1) (the mean
        //     DN of each kind of site is that of its channel, within 4 standard errors), the demosaiced picture keeps the colour (0.5 %),
        //     and without the mosaic there is no recorded mosaic.
        {
            const int side = 128;
            Film f(side, side), g(side, side);
            const Vector3 colour(0.3, 0.5, 0.2);   // raw values (no camera response): 15000, 25000 and 10000 electrons
            for (int y = 0; y < side; ++y) for (int x = 0; x < side; ++x) { f.setPixel(x, y, colour); g.setPixel(x, y, colour); }
            const SensorSetting s = MakeSensorSetting(model, 0);
            std::vector<std::uint16_t> dn, none;
            ApplySensor(f, model, s, 1, 9, &dn);
            SensorModel plain = model;
            plain.mosaic = false;
            ApplySensor(g, plain, s, 1, 9, &none);
            double sum[4] = { 0, 0, 0, 0 }, n[4] = { 0, 0, 0, 0 };
            Vector3 mean(0.0);
            for (int y = 0; y < side; ++y)
                for (int x = 0; x < side; ++x) {
                    const int site = (y & 1) * 2 + (x & 1);   // 0 R, 1 G of the R row, 2 G of the B row, 3 B
                    sum[site] += dn.size() == size_t(side * side) ? double(dn[size_t(y * side + x)]) : 0.0;
                    n[site] += 1;
                    mean += f.getData()[y * side + x] / Real(side * side);
                }
            bool ok = dn.size() == size_t(side * side) && none.empty();
            const Real want[4] = { colour.x, colour.y, colour.y, colour.z };
            for (int k = 0; k < 4; ++k) {
                const double electrons = double(want[k] * s.whiteElectrons), expected = s.black + electrons * double(s.dnPerElectron);
                const double se = std::sqrt(electrons + double(s.readElectrons * s.readElectrons)) * double(s.dnPerElectron) / std::sqrt(n[k]);
                ok = ok && std::fabs(sum[k] / n[k] - expected) < 4 * se;
            }
            const Vector3 rel = mean / colour - Vector3(1.0);
            ok = ok && std::fabs(rel.x) < 0.005 && std::fabs(rel.y) < 0.005 && std::fabs(rel.z) < 0.005;
            report("w  the mosaic: R G / G B as recorded (4 standard errors), the demosaiced colour kept (0.5 %), none without it", ok,
                "mean DN of the sites " + std::to_string(sum[0] / n[0]) + " " + std::to_string(sum[1] / n[1]) + " " + std::to_string(sum[2] / n[2]) + " "
                + std::to_string(sum[3] / n[3]));
        }

        // (x) The DNG file, read back: a little-endian TIFF ("II", 42) whose directory's tags rise, a CFA image (32803) of R G / G B, its
        //     black and white levels and colour matrix as given, and the mosaic itself where StripOffsets points.
        {
            const int w = 4, h = 2;
            const std::vector<std::uint16_t> mosaic = { 512, 1000, 2000, 3000, 4000, 5000, 6000, 16383 };
            io::DngInfo info;
            info.black = 512;
            info.white = 16383;
            info.colorMatrix[0] = 0.75;
            const std::string path = (std::filesystem::temp_directory_path() / "gold_rayt_film_check.dng").string();
            bool ok = io::WriteDng(path, w, h, mosaic, info);
            std::ifstream in(path, std::ios::binary);
            const std::vector<unsigned char> b((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
            in.close();
            std::filesystem::remove(path);
            auto u16 = [&](size_t o) { return o + 2 <= b.size() ? unsigned(b[o] | (b[o + 1] << 8)) : 0u; };
            auto u32 = [&](size_t o) { return o + 4 <= b.size() ? unsigned(b[o] | (b[o + 1] << 8) | (b[o + 2] << 16) | (unsigned(b[o + 3]) << 24)) : 0u; };
            ok = ok && b.size() > 8 && b[0] == 'I' && b[1] == 'I' && u16(2) == 42;
            const size_t ifd = u32(4);
            const unsigned count = u16(ifd);
            unsigned last = 0, photometric = 0, black = 0, white = 0, stripOffset = 0;
            bool rising = true, cfa = false, matrix = false;
            for (unsigned k = 0; k < count && ok; ++k) {
                const size_t e = ifd + 2 + 12 * size_t(k);
                const unsigned tag = u16(e), value = u32(e + 8);
                rising = rising && tag > last;
                last = tag;
                if (tag == 262) photometric = u16(e + 8);
                if (tag == 273) stripOffset = value;
                if (tag == 33422) cfa = b[e + 8] == 0 && b[e + 9] == 1 && b[e + 10] == 1 && b[e + 11] == 2;
                if (tag == 50714) black = value;
                if (tag == 50717) white = value;
                if (tag == 50721) matrix = int(u32(value)) == 7500 && u32(value + 4) == 10000;   // 0.75 as 7500 / 10000
            }
            bool data = stripOffset > 0;
            for (int i = 0; i < w * h && data; ++i) data = u16(stripOffset + 2 * size_t(i)) == mosaic[size_t(i)];
            ok = ok && rising && photometric == 32803 && cfa && black == 512 && white == 16383 && matrix && data;
            report("x  DNG read back: TIFF header, rising tags, CFA R G / G B, black / white levels, colour matrix, the mosaic", ok,
                std::to_string(count) + " tags, " + std::to_string(b.size()) + " bytes");
        }
    }

    std::cout << "[Film] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
