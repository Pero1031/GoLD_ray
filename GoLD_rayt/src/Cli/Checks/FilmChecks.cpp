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
#include "Scene/Aggregate.hpp"
#include "Scene/Scene.hpp"

#include "stb_image.h"

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
            report("p  camera response of Gaussian curves: D65 -> (1, 1, 1), grey stays grey, the raw of D65 as integrated",
                set && glm::length(white - Vector3(1.0)) < 1e-9 && glm::length(grey - Vector3(0.3)) < 1e-9
                && glm::length(raw - rawDirect / rawDirect.y) < 1e-9,
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

    std::cout << "[Film] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
