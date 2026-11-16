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
#include "Render/Camera.hpp"
#include "Render/Film.hpp"
#include "Render/MaskDepth.hpp"
#include "Scene/Aggregate.hpp"
#include "Scene/Scene.hpp"

#include "stb_image.h"

namespace rayt::cli {

int RunFilmTests() {
    bool allPassed = true;
    auto report = [&](const std::string& name, bool ok) {
        std::cout << "[Film]   " << name << "   " << (ok ? "OK" : "FAIL") << "\n";
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

    std::cout << "[Film] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
