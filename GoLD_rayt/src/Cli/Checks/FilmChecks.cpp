// Display pipeline checks (--scene film): the exposure and the tone map of the written pictures (--exposure, --tonemap).
#include "pch.h"
#include "Cli/Checks.hpp"
#include "Color/ColorSpace.hpp"
#include "Color/ColorTransform.hpp"
#include "Render/Film.hpp"

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

    std::cout << "[Film] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed ? 0 : 2;
}

} // namespace rayt::cli
