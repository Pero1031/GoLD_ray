// The camera's sensor (--sensor); see Render/Sensor.hpp.
#include "pch.h"
#include "Render/Sensor.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "Color/CameraResponse.hpp"
#include "Render/Film.hpp"

namespace rayt {

    namespace {

        /// splitmix64: a well-mixed 64-bit number from any 64-bit input (the seed of a row's random numbers).
        std::uint64_t Mix(std::uint64_t z) {
            z += 0x9e3779b97f4a7c15ull;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
            return z ^ (z >> 31);
        }

        /// An index mirrored about the edge pixel: -1 -> 1, n -> n - 2 (keeps the parity, so the Bayer pattern goes on).
        int Mirror(int i, int n) {
            while (i < 0 || i >= n) i = i < 0 ? -i : 2 * (n - 1) - i;
            return i;
        }

    } // namespace

    SensorSetting MakeSensorSetting(const SensorModel& m, Real iso) {
        SensorSetting s;
        s.iso = iso > 0 ? iso : m.baseIso;
        s.whiteElectrons = m.fullWell * m.baseIso / s.iso;
        const Real after = m.adcNoise * m.baseIso / s.iso;   // the noise after the amplifier, counted in electrons
        s.readElectrons = std::sqrt(m.readNoise * m.readNoise + after * after);
        s.maxDn = (1 << m.bits) - 1;
        s.black = 1 << std::max(m.bits - 5, 0);   // 1/32 of the range, as cameras keep (512 of 16383 at 14 bits)
        s.dnPerElectron = Real(s.maxDn - s.black) / std::min(s.whiteElectrons, m.fullWell);
        return s;
    }

    int SenseDn(Real raw, const SensorSetting& s, std::mt19937_64& rng) {
        const Real mean = std::max(raw, Real(0)) * s.whiteElectrons;
        Real electrons = 0;
        if (mean > 0) electrons = Real(std::poisson_distribution<long long>(static_cast<double>(mean))(rng));
        // The full well is where the ADC's range ends (as a camera sets its white level): a full pixel reads the top DN, flat, without
        // the read noise that a clip of the charge before the readout would leave on it.
        electrons += Real(std::normal_distribution<double>(0.0, static_cast<double>(s.readElectrons))(rng));
        const Real dn = std::round(Real(s.black) + electrons * s.dnPerElectron);
        return static_cast<int>(std::clamp(dn, Real(0), Real(s.maxDn)));
    }

    std::vector<Vector3> DemosaicBayer(const std::vector<Real>& mosaic, int w, int h, DemosaicMethod method) {
        std::vector<Vector3> out(static_cast<size_t>(w) * static_cast<size_t>(h), Vector3(0.0));
        if (w < 3 || h < 3 || mosaic.size() != out.size()) return out;
#pragma omp parallel for schedule(static)
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x) {
                // v(dx, dy): the mosaic around the pixel; the sums of the four at +-1 / +-2 on the axes and of the four diagonals
                auto v = [&](int dx, int dy) { return mosaic[static_cast<size_t>(Mirror(y + dy, h)) * static_cast<size_t>(w) + static_cast<size_t>(Mirror(x + dx, w))]; };
                const Real c = v(0, 0);
                const Real h1 = v(-1, 0) + v(1, 0), v1 = v(0, -1) + v(0, 1);   // nearest left + right, up + down
                const Real h2 = v(-2, 0) + v(2, 0), v2 = v(0, -2) + v(0, 2);   // two away
                const Real d1 = v(-1, -1) + v(1, -1) + v(-1, 1) + v(1, 1);     // the diagonals
                const int site = BayerChannel(x, y);
                Vector3 rgb(0.0);
                if (method == DemosaicMethod::Bilinear) {
                    if (site == 1) {   // green: the other two from the two nearest of each (left/right on an R row are R)
                        const bool rRow = (y & 1) == 0;
                        rgb = Vector3((rRow ? h1 : v1) / 2, c, (rRow ? v1 : h1) / 2);
                    }
                    else {             // red or blue: green from the four nearest, the opposite colour from the four diagonals
                        const Real g = (h1 + v1) / 4, o = d1 / 4;
                        rgb = site == 0 ? Vector3(c, g, o) : Vector3(o, g, c);
                    }
                }
                else {   // Malvar, He and Cutler 2004: the bilinear estimate corrected by the pixel's own colour's Laplacian (weights / 8)
                    if (site == 1) {
                        // the colour of the row's neighbours: 5 c + 4 (row's nearest) - (diagonals) - (row's two away) + 1/2 (column's two away)
                        const Real alongRow = (5 * c + 4 * h1 - d1 - h2 + Real(0.5) * v2) / 8;
                        const Real alongCol = (5 * c + 4 * v1 - d1 - v2 + Real(0.5) * h2) / 8;
                        const bool rRow = (y & 1) == 0;   // on an R row the row's neighbours are R, the column's B
                        rgb = Vector3(rRow ? alongRow : alongCol, c, rRow ? alongCol : alongRow);
                    }
                    else {
                        const Real g = (4 * c + 2 * (h1 + v1) - (h2 + v2)) / 8;
                        const Real o = (6 * c + 2 * d1 - Real(1.5) * (h2 + v2)) / 8;
                        rgb = site == 0 ? Vector3(c, g, o) : Vector3(o, g, c);
                    }
                }
                out[static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x)] = rgb;
            }
        return out;
    }

    Real ApplySensor(Film& film, const SensorModel& m, const SensorSetting& s, Real gain, std::uint64_t seed, std::vector<std::uint16_t>* mosaicDn) {
        const int w = film.width(), h = film.height();
        const size_t n = static_cast<size_t>(w) * static_cast<size_t>(h);
        auto at = [w](int x, int y) { return static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x); };
        long long clipped = 0;
        if (!m.mosaic || w < 3 || h < 3) {   // every pixel records all three channels
#pragma omp parallel for schedule(static) reduction(+ : clipped)
            for (int y = 0; y < h; ++y) {
                std::mt19937_64 rng(Mix(seed ^ Mix(0x5e4507ull + static_cast<std::uint64_t>(y))));   // a row draws the same numbers on any thread
                for (int x = 0; x < w; ++x) {
                    const Vector3 raw = color::filmToCameraRaw(film.getData()[at(x, y)] * gain);
                    Vector3 out;
                    bool full = false;
                    for (int c = 0; c < 3; ++c) {
                        const int dn = SenseDn(raw[c], s, rng);
                        full = full || dn == s.maxDn;
                        out[c] = DnToRaw(dn, s);
                    }
                    if (full) ++clipped;
                    film.setPixel(x, y, color::cameraRawToFilm(out) / gain);
                }
            }
            if (mosaicDn) mosaicDn->clear();
            return n > 0 ? Real(clipped) / Real(static_cast<long long>(n)) : Real(0);
        }

        // The mosaic: each pixel records its own channel. The demosaicing works on white-balanced values (as a raw developer does: its
        // corrections assume that R, G and B of a grey are equal), then the balance is undone and the camera's matrix applied.
        const Vector3 balance = color::cameraWhiteBalanceGains();
        std::vector<Real> balanced(n);
        std::vector<std::uint16_t> dns(n);
#pragma omp parallel for schedule(static) reduction(+ : clipped)
        for (int y = 0; y < h; ++y) {
            std::mt19937_64 rng(Mix(seed ^ Mix(0x5e4507ull + static_cast<std::uint64_t>(y))));
            for (int x = 0; x < w; ++x) {
                const int c = BayerChannel(x, y);
                const Vector3 raw = color::filmToCameraRaw(film.getData()[at(x, y)] * gain);
                const int dn = SenseDn(raw[c], s, rng);
                if (dn == s.maxDn) ++clipped;
                dns[at(x, y)] = static_cast<std::uint16_t>(dn);
                balanced[at(x, y)] = DnToRaw(dn, s) * balance[c];
            }
        }
        const std::vector<Vector3> rgb = DemosaicBayer(balanced, w, h, m.demosaic);
#pragma omp parallel for schedule(static)
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x) film.setPixel(x, y, color::cameraRawToFilm(rgb[at(x, y)] / balance) / gain);
        if (mosaicDn) *mosaicDn = std::move(dns);
        return n > 0 ? Real(clipped) / Real(static_cast<long long>(n)) : Real(0);
    }

    bool ParseSensorSpec(const std::string& spec, SensorModel& m, std::string& error) {
        m = SensorModel();
        if (spec == "default" || spec.empty()) return true;
        const std::string keys = " (keys: well, base, read, adc, bits, cfa=rggb|none, demosaic=malvar|bilinear)";
        std::stringstream list(spec);
        for (std::string item; std::getline(list, item, ',');) {
            const size_t eq = item.find('=');
            if (eq == std::string::npos) { error = "--sensor: '" + item + "' is not key=value" + keys; return false; }
            const std::string key = item.substr(0, eq), text = item.substr(eq + 1);
            if (key == "cfa") {
                if (text == "rggb") m.mosaic = true;
                else if (text == "none") m.mosaic = false;
                else { error = "--sensor: cfa must be rggb or none"; return false; }
                continue;
            }
            if (key == "demosaic") {
                if (text == "malvar") m.demosaic = DemosaicMethod::Malvar;
                else if (text == "bilinear") m.demosaic = DemosaicMethod::Bilinear;
                else { error = "--sensor: demosaic must be malvar or bilinear"; return false; }
                continue;
            }
            double v = 0;
            try { v = std::stod(text); }
            catch (...) { error = "--sensor: '" + item + "' needs a number"; return false; }
            if (key == "well") m.fullWell = Real(v);
            else if (key == "base") m.baseIso = Real(v);
            else if (key == "read") m.readNoise = Real(v);
            else if (key == "adc") m.adcNoise = Real(v);
            else if (key == "bits") m.bits = static_cast<int>(v);
            else { error = "--sensor: unknown key '" + key + "'" + keys; return false; }
        }
        if (!(m.fullWell >= 1) || !(m.baseIso > 0) || !(m.readNoise >= 0) || !(m.adcNoise >= 0) || m.bits < 8 || m.bits > 16) {
            error = "--sensor: needs well >= 1 electron, base ISO > 0, read and adc noise >= 0, 8 to 16 bits";
            return false;
        }
        return true;
    }

} // namespace rayt
