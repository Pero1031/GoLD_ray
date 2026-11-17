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

    Real ApplySensor(Film& film, const SensorSetting& s, Real gain, std::uint64_t seed) {
        const int w = film.width(), h = film.height();
        long long clipped = 0;
#pragma omp parallel for schedule(static) reduction(+ : clipped)
        for (int y = 0; y < h; ++y) {
            std::mt19937_64 rng(Mix(seed ^ Mix(0x5e4507ull + static_cast<std::uint64_t>(y))));   // a row draws the same numbers on any thread
            for (int x = 0; x < w; ++x) {
                const Vector3 raw = color::filmToCameraRaw(film.getData()[static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x)] * gain);
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
        return w * h > 0 ? Real(clipped) / Real(static_cast<long long>(w) * h) : Real(0);
    }

    bool ParseSensorSpec(const std::string& spec, SensorModel& m, std::string& error) {
        m = SensorModel();
        if (spec == "default" || spec.empty()) return true;
        std::stringstream list(spec);
        for (std::string item; std::getline(list, item, ',');) {
            const size_t eq = item.find('=');
            if (eq == std::string::npos) { error = "--sensor: '" + item + "' is not key=value (keys: well, base, read, adc, bits)"; return false; }
            const std::string key = item.substr(0, eq);
            double v = 0;
            try { v = std::stod(item.substr(eq + 1)); }
            catch (...) { error = "--sensor: '" + item + "' needs a number"; return false; }
            if (key == "well") m.fullWell = Real(v);
            else if (key == "base") m.baseIso = Real(v);
            else if (key == "read") m.readNoise = Real(v);
            else if (key == "adc") m.adcNoise = Real(v);
            else if (key == "bits") m.bits = static_cast<int>(v);
            else { error = "--sensor: unknown key '" + key + "' (keys: well, base, read, adc, bits)"; return false; }
        }
        if (!(m.fullWell >= 1) || !(m.baseIso > 0) || !(m.readNoise >= 0) || !(m.adcNoise >= 0) || m.bits < 8 || m.bits > 16) {
            error = "--sensor: needs well >= 1 electron, base ISO > 0, read and adc noise >= 0, 8 to 16 bits";
            return false;
        }
        return true;
    }

} // namespace rayt
