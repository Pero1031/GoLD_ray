// The light spectra selectable with --light.
#include "pch.h"
#include "Cli/LightPresets.hpp"
#include "IO/SpectralCurveLoader.hpp"

namespace rayt::cli {

bool MakeLightEmission(const std::string& light, EmissionSpectrum& out) {
    if (light == "d65")           out = EmissionSpectrum(RgbColor(1.0));
    else if (light == "d65uv") {   // D65 with its ultraviolet part (300-380 nm): the light that excites fluorescence (Materials/Fluorescent.hpp)
        color::SpectralCurve curve;
        for (int l = 300; l <= 780; l += 5) curve.add(Real(l), color::evaluateD65Full(Real(l)));
        curve.finalize();
        out = EmissionSpectrum::FromCurve(curve);   // the same luminance and visible shape as "d65"
    }
    else if (light == "tungsten") out = EmissionSpectrum::IlluminantA();
    else if (light == "candle")   out = EmissionSpectrum::Blackbody(1850.0);
    else {
        color::SpectralCurve curve;
        std::string err;
        if (!io::loadSpectralCurveCsv("assets/Illuminants/led_synthetic_coolwhite.csv", curve, &err)) {
            std::cerr << "[SPD] " << err << "\n";
            return false;
        }
        out = EmissionSpectrum::FromCurve(curve);
    }
    return true;
}

} // namespace rayt::cli
