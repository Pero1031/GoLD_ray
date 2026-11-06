#pragma once
// The light spectra selectable with --light (d65 | tungsten | candle | led).

#include <string>
#include "Color/EmissionSpectrum.hpp"

namespace rayt::cli {

// Light spectrum chosen with --light (spd and metalprobe scenes). Returns false when the SPD file cannot be read.
bool MakeLightEmission(const std::string& light, EmissionSpectrum& out);
} // namespace rayt::cli
