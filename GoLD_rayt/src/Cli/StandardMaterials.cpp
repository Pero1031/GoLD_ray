// The materials shared by the scenes: floor, Johnson & Christy gold at three roughnesses, glass, a white reference.
#include "pch.h"
#include "Cli/Scenes.hpp"
#include "IO/ComplexIorLoader.hpp"

namespace rayt::cli {

StdMaterials MakeStandardMaterials(const CliOptions& opt) {
    StdMaterials m;
    // For the floor
    m.floor = std::make_shared<Lambertian>(RgbColor(0.5, 0.5, 0.5)); // A little darker, so that the reflections stand out

    // The optical constants of gold (Au): the table is handed to the materials as it is (spectral mode evaluates it per
    // wavelength, RGB mode reduces it to RGB eta/k).
    m.auTable = rayt::io::loadGoldIorTable(AU_IOR_CSV_PATH);

    // For comparison: three levels of roughness
    // 0.01: almost a mirror (to compare with MirrorConductor)
    // 0.20: slightly blurred metal
    // 0.50: matte metal (like a blasted finish)
    m.goldSmooth = std::make_shared<RoughConductor>(m.auTable, 0.01);
    m.goldMedium = std::make_shared<RoughConductor>(m.auTable, 0.20);
    m.goldRough = std::make_shared<RoughConductor>(m.auTable, 0.50);

    // Glass and reference materials for the "glass" and "furnace" scenes.
    m.glass = std::make_shared<Dielectric>(1.5, 0.0);                      // roughness 0 = fully transparent (no dispersion)
    m.glassDisp = std::make_shared<Dielectric>(1.5, 0.0, 0.0, opt.abbe);   // clear glass with dispersion (spectral mode)
    m.frosted = std::make_shared<Dielectric>(1.5, 0.2);                    // roughness 0.2 = frosted glass
    m.white = std::make_shared<Lambertian>(RgbColor(1.0));                 // albedo 1 (furnace reference)

    // -------------------------------------------------------------------------
    // 2. Place the objects (Scene)
    return m;
}

} // namespace rayt::cli
