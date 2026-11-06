// The materials shared by the scenes: floor, Johnson & Christy gold at three roughnesses, glass, a white reference.
#include "pch.h"
#include "Cli/Scenes.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

StdMaterials MakeStandardMaterials(const CliOptions& opt) {
    StdMaterials m;
    // 床用
    m.floor = std::make_shared<Lambertian>(RgbColor(0.5, 0.5, 0.5)); // 少し暗くして反射を目立たせる

    // --------------------------------------------------
    // 金の光学定数 (Au) をCSVからロード
    // The table is handed to the material as-is: spectral mode evaluates it per wavelength,
    // RGB mode reduces it to RGB eta/k by CMF integration.
    // --------------------------------------------------
    m.auTable = std::make_shared<rayt::color::ComplexIorTable>();
    std::string auErr;

    if (!rayt::io::loadComplexIor_RefractiveIndexInfoCsv(AU_IOR_CSV_PATH, *m.auTable, &auErr)) {
        std::cerr << "[Au IOR] Failed: " << auErr << "\n";
        std::cerr << "[Au IOR] Fallback to hard-coded eta/k.\n";

        // Fallback: a tiny table built from the old hard-coded RGB constants (450/550/650 nm).
        m.auTable->assign({
            { 450.0, 1.45, 1.77 },
            { 550.0, 0.42, 2.45 },
            { 650.0, 0.16, 3.48 },
        });
    }
    else {
        std::cout << "[Au IOR] Loaded: " << AU_IOR_CSV_PATH << "\n";
    }

    // Log the RGB reduction of the table (what RGB mode renders with).
    {
        Vector3 eta_rgb(1.0), k_rgb(0.0);
        rayt::color::integrateIorToRgb_CMF(*m.auTable, eta_rgb, k_rgb);
        std::cout << "[Au IOR] RGB (CMF-integrated) eta=("
            << eta_rgb.x << "," << eta_rgb.y << "," << eta_rgb.z << ") "
            << "k=(" << k_rgb.x << "," << k_rgb.y << "," << k_rgb.z << ")\n";
        std::cout << "[Au IOR] F0 (normal incidence) = ("
            << rayt::color::conductorNormalIncidenceReflectance(eta_rgb.x, k_rgb.x) << ","
            << rayt::color::conductorNormalIncidenceReflectance(eta_rgb.y, k_rgb.y) << ","
            << rayt::color::conductorNormalIncidenceReflectance(eta_rgb.z, k_rgb.z) << ")\n";
    }

    // ★比較用: 3段階の粗さを作成
    // 0.01: ほぼ鏡 (MirrorConductorと比較用)
    // 0.20: 少しぼやけた金属
    // 0.50: マットな金属（ブラスト仕上げ風）
    m.goldSmooth = std::make_shared<RoughConductor>(m.auTable, 0.01);
    m.goldMedium = std::make_shared<RoughConductor>(m.auTable, 0.20);
    m.goldRough = std::make_shared<RoughConductor>(m.auTable, 0.50);

    // Glass and reference materials for the "glass" and "furnace" scenes.
    m.glass = std::make_shared<Dielectric>(1.5, 0.0);                      // 粗さ0 = 完全透明 (no dispersion)
    m.glassDisp = std::make_shared<Dielectric>(1.5, 0.0, 0.0, opt.abbe);   // clear glass with dispersion (spectral mode)
    m.frosted = std::make_shared<Dielectric>(1.5, 0.2);                    // 粗さ0.2 = すりガラス
    m.white = std::make_shared<Lambertian>(RgbColor(1.0));                 // albedo 1 (furnace reference)

    // -------------------------------------------------------------------------
    // 2. 物体の配置 (Scene)
    return m;
}

} // namespace rayt::cli
