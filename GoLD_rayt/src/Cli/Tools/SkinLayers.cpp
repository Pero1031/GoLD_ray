// Monte Carlo radial profile of a layered skin (--scene skinlayers).
#include "pch.h"
#include "Cli/Tools.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

int RunSkinLayers(const CliOptions& opt) {
    using namespace rayt::skin;
    constexpr int kBins = 128;
    constexpr double kRMin = 1e-5, kRMax = 0.03;
    std::vector<double> wavelengths;
    for (double nm = 420; nm <= 700; nm += 20) wavelengths.push_back(nm);

    // ---- check: a homogeneous semi-infinite medium against diffusion theory (index matched) ----
    {
        const double mua = 100, musp = 4000;
        const auto r = simulateRadialProfile({ Layer{ 0, mua, musp } }, kBins, kRMin, kRMax, std::max(200000ll, opt.photons / 10), 1);
        const double a = musp / (musp + mua), q = std::sqrt(3 * (1 - a));
        const double rdDiffusion = a / 2 * (1 + std::exp(-4.0 / 3.0 * q)) * std::exp(-q);
        std::cout << "[SkinLayers] check, homogeneous mua=100/m mu_s'=4000/m: Monte Carlo Rd = " << r.reflected
            << " (+ " << r.lostOutside << " beyond 3 cm), diffusion theory (approximate) " << rdDiffusion
            << ", energy " << r.reflected + r.lostOutside + r.absorbed << " (must be 1)\n";
    }

    SkinModel model;
    model.melanin = opt.melanin;
    model.blood = opt.bloodFraction;
    model.oxygen = opt.oxygen;
    model.epidermis = opt.epidermisMm * 1e-3;
    model.singleLayer = opt.singleLayer;
    std::cout << "[SkinLayers] " << (model.singleLayer ? "dermis only" : "epidermis + dermis") << ": melanin " << model.melanin
        << ", blood " << model.blood << ", oxygen " << model.oxygen << ", epidermis " << opt.epidermisMm << " mm, "
        << opt.photons << " photons per wavelength\n";

    std::vector<std::vector<double>> density;
    std::vector<double> edges;
    for (double nm : wavelengths) {
        const auto layers = layersAt(model, nm);
        const auto res = simulateRadialProfile(layers, kBins, kRMin, kRMax, opt.photons, 12345 + static_cast<uint64_t>(nm));
        edges = res.edges;
        density.push_back(res.density);
        std::cout << "[SkinLayers]   " << nm << " nm: dermis mu_a " << layers.back().mua / 100 << " /cm, Rd " << res.reflected
            << ", beyond 3 cm " << res.lostOutside << ", energy " << res.reflected + res.lostOutside + res.absorbed << "\n";
    }

    std::filesystem::create_directories(std::filesystem::path(opt.skinOut).parent_path());
    std::ostringstream note;
    note << "melanin " << model.melanin << " blood " << model.blood << " oxygen " << model.oxygen << " epidermis_mm " << opt.epidermisMm
        << (model.singleLayer ? " single-layer" : "") << " photons " << opt.photons;
    if (!saveTableCsv(opt.skinOut, wavelengths, edges, density, note.str())) { std::cerr << "[SkinLayers] cannot write " << opt.skinOut << "\n"; return 1; }
    std::cout << "[SkinLayers] table written: " << opt.skinOut << "\n";

    const auto table = buildLayeredTable(wavelengths, edges, density);
    const char* names[3] = { "R", "G", "B" };
    std::cout << "[SkinLayers] RGB profiles (D65 projection)         albedo  mean r   r50     r90    Burley l (matched)\n";
    for (int c = 0; c < 3; ++c) {
        const auto& t = *table->rgb[c];
        std::cout << "[SkinLayers]   " << names[c] << "                                          " << t.albedo() << "  "
            << t.meanRadius() * 1000 << " mm  " << t.sampleRadius(0.5) * 1000 << " mm  " << t.sampleRadius(0.9) * 1000 << " mm  "
            << matchedScatterDistance(t) * 1000 << " mm\n";
    }
    std::cout << "[SkinLayers] share of the light that leaves within r of the entry point: table vs Burley with the same albedo and mean radius\n";
    for (int c = 0; c < 3; ++c) {
        const auto& t = *table->rgb[c];
        const Real l = matchedScatterDistance(t), d = l / sss::shape(t.albedo());
        std::cout << "[SkinLayers]   " << names[c] << ":";
        for (double rmm : { 0.25, 0.5, 1.0, 2.0, 4.0, 8.0 })
            std::cout << "  " << rmm << "mm " << std::fixed << std::setprecision(2) << t.cdf(rmm * 1e-3) << "/" << sss::cdf(rmm * 1e-3, d);
        std::cout << std::defaultfloat << "\n";
    }
    return 0;
}

} // namespace rayt::cli
