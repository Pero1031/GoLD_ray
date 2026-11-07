// Command line parsing of the CLI renderer.
#include "pch.h"
#include "Cli/Options.hpp"
#include "Cli/Commands.hpp"
#include "Cli/Help.hpp"
#include "Cli/Scene.hpp"
#include "Cli/SceneDescription.hpp"
#include "Cli/SceneFile.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

namespace {

// "r,g,b" -> a color (a single number is a gray).
Vector3 ParseRgb(const std::string& text, const char* name) {
    std::vector<double> v;
    std::stringstream ss(text);
    for (std::string item; std::getline(ss, item, ',');) v.push_back(std::stod(item));
    if (v.size() == 1) return Vector3(v[0]);
    if (v.size() != 3) throw std::runtime_error(std::string(name) + " needs r,g,b (three numbers separated by commas)");
    return Vector3(v[0], v[1], v[2]);
}

// Applies the options of one command line; throws std::runtime_error with a message for a bad option or value.
void ApplyArgs(const std::vector<std::string>& args, CliOptions& o) {
    {
        for (size_t i = 0; i < args.size(); ++i) {
            const std::string& a = args[i];
            auto value = [&](const char* name) -> std::string {
                if (i + 1 >= args.size()) throw std::runtime_error(std::string("missing value for ") + name);
                return args[++i];
            };
            if (a == "--scene")      o.scene = value("--scene");
            else if (a == "--width") o.width = std::stoi(value("--width"));
            else if (a == "--height") o.height = std::stoi(value("--height"));
            else if (a == "--spp")   o.spp = std::stoi(value("--spp"));
            else if (a == "--depth") o.maxDepth = std::stoi(value("--depth"));
            else if (a == "--env")   o.env = value("--env");
            else if (a == "--abbe")  o.abbe = std::stod(value("--abbe"));
            else if (a == "--out")   o.out = value("--out");
            else if (a == "--light") o.light = value("--light");
            else if (a == "--metal") o.metal = value("--metal");
            else if (a == "--roughness") o.probeRoughness = std::stod(value("--roughness"));
            else if (a == "--melanin") o.melanin = std::stod(value("--melanin"));
            else if (a == "--blood") o.bloodFraction = std::stod(value("--blood"));
            else if (a == "--oxygen") o.oxygen = std::stod(value("--oxygen"));
            else if (a == "--epidermis-mm") o.epidermisMm = std::stod(value("--epidermis-mm"));
            else if (a == "--photons") o.photons = std::stoll(value("--photons"));
            else if (a == "--single-layer") o.singleLayer = true;
            else if (a == "--skin-out") o.skinOut = value("--skin-out");
            else if (a == "--skin-table") o.skinTable = value("--skin-table");
            else if (a == "--skin-matched") o.skinMatched = true;
            else if (a == "--sigma-t") o.sigmaT = std::stod(value("--sigma-t"));
            else if (a == "--albedo") o.mediumAlbedo = std::stod(value("--albedo"));
            else if (a == "--sigma-t-rgb") o.sigmaTRgb = ParseRgb(value("--sigma-t-rgb"), "--sigma-t-rgb");
            else if (a == "--albedo-rgb") o.albedoRgb = ParseRgb(value("--albedo-rgb"), "--albedo-rgb");
            else if (a == "--phase-g") o.phaseG = std::stod(value("--phase-g"));
            else if (a == "--glass-roughness") o.glassRoughness = std::stod(value("--glass-roughness"));
            else if (a == "--mesh") o.mesh = value("--mesh");
            else if (a == "--bg-gray") o.bgGray = std::stod(value("--bg-gray"));
            else if (a == "--volume-cap") o.volumeCap = std::stoi(value("--volume-cap"));
            else if (a == "--absorb-scale") o.absorbScale = std::stod(value("--absorb-scale"));
            else if (a == "--sea-depth") o.seaDepth = std::stod(value("--sea-depth"));
            else if (a == "--sea-steep") o.seaSteep = std::stod(value("--sea-steep"));
            else if (a == "--pool-depth") o.poolDepth = std::stod(value("--pool-depth"));
            else if (a == "--pool-steep") o.poolSteep = std::stod(value("--pool-steep"));
            else if (a == "--pool-sun-elev") o.poolSunElev = std::stod(value("--pool-sun-elev"));
            else if (a == "--pool-sun-az") o.poolSunAz = std::stod(value("--pool-sun-az"));
            else if (a == "--sea-hdr") o.seaHdr = true;
            else if (a == "--sea-sun-elev") o.seaSunElev = std::stod(value("--sea-sun-elev"));
            else if (a == "--sea-seed") o.seaSeed = static_cast<unsigned>(std::stoul(value("--sea-seed")));
            else if (a == "--water-sky") o.waterSky = true;
            else if (a == "--water-lamp") o.waterLamp = true;
            else if (a == "--crystal-inside") o.crystalInside = true;
            else if (a == "--caustic-method") o.causticMethod = value("--caustic-method");
            else if (a == "--caustics") o.causticPhotons = std::stoll(value("--caustics"));
            else if (a == "--caustic-passes") o.causticPasses = std::stoi(value("--caustic-passes"));
            else if (a == "--caustic-radius")o.causticRadius = std::stod(value("--caustic-radius"));
            else if (a == "--cam-elev") o.camElev = std::stod(value("--cam-elev"));
            else if (a == "--cam-yaw") o.camYaw = std::stod(value("--cam-yaw"));
            else if (a == "--light-sampler") o.lightSampler = value("--light-sampler");
            else if (a == "--accel") o.accel = value("--accel");
            else if (a == "--no-ms") o.noMs = true;
            else if (a == "--pmx-index") o.pmxIndex = std::stoi(value("--pmx-index"));
            else if (a == "--pmx-list") o.pmxList = true;
            else if (a == "--pmx-info") o.pmxInfoOnly = true;
            else if (a == "--pmx-all") o.pmxAll = true;
            else if (a == "--pmx-matte") o.pmxMatte = true;
            else if (a == "--sss") o.sss = true;
            else if (a == "--sss-scale") { o.sssScale = std::stod(value("--sss-scale")); o.sss = true; }
            else if (a == "--pmx-light") o.pmxLight = true;
            else if (a == "--pmx-studio") o.pmxStudio = true;
            else if (a == "--studio-look") o.studioLook = value("--studio-look");
            else if (a == "--aov") o.aov = true;
            else if (a == "--denoise") o.denoise = true;
            else if (a == "--clean-aux") o.cleanAux = true;
            else if (a == "--pose-index") o.poseIndex = std::stoi(value("--pose-index"));
            else if (a == "--pose-list") o.poseList = true;
            else if (a == "--pmx-check") o.pmxCheck = true;
            else if (a == "--pmx-bench") o.pmxBench = true;
            else if (a == "--morph-list") o.morphList = true;
            else if (a == "--morph") o.morphSpec = value("--morph");
            else if (a == "--expr") o.exprFile = value("--expr");
            else if (a == "--pmx-zoom") o.pmxZoom = std::stod(value("--pmx-zoom"));
            else if (a == "--pmx-focus") o.pmxFocus = std::stod(value("--pmx-focus"));
            else if (a == "--pmx-yaw") o.pmxYaw = std::stod(value("--pmx-yaw"));
            else if (a == "--wb")    o.whiteBalance = true;
            else if (a == "--trace") {
                const std::string xy = value("--trace");
                const size_t comma = xy.find(',');
                if (comma == std::string::npos) throw std::runtime_error("--trace needs x,y");
                o.traceX = std::stoi(xy.substr(0, comma));
                o.traceY = std::stoi(xy.substr(comma + 1));
            }
            else throw std::runtime_error("unknown option: " + a);
        }
        if (!FindScene(o.scene) && !FindCommand(o.scene))
            throw std::runtime_error("unknown scene or command '" + o.scene + "' (see --list-scenes)");
        if (!o.causticMethod.empty() && o.causticMethod != "manifold" && o.causticMethod != "photons" && o.causticMethod != "none")
            throw std::runtime_error("--caustic-method must be manifold, photons or none");
        if (o.metal != "Au" && o.metal != "Ag" && o.metal != "Cu")
            throw std::runtime_error("--metal must be Au, Ag or Cu");
        if (o.accel != "own" && o.accel != "embree")
            throw std::runtime_error("--accel must be own or embree");
        if (o.studioLook != "default" && o.studioLook != "warm" && o.studioLook != "cool" && o.studioLook != "dramatic" && o.studioLook != "soft")
            throw std::runtime_error("--studio-look must be default, warm, cool, dramatic or soft");
        if (o.lightSampler != "uniform" && o.lightSampler != "power" && o.lightSampler != "spatial" && o.lightSampler != "all")
            throw std::runtime_error("--light-sampler must be uniform, power, spatial or all");
        if (o.light != "d65" && o.light != "tungsten" && o.light != "candle" && o.light != "led")
            throw std::runtime_error("--light must be d65, tungsten, candle or led");
        if (o.width < 16 || o.spp < 1 || o.maxDepth < 1)
            throw std::runtime_error("--width must be >= 16, --spp and --depth >= 1");
        if (o.height < 0) throw std::runtime_error("--height must be >= 0");
    }
}

} // namespace

std::optional<int> ParseArgs(int argc, char** argv, CliOptions& o) {
    std::vector<std::string> args(argv + (argc > 0 ? 1 : 0), argv + argc);
    for (const std::string& a : args) {
        if (a == "--help" || a == "-h") { PrintHelp(std::cout); return 0; }
        if (a == "--list-scenes") { PrintSceneList(std::cout); return 0; }
    }
    try {
        // A scene file comes first; what it stands for is applied before the command line, so the command line wins.
        std::vector<std::string> all;
        if (!args.empty() && args[0].rfind("--", 0) != 0 && args[0].size() > 5 && args[0].substr(args[0].size() - 5) == ".json") {
            std::string error;
            if (!LoadSceneFile(args[0], all, o.sceneDescription, error)) throw std::runtime_error(error);
            o.sceneFile = args[0];
            args.erase(args.begin());
        }
        all.insert(all.end(), args.begin(), args.end());
        ApplyArgs(all, o);
    }
    catch (const std::exception& e) {
        std::cerr << "[Args] " << e.what() << (o.sceneFile.empty() ? "" : "  (the options come from " + o.sceneFile + " and the command line)")
            << "\n(see --help and --list-scenes)\n";
        return 1;
    }
    if (o.out.empty()) {
        o.out = (o.scene == "gold") ? "result_gold_pbr.png" : "result_" + o.scene + ".png";
        if (o.scene == "spd") o.out = "result_spd_" + o.light + (o.whiteBalance ? "_wb" : "") + ".png";
    }
    return std::nullopt;
}

} // namespace rayt::cli
