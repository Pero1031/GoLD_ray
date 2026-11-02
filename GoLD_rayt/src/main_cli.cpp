// GoLD_rayt.cpp 
//
/// @brief 金属をJohnsonの物理データを基にレンダリング
///
/// 

// C++20を仮定している
// Precompiled Header (Must be first)
#include "pch.h"    

#define GLFW_INCLUDE_NONE
#include <glad/glad.h>
#include <GLFW/glfw3.h>

//#include <imgui.h>
//#include <backends/imgui_impl_glfw.h>
//#include <backends/imgui_impl_opengl3.h>

#include "UI/ImGuiLayer.hpp"
#include "Core/RenderSettings.hpp"

// 自作ヘッダー
// Project Headers
// Core
#include "Core/Core.hpp"
#include "Core/Interaction.hpp"
#include "Core/Sampling.hpp"
#include "Core/Image.hpp"

// accelerates
#include "Accelerators/BVH.hpp" 

// Geometry
#include "Scene/Hittable.hpp"
#include "Scene/Aggregate.hpp"
#include "Geometry/Sphere.hpp"
#include "Geometry/Frame.hpp"
#include "Scene/Primitive.hpp"
#include "Scene/Scene.hpp"

// Renderer
#include "Render/Film.hpp"
#include "Render/Camera.hpp"
#include "Scene/Scene.hpp"
#include "Render/Integrator.hpp"

// Materials
#include "Materials/Material.hpp"
#include "Materials/Lambertian.hpp"
#include "Materials/DiffuseLight.hpp"
#include "Materials/RoughConductor.hpp"
#include "Materials/Dielectric.hpp"

// IO
#include "IO/ImageLoader.hpp"
#include "IO/EnvMap.hpp"
#include "IO/ComplexIorLoader.hpp"
#include "Color/IorRgbApprox.hpp"

#include <filesystem>
#include <map>
#include <stdexcept>
#include <vector>

#include "DebugTools/FrameDebug.hpp"


// 画像生成のためのヘッダー
// マクロを書く必要はなし
#include "stb_image_write.h"

using namespace rayt;

// -----------------------------------------------------------------------------
// Scene Configuration
// -----------------------------------------------------------------------------
const int IMAGE_WIDTH = 800;
const int IMAGE_HEIGHT = 450;      // 16:9 Aspect Ratio
const int SAMPLES_PER_PIXEL = 256; // Higher = less noise, slower
const int MAX_DEPTH = 16;          // Max recursion depth for rays

// env path
const std::string ENV_HDR_PATH = "assets/env/bryanston_park_sunrise_2k.hdr";

const std::string AU_IOR_CSV_PATH = "assets/Metal_data/Johnson_Au.csv";

// -----------------------------------------------------------------------------
// Command line
//   main_cli [--scene gold|glass|furnace] [--width N] [--spp N] [--depth N]
//            [--env file.hdr] [--abbe V] [--out file.png] [--trace x,y]
//
//   --trace x,y : debug aid. Prints a few light paths of the pixel (x,y) (image coordinates,
//                 origin at the top-left) with the object hit at every bounce, then exits.
//
//   gold    : three gold spheres (default, as before)
//   glass   : clear / dispersive / frosted glass spheres on a floor
//   furnace : objects in a uniform white environment. The mean color of every object must match
//             its expected value (1.0 for lossless materials; a value below/above means energy is
//             lost/created, i.e. a bug). Prints OK/FAIL per object and exits with code 2 on FAIL.
//             Debug builds also assert that all ray directions have unit length.
// Without options the behavior is the same as before (gold scene, constants above).
// -----------------------------------------------------------------------------
struct CliOptions {
    std::string scene = "gold";
    int width = IMAGE_WIDTH;           // height follows the 16:9 aspect ratio of the constants
    int spp = SAMPLES_PER_PIXEL;
    int maxDepth = MAX_DEPTH;
    std::string env = ENV_HDR_PATH;
    Real abbe = 20.0;                  // Abbe number of the dispersive glass (smaller = stronger)
    std::string out;                   // default: result_<scene>.png ("gold" keeps result_gold_pbr.png)
    int traceX = -1, traceY = -1;      // --trace x,y : print a few light paths of one pixel and exit
};

static bool ParseArgs(int argc, char** argv, CliOptions& o) {
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            auto value = [&](const char* name) -> std::string {
                if (i + 1 >= argc) throw std::runtime_error(std::string("missing value for ") + name);
                return argv[++i];
            };
            if (a == "--scene")      o.scene = value("--scene");
            else if (a == "--width") o.width = std::stoi(value("--width"));
            else if (a == "--spp")   o.spp = std::stoi(value("--spp"));
            else if (a == "--depth") o.maxDepth = std::stoi(value("--depth"));
            else if (a == "--env")   o.env = value("--env");
            else if (a == "--abbe")  o.abbe = std::stod(value("--abbe"));
            else if (a == "--out")   o.out = value("--out");
            else if (a == "--trace") {
                const std::string xy = value("--trace");
                const size_t comma = xy.find(',');
                if (comma == std::string::npos) throw std::runtime_error("--trace needs x,y");
                o.traceX = std::stoi(xy.substr(0, comma));
                o.traceY = std::stoi(xy.substr(comma + 1));
            }
            else throw std::runtime_error("unknown option: " + a);
        }
        if (o.scene != "gold" && o.scene != "glass" && o.scene != "furnace")
            throw std::runtime_error("--scene must be gold, glass or furnace");
        if (o.width < 16 || o.spp < 1 || o.maxDepth < 1)
            throw std::runtime_error("--width must be >= 16, --spp and --depth >= 1");
    }
    catch (const std::exception& e) {
        std::cerr << "[Args] " << e.what() << "\n"
            << "usage: main_cli [--scene gold|glass|furnace] [--width N] [--spp N] [--depth N]"
            << " [--env file.hdr] [--abbe V] [--out file.png] [--trace x,y]\n";
        return false;
    }
    if (o.out.empty()) o.out = (o.scene == "gold") ? "result_gold_pbr.png" : "result_" + o.scene + ".png";
    return true;
}

// One object of the furnace test with the mean color it must reach.
//   expected: 1.0 for lossless materials. Models that lose energy by construction use the value
//   of an independent numerical integration (conductors) or a measured regression value (frosted
//   glass, which has no closed form here).
struct FurnaceObject {
    std::string label;
    const Material* material;   // nullptr = the background (environment)
    Real expected;
    Real tolerance;             // allowed |mean - expected| per color channel
};

// Mean film color per object seen through the pixel centers (furnace test report).
// Returns false when any object is outside its tolerance (the caller turns that into an exit code).
static bool CheckFurnace(const Scene& scene, const Camera& camera, const Film& film,
    const std::vector<FurnaceObject>& objects)
{
    struct Acc { Vector3 sum = Vector3(0.0); long long n = 0; };
    std::map<const Material*, Acc> acc;   // key nullptr = background (environment)

    const int w = film.width();
    const int h = film.height();
    for (int j = 0; j < h; ++j) {
        for (int i = 0; i < w; ++i) {
            // Same pixel convention as PathIntegrator::render (film row = h - 1 - j).
            const Ray r = camera.getRay((Real(i) + 0.5) / Real(w), (Real(j) + 0.5) / Real(h), Point2(0.5, 0.5));
            SurfaceInteraction rec;
            const Material* m = scene.hit(r, rec) ? rec.matPtr : nullptr;
            Acc& a = acc[m];
            a.sum += film.getData()[static_cast<size_t>(h - 1 - j) * static_cast<size_t>(w) + static_cast<size_t>(i)];
            a.n++;
        }
    }

    std::cout << "\n[Furnace] mean linear RGB per object vs the expected value (+- tolerance, per channel)\n";
    bool allPassed = true;
    for (const FurnaceObject& obj : objects) {
        const auto it = acc.find(obj.material);
        if (it == acc.end() || it->second.n == 0) {
            std::cout << "[Furnace]   " << obj.label << ": no pixels found -> FAIL\n";
            allPassed = false;
            continue;
        }
        const Vector3 mean = it->second.sum / Real(it->second.n);
        const bool ok = std::abs(mean.x - obj.expected) <= obj.tolerance
            && std::abs(mean.y - obj.expected) <= obj.tolerance
            && std::abs(mean.z - obj.expected) <= obj.tolerance;
        allPassed = allPassed && ok;
        std::cout << "[Furnace]   " << obj.label << ": (" << mean.x << ", " << mean.y << ", " << mean.z
            << ")  expected " << obj.expected << " +- " << obj.tolerance
            << "  [" << it->second.n << " px]  " << (ok ? "OK" : "FAIL") << "\n";
    }
    std::cout << "[Furnace] RESULT: " << (allPassed ? "PASS" : "FAIL") << "\n";
    return allPassed;
}

// Debug aid: follow a few light paths of one pixel and print which object each bounce hits.
//   [in]/[out] = the ray hits the surface from outside / from inside, (R)/(T) = reflected / transmitted.
static void TracePixelDebug(const Scene& scene, const Camera& camera, int px, int py, int w, int h,
    int maxDepth, int numPaths,
    const std::vector<std::pair<std::string, const Material*>>& names)
{
    auto label = [&](const Material* m) {
        for (const auto& [n, p] : names) if (p == m) return n;
        return std::string("?");
    };
    const SampledWavelengths lambda;   // RGB-mode default (wavelength-independent materials)
    std::cout << "[Trace] pixel (" << px << ", " << py << ") of " << w << "x" << h << ", " << numPaths << " paths\n";
    for (int k = 0; k < numPaths; ++k) {
        // Same convention as PathIntegrator: image row 0 is the top, camera t = 1 at the top.
        Ray r = camera.getRay((Real(px) + 0.5) / Real(w), Real(1) - (Real(py) + 0.5) / Real(h), Point2(0.5, 0.5));
        std::cout << "  path " << k << ":";
        for (int depth = 0; depth < maxDepth; ++depth) {
            SurfaceInteraction rec;
            if (!scene.hit(r, rec)) {
                std::cout << " -> ENV [ray o=(" << r.o.x << "," << r.o.y << "," << r.o.z << ") d=("
                    << r.d.x << "," << r.d.y << "," << r.d.z << ") tMin=" << r.tMin << " tMax=" << r.tMax << "]";
                break;
            }
            std::cout << " -> " << label(rec.matPtr) << (rec.frontFace ? "[in]" : "[out]")
                << "@(" << rec.p.x << "," << rec.p.y << "," << rec.p.z << ")";
            const auto s = rec.matPtr->sample(rec, -r.d, Point2(sampling::Random(), sampling::Random()), lambda);
            if (!s) { std::cout << " (no sample: " << (rec.frontFace ? "?" : "TIR") << ")"; break; }
            std::cout << (s->isTransmission() ? "(T)" : "(R)");
            r = rec.SpawnRay(s->wi);
        }
        std::cout << "\n";
    }
}

// -----------------------------------------------------------------------------
// Main Entry Point
// -----------------------------------------------------------------------------
int main(int argc, char** argv) {

    CliOptions opt;
    if (!ParseArgs(argc, argv, opt)) return 1;
    const int imageWidth = opt.width;
    const int imageHeight = std::max(1, static_cast<int>(std::lround(double(opt.width) * IMAGE_HEIGHT / IMAGE_WIDTH)));
    std::cout << "[Render] mode = " << kRenderModeName << "  (switch: Core/RenderMode.hpp)\n";
    std::cout << "[Options] scene=" << opt.scene << " " << imageWidth << "x" << imageHeight
        << " spp=" << opt.spp << " maxDepth=" << opt.maxDepth << " out=" << opt.out << "\n";

    // --------------------------------------------------
    // 2. 環境マップ (EnvMap) 読み込み
    // --------------------------------------------------
    std::cout << "CWD = " << std::filesystem::current_path() << std::endl;

    std::shared_ptr<EnvMap> envMap = nullptr;  // 環境マップクラスをインスタンス化
    std::shared_ptr<EnvLight> envLight = nullptr;

    if (opt.scene == "furnace") {
        // Uniform white environment: radiance (1,1,1) in every direction.
        envMap = std::make_shared<EnvMap>(Image(64, 32, std::vector<Vector3>(64 * 32, Vector3(1.0))));
        envLight = std::make_shared<EnvLight>(envMap, 1.0);
        std::cout << "[EnvMap] Uniform white environment (furnace test).\n";
    }
    else try {
        std::cout << "[EnvMap] Loading: " << opt.env << std::endl;
        auto envImg = io::loadHDR(opt.env);   // HDRでダウンロード
        envMap = std::make_shared<EnvMap>(std::move(envImg));
        std::cout << "[EnvMap] Success.\n";

        // EnvMapからEnvLightを作成 (スケール1.0)
        envLight = std::make_shared<EnvLight>(envMap, 1.0);
        std::cout << "[EnvLight] Created from environment map.\n";
    }
    catch (const std::exception& e) {
        std::cerr << "[EnvMap] Failed: " << e.what() << "\n";
        std::cerr << "[EnvMap] Fallback to black background.\n";
        envMap = nullptr;
        envLight = nullptr;
    }

    //-------------------------------------------------------------------------------------------
    // ラフのテスト
    // -------------------------------------------------------------------------
    // 1. マテリアルの作成 (Roughness Test)
    // -------------------------------------------------------------------------

    // 床用
    auto matFloor = std::make_shared<Lambertian>(RgbColor(0.5, 0.5, 0.5)); // 少し暗くして反射を目立たせる

    // --------------------------------------------------
    // 金の光学定数 (Au) をCSVからロード
    // The table is handed to the material as-is: spectral mode evaluates it per wavelength,
    // RGB mode reduces it to RGB eta/k by CMF integration.
    // --------------------------------------------------
    auto auTable = std::make_shared<rayt::color::ComplexIorTable>();
    std::string auErr;

    if (!rayt::io::loadComplexIor_RefractiveIndexInfoCsv(AU_IOR_CSV_PATH, *auTable, &auErr)) {
        std::cerr << "[Au IOR] Failed: " << auErr << "\n";
        std::cerr << "[Au IOR] Fallback to hard-coded eta/k.\n";

        // Fallback: a tiny table built from the old hard-coded RGB constants (450/550/650 nm).
        auTable->samples() = {
            { 450.0, 1.45, 1.77 },
            { 550.0, 0.42, 2.45 },
            { 650.0, 0.16, 3.48 },
        };
    }
    else {
        std::cout << "[Au IOR] Loaded: " << AU_IOR_CSV_PATH << "\n";
    }

    // Log the RGB reduction of the table (what RGB mode renders with).
    {
        Vector3 eta_rgb(1.0), k_rgb(0.0);
        rayt::color::integrateIorToRgb_CMF(*auTable, eta_rgb, k_rgb);
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
    auto matGoldSmooth = std::make_shared<RoughConductor>(auTable, 0.01);
    auto matGoldMedium = std::make_shared<RoughConductor>(auTable, 0.20);
    auto matGoldRough = std::make_shared<RoughConductor>(auTable, 0.50);

    // Glass and reference materials for the "glass" and "furnace" scenes.
    auto matGlass = std::make_shared<Dielectric>(1.5, 0.0);                      // 粗さ0 = 完全透明 (no dispersion)
    auto matGlassDisp = std::make_shared<Dielectric>(1.5, 0.0, 0.0, opt.abbe);   // clear glass with dispersion (spectral mode)
    auto matFrosted = std::make_shared<Dielectric>(1.5, 0.2);                    // 粗さ0.2 = すりガラス
    auto matWhite = std::make_shared<Lambertian>(RgbColor(1.0));                 // albedo 1 (furnace reference)

    // -------------------------------------------------------------------------
    // 2. 物体の配置 (Scene)
    // -------------------------------------------------------------------------
    auto worldList = std::make_shared<Aggregate>();

    std::vector<FurnaceObject> furnaceObjects;                    // what the furnace test checks
    std::vector<std::shared_ptr<Material>> furnaceMaterials;      // owns materials only used there

    if (opt.scene == "gold") {
        // リストに対して add する
        worldList->add(std::make_shared<Sphere>(Point3(-1.2, 0, -1), 0.5, matGoldSmooth));
        worldList->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), 0.5, matGoldMedium));
        worldList->add(std::make_shared<Sphere>(Point3(1.2, 0, -1), 0.5, matGoldRough));

        // 大きな球で床を作製
        worldList->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));
    }
    else if (opt.scene == "glass") {
        // left: clear glass, center: dispersive glass, right: frosted glass
        worldList->add(std::make_shared<Sphere>(Point3(-1.2, 0, -1), 0.5, matGlass));
        worldList->add(std::make_shared<Sphere>(Point3(0.0, 0, -1), 0.5, matGlassDisp));
        worldList->add(std::make_shared<Sphere>(Point3(1.2, 0, -1), 0.5, matFrosted));
        worldList->add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100.0, matFloor));
    }
    else {
        // furnace: no floor, so the only light path is environment <-> object.
        // The two metals use a near-perfect conductor (k = 100, F ~ 0.9996 at every angle), so
        // only the microfacet model's own single-scattering energy loss is left.
        auto perfectMetal = std::make_shared<rayt::color::ComplexIorTable>();
        perfectMetal->samples() = { { 360.0, 1.0, 100.0 }, { 830.0, 1.0, 100.0 } };
        auto matMetalMid = std::make_shared<RoughConductor>(perfectMetal, 0.2);
        auto matMetalRough = std::make_shared<RoughConductor>(perfectMetal, 0.5);
        furnaceMaterials = { matMetalMid, matMetalRough };   // keep them alive until the report

        const Real r = 0.28;
        worldList->add(std::make_shared<Sphere>(Point3(-1.55, 0, -1), r, matGlass));
        worldList->add(std::make_shared<Sphere>(Point3(-0.93, 0, -1), r, matGlassDisp));
        worldList->add(std::make_shared<Sphere>(Point3(-0.31, 0, -1), r, matFrosted));
        worldList->add(std::make_shared<Sphere>(Point3(0.31, 0, -1), r, matWhite));
        worldList->add(std::make_shared<Sphere>(Point3(0.93, 0, -1), r, matMetalMid));
        worldList->add(std::make_shared<Sphere>(Point3(1.55, 0, -1), r, matMetalRough));
        // expected / tolerance:  1.0 for lossless materials;  the near-perfect conductors use an
        // independent numerical integration of single-scattering GGX (0.994 / 0.882);  frosted
        // glass uses the measured value 0.986 (no closed form).
        furnaceObjects = {
            { "background      ", nullptr,              1.000, 0.010 },
            { "clear glass     ", matGlass.get(),       1.000, 0.030 },
            { "dispersive glass", matGlassDisp.get(),   1.000, 0.030 },
            { "white diffuse   ", matWhite.get(),       1.000, 0.030 },
            { "frosted glass   ", matFrosted.get(),     0.986, 0.030 },
            { "metal rough 0.2 ", matMetalMid.get(),    0.994, 0.030 },
            { "metal rough 0.5 ", matMetalRough.get(),  0.882, 0.030 },
        };
    }

    // Area Light (メッシュに紐付いた光源) を作成してシーンに配置
    /*auto lightMat = std::make_shared<DiffuseLight>(Spectrum(50.0, 0.0, 0.0));
    auto lightShape = std::make_shared<Sphere>(Point3(5.0, 10.0, 6.0), 4.0, lightMat);
    auto areaLight = std::make_shared<AreaLight>(lightShape, lightMat, Spectrum(1.0));
    auto lightPrim = std::make_shared<Primitive>(lightShape, lightMat, areaLight);
    worldList->add(lightPrim);*/

    // --------------------------------------------------
    // 5. Scene作成と光源設定
    // --------------------------------------------------
    auto scene = std::make_shared<Scene>(worldList);

    // 環境光を設定（存在する場合）
    if (envLight) {
        scene->setEnvLight(envLight);
        std::cout << "[Scene] Environment light set.\n";
    }

    // オプション: ポイントライトや他の光源を追加
    // 例: シーン内にポイントライトを追加する場合

    /*auto pointLight = std::make_shared<PointLight>(
        Point3(-1.2, 3, 0),           // 位置
        Spectrum(0, 0, 100)       // 強度
    );
    scene->addLight(pointLight);
    std::cout << "[Scene] Point light added.\n";*/

    // エリアライトの設定
    //scene->addLight(areaLight);
    //std::cout << "[Scene] Area light added.\n";

    // IMPORTANT: シーンをファイナライズ（LightSamplerを構築）
    scene->finalize("uniform");
    std::cout << "[Scene] Finalized with "
        << scene->numLights() << " scene lights";
    if (scene->hasEnvLight()) {
        std::cout << " + environment light";
    }
    std::cout << ".\n";

    // -------------------------------------------------------------------------
    // 3. カメラ設定
    // -------------------------------------------------------------------------
    Point3 lookFrom(0, 0.5, 2.5); // 少し高い位置から見下ろす
    Point3 lookAt(0, 0, -1);
    Vector3 vUp(0, 1, 0);

    Real distToFocus = glm::length(lookFrom - lookAt);
    Real aperture = 0.0; // ピンホールカメラ（ボケなし）でテスト

    auto camera = std::make_shared<Camera>(
        lookFrom, lookAt, vUp,
        35.0, // FOV
        double(imageWidth) / imageHeight,
        aperture,
        distToFocus
    );

    if (opt.traceX >= 0) {
        TracePixelDebug(*scene, *camera, opt.traceX, opt.traceY, imageWidth, imageHeight, 12, 24,
            { { "clear", matGlass.get() }, { "disp", matGlassDisp.get() }, { "frosted", matFrosted.get() },
              { "white", matWhite.get() }, { "floor", matFloor.get() },
              { "goldS", matGoldSmooth.get() }, { "goldM", matGoldMedium.get() }, { "goldR", matGoldRough.get() } });
        return 0;
    }

    // --------------------------------------------------
    // 7. レンダラー生成
    // --------------------------------------------------
    Film film(imageWidth, imageHeight);

    // 新しいPathIntegratorのコンストラクタ: (camera, maxDepth, spp)
    auto integrator = std::make_unique<PathIntegrator>(
        camera,
        opt.maxDepth,  // 最大深度
        opt.spp        // プログレッシブ用のspp
    );

    std::cout << "[Integrator] PathIntegrator created (maxDepth="
        << opt.maxDepth << ", spp=" << opt.spp << ").\n";

    // -------------------------------------------------------------------------
    // 6. レンダリング実行 & 保存
    // -------------------------------------------------------------------------
    std::cout << "[Render] Start PBR rendering..." << std::endl;
    integrator->render(*scene, film);

    std::cout << "[Output] Saving images..." << std::endl;
    film.save(opt.out);
    // film.save("result_gold_pbr.hdr");

    bool furnacePassed = true;
    if (opt.scene == "furnace") {
        furnacePassed = CheckFurnace(*scene, *camera, film, furnaceObjects);
    }

    std::cout << "[System] Finished." << std::endl;

    // Furnace test: exit code 2 when an object is outside its expected value, so a script can check it.
    return furnacePassed ? 0 : 2;
}
