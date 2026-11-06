// MMD (PMX) helpers of the CLI.
#include "pch.h"
#include "Cli/PmxSupport.hpp"
#include "Cli/EngineHeaders.hpp"

namespace rayt::cli {

std::string PathToUtf8(const std::filesystem::path& p) {
    const std::u8string s = p.u8string();
    return std::string(s.begin(), s.end());
}

std::filesystem::path PathFromUtf8(const std::string& s) {
    return std::filesystem::path(std::u8string(s.begin(), s.end()));
}

std::vector<std::filesystem::path> FindMmdFiles(const std::string& subdir, const std::string& extension) {
    namespace fs = std::filesystem;
    std::vector<fs::path> found;
    for (const std::string base : { "../MMD_model", "MMD_model" }) {
        const fs::path dir = subdir.empty() ? fs::path(base) : fs::path(base) / subdir;
        std::error_code ec;
        if (!fs::exists(dir, ec)) continue;
        for (const auto& e : fs::recursive_directory_iterator(dir, ec)) {
            if (e.is_regular_file() && e.path().extension() == extension) found.push_back(e.path());
        }
        if (!found.empty()) break;
    }
    std::sort(found.begin(), found.end());
    return found;
}

std::vector<std::filesystem::path> FindPmxFiles() { return FindMmdFiles("", ".pmx"); }

Real SrgbToLinear(Real c) { return std::pow(std::clamp(c, Real(0), Real(1)), Real(2.2)); }

void PrintPmxInfo(const io::PmxModel& m, const std::filesystem::path& file) {
    namespace fs = std::filesystem;
    std::cout << "[PMX] " << PathToUtf8(file.filename()) << "   model name: " << m.name << "  (PMX " << m.version
        << ", " << m.vertices.size() << " vertices, " << m.indices.size() / 3 << " triangles, " << m.textures.size()
        << " textures, " << m.materials.size() << " materials, " << m.additionalUv << " extra UV sets)\n";

    int skin[5] = {};
    for (const auto& v : m.vertices) ++skin[std::clamp(v.weightType, 0, 4)];
    std::cout << "[PMX]   skinning: BDEF1 " << skin[0] << ", BDEF2 " << skin[1] << ", BDEF4 " << skin[2] << ", SDEF " << skin[3] << ", QDEF " << skin[4]
        << "   (bones/morphs are in the file but not read: this is a static render)\n";

    int existing = 0;
    for (const auto& t : m.textures) {
        std::string rel = t;
        std::replace(rel.begin(), rel.end(), '\\', '/');
        std::error_code ec;
        if (fs::exists(file.parent_path() / PathFromUtf8(rel), ec)) ++existing;
    }
    std::cout << "[PMX]   texture files found next to the model: " << existing << " of " << m.textures.size() << "\n";

    for (size_t i = 0; i < m.materials.size(); ++i) {
        const auto& mat = m.materials[i];
        std::cout << "[PMX]   material " << i << ": " << mat.name << "  tris " << mat.triangleCount()
            << "  diffuse (" << mat.diffuse.x << ", " << mat.diffuse.y << ", " << mat.diffuse.z << ") alpha " << mat.diffuse.w
            << (mat.doubleSided() ? "  double-sided" : "")
            << (mat.texture >= 0 ? "  tex " + m.textures[static_cast<size_t>(mat.texture)] : std::string("  no texture"))
            << "  spec (" << mat.specular.x << "," << mat.specular.y << "," << mat.specular.z << ") shin " << mat.shininess
            << (mat.sphereTexture >= 0 ? "  sphere" + std::to_string(mat.sphereMode) + " " + m.textures[static_cast<size_t>(mat.sphereTexture)] : std::string())
            << (mat.sharedToon ? "  toon#" + std::to_string(mat.toon) : std::string("  own toon"))
            << "\n";
    }
}

std::shared_ptr<Material> MakePmxMaterial(const io::PmxMaterial& pm, const RgbColor& base,
    const std::shared_ptr<Texture2D>& tex, bool matte, const std::shared_ptr<const SubsurfaceProfile>& skinProfile)
{
    if (matte) return std::make_shared<Lambertian>(base, tex);

    const Real specLum = Real(0.2126) * pm.specular.x + Real(0.7152) * pm.specular.y + Real(0.0722) * pm.specular.z;
    if (specLum > Real(0.01)) {
        const Real alpha = std::sqrt(Real(2) / (std::max<Real>(pm.shininess, Real(1)) + Real(2)));
        return std::make_shared<CoatedDiffuse>(base, std::sqrt(alpha), std::min(Real(1), Real(2) * specLum), Real(0), tex);
    }

    const auto has = [&](const char* s) { return pm.name.find(s) != std::string::npos; };
    if (has("髮") || has("髪")) return std::make_shared<CoatedDiffuse>(base, Real(0.30), Real(1.0), Real(0), tex);   // hair: soft shine
    if (has("目")) return std::make_shared<CoatedDiffuse>(base, Real(0.08), Real(1.0), Real(0), tex);                 // eyes: wet cornea
    if (pm.name == "饰") return std::make_shared<CoatedDiffuse>(base, Real(0.30), Real(1.0), Real(0.85), tex);        // gold ornaments
    if (has("体") || has("肌") || has("颜") || has("首")) {   // skin: faint sheen, and subsurface scattering with --sss (one profile for all skin)
        auto skin = std::make_shared<CoatedDiffuse>(base, Real(0.50), Real(0.4), Real(0), tex);
        if (skinProfile) skin->setSubsurface(skinProfile);
        return skin;
    }
    return std::make_shared<Lambertian>(base, tex);
}

std::shared_ptr<TriangleMesh> PmxToMesh(const io::PmxModel& m, Real yawDegrees, Real* outHeight, Vector3* outSize) {
    constexpr Real kScale = 0.08;
    auto mesh = std::make_shared<TriangleMesh>();
    const Real yaw = glm::radians(yawDegrees), cy = std::cos(yaw), sy = std::sin(yaw);
    auto rotate = [&](const Vector3& v) { return Vector3(cy * v.x + sy * v.z, v.y, -sy * v.x + cy * v.z); };

    Vector3 lo(1e30), hi(-1e30);
    for (const auto& v : m.vertices) {
        const Vector3 p = rotate(Vector3(v.position.x, v.position.y, -v.position.z) * kScale);
        mesh->p.push_back(Point3(p));
        mesh->n.push_back(glm::normalize(rotate(Vector3(v.normal.x, v.normal.y, -v.normal.z))));
        mesh->uv.push_back(v.uv);
        lo = glm::min(lo, p);
        hi = glm::max(hi, p);
    }
    const Vector3 shift(-(lo.x + hi.x) / 2, -lo.y, -(lo.z + hi.z) / 2);
    for (auto& p : mesh->p) p += shift;
    mesh->indices = m.indices;
    if (outHeight) *outHeight = hi.y - lo.y;
    if (outSize) *outSize = hi - lo;
    return mesh;
}

} // namespace rayt::cli
