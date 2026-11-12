// Reading scene files (JSON): see Cli/SceneFile.hpp. The materials are read in SceneFileMaterials.cpp.
#include "pch.h"
#include "Cli/SceneFile.hpp"
#include "Cli/SceneDescription.hpp"
#include "Cli/SceneFileReaders.hpp"

#include <fstream>
#include <sstream>

namespace rayt::cli {

using namespace scenefile;   // the readers of JSON values, and ParseMaterial

namespace {

// ---- the parts of a scene -----------------------------------------------------------------------------------------------------------

// The profile of a lathe: a list of [radius, height] pairs, or {"shell": {...}} (the wall of a bowl: two arcs of a sphere).
std::vector<Point2> ParseProfile(const json& v, const std::string& where) {
    std::vector<Point2> profile;
    if (v.is_array()) {
        for (size_t i = 0; i < v.size(); ++i) {
            const std::string w = where + "[" + std::to_string(i) + "]";
            if (!v[i].is_array() || v[i].size() != 2 || !v[i][0].is_number() || !v[i][1].is_number()) Fail(w, "must be [radius, height]");
            profile.push_back(Point2(v[i][0].get<double>(), v[i][1].get<double>()));
        }
    }
    else if (v.is_object() && v.contains("shell")) {
        CheckKeys(v, { "shell" }, where);
        const json& s = v.at("shell");
        const std::string w = where + ".shell";
        CheckKeys(s, { "radius", "wall", "endDegrees", "steps" }, w);
        Real R = 0.08, wall = 0.003, endDeg = 78;
        int steps = 40;
        Optional<double>(s, "radius", w, R, ReadNumber);
        Optional<double>(s, "wall", w, wall, ReadNumber);
        Optional<double>(s, "endDegrees", w, endDeg, ReadNumber);
        if (s.contains("steps")) steps = static_cast<int>(ReadNumber(s.at("steps"), w + ".steps"));
        if (!(R > 0) || !(wall > 0) || wall >= R) Fail(w, "needs 0 < wall < radius");
        if (steps < 2) Fail(w + ".steps", "must be at least 2");
        const Real tEnd = glm::radians(endDeg);
        for (int i = 0; i <= steps; ++i) {   // the outer surface from the bottom up to the rim ...
            const Real t = tEnd * Real(i) / Real(steps);
            profile.push_back(Point2(R * std::sin(t), R - R * std::cos(t)));
        }
        for (int i = steps; i >= 0; --i) {   // ... and the inner surface back down to the axis
            const Real t = tEnd * Real(i) / Real(steps);
            profile.push_back(Point2((R - wall) * std::sin(t), R - (R - wall) * std::cos(t)));
        }
    }
    else {
        Fail(where, "must be a list of [radius, height] pairs, or {\"shell\": {\"radius\", \"wall\", \"endDegrees\", \"steps\"}}");
    }
    if (profile.size() < 2) Fail(where, "needs at least two points");
    return profile;
}

ObjectSpec ParseObject(const json& o, const std::string& where) {
    if (!o.is_object() || !o.contains("shape")) Fail(where, "must be an object with a \"shape\" (sphere, quad, mesh or lathe)");
    const std::string shape = ReadString(o.at("shape"), where + ".shape");
    ObjectSpec s;
    if (shape == "splats") {   // a captured cloud of Gaussians (3DGS point_cloud.ply): no material, its color is the capture's (docs/CAPTURE.md)
        CheckKeys(o, { "shape", "file", "scale", "rotate", "translate", "alphaMin", "intensity", "linearize", "degree", "keep", "remove" }, where);
        s.shape = ObjectSpec::Shape::Splats;
        if (!o.contains("file")) Fail(where, "splats need a \"file\" (point_cloud.ply of 3D Gaussian Splatting)");
        s.file = ReadString(o.at("file"), where + ".file");
        Optional<double>(o, "scale", where, s.scale, ReadNumber);
        if (o.contains("translate")) s.translate = ReadVec3(o.at("translate"), where + ".translate");
        if (o.contains("rotate")) {
            const json& r = o.at("rotate");
            CheckKeys(r, { "axis", "degrees" }, where + ".rotate");
            if (r.contains("axis")) s.rotateAxis = ReadVec3(r.at("axis"), where + ".rotate.axis");
            if (r.contains("degrees")) s.rotateDegrees = ReadNumber(r.at("degrees"), where + ".rotate.degrees");
            if (glm::length(s.rotateAxis) < 1e-12) Fail(where + ".rotate.axis", "must not be zero");
        }
        Optional<double>(o, "alphaMin", where, s.splatAlphaMin, ReadNumber);
        Optional<double>(o, "intensity", where, s.splatIntensity, ReadNumber);
        if (o.contains("linearize")) s.splatLinearize = ReadBool(o.at("linearize"), where + ".linearize");
        if (o.contains("degree")) s.splatDegree = static_cast<int>(ReadNumber(o.at("degree"), where + ".degree"));
        if (s.splatAlphaMin <= 0 || s.splatAlphaMin >= 1) Fail(where + ".alphaMin", "must be between 0 and 1 (default 1/255)");
        if (s.splatDegree < 0 || s.splatDegree > 3) Fail(where + ".degree", "must be 0 to 3");
        for (const char* key : { "keep", "remove" }) {
            if (!o.contains(key)) continue;
            const json& b = o.at(key);
            const std::string bw = where + "." + key;
            CheckKeys(b, { "min", "max" }, bw);
            if (!b.contains("min") || !b.contains("max")) Fail(bw, "needs \"min\" and \"max\" (a box in world space)");
            const Vector3 lo = ReadVec3(b.at("min"), bw + ".min"), hi = ReadVec3(b.at("max"), bw + ".max");
            if (std::string(key) == "keep") { s.splatKeep = true; s.splatKeepMin = lo; s.splatKeepMax = hi; }
            else { s.splatRemove = true; s.splatRemoveMin = lo; s.splatRemoveMax = hi; }
        }
        return s;
    }
    if (!o.contains("material")) Fail(where, "needs a \"material\" (the name of one of the materials)");
    s.material = ReadString(o.at("material"), where + ".material");
    if (o.contains("light")) s.light = ReadBool(o.at("light"), where + ".light") ? 1 : 0;
    if (shape == "sphere") {
        CheckKeys(o, { "shape", "material", "light", "center", "radius" }, where);
        s.shape = ObjectSpec::Shape::Sphere;
        if (!o.contains("center") || !o.contains("radius")) Fail(where, "a sphere needs \"center\" and \"radius\"");
        s.center = ReadVec3(o.at("center"), where + ".center");
        s.radius = ReadNumber(o.at("radius"), where + ".radius");
        if (!(s.radius > 0)) Fail(where + ".radius", "must be above 0");
    }
    else if (shape == "quad") {
        CheckKeys(o, { "shape", "material", "light", "corner", "edge1", "edge2" }, where);
        s.shape = ObjectSpec::Shape::Quad;
        if (!o.contains("corner") || !o.contains("edge1") || !o.contains("edge2")) Fail(where, "a quad needs \"corner\", \"edge1\" and \"edge2\"");
        s.corner = ReadVec3(o.at("corner"), where + ".corner");
        s.edge1 = ReadVec3(o.at("edge1"), where + ".edge1");
        s.edge2 = ReadVec3(o.at("edge2"), where + ".edge2");
    }
    else if (shape == "mesh") {
        CheckKeys(o, { "shape", "material", "light", "file", "cutout", "scale", "rotate", "translate" }, where);
        s.shape = ObjectSpec::Shape::Mesh;
        if (!o.contains("file")) Fail(where, "a mesh needs a \"file\" (an OBJ or a PLY)");
        s.file = ReadString(o.at("file"), where + ".file");
        Optional<std::string>(o, "cutout", where, s.cutout, ReadString);   // an RGBA image over the uv: alpha < 1/2 = a hole (the glass under a water drop)
        Optional<double>(o, "scale", where, s.scale, ReadNumber);
        if (o.contains("translate")) s.translate = ReadVec3(o.at("translate"), where + ".translate");
        if (o.contains("rotate")) {
            const json& r = o.at("rotate");
            CheckKeys(r, { "axis", "degrees" }, where + ".rotate");
            if (r.contains("axis")) s.rotateAxis = ReadVec3(r.at("axis"), where + ".rotate.axis");
            if (r.contains("degrees")) s.rotateDegrees = ReadNumber(r.at("degrees"), where + ".rotate.degrees");
            if (glm::length(s.rotateAxis) < 1e-12) Fail(where + ".rotate.axis", "must not be zero");
        }
    }
    else if (shape == "lathe") {
        CheckKeys(o, { "shape", "material", "light", "profile", "segments", "translate" }, where);
        s.shape = ObjectSpec::Shape::Lathe;
        if (!o.contains("profile")) Fail(where, "a lathe needs a \"profile\"");
        s.profile = ParseProfile(o.at("profile"), where + ".profile");
        if (o.contains("translate")) s.translate = ReadVec3(o.at("translate"), where + ".translate");
        if (o.contains("segments")) s.segments = static_cast<int>(ReadNumber(o.at("segments"), where + ".segments"));
        if (s.segments < 3) Fail(where + ".segments", "must be at least 3");
    }
    else if (shape == "ocean") {   // a sea surface from a wave spectrum
        CheckKeys(o, { "shape", "material", "light", "windSpeed", "windDirection", "spread", "peakEnhance", "chop", "time", "minWavelength",
            "seed", "rings", "sectors", "rInner", "rOuter" }, where);
        s.shape = ObjectSpec::Shape::Ocean;
        Optional<double>(o, "windSpeed", where, s.windSpeed, ReadNumber);
        Optional<double>(o, "windDirection", where, s.windDirection, ReadNumber);
        Optional<double>(o, "spread", where, s.spread, ReadNumber);
        Optional<double>(o, "peakEnhance", where, s.peakEnhance, ReadNumber);
        Optional<double>(o, "chop", where, s.chop, ReadNumber);
        Optional<double>(o, "time", where, s.time, ReadNumber);
        Optional<double>(o, "minWavelength", where, s.minWavelength, ReadNumber);
        Optional<double>(o, "rInner", where, s.rInner, ReadNumber);
        Optional<double>(o, "rOuter", where, s.rOuter, ReadNumber);
        double v = s.oceanSeed;
        Optional<double>(o, "seed", where, v, ReadNumber); s.oceanSeed = int(v);
        v = s.rings; Optional<double>(o, "rings", where, v, ReadNumber); s.rings = int(v);
        v = s.sectors; Optional<double>(o, "sectors", where, v, ReadNumber); s.sectors = int(v);
        if (!(s.windSpeed >= 1 && s.windSpeed <= 40)) Fail(where + ".windSpeed", "must be between 1 and 40 m/s");
        if (s.chop < 0 || s.chop > 1) Fail(where + ".chop", "must be between 0 and 1");
        if (!(s.rInner > 0 && s.rOuter > s.rInner * 4)) Fail(where, "needs 0 < rInner and rOuter > 4 rInner");
        if (s.rings < 8 || s.rings > 2000 || s.sectors < 16 || s.sectors > 4096) Fail(where, "rings must be 8-2000 and sectors 16-4096");
        if (!(s.minWavelength >= 0.005)) Fail(where + ".minWavelength", "must be at least 0.005 m");
    }
    else if (shape == "hair") {   // procedural strands: "style": "lock" | "head"
        CheckKeys(o, { "shape", "material", "light", "style", "strands", "segments", "length", "rootRadius", "tipRadius", "seed", "root", "width", "tipWidth", "waveAmp",
            "waveLength", "twist", "frizz", "flyaway", "clumps", "clumpStrength", "center", "radius", "droop", "outward", "sway", "swayFreq", "headFrizz", "maxPolar", "hairline" }, where);
        s.shape = ObjectSpec::Shape::Hair;
        HairSpec& h = s.hair;
        if (o.contains("style")) h.style = ReadString(o.at("style"), where + ".style");
        if (h.style != "lock" && h.style != "head") Fail(where + ".style", "must be \"lock\" or \"head\"");
        if (h.style == "head") { h.strands = 40000; h.length = 0.28; h.segments = 26; }
        double v = h.strands; Optional<double>(o, "strands", where, v, ReadNumber); h.strands = int(v);
        v = h.segments; Optional<double>(o, "segments", where, v, ReadNumber); h.segments = int(v);
        v = h.seed; Optional<double>(o, "seed", where, v, ReadNumber); h.seed = int(v);
        Optional<double>(o, "length", where, h.length, ReadNumber);
        Optional<double>(o, "rootRadius", where, h.rootRadius, ReadNumber);
        Optional<double>(o, "tipRadius", where, h.tipRadius, ReadNumber);
        if (o.contains("root")) h.root = ReadVec3(o.at("root"), where + ".root");
        Optional<double>(o, "width", where, h.width, ReadNumber);
        Optional<double>(o, "tipWidth", where, h.tipWidth, ReadNumber);
        Optional<double>(o, "waveAmp", where, h.waveAmp, ReadNumber);
        Optional<double>(o, "waveLength", where, h.waveLength, ReadNumber);
        Optional<double>(o, "twist", where, h.twist, ReadNumber);
        Optional<double>(o, "frizz", where, h.frizz, ReadNumber);
        Optional<double>(o, "flyaway", where, h.flyaway, ReadNumber);
        Optional<double>(o, "clumpStrength", where, h.clumpStrength, ReadNumber);
        v = h.clumps; Optional<double>(o, "clumps", where, v, ReadNumber); h.clumps = int(v);
        if (o.contains("center")) h.center = ReadVec3(o.at("center"), where + ".center");
        Optional<double>(o, "radius", where, h.headRadius, ReadNumber);
        Optional<double>(o, "droop", where, h.droop, ReadNumber);
        Optional<double>(o, "outward", where, h.outward, ReadNumber);
        Optional<double>(o, "sway", where, h.sway, ReadNumber);
        Optional<double>(o, "swayFreq", where, h.swayFreq, ReadNumber);
        Optional<double>(o, "headFrizz", where, h.headFrizz, ReadNumber);
        Optional<double>(o, "maxPolar", where, h.maxPolarDeg, ReadNumber);
        Optional<double>(o, "hairline", where, h.hairline, ReadNumber);
        if (h.strands < 1 || h.strands > 2000000) Fail(where + ".strands", "must be between 1 and 2000000");
        if (h.segments < 2 || h.segments > 200) Fail(where + ".segments", "must be between 2 and 200");
        if (!(h.length > 0)) Fail(where + ".length", "must be above 0");
        if (!(h.rootRadius > 0 && h.tipRadius > 0)) Fail(where, "needs rootRadius and tipRadius above 0 (metres; a real hair is about 0.00005)");
        if (h.style == "lock" && !(h.width > 0 && h.tipWidth >= 0 && h.waveLength > 0)) Fail(where, "a lock needs width > 0, tipWidth >= 0 and waveLength > 0");
        if (h.flyaway < 0 || h.flyaway > 1 || h.clumpStrength < 0 || h.clumpStrength > 1) Fail(where, "flyaway and clumpStrength must be between 0 and 1");
        if (h.clumps == 0 || h.clumps > 100000) Fail(where + ".clumps", "must be between 1 and 100000");
        if (h.style == "head" && !(h.headRadius > 0)) Fail(where + ".radius", "must be above 0");
    }
    else {
        Fail(where + ".shape", "is \"" + shape + "\", must be sphere, quad, mesh, lathe, ocean or hair");
    }
    return s;
}

// The scene of form 2. `args` gets "--env FILE" when the environment names a map.
std::shared_ptr<const SceneDescription> ParseDescription(const json& doc, std::vector<std::string>& args) {
    auto d = std::make_shared<SceneDescription>();

    if (doc.contains("camera")) {
        const json& c = doc.at("camera");
        CheckKeys(c, { "from", "at", "up", "fov", "aperture", "focus" }, "camera");
        if (c.contains("from")) d->camera.from = ReadVec3(c.at("from"), "camera.from");
        if (c.contains("at")) d->camera.at = ReadVec3(c.at("at"), "camera.at");
        if (c.contains("up")) d->camera.up = ReadVec3(c.at("up"), "camera.up");
        if (c.contains("fov")) d->camera.fovDeg = ReadNumber(c.at("fov"), "camera.fov");
        if (!(d->camera.fovDeg > 0 && d->camera.fovDeg < 180)) Fail("camera.fov", "must be between 0 and 180 degrees");
        if (c.contains("aperture")) d->camera.aperture = ReadNumber(c.at("aperture"), "camera.aperture");   // lens diameter in metres (depth of field)
        if (c.contains("focus")) d->camera.focus = ReadNumber(c.at("focus"), "camera.focus");             // focus distance in metres (default: to "at")
        if (d->camera.aperture < 0) Fail("camera.aperture", "must be 0 (pinhole) or the lens diameter in metres");
        if (d->camera.focus < 0) Fail("camera.focus", "must be 0 (the distance to \"at\") or a distance in metres");
    }

    if (doc.contains("environment")) {
        const json& e = doc.at("environment");
        std::string type = "none";
        if (e.is_string()) type = e.get<std::string>();
        else {
            CheckKeys(e, { "type", "file", "scale" }, "environment");
            if (e.contains("type")) type = ReadString(e.at("type"), "environment.type");
            if (e.contains("scale")) d->environment.scale = ReadNumber(e.at("scale"), "environment.scale");
            if (e.contains("file")) { args.push_back("--env"); args.push_back(ReadString(e.at("file"), "environment.file")); }
        }
        if (type == "none") d->environment.type = EnvironmentSpec::Type::None;
        else if (type == "hdr") d->environment.type = EnvironmentSpec::Type::Hdr;
        else if (type == "uniform") d->environment.type = EnvironmentSpec::Type::Uniform;
        else Fail("environment.type", "is \"" + type + "\", must be none, hdr or uniform");
    }

    if (doc.contains("materials")) {
        const json& ms = doc.at("materials");
        if (!ms.is_object()) Fail("materials", "must be an object: material name -> material");
        for (const auto& [name, m] : ms.items()) d->materials[name] = ParseMaterial(m, "materials." + name);
    }

    if (doc.contains("lights")) {
        const json& ls = doc.at("lights");
        if (!ls.is_array()) Fail("lights", "must be a list");
        for (size_t i = 0; i < ls.size(); ++i) {
            const std::string where = "lights[" + std::to_string(i) + "]";
            CheckKeys(ls[i], { "type", "position", "color", "intensity" }, where);
            if (!ls[i].contains("type") || ReadString(ls[i].at("type"), where + ".type") != "point") Fail(where + ".type", "must be \"point\" (the only light that is not an object)");
            if (!ls[i].contains("position")) Fail(where, "needs a \"position\"");
            PointLightSpec p;
            p.position = ReadVec3(ls[i].at("position"), where + ".position");
            if (ls[i].contains("color")) p.color = ReadColor(ls[i].at("color"), where + ".color");
            if (ls[i].contains("intensity")) p.intensity = ReadNumber(ls[i].at("intensity"), where + ".intensity");
            d->lights.push_back(p);
        }
    }

    if (!doc.contains("objects") || !doc.at("objects").is_array() || doc.at("objects").empty()) Fail("objects", "must be a non-empty list");
    const json& os = doc.at("objects");
    for (size_t i = 0; i < os.size(); ++i) {
        const std::string where = "objects[" + std::to_string(i) + "]";
        ObjectSpec o = ParseObject(os[i], where);
        if (o.shape == ObjectSpec::Shape::Splats) { d->objects.push_back(o); continue; }   // no material
        const auto it = d->materials.find(o.material);
        if (it == d->materials.end()) {
            std::string names;
            for (const auto& [n, m] : d->materials) names += std::string(names.empty() ? "" : ", ") + n;
            Fail(where + ".material", "names \"" + o.material + "\", which is not in \"materials\" (" + (names.empty() ? "none defined" : names) + ")");
        }
        const bool emission = it->second.type == MaterialSpec::Type::Emission;
        if (emission && o.shape != ObjectSpec::Shape::Sphere && o.shape != ObjectSpec::Shape::Quad)
            Fail(where, "an emission material can be used on a sphere or a quad only");
        if (!emission && o.light != -1) Fail(where + ".light", "only applies to an object with an emission material");
        d->objects.push_back(std::move(o));
    }
    return d;
}

// "options": every entry becomes "--name value" (true: "--name" alone, false: left out).
void OptionsToArgs(const json& options, std::vector<std::string>& args) {
    if (!options.is_object()) throw std::runtime_error("\"options\" must be an object");
    for (const auto& [key, value] : options.items()) {
        if (value.is_boolean()) {
            if (value.get<bool>()) args.push_back("--" + key);
        }
        else if (value.is_string()) {
            args.push_back("--" + key);
            args.push_back(value.get<std::string>());
        }
        else if (value.is_number()) {
            args.push_back("--" + key);
            args.push_back(value.dump());
        }
        else {
            throw std::runtime_error("option \"" + key + "\" must be a string, a number or true / false");
        }
    }
}

} // namespace

bool LoadSceneFile(const std::string& path, std::vector<std::string>& args, std::shared_ptr<const SceneDescription>& description,
                   std::string& error) {
    std::ifstream in(path, std::ios::binary);
    if (!in) { error = "cannot open scene file " + path; return false; }
    std::stringstream buffer;
    buffer << in.rdbuf();
    const std::string text = buffer.str();
    try {
        const json doc = json::parse(text, nullptr, /*allow_exceptions=*/true, /*ignore_comments=*/true);
        CheckKeys(doc, { "comment", "scene", "options", "camera", "environment", "materials", "lights", "objects" }, "scene file");
        const bool explicitScene = doc.contains("objects");
        if (explicitScene == doc.contains("scene"))
            throw std::runtime_error(explicitScene ? "give either \"scene\" (a built-in scene) or \"objects\" (a scene written as data), not both"
                                                   : "missing \"scene\" (the name of a built-in scene, see --list-scenes) or \"objects\"");
        if (!explicitScene) {
            for (const char* k : { "camera", "environment", "materials", "lights" })
                if (doc.contains(k)) throw std::runtime_error(std::string("\"") + k + "\" belongs to a scene written as data (with \"objects\")");
            if (!doc["scene"].is_string()) throw std::runtime_error("\"scene\" must be a string");
            args.push_back("--scene");
            args.push_back(doc["scene"].get<std::string>());
        }
        else {
            args.push_back("--scene");
            args.push_back("file");
            description = ParseDescription(doc, args);
        }
        if (doc.contains("options")) OptionsToArgs(doc["options"], args);
    }
    catch (const json::parse_error& e) {
        error = path + ": " + e.what();
        return false;
    }
    catch (const std::exception& e) {
        error = path + ": " + e.what();
        return false;
    }
    return true;
}

} // namespace rayt::cli
