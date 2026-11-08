// Reading scene files (JSON): see Cli/SceneFile.hpp.
#include "pch.h"
#include "Cli/SceneFile.hpp"
#include "Cli/SceneDescription.hpp"

#include <fstream>
#include <initializer_list>
#include <sstream>

#include <nlohmann/json.hpp>
#include "Color/Sellmeier.hpp"

namespace rayt::cli {

namespace {

using json = nlohmann::json;

// ---- small readers: every one names the place it is reading in its error message --------------------------------------------------

[[noreturn]] void Fail(const std::string& where, const std::string& what) { throw std::runtime_error(where + " " + what); }

// Throws unless every key of `obj` is one of `allowed`; the message lists the allowed ones.
void CheckKeys(const json& obj, std::initializer_list<const char*> allowed, const std::string& where) {
    if (!obj.is_object()) Fail(where, "must be an object");
    for (const auto& [key, value] : obj.items()) {
        bool known = false;
        for (const char* k : allowed) known = known || key == k;
        if (!known) {
            std::string list;
            for (const char* k : allowed) list += std::string(list.empty() ? "" : ", ") + k;
            throw std::runtime_error(where + ": unknown key \"" + key + "\" (allowed: " + list + ")");
        }
    }
}

double ReadNumber(const json& v, const std::string& where) {
    if (!v.is_number()) Fail(where, "must be a number");
    return v.get<double>();
}

bool ReadBool(const json& v, const std::string& where) {
    if (!v.is_boolean()) Fail(where, "must be true or false");
    return v.get<bool>();
}

std::string ReadString(const json& v, const std::string& where) {
    if (!v.is_string()) Fail(where, "must be a string");
    return v.get<std::string>();
}

Vector3 ReadVec3(const json& v, const std::string& where) {
    if (!v.is_array() || v.size() != 3 || !v[0].is_number() || !v[1].is_number() || !v[2].is_number())
        Fail(where, "must be an array of three numbers, like [0, 1, 0]");
    return Vector3(v[0].get<double>(), v[1].get<double>(), v[2].get<double>());
}

// An RGB color: three numbers, or one number for a gray.
Vector3 ReadColor(const json& v, const std::string& where) {
    if (v.is_number()) return Vector3(v.get<double>());
    return ReadVec3(v, where);
}

template <class T>
void Optional(const json& obj, const char* key, const std::string& where, T& target, T (*read)(const json&, const std::string&)) {
    if (obj.contains(key)) target = read(obj.at(key), where + "." + key);
}

// ---- the parts of a scene -----------------------------------------------------------------------------------------------------------

SpectrumSpec ParseSpectrum(const json& v, const std::string& where) {
    SpectrumSpec s;
    if (v.is_string()) {
        s.kind = SpectrumSpec::Kind::Preset;
        s.preset = v.get<std::string>();
    }
    else if (v.is_object() && v.size() == 1) {
        CheckKeys(v, { "preset", "blackbody", "rgb", "file" }, where);
        if (v.contains("preset")) { s.kind = SpectrumSpec::Kind::Preset; s.preset = ReadString(v.at("preset"), where + ".preset"); }
        else if (v.contains("blackbody")) { s.kind = SpectrumSpec::Kind::Blackbody; s.kelvin = ReadNumber(v.at("blackbody"), where + ".blackbody"); }
        else if (v.contains("rgb")) { s.kind = SpectrumSpec::Kind::Rgb; s.rgb = ReadColor(v.at("rgb"), where + ".rgb"); }
        else { s.kind = SpectrumSpec::Kind::Curve; s.file = ReadString(v.at("file"), where + ".file"); }
    }
    else {
        Fail(where, "must be a preset name (d65, tungsten, candle, led) or an object with one of: preset, blackbody (kelvin), rgb, file (csv)");
    }
    if (s.kind == SpectrumSpec::Kind::Preset && s.preset != "d65" && s.preset != "tungsten" && s.preset != "candle" && s.preset != "led")
        Fail(where, "has the unknown preset \"" + s.preset + "\" (d65, tungsten, candle or led)");
    if (s.kind == SpectrumSpec::Kind::Blackbody && !(s.kelvin > 0)) Fail(where + ".blackbody", "must be a temperature in kelvin above 0");
    return s;
}

MaterialSpec ParseMaterial(const json& m, const std::string& where) {
    if (!m.is_object() || !m.contains("type")) Fail(where, "must be an object with a \"type\" (diffuse, conductor, layered, dielectric, thinFilm or emission)");
    const std::string type = ReadString(m.at("type"), where + ".type");
    MaterialSpec s;
    if (type == "diffuse") {
        CheckKeys(m, { "type", "albedo", "reflectance" }, where);
        s.type = MaterialSpec::Type::Diffuse;
        if (m.contains("albedo")) s.albedo = ReadColor(m.at("albedo"), where + ".albedo");
        Optional<std::string>(m, "reflectance", where, s.reflectance, ReadString);
        if (m.contains("albedo") && !s.reflectance.empty()) Fail(where, "takes either \"albedo\" (RGB) or \"reflectance\" (a CSV of R(wavelength)), not both");
    }
    else if (type == "conductor") {
        CheckKeys(m, { "type", "metal", "table", "reflectance", "roughness", "anisotropy", "rotation", "multipleScattering" }, where);
        s.type = MaterialSpec::Type::Conductor;
        Optional<std::string>(m, "metal", where, s.metal, ReadString);
        Optional<std::string>(m, "table", where, s.table, ReadString);
        Optional<std::string>(m, "reflectance", where, s.reflectance, ReadString);
        if (int(!s.metal.empty()) + int(!s.table.empty()) + int(!s.reflectance.empty()) != 1)
            Fail(where, "needs exactly one of \"metal\" (Au, Ag or Cu), \"table\" (a CSV of n and k) or \"reflectance\" (a CSV of R(wavelength))");
        if (!s.metal.empty() && s.metal != "Au" && s.metal != "Ag" && s.metal != "Cu") Fail(where + ".metal", "must be Au, Ag or Cu");
        Optional<double>(m, "roughness", where, s.roughness, ReadNumber);
        Optional<double>(m, "anisotropy", where, s.anisotropy, ReadNumber);
        Optional<double>(m, "rotation", where, s.rotation, ReadNumber);
        Optional<bool>(m, "multipleScattering", where, s.multipleScattering, ReadBool);
        if (s.roughness < 0 || s.roughness > 1) Fail(where + ".roughness", "must be between 0 and 1");
        if (s.anisotropy < -1 || s.anisotropy > 1) Fail(where + ".anisotropy", "must be between -1 and 1");
    }
    else if (type == "layered") {   // "coat": { ior, roughness, tint, samples }, "base": { albedo | metal | table | reflectance, roughness, multipleScattering }
        CheckKeys(m, { "type", "coat", "base" }, where);
        s.type = MaterialSpec::Type::Layered;
        if (m.contains("coat")) {
            const json& c = m.at("coat");
            const std::string cw = where + ".coat";
            if (!c.is_object()) Fail(cw, "must be an object (ior, roughness, tint, samples)");
            CheckKeys(c, { "ior", "roughness", "tint", "samples" }, cw);
            Optional<double>(c, "ior", cw, s.coatIor, ReadNumber);
            Optional<double>(c, "roughness", cw, s.coatRoughness, ReadNumber);
            if (c.contains("tint")) s.coatTint = ReadColor(c.at("tint"), cw + ".tint");
            double smp = s.coatSamples;
            Optional<double>(c, "samples", cw, smp, ReadNumber);
            s.coatSamples = int(smp);
            if (!(s.coatIor >= 1.0 && s.coatIor <= 3.0)) Fail(cw + ".ior", "must be between 1 and 3");
            if (s.coatRoughness < 0 || s.coatRoughness > 1) Fail(cw + ".roughness", "must be between 0 and 1");
            for (int k = 0; k < 3; ++k) if (s.coatTint[k] < 0 || s.coatTint[k] > 1) Fail(cw + ".tint", "must be between 0 and 1");
            if (s.coatSamples < 1 || s.coatSamples > 64) Fail(cw + ".samples", "must be between 1 and 64");
        }
        if (!m.contains("base")) Fail(where, "needs a \"base\" (albedo, or metal / table / reflectance)");
        const json& b = m.at("base");
        const std::string bw = where + ".base";
        if (!b.is_object()) Fail(bw, "must be an object (albedo, or metal / table / reflectance, with roughness)");
        CheckKeys(b, { "albedo", "metal", "table", "reflectance", "roughness", "anisotropy", "rotation", "multipleScattering" }, bw);
        if (b.contains("albedo")) s.albedo = ReadColor(b.at("albedo"), bw + ".albedo");
        Optional<std::string>(b, "metal", bw, s.metal, ReadString);
        Optional<std::string>(b, "table", bw, s.table, ReadString);
        Optional<std::string>(b, "reflectance", bw, s.reflectance, ReadString);
        const int kinds = int(b.contains("albedo")) + int(!s.metal.empty()) + int(!s.table.empty()) + int(!s.reflectance.empty());
        if (kinds != 1) Fail(bw, "needs exactly one of \"albedo\" (diffuse), \"metal\" (Au, Ag or Cu), \"table\" or \"reflectance\"");
        s.baseIsConductor = !b.contains("albedo");
        if (!s.metal.empty() && s.metal != "Au" && s.metal != "Ag" && s.metal != "Cu") Fail(bw + ".metal", "must be Au, Ag or Cu");
        Optional<double>(b, "roughness", bw, s.roughness, ReadNumber);
        Optional<double>(b, "anisotropy", bw, s.anisotropy, ReadNumber);
        Optional<double>(b, "rotation", bw, s.rotation, ReadNumber);
        Optional<bool>(b, "multipleScattering", bw, s.multipleScattering, ReadBool);
        if (s.roughness < 0 || s.roughness > 1) Fail(bw + ".roughness", "must be between 0 and 1");
        if (s.anisotropy < -1 || s.anisotropy > 1) Fail(bw + ".anisotropy", "must be between -1 and 1");
    }
    else if (type == "dielectric") {
        CheckKeys(m, { "type", "ior", "roughness", "abbe", "sellmeier", "medium" }, where);
        s.type = MaterialSpec::Type::Dielectric;
        if (m.contains("medium")) {
            const json& md = m.at("medium");
            const std::string mw = where + ".medium";
            if (!md.is_object()) Fail(mw, "must be an object with \"sigmaT\" (1/m), \"albedo\" and \"g\"");
            CheckKeys(md, { "sigmaT", "albedo", "g", "g2", "weight", "phaseTable", "noise" }, mw);
            s.hasMedium = true;
            if (md.contains("sigmaT")) s.mediumSigmaT = ReadColor(md.at("sigmaT"), mw + ".sigmaT");
            if (md.contains("albedo")) s.mediumAlbedo = ReadColor(md.at("albedo"), mw + ".albedo");
            Optional<double>(md, "g", mw, s.mediumG, ReadNumber);
            Optional<double>(md, "g2", mw, s.mediumG2, ReadNumber);
            Optional<double>(md, "weight", mw, s.mediumWeight, ReadNumber);
            Optional<std::string>(md, "phaseTable", mw, s.mediumPhaseTable, ReadString);
            if (!s.mediumPhaseTable.empty() && (md.contains("g") || md.contains("g2") || md.contains("weight")))
                Fail(mw, "\"phaseTable\" replaces \"g\", \"g2\" and \"weight\": give one or the other");
            if (md.contains("noise")) {   // "noise": { "frequency", "octaves", "lo", "hi", "seed" } (all optional)
                const json& nz = md.at("noise");
                const std::string nw = mw + ".noise";
                if (!nz.is_object()) Fail(nw, "must be an object (frequency, octaves, lo, hi, seed)");
                CheckKeys(nz, { "frequency", "octaves", "lo", "hi", "seed", "maskCenter", "maskRadius", "maskStart" }, nw);
                s.mediumNoise = true;
                Optional<double>(nz, "frequency", nw, s.noiseFrequency, ReadNumber);
                Optional<double>(nz, "lo", nw, s.noiseLo, ReadNumber);
                Optional<double>(nz, "hi", nw, s.noiseHi, ReadNumber);
                Optional<double>(nz, "maskRadius", nw, s.maskRadius, ReadNumber);
                Optional<double>(nz, "maskStart", nw, s.maskStart, ReadNumber);
                if (nz.contains("maskCenter")) s.maskCenter = ReadVec3(nz.at("maskCenter"), nw + ".maskCenter");
                double oct = s.noiseOctaves, seed = s.noiseSeed;
                Optional<double>(nz, "octaves", nw, oct, ReadNumber);
                Optional<double>(nz, "seed", nw, seed, ReadNumber);
                s.noiseOctaves = int(oct); s.noiseSeed = int(seed);
                if (!(s.noiseFrequency > 0)) Fail(nw + ".frequency", "must be above 0");
                if (s.noiseOctaves < 1 || s.noiseOctaves > 12) Fail(nw + ".octaves", "must be between 1 and 12");
                if (!(s.noiseHi > s.noiseLo)) Fail(nw, "needs hi above lo");
            }
            for (int c = 0; c < 3; ++c) {
                if (!(s.mediumSigmaT[c] > 0)) Fail(mw + ".sigmaT", "must be above 0");
                if (s.mediumAlbedo[c] < 0 || s.mediumAlbedo[c] > 1) Fail(mw + ".albedo", "must be between 0 and 1");
            }
            if (s.mediumG <= -1 || s.mediumG >= 1) Fail(mw + ".g", "must be between -1 and 1");
            if (s.mediumG2 <= -1 || s.mediumG2 >= 1) Fail(mw + ".g2", "must be between -1 and 1");
            if (s.mediumWeight < 0 || s.mediumWeight > 1) Fail(mw + ".weight", "must be between 0 and 1");
        }
        Optional<double>(m, "ior", where, s.ior, ReadNumber);
        Optional<double>(m, "roughness", where, s.roughness, ReadNumber);
        Optional<double>(m, "abbe", where, s.abbe, ReadNumber);
        if (m.contains("sellmeier")) {   // "sellmeier": "diamond" | { "b": [b1, b2, b3], "c": [c1, c2, c3] }  (c in um^2)
            const json& sm = m.at("sellmeier");
            const std::string sw = where + ".sellmeier";
            if (s.abbe != 0) Fail(sw, "and \"abbe\" both give the dispersion: use one of them");
            if (sm.is_string()) {
                Sellmeier test;
                s.sellmeierPreset = sm.get<std::string>();
                if (!sellmeierPreset(s.sellmeierPreset, test)) Fail(sw, "is the unknown glass \"" + s.sellmeierPreset + "\" (diamond, sapphire, bk7 or silica)");
            }
            else if (sm.is_object()) {
                CheckKeys(sm, { "b", "c" }, sw);
                if (!sm.contains("b") || !sm.contains("c") || !sm.at("b").is_array() || !sm.at("c").is_array() || sm.at("b").size() < 1 || sm.at("b").size() > 3 || sm.at("b").size() != sm.at("c").size())
                    Fail(sw, "needs \"b\" and \"c\": lists of 1 to 3 numbers each, the same length");
                s.hasSellmeier = true;
                for (size_t i = 0; i < sm.at("b").size(); ++i) { s.sellmeierB[i] = ReadNumber(sm.at("b")[i], sw + ".b"); s.sellmeierC[i] = ReadNumber(sm.at("c")[i], sw + ".c"); }
            }
            else Fail(sw, "must be a glass name (diamond, sapphire, bk7, silica) or { \"b\": [...], \"c\": [...] }");
        }
        if (!(s.ior > 0)) Fail(where + ".ior", "must be above 0");
        if (s.roughness < 0 || s.roughness > 1) Fail(where + ".roughness", "must be between 0 and 1");
    }
    else if (type == "thinFilm") {
        // { "type": "thinFilm", "ior": 1.33, "thickness": 450 | { "nm": 450 } | { "from": 150, "to": 700, "axis": [0,1,0], "range": [y0, y1],
        //   "noise": { "amplitude", "frequency", "octaves", "seed" } },  "substrate": "free" | { "dielectric": 1.33 } | { "metal": "Au" } | { "table": "file.csv" } }
        CheckKeys(m, { "type", "ior", "thickness", "substrate", "coverage" }, where);
        s.type = MaterialSpec::Type::ThinFilm;
        Optional<double>(m, "ior", where, s.filmIor, ReadNumber);
        if (!(s.filmIor >= 1)) Fail(where + ".ior", "must be at least 1");
        if (m.contains("thickness")) {
            const json& t = m.at("thickness");
            const std::string tw = where + ".thickness";
            if (t.is_number()) s.filmThickness = t.get<double>();
            else if (t.is_object()) {
                CheckKeys(t, { "nm", "from", "to", "axis", "range", "exponent", "warp", "noise" }, tw);
                Optional<double>(t, "nm", tw, s.filmThickness, ReadNumber);
                if (t.contains("from") || t.contains("to")) {
                    s.filmGradient = true;
                    if (t.contains("from")) s.filmThickness = ReadNumber(t.at("from"), tw + ".from");
                    s.filmTop = t.contains("to") ? ReadNumber(t.at("to"), tw + ".to") : s.filmThickness;
                    if (t.contains("axis")) s.filmAxis = glm::normalize(ReadVec3(t.at("axis"), tw + ".axis"));
                    if (t.contains("range")) {
                        const json& r = t.at("range");
                        if (!r.is_array() || r.size() != 2 || !r[0].is_number() || !r[1].is_number()) Fail(tw + ".range", "must be [start, end] along the axis (metres)");
                        s.filmRange0 = r[0].get<double>(); s.filmRange1 = r[1].get<double>();
                        if (s.filmRange0 == s.filmRange1) Fail(tw + ".range", "needs two different values");
                    }
                }
                Optional<double>(t, "exponent", tw, s.filmExponent, ReadNumber);
                if (!(s.filmExponent > 0)) Fail(tw + ".exponent", "must be above 0");
                if (t.contains("warp")) {   // "warp": { "amplitude" (metres), "frequency", "seed" }
                    const json& wp = t.at("warp");
                    const std::string ww = tw + ".warp";
                    if (!wp.is_object()) Fail(ww, "must be an object (amplitude, frequency, seed)");
                    CheckKeys(wp, { "amplitude", "frequency", "seed" }, ww);
                    Optional<double>(wp, "amplitude", ww, s.filmWarpAmp, ReadNumber);
                    Optional<double>(wp, "frequency", ww, s.filmWarpFreq, ReadNumber);
                    double wseed = s.filmWarpSeed;
                    Optional<double>(wp, "seed", ww, wseed, ReadNumber);
                    s.filmWarpSeed = int(wseed);
                    if (!(s.filmWarpFreq > 0)) Fail(ww + ".frequency", "must be above 0");
                }
                if (t.contains("noise")) {
                    const json& nz = t.at("noise");
                    const std::string nw = tw + ".noise";
                    if (!nz.is_object()) Fail(nw, "must be an object (amplitude, frequency, octaves, seed)");
                    CheckKeys(nz, { "amplitude", "frequency", "octaves", "seed", "stretch" }, nw);
                    if (nz.contains("stretch")) s.filmNoiseStretch = ReadVec3(nz.at("stretch"), nw + ".stretch");
                    Optional<double>(nz, "amplitude", nw, s.filmNoiseAmp, ReadNumber);
                    Optional<double>(nz, "frequency", nw, s.filmNoiseFreq, ReadNumber);
                    double oct = s.filmNoiseOctaves, seed = s.filmNoiseSeed;
                    Optional<double>(nz, "octaves", nw, oct, ReadNumber);
                    Optional<double>(nz, "seed", nw, seed, ReadNumber);
                    s.filmNoiseOctaves = int(oct); s.filmNoiseSeed = int(seed);
                    if (!(s.filmNoiseFreq > 0)) Fail(nw + ".frequency", "must be above 0");
                }
            }
            else Fail(tw, "must be a thickness in nm or an object (nm, from / to / axis / range, noise)");
        }
        if (m.contains("coverage")) {   // "coverage": { "frequency", "octaves", "lo", "hi", "seed" }: patches of film where the noise is above lo (sharp edge when hi is close to lo)
            const json& cv = m.at("coverage");
            const std::string cw = where + ".coverage";
            if (!cv.is_object()) Fail(cw, "must be an object (frequency, octaves, lo, hi, seed)");
            CheckKeys(cv, { "frequency", "octaves", "lo", "hi", "seed" }, cw);
            s.filmCoverage = true;
            Optional<double>(cv, "frequency", cw, s.covFrequency, ReadNumber);
            Optional<double>(cv, "lo", cw, s.covLo, ReadNumber);
            Optional<double>(cv, "hi", cw, s.covHi, ReadNumber);
            double oct = s.covOctaves, seed = s.covSeed;
            Optional<double>(cv, "octaves", cw, oct, ReadNumber);
            Optional<double>(cv, "seed", cw, seed, ReadNumber);
            s.covOctaves = int(oct); s.covSeed = int(seed);
            if (!(s.covFrequency > 0)) Fail(cw + ".frequency", "must be above 0");
            if (!(s.covHi > s.covLo)) Fail(cw, "needs hi above lo");
        }
        if (s.filmThickness < 0 || s.filmTop < 0) Fail(where + ".thickness", "must not be negative");
        if (m.contains("substrate")) {
            const json& sub = m.at("substrate");
            const std::string sw = where + ".substrate";
            if (sub.is_string()) {
                if (sub.get<std::string>() != "free") Fail(sw, "is a string only for \"free\"; use { \"dielectric\": ior } or { \"metal\": \"Au\" } otherwise");
            }
            else if (sub.is_object() && sub.size() == 1) {
                CheckKeys(sub, { "dielectric", "metal", "table" }, sw);
                if (sub.contains("dielectric")) { s.substrateKind = "dielectric"; s.substrateIor = ReadNumber(sub.at("dielectric"), sw + ".dielectric"); }
                else if (sub.contains("metal")) {
                    s.substrateKind = "conductor"; s.substrateMetal = ReadString(sub.at("metal"), sw + ".metal");
                    if (s.substrateMetal != "Au" && s.substrateMetal != "Ag" && s.substrateMetal != "Cu") Fail(sw + ".metal", "must be Au, Ag or Cu");
                }
                else { s.substrateKind = "conductor"; s.substrateTable = ReadString(sub.at("table"), sw + ".table"); }
                if (s.substrateKind == "dielectric" && !(s.substrateIor >= 1)) Fail(sw + ".dielectric", "must be at least 1");
            }
            else Fail(sw, "must be \"free\" or an object with one of: dielectric (index), metal (Au, Ag, Cu), table (a CSV of n and k)");
        }
    }
    else if (type == "emission") {
        CheckKeys(m, { "type", "spectrum", "radiance", "twoSided" }, where);
        s.type = MaterialSpec::Type::Emission;
        if (!m.contains("spectrum")) Fail(where, "needs a \"spectrum\"");
        s.spectrum = ParseSpectrum(m.at("spectrum"), where + ".spectrum");
        Optional<double>(m, "radiance", where, s.radiance, ReadNumber);
        Optional<bool>(m, "twoSided", where, s.twoSided, ReadBool);
    }
    else {
        Fail(where + ".type", "is \"" + type + "\", must be diffuse, conductor, layered, dielectric, thinFilm or emission");
    }
    return s;
}

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
        CheckKeys(o, { "shape", "material", "light", "file", "scale", "rotate", "translate" }, where);
        s.shape = ObjectSpec::Shape::Mesh;
        if (!o.contains("file")) Fail(where, "a mesh needs a \"file\" (an OBJ)");
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
    }
    else if (shape == "lathe") {
        CheckKeys(o, { "shape", "material", "light", "profile", "segments" }, where);
        s.shape = ObjectSpec::Shape::Lathe;
        if (!o.contains("profile")) Fail(where, "a lathe needs a \"profile\"");
        s.profile = ParseProfile(o.at("profile"), where + ".profile");
        if (o.contains("segments")) s.segments = static_cast<int>(ReadNumber(o.at("segments"), where + ".segments"));
        if (s.segments < 3) Fail(where + ".segments", "must be at least 3");
    }
    else {
        Fail(where + ".shape", "is \"" + shape + "\", must be sphere, quad, mesh or lathe");
    }
    return s;
}

// The scene of form 2. `args` gets "--env FILE" when the environment names a map.
std::shared_ptr<const SceneDescription> ParseDescription(const json& doc, std::vector<std::string>& args) {
    auto d = std::make_shared<SceneDescription>();

    if (doc.contains("camera")) {
        const json& c = doc.at("camera");
        CheckKeys(c, { "from", "at", "up", "fov" }, "camera");
        if (c.contains("from")) d->camera.from = ReadVec3(c.at("from"), "camera.from");
        if (c.contains("at")) d->camera.at = ReadVec3(c.at("at"), "camera.at");
        if (c.contains("up")) d->camera.up = ReadVec3(c.at("up"), "camera.up");
        if (c.contains("fov")) d->camera.fovDeg = ReadNumber(c.at("fov"), "camera.fov");
        if (!(d->camera.fovDeg > 0 && d->camera.fovDeg < 180)) Fail("camera.fov", "must be between 0 and 180 degrees");
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
