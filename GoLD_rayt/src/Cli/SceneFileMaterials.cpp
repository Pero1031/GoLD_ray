// Reading the materials of a scene file (JSON): one branch per material type. See Cli/SceneFile.hpp for the format.
#include "pch.h"
#include "Cli/SceneFileReaders.hpp"

#include "Color/Sellmeier.hpp"
#include "Materials/CosmeticLayer.hpp"
#include "Materials/Hair.hpp"
#include "Materials/SkinLotion.hpp"

namespace rayt::cli::scenefile {

namespace {

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

} // namespace

MaterialSpec ParseMaterial(const json& m, const std::string& where) {
    if (!m.is_object() || !m.contains("type")) Fail(where, "must be an object with a \"type\" (diffuse, conductor, layered, dielectric, thinFilm, hair, skin or emission)");
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
        CheckKeys(m, { "type", "ior", "roughness", "distanceRoughness", "abbe", "sellmeier", "medium" }, where);
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
        if (m.contains("distanceRoughness")) {   // "distanceRoughness": { "roughness", "start", "end" } (metres from the y axis)
            const json& dr = m.at("distanceRoughness");
            const std::string dw = where + ".distanceRoughness";
            if (!dr.is_object()) Fail(dw, "must be an object (roughness, start, end)");
            CheckKeys(dr, { "roughness", "start", "end" }, dw);
            s.farRoughness = 0.3;
            Optional<double>(dr, "roughness", dw, s.farRoughness, ReadNumber);
            Optional<double>(dr, "start", dw, s.farStart, ReadNumber);
            Optional<double>(dr, "end", dw, s.farEnd, ReadNumber);
            if (s.farRoughness < 0 || s.farRoughness > 1) Fail(dw + ".roughness", "must be between 0 and 1");
            if (!(s.farEnd > s.farStart)) Fail(dw, "needs end above start");
        }
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
    else if (type == "hair") {   // "melanin": [eumelanin, pheomelanin] | "color": [r, g, b] | "sigmaA": [r, g, b], "eta", "betaM", "betaN", "tilt"
        CheckKeys(m, { "type", "melanin", "color", "sigmaA", "eta", "betaM", "betaN", "tilt", "lobes" }, where);
        s.type = MaterialSpec::Type::Hair;
        Optional<double>(m, "eta", where, s.hairEta, ReadNumber);
        Optional<double>(m, "betaM", where, s.hairBetaM, ReadNumber);
        Optional<double>(m, "betaN", where, s.hairBetaN, ReadNumber);
        Optional<double>(m, "tilt", where, s.hairTiltDeg, ReadNumber);
        if (m.contains("lobes")) {   // "lobes": ["R", "TT", "TRT", "rest"]: the lobes that are on (to show what each one does)
            const json& l = m.at("lobes");
            if (!l.is_array()) Fail(where + ".lobes", "must be a list of lobe names (R, TT, TRT, rest)");
            s.hairLobes = 0;
            for (const json& e : l) {
                const std::string name = ReadString(e, where + ".lobes");
                if (name == "R") s.hairLobes |= 1u; else if (name == "TT") s.hairLobes |= 2u; else if (name == "TRT") s.hairLobes |= 4u; else if (name == "rest") s.hairLobes |= 8u;
                else Fail(where + ".lobes", "has the unknown lobe \"" + name + "\" (R, TT, TRT or rest)");
            }
        }
        if (!(s.hairEta >= 1.0 && s.hairEta <= 3.0)) Fail(where + ".eta", "must be between 1 and 3");
        if (s.hairBetaM < 0 || s.hairBetaM > 1) Fail(where + ".betaM", "must be between 0 and 1");
        if (s.hairBetaN < 0 || s.hairBetaN > 1) Fail(where + ".betaN", "must be between 0 and 1");
        if (int(m.contains("melanin")) + int(m.contains("color")) + int(m.contains("sigmaA")) > 1)
            Fail(where, "takes one of \"melanin\", \"color\" and \"sigmaA\"");
        if (m.contains("color")) {   // the apparent color of the fiber (d'Eon et al.): the absorption that gives it
            const Vector3 c = ReadColor(m.at("color"), where + ".color");
            for (int i = 0; i < 3; ++i) if (!(c[i] > 0 && c[i] <= 1)) Fail(where + ".color", "needs values in (0, 1]");
            s.hairSigmaA = Hair::sigmaFromColor(c, Real(s.hairBetaN));
        }
        else if (m.contains("sigmaA")) {
            s.hairSigmaA = ReadVec3(m.at("sigmaA"), where + ".sigmaA");
            for (int i = 0; i < 3; ++i) if (s.hairSigmaA[i] < 0) Fail(where + ".sigmaA", "must not be negative");
        }
        else {   // melanin: brown by default; 0 / 0 is white (gray) hair
            double eu = 1.3, ph = 0.0;
            if (m.contains("melanin")) {
                const json& mel = m.at("melanin");
                if (!mel.is_array() || mel.size() != 2 || !mel[0].is_number() || !mel[1].is_number()) Fail(where + ".melanin", "must be [eumelanin, pheomelanin], like [1.3, 0]");
                eu = mel[0].get<double>(); ph = mel[1].get<double>();
                if (eu < 0 || ph < 0) Fail(where + ".melanin", "must not be negative");
            }
            s.hairSigmaA = Hair::sigmaFromMelanin(Real(eu), Real(ph));
        }
    }
    else if (type == "skin") {   // "albedo", "roughness", "specular", "scatterDistanceMm": [r, g, b] (turns subsurface scattering on), "sssScale", "relief": { ... }
        CheckKeys(m, { "type", "albedo", "roughness", "specular", "scatterDistanceMm", "sssScale", "relief", "pigment", "gloss", "film", "scatter", "lotion", "foundation" }, where);
        s.type = MaterialSpec::Type::Skin;
        s.albedo = Vector3(0.80, 0.52, 0.42);
        s.roughness = 0.5;
        if (m.contains("albedo")) s.albedo = ReadColor(m.at("albedo"), where + ".albedo");
        Optional<double>(m, "roughness", where, s.roughness, ReadNumber);
        Optional<double>(m, "specular", where, s.skinSpecular, ReadNumber);
        Optional<double>(m, "sssScale", where, s.skinSssScale, ReadNumber);
        if (m.contains("scatterDistanceMm")) { s.skinSss = true; s.skinScatterMm = ReadColor(m.at("scatterDistanceMm"), where + ".scatterDistanceMm"); }
        if (s.roughness <= 0 || s.roughness > 1) Fail(where + ".roughness", "must be above 0 and at most 1");
        if (s.skinSpecular < 0 || s.skinSpecular > 1) Fail(where + ".specular", "must be between 0 and 1");
        if (s.skinSssScale <= 0) Fail(where + ".sssScale", "must be above 0");
        for (int i = 0; i < 3; ++i) if (s.skinSss && !(s.skinScatterMm[i] > 0)) Fail(where + ".scatterDistanceMm", "must be above 0 in every channel");
        if (m.contains("relief")) {   // "relief": { "cellMm", "furrowDepthMm", "furrowWidth", "poreDepthMm", "poreRadiusMm", "poreFraction", "seed" } (all optional)
            const json& r = m.at("relief");
            const std::string rw = where + ".relief";
            if (!r.is_object()) Fail(rw, "must be an object (cellMm, furrowDepthMm, furrowWidth, poreDepthMm, poreRadiusMm, poreFraction, seed)");
            CheckKeys(r, { "cellMm", "furrowDepthMm", "furrowWidth", "stretch", "coarseScale", "coarseWeight", "poreDepthMm", "poreRadiusMm", "poreFraction", "poreSizeVariation", "poreDensityVariation", "poreDensityFreq", "poreDarkness", "wrinkles", "seed" }, rw);
            s.skinRelief = true;
            Optional<double>(r, "cellMm", rw, s.reliefCellMm, ReadNumber);
            Optional<double>(r, "furrowDepthMm", rw, s.reliefFurrowDepthMm, ReadNumber);
            Optional<double>(r, "furrowWidth", rw, s.reliefFurrowWidth, ReadNumber);
            Optional<double>(r, "stretch", rw, s.reliefStretch, ReadNumber);
            Optional<double>(r, "coarseScale", rw, s.reliefCoarseScale, ReadNumber);
            Optional<double>(r, "coarseWeight", rw, s.reliefCoarseWeight, ReadNumber);
            if (r.contains("wrinkles")) {   // "wrinkles": { "tileMm", "resolution", "steps", "depthScale", "seed", "levels": [ { "spacingMm", "angleDeg", ... } ] }
                const json& wk = r.at("wrinkles");
                const std::string ww = rw + ".wrinkles";
                if (!wk.is_object()) Fail(ww, "must be an object (tileMm, resolution, steps, depthScale, seed, levels)");
                CheckKeys(wk, { "tileMm", "resolution", "steps", "depthScale", "seed", "levels" }, ww);
                s.reliefWrinkles = true;
                Optional<double>(wk, "tileMm", ww, s.wrTileMm, ReadNumber);
                Optional<double>(wk, "depthScale", ww, s.wrDepthScale, ReadNumber);
                double v1 = s.wrResolution; Optional<double>(wk, "resolution", ww, v1, ReadNumber); s.wrResolution = int(v1);
                v1 = s.wrSteps; Optional<double>(wk, "steps", ww, v1, ReadNumber); s.wrSteps = int(v1);
                v1 = s.wrSeed; Optional<double>(wk, "seed", ww, v1, ReadNumber); s.wrSeed = int(v1);
                if (!(s.wrTileMm > 1) || s.wrResolution < 64 || s.wrResolution > 8192 || s.wrSteps < 1 || s.wrSteps > 5000 || !(s.wrDepthScale > 0)) Fail(ww, "needs tileMm > 1, resolution 64-8192, steps 1-5000 and depthScale > 0");
                if (wk.contains("levels")) {
                    const json& lv = wk.at("levels");
                    if (!lv.is_array() || lv.empty() || lv.size() > 8) Fail(ww + ".levels", "must be a list of 1 to 8 levels");
                    for (size_t li = 0; li < lv.size(); ++li) {
                        const std::string lw = ww + ".levels[" + std::to_string(li) + "]";
                        const json& l = lv[li];
                        CheckKeys(l, { "spacingMm", "angleDeg", "base", "distance", "continuation", "similarity", "crossing", "deposit", "width", "poreWidth", "poreMin", "poreBlend", "poreSkew", "wobble", "wobbleFreq", "depthUm", "bulge", "darkPores" }, lw);
                        WrinkleLevelSpec ls;
                        Optional<double>(l, "spacingMm", lw, ls.spacingMm, ReadNumber);
                        Optional<double>(l, "angleDeg", lw, ls.angleDeg, ReadNumber);
                        Optional<double>(l, "base", lw, ls.base, ReadNumber);
                        Optional<double>(l, "distance", lw, ls.distance, ReadNumber);
                        Optional<double>(l, "continuation", lw, ls.continuation, ReadNumber);
                        Optional<double>(l, "similarity", lw, ls.similarity, ReadNumber);
                        Optional<double>(l, "crossing", lw, ls.crossing, ReadNumber);
                        Optional<double>(l, "deposit", lw, ls.deposit, ReadNumber);
                        Optional<double>(l, "width", lw, ls.width, ReadNumber);
                        Optional<double>(l, "poreWidth", lw, ls.poreWidth, ReadNumber);
                        Optional<double>(l, "poreMin", lw, ls.poreMin, ReadNumber);
                        Optional<double>(l, "poreBlend", lw, ls.poreBlend, ReadNumber);
                        Optional<double>(l, "poreSkew", lw, ls.poreSkew, ReadNumber);
                        Optional<double>(l, "wobble", lw, ls.wobble, ReadNumber);
                        Optional<double>(l, "wobbleFreq", lw, ls.wobbleFreq, ReadNumber);
                        Optional<double>(l, "depthUm", lw, ls.depthUm, ReadNumber);
                        Optional<bool>(l, "bulge", lw, ls.bulge, ReadBool);
                        Optional<bool>(l, "darkPores", lw, ls.darkPores, ReadBool);
                        if (!(ls.spacingMm >= 0.2 && ls.spacingMm <= 20)) Fail(lw + ".spacingMm", "must be between 0.2 and 20");
                        if (ls.deposit <= 0 || ls.deposit > 1 || ls.base < 0 || ls.width <= 0 || ls.depthUm < 0) Fail(lw, "needs 0 < deposit <= 1, base >= 0, width > 0 and depthUm >= 0");
                        s.wrLevels.push_back(ls);
                    }
                }
            }
            Optional<double>(r, "poreSizeVariation", rw, s.reliefPoreSizeVariation, ReadNumber);
            Optional<double>(r, "poreDensityVariation", rw, s.reliefPoreDensityVariation, ReadNumber);
            Optional<double>(r, "poreDensityFreq", rw, s.reliefPoreDensityFreq, ReadNumber);
            Optional<double>(r, "poreDarkness", rw, s.reliefPoreDarkness, ReadNumber);
            if (s.reliefPoreSizeVariation < 0 || s.reliefPoreSizeVariation >= 1 || s.reliefPoreDensityVariation < 0 || s.reliefPoreDensityVariation > 1 || s.reliefPoreDarkness < 0 || s.reliefPoreDarkness > 1)
                Fail(rw, "poreSizeVariation must be in [0, 1), poreDensityVariation and poreDarkness in [0, 1]");
            if (!(s.reliefPoreDensityFreq > 0)) Fail(rw + ".poreDensityFreq", "must be above 0");
            Optional<double>(r, "poreDepthMm", rw, s.reliefPoreDepthMm, ReadNumber);
            Optional<double>(r, "poreRadiusMm", rw, s.reliefPoreRadiusMm, ReadNumber);
            Optional<double>(r, "poreFraction", rw, s.reliefPoreFraction, ReadNumber);
            double sd = s.reliefSeed; Optional<double>(r, "seed", rw, sd, ReadNumber); s.reliefSeed = int(sd);
            if (!(s.reliefCellMm > 0 && s.reliefFurrowWidth > 0 && s.reliefPoreRadiusMm > 0 && s.reliefStretch > 0 && s.reliefCoarseScale >= 1 && s.reliefCoarseWeight >= 0)) Fail(rw, "needs cellMm, furrowWidth, poreRadiusMm and stretch above 0, coarseScale >= 1 and coarseWeight >= 0");
            if (s.reliefFurrowDepthMm < 0 || s.reliefPoreDepthMm < 0) Fail(rw, "depths must not be negative");
            if (s.reliefPoreFraction < 0 || s.reliefPoreFraction > 1) Fail(rw + ".poreFraction", "must be between 0 and 1");
        }
        if (m.contains("gloss")) {   // "gloss": { "variation": 0..0.95, "frequency": patches per metre, "seed" }: oily and dry patches
            const json& gl = m.at("gloss");
            const std::string glw = where + ".gloss";
            if (!gl.is_object()) Fail(glw, "must be an object (variation, frequency, seed)");
            CheckKeys(gl, { "variation", "frequency", "seed" }, glw);
            Optional<double>(gl, "variation", glw, s.skinGlossVariation, ReadNumber);
            Optional<double>(gl, "frequency", glw, s.skinGlossFreq, ReadNumber);
            double gs = s.skinGlossSeed; Optional<double>(gl, "seed", glw, gs, ReadNumber); s.skinGlossSeed = int(gs);
            if (s.skinGlossVariation < 0 || s.skinGlossVariation > 0.95 || !(s.skinGlossFreq > 0)) Fail(glw, "variation must be between 0 and 0.95 and frequency above 0");
        }
        if (m.contains("film")) {   // "film": { "thicknessUm": mean thickness, "ior", "roughness", "softnessUm" }: a liquid pooling in the furrows (needs "relief")
            const json& fm = m.at("film");
            const std::string fw = where + ".film";
            if (!fm.is_object()) Fail(fw, "must be an object (thicknessUm, ior, roughness, softnessUm)");
            CheckKeys(fm, { "thicknessUm", "ior", "roughness", "softnessUm" }, fw);
            if (!s.skinRelief) Fail(fw, "needs a \"relief\" (the film fills its furrows)");
            Optional<double>(fm, "thicknessUm", fw, s.skinFilmThicknessUm, ReadNumber);
            Optional<double>(fm, "ior", fw, s.skinFilmIor, ReadNumber);
            Optional<double>(fm, "roughness", fw, s.skinFilmRoughness, ReadNumber);
            Optional<double>(fm, "softnessUm", fw, s.skinFilmSoftnessUm, ReadNumber);
            if (!(s.skinFilmThicknessUm >= 0) || !(s.skinFilmIor >= 1 && s.skinFilmIor <= 2.5) || !(s.skinFilmRoughness >= 0 && s.skinFilmRoughness <= 1) || !(s.skinFilmSoftnessUm > 0))
                Fail(fw, "needs thicknessUm >= 0, ior between 1 and 2.5, roughness between 0 and 1 and softnessUm above 0");
            s.skinFilmOn = s.skinFilmThicknessUm > 0;   // 0 = no film
        }
        if (m.contains("scatter")) {   // "scatter": { "rho0", "frequency", "threshold", "softness", "hydration", "kappaH", "gapUm", "seed" }: the white scatter of flaking (dry) skin
            const json& sc = m.at("scatter");
            const std::string sw = where + ".scatter";
            if (!sc.is_object()) Fail(sw, "must be an object (rho0, frequency, threshold, softness, hydration, kappaH, gapUm, seed)");
            CheckKeys(sc, { "rho0", "frequency", "threshold", "softness", "hydration", "kappaH", "gapUm", "seed" }, sw);
            Optional<double>(sc, "rho0", sw, s.scatRho0, ReadNumber);
            Optional<double>(sc, "frequency", sw, s.scatFrequency, ReadNumber);
            Optional<double>(sc, "threshold", sw, s.scatThreshold, ReadNumber);
            Optional<double>(sc, "softness", sw, s.scatSoftness, ReadNumber);
            Optional<double>(sc, "hydration", sw, s.scatHydration, ReadNumber);
            Optional<double>(sc, "kappaH", sw, s.scatKappaH, ReadNumber);
            Optional<double>(sc, "gapUm", sw, s.scatGapUm, ReadNumber);
            double sd = s.scatSeed; Optional<double>(sc, "seed", sw, sd, ReadNumber); s.scatSeed = int(sd);
            if (!(s.scatRho0 >= 0 && s.scatRho0 <= 0.9) || !(s.scatFrequency > 0) || !(s.scatSoftness > 0) || !(s.scatHydration >= 0 && s.scatHydration <= 1) || !(s.scatKappaH >= 0 && s.scatKappaH <= 1) || !(s.scatGapUm > 0))
                Fail(sw, "needs rho0 between 0 and 0.9, frequency, softness and gapUm above 0, hydration and kappaH between 0 and 1");
            s.skinScatterOn = s.scatRho0 > 0;
        }
        if (m.contains("foundation")) {   // "foundation": { "preset", "thickness", "diffuserFraction", "diffuserAlbedo", "plateletAlbedo", "g1", "g2", "lobeWeight", "plateletRoughness", "plateletTiltDeg" }: a layer of foundation make-up on top
            const json& fd = m.at("foundation");
            const std::string fw = where + ".foundation";
            if (!fd.is_object()) Fail(fw, "must be an object (preset, thickness, diffuserFraction, diffuserAlbedo, plateletAlbedo, g1, g2, lobeWeight, plateletRoughness, plateletTiltDeg)");
            CheckKeys(fd, { "preset", "thickness", "diffuserFraction", "diffuserAlbedo", "plateletAlbedo", "g1", "g2", "lobeWeight", "plateletRoughness", "plateletTiltDeg", "leveling" }, fw);
            Optional<double>(fd, "leveling", fw, s.fndLeveling, ReadNumber);
            if (!(s.fndLeveling >= 0 && s.fndLeveling <= 1)) Fail(fw + ".leveling", "must be between 0 (the layer follows the skin relief) and 1 (it levels it)");
            CosmeticLayer::Params cp;
            if (fd.contains("preset")) {   // one of the four foundations fitted in the paper; the keys after it change single values
                const std::string preset = ReadString(fd.at("preset"), fw + ".preset");
                if (!CosmeticLayer::preset(preset, cp)) Fail(fw + ".preset", "is \"" + preset + "\", must be dewy1, dewy2, matte1 or matte2");
            }
            auto number = [&](const char* key, Real& target) { double v = double(target); Optional<double>(fd, key, fw, v, ReadNumber); target = Real(v); };
            number("thickness", cp.thickness);
            number("diffuserFraction", cp.diffuserFraction);
            number("g1", cp.g1);
            number("g2", cp.g2);
            number("lobeWeight", cp.lobeWeight);
            number("plateletRoughness", cp.plateletRoughness);
            number("plateletTiltDeg", cp.plateletTiltDeg);
            if (fd.contains("diffuserAlbedo")) cp.diffuserAlbedo = ReadColor(fd.at("diffuserAlbedo"), fw + ".diffuserAlbedo");
            if (fd.contains("plateletAlbedo")) cp.plateletAlbedo = ReadColor(fd.at("plateletAlbedo"), fw + ".plateletAlbedo");
            if (const char* problem = CosmeticLayer::validate(cp)) Fail(fw, std::string("is not usable: ") + problem);
            s.fndThickness = cp.thickness; s.fndDiffuserFraction = cp.diffuserFraction; s.fndG1 = cp.g1; s.fndG2 = cp.g2; s.fndLobeWeight = cp.lobeWeight;
            s.fndPlateletRoughness = cp.plateletRoughness; s.fndPlateletTiltDeg = cp.plateletTiltDeg;
            s.fndDiffuserAlbedo = cp.diffuserAlbedo; s.fndPlateletAlbedo = cp.plateletAlbedo;
            s.skinFoundationOn = cp.thickness > 0;   // 0 = no layer
        }
        if (m.contains("lotion")) {   // "lotion": { "doseUm", "glycerin", "timeS", "waterTauS", ... }: the film and the hydration `timeS` seconds after a lotion was applied (replaces "film")
            const json& lo = m.at("lotion");
            const std::string lw = where + ".lotion";
            if (!lo.is_object()) Fail(lw, "must be an object (doseUm, glycerin, timeS, waterTauS, hydrationPeak, hydrationRiseS, hydrationDecayS, waterIor, humectantIor, waterRoughness, residueRoughness)");
            CheckKeys(lo, { "doseUm", "glycerin", "timeS", "waterTauS", "hydrationPeak", "hydrationRiseS", "hydrationDecayS", "waterIor", "humectantIor", "waterRoughness", "residueRoughness", "bulk", "oil", "oilIor", "dropletUm", "coalesceTauS", "dropletG", "occlusion", "occlusionOilUm", "layers" }, lw);
            if (!s.skinRelief) Fail(lw, "needs a \"relief\" (the liquid pools in its furrows)");
            if (m.contains("film")) Fail(lw, "sets the film itself: write either \"film\" or \"lotion\"");
            SkinLotion::Params lp;
            double doseUm = double(lp.dose) * 1e6, timeS = 0, tau = lp.waterTau, hp = lp.hydrationPeak, hr = lp.hydrationRise, hd = lp.hydrationDecay, wi = lp.waterIor, hi = lp.humectantIor, wr = lp.waterRoughness, rr = lp.residueRoughness, gl = lp.glycerin;
            Optional<double>(lo, "doseUm", lw, doseUm, ReadNumber);
            Optional<double>(lo, "glycerin", lw, gl, ReadNumber);
            Optional<double>(lo, "timeS", lw, timeS, ReadNumber);
            Optional<double>(lo, "waterTauS", lw, tau, ReadNumber);
            Optional<double>(lo, "hydrationPeak", lw, hp, ReadNumber);
            Optional<double>(lo, "hydrationRiseS", lw, hr, ReadNumber);
            Optional<double>(lo, "hydrationDecayS", lw, hd, ReadNumber);
            Optional<double>(lo, "waterIor", lw, wi, ReadNumber);
            Optional<double>(lo, "humectantIor", lw, hi, ReadNumber);
            Optional<double>(lo, "waterRoughness", lw, wr, ReadNumber);
            Optional<double>(lo, "residueRoughness", lw, rr, ReadNumber);
            double bulk = 0; Optional<double>(lo, "bulk", lw, bulk, ReadNumber);   // the scatter distance of the skin grows by (1 + bulk * hydration): a hydrated skin scatters less inside
            if (!(bulk >= 0)) Fail(lw + ".bulk", "must not be negative");
            double oil = lp.oil, oilIor = lp.oilIor, dropletUm = double(lp.dropletDiameter) * 1e6, coalesce = lp.coalesceTau, dropletG = lp.dropletG;   // an emulsion: oil droplets in the water
            Optional<double>(lo, "oil", lw, oil, ReadNumber);
            Optional<double>(lo, "oilIor", lw, oilIor, ReadNumber);
            Optional<double>(lo, "dropletUm", lw, dropletUm, ReadNumber);
            Optional<double>(lo, "coalesceTauS", lw, coalesce, ReadNumber);
            Optional<double>(lo, "dropletG", lw, dropletG, ReadNumber);
            lp.oil = Real(oil); lp.oilIor = Real(oilIor); lp.dropletDiameter = Real(dropletUm * 1e-6); lp.coalesceTau = Real(coalesce); lp.dropletG = Real(dropletG);
            lp.dose = Real(doseUm * 1e-6); lp.glycerin = Real(gl); lp.waterTau = Real(tau); lp.hydrationPeak = Real(hp); lp.hydrationRise = Real(hr); lp.hydrationDecay = Real(hd);
            lp.waterIor = Real(wi); lp.humectantIor = Real(hi); lp.waterRoughness = Real(wr); lp.residueRoughness = Real(rr);
            double occlusion = lp.occlusion, occlusionOilUm = double(lp.occlusionOil) * 1e6;   // the merged oil film slows the loss of water from the stratum corneum
            Optional<double>(lo, "occlusion", lw, occlusion, ReadNumber);
            Optional<double>(lo, "occlusionOilUm", lw, occlusionOilUm, ReadNumber);
            lp.occlusion = Real(occlusion); lp.occlusionOil = Real(occlusionOilUm * 1e-6);
            const std::string why = SkinLotion::validate(lp);
            if (!why.empty()) Fail(lw, why);
            if (!(timeS >= 0)) Fail(lw + ".timeS", "must not be negative");
            // "layers": [ { "atS", "doseUm", "glycerin", "hydrationPeak", "oil", "oilIor", "dropletUm", "coalesceTauS", "dropletG" }, ... ]: several products, each applied at its time
            // `atS` (on the clock of `timeS`). Without it, the one product is the one written in "lotion" itself, applied at 0.
            std::vector<SkinLotion::Layer> layers;
            if (lo.contains("layers")) {
                const json& arr = lo.at("layers");
                if (!arr.is_array() || arr.empty()) Fail(lw + ".layers", "must be a list of products with at least one entry");
                for (const char* key : { "doseUm", "glycerin", "hydrationPeak", "oil", "oilIor", "dropletUm", "coalesceTauS", "dropletG" })
                    if (lo.contains(key)) Fail(lw + "." + key, "belongs to a product: with \"layers\", write it inside the layer");
                const SkinLotion::Layer defaults;
                for (size_t i = 0; i < arr.size(); ++i) {
                    const json& la = arr[i];
                    const std::string aw = lw + ".layers[" + std::to_string(i) + "]";
                    if (!la.is_object()) Fail(aw, "must be an object (atS, doseUm, glycerin, hydrationPeak, oil, oilIor, dropletUm, coalesceTauS, dropletG)");
                    CheckKeys(la, { "atS", "doseUm", "glycerin", "hydrationPeak", "oil", "oilIor", "dropletUm", "coalesceTauS", "dropletG" }, aw);
                    double atS = 0, dose = double(defaults.dose) * 1e6, glyc = defaults.glycerin, peak = defaults.hydrationPeak, loil = defaults.oil, loilIor = defaults.oilIor,
                        drop = double(defaults.dropletDiameter) * 1e6, coal = defaults.coalesceTau, g = defaults.dropletG;
                    Optional<double>(la, "atS", aw, atS, ReadNumber);
                    Optional<double>(la, "doseUm", aw, dose, ReadNumber);
                    Optional<double>(la, "glycerin", aw, glyc, ReadNumber);
                    Optional<double>(la, "hydrationPeak", aw, peak, ReadNumber);
                    Optional<double>(la, "oil", aw, loil, ReadNumber);
                    Optional<double>(la, "oilIor", aw, loilIor, ReadNumber);
                    Optional<double>(la, "dropletUm", aw, drop, ReadNumber);
                    Optional<double>(la, "coalesceTauS", aw, coal, ReadNumber);
                    Optional<double>(la, "dropletG", aw, g, ReadNumber);
                    if (!(atS >= 0)) Fail(aw + ".atS", "must not be negative");
                    SkinLotion::Layer l;
                    l.at = Real(atS); l.dose = Real(dose * 1e-6); l.glycerin = Real(glyc); l.hydrationPeak = Real(peak);
                    l.oil = Real(loil); l.oilIor = Real(loilIor); l.dropletDiameter = Real(drop * 1e-6); l.coalesceTau = Real(coal); l.dropletG = Real(g);
                    const std::string whyLayer = SkinLotion::validate(SkinLotion::withLayer(lp, l));
                    if (!whyLayer.empty()) Fail(aw, whyLayer);
                    layers.push_back(l);
                }
            }
            else layers.push_back(SkinLotion::layerOf(lp));
            int applied = 0;
            for (const SkinLotion::Layer& l : layers) applied += double(l.at) <= timeS ? 1 : 0;
            oil = 0;
            for (const SkinLotion::Layer& l : layers) oil = std::max(oil, double(l.at) <= timeS ? double(l.oil) : 0.0);
            const SkinLotion::State st = SkinLotion::at(lp, layers, Real(timeS));
            s.skinFilmOn = st.thickness > 0;
            s.skinFilmThicknessUm = double(st.thickness) * 1e6; s.skinFilmIor = double(st.ior); s.skinFilmRoughness = double(st.roughness);
            if (s.skinScatterOn) s.scatHydration = double(st.hydration);
            if (bulk > 0) s.skinSssScale *= 1.0 + bulk * double(st.hydration);
            s.skinEmulsionTau = double(st.emulsionTau);
            std::ostringstream note;
            note << "t = " << timeS << " s" << (lo.contains("layers") ? " (" + std::to_string(applied) + " of " + std::to_string(layers.size()) + " layers applied)" : std::string())
                << ": film " << s.skinFilmThicknessUm << " um (water " << double(st.water) * 1e6 << " + humectant " << double(st.residue) * 1e6 << "), index " << s.skinFilmIor
                << ", roughness " << s.skinFilmRoughness << ", hydration " << double(st.hydration) << (bulk > 0 ? ", scatter distance x " + std::to_string(1.0 + bulk * double(st.hydration)) : std::string())
                << (oil > 0 ? ", oil " + std::to_string(double(st.oil) * 1e6) + " um, droplets: transport optical depth " + std::to_string(double(st.emulsionTau))
                    + ", occlusion " + std::to_string(double(st.occlusion)) : std::string());
            s.skinLotionNote = note.str();
        }
        if (m.contains("pigment")) {   // "pigment": { "freckleCellMm", "freckleFraction", "freckleRadiusMm", "freckleStrength", "clusterFreq", "blotchFreq", "blotchStrength", "rednessFreq", "rednessStrength", "seed" } (all optional)
            const json& g = m.at("pigment");
            const std::string gw = where + ".pigment";
            if (!g.is_object()) Fail(gw, "must be an object (freckleCellMm, freckleFraction, freckleRadiusMm, freckleStrength, clusterFreq, blotchFreq, blotchStrength, rednessFreq, rednessStrength, seed)");
            CheckKeys(g, { "freckleCellMm", "freckleFraction", "freckleRadiusMm", "freckleStrength", "clusterFreq", "blotchFreq", "blotchStrength", "rednessFreq", "rednessStrength", "seed" }, gw);
            s.skinPigment = true;
            Optional<double>(g, "freckleCellMm", gw, s.pigFreckleCellMm, ReadNumber);
            Optional<double>(g, "freckleFraction", gw, s.pigFreckleFraction, ReadNumber);
            Optional<double>(g, "freckleRadiusMm", gw, s.pigFreckleRadiusMm, ReadNumber);
            Optional<double>(g, "freckleStrength", gw, s.pigFreckleStrength, ReadNumber);
            Optional<double>(g, "clusterFreq", gw, s.pigClusterFreq, ReadNumber);
            Optional<double>(g, "blotchFreq", gw, s.pigBlotchFreq, ReadNumber);
            Optional<double>(g, "blotchStrength", gw, s.pigBlotchStrength, ReadNumber);
            Optional<double>(g, "rednessFreq", gw, s.pigRednessFreq, ReadNumber);
            Optional<double>(g, "rednessStrength", gw, s.pigRednessStrength, ReadNumber);
            double sd = s.pigSeed; Optional<double>(g, "seed", gw, sd, ReadNumber); s.pigSeed = int(sd);
            if (!(s.pigFreckleCellMm > 0 && s.pigFreckleRadiusMm > 0 && s.pigClusterFreq > 0 && s.pigBlotchFreq > 0 && s.pigRednessFreq > 0)) Fail(gw, "cell sizes, radius and frequencies must be above 0");
            for (double v : { s.pigFreckleFraction, s.pigFreckleStrength, s.pigBlotchStrength, s.pigRednessStrength })
                if (v < 0 || v > 1) Fail(gw, "fractions and strengths must be between 0 and 1");
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
        Fail(where + ".type", "is \"" + type + "\", must be diffuse, conductor, layered, dielectric, thinFilm, hair, skin or emission");
    }
    return s;
}

} // namespace rayt::cli::scenefile
