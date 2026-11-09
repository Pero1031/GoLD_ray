#pragma once
// What the parts of the scene-file reader share: the small readers of JSON values. Internal to the reader
// (SceneFile.cpp: the scene and its objects; SceneFileMaterials.cpp: the materials); other code uses Cli/SceneFile.hpp.

#include <initializer_list>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include "Cli/SceneDescription.hpp"

namespace rayt::cli::scenefile {

using json = nlohmann::json;

// ---- small readers: every one names the place it is reading in its error message --------------------------------------------------

[[noreturn]] inline void Fail(const std::string& where, const std::string& what) { throw std::runtime_error(where + " " + what); }

// Throws unless every key of `obj` is one of `allowed`; the message lists the allowed ones.
inline void CheckKeys(const json& obj, std::initializer_list<const char*> allowed, const std::string& where) {
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

inline double ReadNumber(const json& v, const std::string& where) {
    if (!v.is_number()) Fail(where, "must be a number");
    return v.get<double>();
}

inline bool ReadBool(const json& v, const std::string& where) {
    if (!v.is_boolean()) Fail(where, "must be true or false");
    return v.get<bool>();
}

inline std::string ReadString(const json& v, const std::string& where) {
    if (!v.is_string()) Fail(where, "must be a string");
    return v.get<std::string>();
}

inline Vector3 ReadVec3(const json& v, const std::string& where) {
    if (!v.is_array() || v.size() != 3 || !v[0].is_number() || !v[1].is_number() || !v[2].is_number())
        Fail(where, "must be an array of three numbers, like [0, 1, 0]");
    return Vector3(v[0].get<double>(), v[1].get<double>(), v[2].get<double>());
}

// An RGB color: three numbers, or one number for a gray.
inline Vector3 ReadColor(const json& v, const std::string& where) {
    if (v.is_number()) return Vector3(v.get<double>());
    return ReadVec3(v, where);
}

template <class T>
void Optional(const json& obj, const char* key, const std::string& where, T& target, T (*read)(const json&, const std::string&)) {
    if (obj.contains(key)) target = read(obj.at(key), where + "." + key);
}

// One material of "materials" (SceneFileMaterials.cpp). `where` names it in the error messages, e.g. "materials.gold".
MaterialSpec ParseMaterial(const json& m, const std::string& where);

} // namespace rayt::cli::scenefile
