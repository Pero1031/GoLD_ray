#pragma once
// Scene files: JSON documents with two forms.
//
//   1. A built-in scene with options (what the command line says, saved in a file):
//        { "scene": "goldbowl", "options": { "light": "tungsten", "width": 800, "spp": 256 } }
//      "options" uses the names of the command line without the leading "--" (a number or string is the value, true switches a flag on).
//
//   2. A scene written as data: "camera", "environment", "materials", "lights" and "objects" (see SceneDescription.hpp and scenes/*.json).
//      A file with "objects" is rendered as the scene "file". "options" may be given as well (width, spp, out ...).
//
// `//` comments are allowed. Relative paths in a scene file are relative to the folder the program runs in, like the paths on the
// command line. Every misspelled key is an error that lists the keys that are allowed there.

#include <memory>
#include <string>
#include <vector>

namespace rayt::cli {

struct SceneDescription;

// Reads a scene file. `args` receives the command line it stands for ("--scene", "goldbowl", "--light", "tungsten", ...), so that the
// normal option parser checks every value; `description` receives the scene of form 2 (null for form 1).
// Returns false and a message (with line and column for a syntax error) when the file cannot be read or is malformed.
bool LoadSceneFile(const std::string& path, std::vector<std::string>& args, std::shared_ptr<const SceneDescription>& description,
                   std::string& error);

} // namespace rayt::cli
