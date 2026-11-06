#pragma once
// The commands of the CLI that do not render a picture: self-checks, measurement tools and the model / pose lists.

#include <optional>
#include <string>
#include <vector>

#include "Cli/Options.hpp"

namespace rayt::cli {

// A command is selected with --scene NAME like a scene, but it prints a result instead of writing a picture.
struct CommandDef {
    const char* name;
    const char* description;   // one line, shown by --list-scenes
};

const CommandDef* FindCommand(const std::string& name);
std::vector<const CommandDef*> AllCommands();

// Runs the command selected by --scene (lights, aov, sss, ..., bench) or by --pmx-list / --pose-list.
// Returns its exit code, or std::nullopt when the options ask for a picture.
std::optional<int> RunNonRenderCommand(const CliOptions& opt);

} // namespace rayt::cli
