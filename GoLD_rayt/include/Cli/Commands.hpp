#pragma once
// The commands of the CLI that do not render a picture: self-checks, measurement tools and the model / pose lists.

#include <optional>

#include "Cli/Options.hpp"

namespace rayt::cli {

// Runs the command selected by --scene (lights, aov, sss, ..., bench) or by --pmx-list / --pose-list.
// Returns its exit code, or std::nullopt when the options ask for a picture.
std::optional<int> RunNonRenderCommand(const CliOptions& opt);

} // namespace rayt::cli
