#pragma once
// The text printed by --help and --list-scenes. The scene and command lists come from the tables (Scene.hpp, Commands.hpp).

#include <iosfwd>

namespace rayt::cli {

void PrintHelp(std::ostream& out);
void PrintSceneList(std::ostream& out);

} // namespace rayt::cli
