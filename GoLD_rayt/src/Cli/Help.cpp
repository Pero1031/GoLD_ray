// The text printed by --help and --list-scenes.
#include "pch.h"
#include "Cli/Help.hpp"
#include "Cli/Commands.hpp"
#include "Cli/Scene.hpp"

#include <iomanip>
#include <iostream>

namespace rayt::cli {

void PrintSceneList(std::ostream& out) {
    out << "scenes (they render a picture):\n";
    for (const SceneDef* s : AllScenes()) out << "  " << std::left << std::setw(13) << s->name << s->description << "\n";
    out << "\ncommands (they print a result; select them with --scene like a scene):\n";
    for (const CommandDef* c : AllCommands()) out << "  " << std::left << std::setw(13) << c->name << c->description << "\n";
}

void PrintHelp(std::ostream& out) {
    out << "usage: GoLD_rayt [scene.json] [options]\n"
        << "\n"
        << "  scene.json        a scene file (see the scenes/ folder); it must come first, and the options after it win over the file\n"
        << "  --scene NAME      a built-in scene or command (--list-scenes)\n"
        << "  --list-scenes     list the scenes and commands\n"
        << "  --help            this text\n"
        << "\n"
        << "common options:\n"
        << "  --width N         picture width (the height follows the format of the scene)   --height N   set the height as well\n"
        << "  --spp N           samples per pixel      --depth N   maximum path depth\n"
        << "  --out FILE        .png / .bmp / .jpg (tone mapped), .hdr or .pfm (linear float; .pfm keeps negative values)\n"
        << "  --gamut srgb|p3|rec2020   primaries of the output (default srgb): the .png gets a cICP tag, .pfm / .hdr hold linear values in them\n"
        << "  --bits 8|16       bits per channel of a .png (default 8; 16 keeps a wide gamut free of banding)\n"
        << "  --env FILE        environment map (.hdr)      --light d65|tungsten|candle|led   light spectrum of some scenes\n"
        << "  --denoise  --aov  Intel OIDN denoised copy; albedo / normal images\n"
        << "  --caustic-trace   with --caustics: light tracing of the caustics on surfaces the camera sees directly (sharp edges, no radius)\n"
        << "  --light-sampler uniform|power|spatial|all     --accel own|embree\n"
        << "\n"
        << "The options of single scenes (--sea-steep, --pmx-yaw, --caustics, ...) are listed in include/Cli/Options.hpp.\n"
        << "\n";
    PrintSceneList(out);
}

} // namespace rayt::cli
