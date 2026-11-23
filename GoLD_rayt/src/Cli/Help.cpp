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
        << "  --seed N          the same picture again, noise included (path tracing; not the photon map or --lvc). Without it every run\n"
        << "                    draws a new seed and prints it\n"
        << "  --out FILE        .png / .bmp / .jpg (tone mapped), .hdr or .pfm (linear float; .pfm keeps negative values)\n"
        << "  --gamut srgb|p3|rec2020   primaries of the output (default srgb): the .png gets a cICP tag, .pfm / .hdr hold linear values in them\n"
        << "  --bits 8|16       bits per channel of a .png (default 8; 16 keeps a wide gamut free of banding)\n"
        << "  --exposure EV     the .png / .bmp / .jpg are 2^EV times brighter (.pfm / .hdr keep the radiance as rendered)\n"
        << "  --tonemap reinhard|clip   reinhard (default): x / (1 + x), highlights compressed; clip: linear up to white, clipped above it\n"
        << "                    (for stimuli: brightness ratios kept below white; it prints how many pixels clipped and how far the brightest is from white)\n"
        << "  --focal-length MM  lens as in photography (35 mm film; the long side of the picture is 36 mm): replaces the field of view of the scene\n"
        << "  --lens FILE       a real lens: its prescription (assets/Lenses/*.txt) traced surface by surface replaces the thin lens; the lens and\n"
        << "                    the 36 mm film decide the field of view; --fstop closes its stop, --blades shapes it, the focus is as for the thin lens\n"
        << "  --ortho H         orthographic view H (scene units) high: parallel rays, sizes independent of the distance (0 = perspective; not with light tracing)\n"
        << "  --fstop N         defocus blur: aperture = focal length / N (scenes are in metres; without it a pinhole)\n"
        << "  --focus-dist D | --focus-at x,y,z   where the plane of focus is (default: the point the camera looks at)\n"
        << "  --blades N [--blade-rotation DEG]   polygonal aperture (3 to 32 blades): out-of-focus highlights take its shape\n"
        << "  --diffraction S   diffraction at the aperture (needs --fstop): the star of bright points, n spikes for n blades (2n if odd), colored\n"
        << "                    by the wavelength. S = 1: the real size (pixels of 36 mm / long side); S > 1 enlarges the pattern S times\n"
        << "  --ghosts coated|uncoated   with --lens: the lens flare of the scene's lights, light reflected twice between the lens surfaces\n"
        << "                    (one ghost per pair of surfaces), behind a single anti-reflection coating or bare glass (spectral build: with dispersion)\n"
        << "  --distortion K1[,K2]   radial distortion r' = r (1 + K1 r^2 + K2 r^4), r in half diagonals: K1 > 0 barrel, K1 < 0 pincushion\n"
        << "  --lateral-ca A    lateral chromatic aberration: the picture scales by A per 100 nm (0.004: red and blue differ by 0.8 %) (spectral build)\n"
        << "  --axial-ca D      axial chromatic aberration: the focal length changes by D per 100 nm, colored defocus with --fstop (0.01: a singlet) (spectral build)\n"
        << "  --vignette N [--vignette-dist MM]   lens barrel opening of f/N (diameter = focal length / N) MM in front of the stop (default: half the focal length):\n"
        << "                    darkens the corners and makes the out-of-focus highlights there cat eyes; gone once the lens is stopped down (needs --fstop)\n"
        << "  --natural-vignette S   the cos^4 law of the film (0 = off, 1 = full; no aperture needed)\n"
        << "  --shutter T       motion blur: the shutter stays open T seconds (0.02 or 1/50); scene objects with a \"velocity\" (m/s) and a moving camera blur over it\n"
        << "  --iso S --nits K  physical exposure: the f-number, the shutter and the ISO speed S set the brightness as in a camera (f/16 is 64 x darker than\n"
        << "                    f/2); K is the luminance (cd/m^2) of a radiance whose Y is 1; mid grey 0.18 = the exposure 10 / S lux s (.png / .bmp / .jpg)\n"
        << "  --cam-velocity x,y,z   the camera moves with this velocity (m/s) during the exposure      --pan YAW[,PITCH]   it turns (degrees per second; right / up)\n"
        << "  --env FILE        environment map (.hdr)     --light d65|d65uv|tungsten|candle|led  light spectrum of some scenes\n"
        << "  --denoise  --aov  Intel OIDN denoised copy; albedo / normal images\n"
        << "  --mask A[,B..]    <out>_mask.png: the share of each pixel that objects of these materials (names in the scene file) cover, as the camera\n"
        << "                    sees them first (edges anti-aliased, through the lens and shutter; linear gray, 16 bit with --bits 16; .pfm if --out is)\n"
        << "  --toon  --toon-line PX   experimental cel shading with outlines (PX: line width, 0 = none) instead of path tracing\n"
        << "  --depth-map       <out>_depth.pfm: distance of the first surface along the view axis through each pixel centre (0 = none)\n"
        << "  --camera-response FILE   the film as this camera sees (spectral build): its R, G, B sensitivities (CSV of wavelength, R, G, B),\n"
        << "                    balanced to D65 (--wb: to the scene's light) and mapped to sRGB by a fitted colour matrix, instead of the CIE\n"
        << "                    observer; also <out>_camraw.pfm\n"
        << "  --awb             auto white balance: the light guessed from the picture (grey world), the picture balanced to it (a camera's\n"
        << "                    gains and matrix with --camera-response, else CAT16 to D65)\n"
        << "  --sensor default|well=E,base=S,read=E,adc=E,bits=B,cfa=rggb|none,demosaic=malvar|bilinear,highlights=clip|keep   the pictures\n"
        << "                    as a sensor records them: shot noise, read noise before and after the ISO amplifier (electrons), the full well\n"
        << "                    (electrons at the base ISO's white), a B-bit ADC, a Bayer mosaic demosaiced, clipped highlights kept white\n"
        << "                    (also <out>.dng, the raw file); with --iso its ISO\n"
        << "                    (default: well=50000,base=100,read=1.5,adc=3,bits=14,cfa=rggb,demosaic=malvar,highlights=clip)\n"
        << "  --spectral-image STEP[,FIRST,LAST]   <out>_spectral.img + .hdr (ENVI): the mean spectral radiance of each pixel in bands of STEP nm\n"
        << "                    (default 380..780 nm; spectral build). Each band gets a share of the 4 wavelengths of a path: use more spp\n"
        << "  --stokes          polarization: Stokes images <out>_S0..S3.pfm, degree / angle of linear polarization (smooth glass and water; docs/POLARIZATION.md)\n"
        << "  --polarizer A[,A2]   a linear polarizer at A degrees in front of the camera (0 = horizontal): <out>_pol<A>.*   (implies --stokes)\n"
        << "  --sky-pol P[,EL,AZ]  the environment is a clear sky polarized by Rayleigh scattering: degree P (0..1, clear sky 0.6-0.8) at 90 deg\n"
        << "                    from the sun at elevation EL, azimuth AZ (degrees, from +z toward +x; default: the brightest spot of the map).\n"
        << "                    Only the Stokes and polarizer images show it (implies --stokes; docs/POLARIZATION.md)\n"
        << "  --caustic-trace   with --caustics: light tracing of the caustics on surfaces the camera sees directly (sharp edges, no radius)\n"
        << "  --fog-sigma S, --fog-light I, --fog-fill F, --fog-prism-deg D   fogbeam scene: extinction of the fog (1/m, default 1.5), intensity of the lamp (60) and of a faint fill light (0.4), rotation of the prism\n"
        << "  --lvc [--lvc-paths N --lvc-depth D --lvc-mask M --lvc-merge R]   fogspot scene: bidirectional path tracing with a cache of light paths (prototype)   --spot-deg A\n"
        << "                    --lvc-merge R: vertex merging with radius R (scene units) in the first pass, shrinking with --lvc-passes (caustics seen through glass; 0 = off)\n"
        << "                    --lvc-mask: 1 = NEE, 2 = connections, 4 = light tracing, 8 = merging, 16 / 32 = manifold NEE / connections through flat-faced glass (scene with a prism)\n"
        << "  --light-sampler uniform|power|spatial|all     --accel own|embree\n"
        << "\n"
        << "The options of single scenes (--sea-steep, --pmx-yaw, --caustics, ...) are listed in include/Cli/Options.hpp.\n"
        << "\n";
    PrintSceneList(out);
}

} // namespace rayt::cli
