/**
 * @file IO/StlLoader.hpp
 * @brief STL reader (binary and ascii) into a TriangleMesh.
 *
 * - Binary STL: an 80-byte header, a uint32 triangle count, then 50 bytes per triangle (normal, 3 vertices, a 2-byte attribute). A file is
 *   taken as binary when its size is exactly 84 + 50 * count: the header of a binary file may itself start with "solid" (Materialise
 *   Magics writes "COLOR=... MATERIAL=..."), so the size is the reliable test. Anything else that starts with "solid" is read as ascii.
 * - STL has no shared vertices (every triangle repeats its three corners), so vertices at exactly the same position are welded.
 *   Without that a mesh would have no smooth normals and 3x the vertices.
 * - Normals are rebuilt from the winding (area weighted, counter-clockwise = outside, the STL convention); the normals stored in the
 *   file are ignored, since writers often leave them zero or flipped. The result is a smooth-shaded mesh.
 * - Triangles that collapse to a line or a point after welding are dropped. Colors (the attribute bytes) and uv are not read: STL has none.
 * - The unit is whatever the file used (usually mm) and Z is usually up: use the scene file's scale and rotate (docs/SCENE_FILE.md, "mesh").
 */

#pragma once

#include <string>

#include "Geometry/TriangleMesh.hpp"

namespace rayt::io {

    /// Loads a binary or ascii STL file. Returns false with a message in outError on failure.
    bool loadStl(const std::string& filename, TriangleMesh& out, std::string* outError = nullptr);

} // namespace rayt::io
