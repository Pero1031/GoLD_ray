"""Cleans a triangle mesh (OBJ) for the renderer: welds the vertices, throws away the file's normals and uv's, recomputes smooth
(area-weighted) vertex normals from the geometry, and optionally closes the open boundaries.

    python tools/mesh_clean.py in.obj out.obj [--cap]

Why: a scan can carry vertex normals that disagree with the triangle winding (the Erato scan: 25% of the triangles), and a metal
or glass material then shades part of the surface black. Closing the boundaries (--cap: a triangle fan from the centroid of each
boundary loop) is needed where the mesh is filled with a participating medium (a closed surface with outward normals).
The output has v and vn only (no texture coordinates), so it is also much smaller and faster to load.

Needs numpy. Triangles only (a polygon is triangulated as a fan). The winding of the input must be consistent (counter-clockwise
seen from outside); the script reports the signed volume, which is positive for outward normals.
"""
import argparse
import collections
import sys

import numpy as np


def read_obj(path):
    verts, faces = [], []
    with open(path, errors="ignore") as f:
        for line in f:
            if line.startswith("v "):
                verts.append(line.split()[1:4])
            elif line.startswith("f "):
                idx = [int(tok.split("/")[0]) for tok in line.split()[1:]]
                idx = [i - 1 if i > 0 else len(verts) + i for i in idx]
                for k in range(1, len(idx) - 1):          # a polygon becomes a fan
                    faces.append((idx[0], idx[k], idx[k + 1]))
    return np.array(verts, dtype=np.float64), np.array(faces, dtype=np.int64)


def weld(verts, faces, digits=5):
    key = np.round(verts * 10 ** digits).astype(np.int64)
    uniq, first, inverse = np.unique(key, axis=0, return_index=True, return_inverse=True)
    faces = inverse[faces]
    ok = (faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2]) & (faces[:, 0] != faces[:, 2])   # drop collapsed triangles
    return verts[first], faces[ok]


def boundary_loops(faces):
    """The loops of boundary edges (edges used by exactly one triangle), as lists of vertex indices in the direction of the edges."""
    directed = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    count = collections.Counter(map(tuple, np.sort(directed, axis=1).tolist()))
    edges = [e for e in map(tuple, directed.tolist()) if count[tuple(sorted(e))] == 1]
    following = collections.defaultdict(list)
    for a, b in edges:
        following[a].append(b)
    used, loops = set(), []
    for a, b in edges:
        if (a, b) in used:
            continue
        loop, cur = [a], (a, b)
        used.add(cur)
        while cur[1] != a:
            loop.append(cur[1])
            nxt = [x for x in following[cur[1]] if (cur[1], x) not in used]
            if not nxt:
                break
            cur = (cur[1], nxt[0])
            used.add(cur)
        loops.append(loop)
    return loops


def cap(verts, faces, loops):
    verts = list(verts)
    new = []
    for loop in loops:
        c = np.mean([verts[i] for i in loop], axis=0)
        ci = len(verts)
        verts.append(c)
        for k, a in enumerate(loop):
            b = loop[(k + 1) % len(loop)]
            new.append((b, a, ci))          # the opposite direction of the boundary edge a -> b
    return np.array(verts), np.concatenate([faces, np.array(new, dtype=np.int64)])


def smooth_normals(verts, faces):
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    fn = np.cross(b - a, c - a)                  # the length is twice the area: an area weight
    n = np.zeros_like(verts)
    for k in range(3):
        np.add.at(n, faces[:, k], fn)
    length = np.linalg.norm(n, axis=1, keepdims=True)
    length[length == 0] = 1
    return n / length


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--cap", action="store_true", help="close the open boundaries with a triangle fan")
    args = ap.parse_args()

    verts, faces = read_obj(args.input)
    print("read: %d vertices, %d triangles" % (len(verts), len(faces)))
    verts, faces = weld(verts, faces)
    print("welded: %d vertices, %d triangles" % (len(verts), len(faces)))
    loops = boundary_loops(faces)
    print("open boundaries: %d loops %s" % (len(loops), sorted(map(len, loops), reverse=True)[:8]))
    if args.cap and loops:
        verts, faces = cap(verts, faces, loops)
        print("closed: %d vertices, %d triangles" % (len(verts), len(faces)))
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    volume = np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6
    print("signed volume %.3f (positive: outward normals)" % volume)
    if volume < 0:
        print("warning: negative volume - the winding is inside out", file=sys.stderr)
    normals = smooth_normals(verts, faces)
    with open(args.output, "w") as f:
        f.write("# cleaned by tools/mesh_clean.py: welded, smooth normals recomputed%s\n" % (", boundaries closed" if args.cap else ""))
        np.savetxt(f, verts, fmt="v %.6f %.6f %.6f")
        np.savetxt(f, normals, fmt="vn %.5f %.5f %.5f")
        np.savetxt(f, np.stack([faces[:, 0] + 1, faces[:, 0] + 1, faces[:, 1] + 1, faces[:, 1] + 1, faces[:, 2] + 1, faces[:, 2] + 1], axis=1),
                   fmt="f %d//%d %d//%d %d//%d")
    print("wrote", args.output)


if __name__ == "__main__":
    main()
