"""Prepares a triangle mesh (OBJ or PLY) for the renderer: welds the vertices, recomputes smooth (area-weighted) vertex normals
from the geometry, and optionally keeps the texture coordinates, drops the small disconnected pieces, turns the up axis, scales the
mesh to a height in metres, puts its base on the floor and closes the open boundaries.

    python tools/data/mesh_clean.py in.obj out.obj [--cap]
    python tools/data/mesh_clean.py mesh.ply out.ply --keep-uv --up=-y --drop-small 0.01 --height 0.3 --ground

Why: a scan can carry vertex normals that disagree with the triangle winding (the Erato scan: 25% of the triangles), and a metal
or glass material then shades part of the surface black. A mesh extracted from 3D Gaussian Splatting (2DGS, GOF, SuGaR) comes in
the camera's units and frame (COLMAP: y down), with floating fragments ("floaters") around the object; --up, --height, --ground
and --drop-small bring it into the renderer's metres with y up and the base at y = 0 (docs/CAPTURE.md). Closing the boundaries
(--cap: a triangle fan from the centroid of each boundary loop) is needed where the mesh is filled with a participating medium
(a closed surface with outward normals).

The file's normals are always thrown away and recomputed. Without --keep-uv the texture coordinates are thrown away too, and
the output is smaller and faster to load; with --keep-uv a vertex is welded only with vertices of the same position AND uv (the
seams of the texture atlas stay), and the output carries vt (OBJ) or s, t (PLY) for a "texture" material.

Needs numpy (and scipy for --drop-small). Triangles only (a polygon is triangulated as a fan). PLY: ascii or binary, the usual
vertex properties (x y z, nx ny nz, s t / u v / texture_u texture_v, colors are ignored), faces as vertex_indices, and MeshLab's
per-face texcoord list. The winding of the input must be consistent (counter-clockwise seen from outside); the script reports the
signed volume, which is positive for outward normals, and --flip inverts the winding if it is negative.
"""
import argparse
import collections
import os
import sys

import numpy as np


# ----------------------------------------------------------------------------------------------------------------------------------
# reading
# ----------------------------------------------------------------------------------------------------------------------------------

def read_obj(path):
    """Returns (verts (N,3), faces (F,3) of vertex indices, corner_uv (F,3,2) or None)."""
    verts, texs, faces, corner_uv = [], [], [], []
    any_uv = False
    with open(path, errors="ignore") as f:
        for line in f:
            if line.startswith("v "):
                verts.append(line.split()[1:4])
            elif line.startswith("vt "):
                texs.append(line.split()[1:3])
            elif line.startswith("f "):
                toks = line.split()[1:]
                idx, uvi = [], []
                for tok in toks:
                    parts = tok.split("/")
                    i = int(parts[0])
                    idx.append(i - 1 if i > 0 else len(verts) + i)
                    if len(parts) > 1 and parts[1]:
                        t = int(parts[1])
                        uvi.append(t - 1 if t > 0 else len(texs) + t)
                        any_uv = True
                    else:
                        uvi.append(-1)
                for k in range(1, len(idx) - 1):          # a polygon becomes a fan
                    faces.append((idx[0], idx[k], idx[k + 1]))
                    corner_uv.append((uvi[0], uvi[k], uvi[k + 1]))
    verts = np.array(verts, dtype=np.float64)
    faces = np.array(faces, dtype=np.int64)
    uv = None
    if any_uv and texs:
        texs = np.array(texs, dtype=np.float64)
        ci = np.array(corner_uv, dtype=np.int64)
        uv = np.zeros((len(faces), 3, 2))
        ok = ci >= 0
        uv[ok] = texs[ci[ok]]
    return verts, faces, uv


_PLY_TYPES = {"char": "i1", "int8": "i1", "uchar": "u1", "uint8": "u1", "short": "i2", "int16": "i2", "ushort": "u2", "uint16": "u2",
              "int": "i4", "int32": "i4", "uint": "u4", "uint32": "u4", "float": "f4", "float32": "f4", "double": "f8", "float64": "f8"}


def _ply_header(f):
    """Parses the header of an open binary file. Returns (format, elements) where elements = [(name, count, [props])] and a
    prop is (name, type, count_type or None)."""
    if f.readline().strip() != b"ply":
        raise ValueError("not a PLY file")
    fmt, elements = None, []
    while True:
        line = f.readline()
        if not line:
            raise ValueError("no end_header")
        toks = line.decode("ascii", errors="ignore").split()
        if not toks:
            continue
        if toks[0] == "format":
            fmt = toks[1]
        elif toks[0] == "element":
            elements.append((toks[1], int(toks[2]), []))
        elif toks[0] == "property":
            if toks[1] == "list":
                elements[-1][2].append((toks[4], _PLY_TYPES[toks[3]], _PLY_TYPES[toks[2]]))
            else:
                elements[-1][2].append((toks[2], _PLY_TYPES[toks[1]], None))
        elif toks[0] == "end_header":
            return fmt, elements


def read_ply(path):
    """Returns (verts (N,3), faces (F,3), corner_uv (F,3,2) or None). Vertex uv (s t / u v / texture_*) or MeshLab's per-face texcoord."""
    with open(path, "rb") as f:
        fmt, elements = _ply_header(f)
        endian = "<" if fmt == "binary_little_endian" else ">"
        verts = faces = None
        vert_uv = corner_uv = None
        for name, count, props in elements:
            if fmt == "ascii":
                rows = [f.readline().split() for _ in range(count)]
                if name == "vertex":
                    cols = {p[0]: i for i, p in enumerate(props)}
                    arr = np.array([[float(x) for x in r] for r in rows], dtype=np.float64) if count else np.zeros((0, len(props)))
                    verts = arr[:, [cols["x"], cols["y"], cols["z"]]]
                    vert_uv = _vertex_uv(arr, cols)
                elif name == "face":
                    faces, corner_uv = _faces_from_rows(rows, props)
                continue
            # binary: scalar-only elements in one go; a face element as fixed-size records when every list has the same length
            if all(p[2] is None for p in props):
                dt = np.dtype([(p[0], endian + p[1]) for p in props])
                arr = np.frombuffer(f.read(dt.itemsize * count), dtype=dt, count=count)
                if name == "vertex":
                    cols = {p[0]: i for i, p in enumerate(props)}
                    arr2 = np.stack([arr[p[0]].astype(np.float64) for p in props], axis=1) if count else np.zeros((0, len(props)))
                    verts = arr2[:, [cols["x"], cols["y"], cols["z"]]]
                    vert_uv = _vertex_uv(arr2, cols)
                continue
            start = f.tell()
            first = _read_record(f, props, endian)                      # the list lengths of the first record
            f.seek(start)
            fields = []
            for p, v in zip(props, first):
                if p[2] is None:
                    fields.append((p[0], endian + p[1]))
                else:
                    fields.append((p[0] + "_n", endian + p[2]))
                    fields.append((p[0], endian + p[1], (len(v),)))
            dt = np.dtype(fields)
            data = f.read(dt.itemsize * count)
            arr = np.frombuffer(data, dtype=dt, count=count) if len(data) == dt.itemsize * count else None
            same = arr is not None and all((arr[p[0] + "_n"] == len(v)).all() for p, v in zip(props, first) if p[2] is not None)
            if not same:                                                # ragged lists: one record at a time
                f.seek(start)
                rows = [_read_record(f, props, endian) for _ in range(count)]
                if name == "face":
                    faces, corner_uv = _faces_from_rows(rows, props)
                continue
            if name == "face":
                key = "vertex_indices" if "vertex_indices" in dt.names else "vertex_index"
                idx = arr[key].astype(np.int64)
                faces = np.concatenate([np.stack([idx[:, 0], idx[:, k], idx[:, k + 1]], axis=1) for k in range(1, idx.shape[1] - 1)])
                if "texcoord" in dt.names:
                    tc = arr["texcoord"].astype(np.float64).reshape(count, -1, 2)
                    corner_uv = np.concatenate([np.stack([tc[:, 0], tc[:, k], tc[:, k + 1]], axis=1) for k in range(1, tc.shape[1] - 1)])
    if verts is None or faces is None:
        raise ValueError("the PLY has no vertex or no face element")
    uv = corner_uv
    if uv is None and vert_uv is not None:
        uv = vert_uv[faces]
    return verts, faces, uv


def _vertex_uv(arr, cols):
    for a, b in (("s", "t"), ("u", "v"), ("texture_u", "texture_v"), ("texture_s", "texture_t")):
        if a in cols and b in cols:
            return arr[:, [cols[a], cols[b]]]
    return None


def _read_record(f, props, endian):
    out = []
    for name, typ, count_type in props:
        if count_type is None:
            out.append(np.frombuffer(f.read(int(typ[1])), dtype=endian + typ)[0])
        else:
            n = int(np.frombuffer(f.read(int(count_type[1])), dtype=endian + count_type)[0])
            out.append(np.frombuffer(f.read(int(typ[1]) * n), dtype=endian + typ))
    return out


def _faces_from_rows(rows, props):
    """Faces and per-corner uv from ascii rows (tokens) or binary records (values)."""
    faces, corner_uv, any_uv = [], [], False
    for r in rows:
        pos, idx, tc = 0, None, None
        for name, typ, count_type in props:
            if count_type is None:
                pos += 1
                continue
            if isinstance(r[pos], (bytes, str)):                  # ascii tokens: count then items
                n = int(r[pos])
                vals = r[pos + 1:pos + 1 + n]
                pos += 1 + n
            else:                                                   # a binary record: the list is one array
                vals = r[pos]
                pos += 1
            if name in ("vertex_indices", "vertex_index"):
                idx = [int(x) for x in vals]
            elif name == "texcoord":
                tc = [float(x) for x in vals]
                any_uv = True
        for k in range(1, len(idx) - 1):
            faces.append((idx[0], idx[k], idx[k + 1]))
            if tc is not None:
                corner_uv.append(((tc[0], tc[1]), (tc[2 * k], tc[2 * k + 1]), (tc[2 * k + 2], tc[2 * k + 3])))
            else:
                corner_uv.append(((0, 0), (0, 0), (0, 0)))
    return np.array(faces, dtype=np.int64), (np.array(corner_uv, dtype=np.float64) if any_uv else None)


def read_mesh(path):
    ext = os.path.splitext(path)[1].lower()
    return read_ply(path) if ext == ".ply" else read_obj(path)


# ----------------------------------------------------------------------------------------------------------------------------------
# processing
# ----------------------------------------------------------------------------------------------------------------------------------

def weld(verts, faces, corner_uv=None, digits=5):
    """Merges vertices of the same position (and the same uv, if corner_uv is given). Returns (verts, faces, uv per vertex or None)."""
    key = np.round(verts * 10 ** digits).astype(np.int64)
    if corner_uv is None:
        uniq, first, inverse = np.unique(key, axis=0, return_index=True, return_inverse=True)
        inverse = inverse.reshape(-1)
        new_faces = inverse[faces]
        new_verts, uv = verts[first], None
    else:
        pos_key = key[faces.reshape(-1)]                                             # a row per corner
        uv_key = np.round(corner_uv.reshape(-1, 2) * 10 ** digits).astype(np.int64)
        corner_key = np.concatenate([pos_key, uv_key], axis=1)
        uniq, first, inverse = np.unique(corner_key, axis=0, return_index=True, return_inverse=True)
        inverse = inverse.reshape(-1)
        new_faces = inverse.reshape(-1, 3)
        new_verts = verts[faces.reshape(-1)[first]]
        uv = corner_uv.reshape(-1, 2)[first]
    ok = (new_faces[:, 0] != new_faces[:, 1]) & (new_faces[:, 1] != new_faces[:, 2]) & (new_faces[:, 0] != new_faces[:, 2])   # drop collapsed triangles
    return new_verts, new_faces[ok], uv


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


def cap(verts, faces, loops, uv=None):
    verts = list(verts)
    uv = None if uv is None else list(uv)
    new = []
    for loop in loops:
        c = np.mean([verts[i] for i in loop], axis=0)
        ci = len(verts)
        verts.append(c)
        if uv is not None:
            uv.append(np.mean([uv[i] for i in loop], axis=0))
        for k, a in enumerate(loop):
            b = loop[(k + 1) % len(loop)]
            new.append((b, a, ci))          # the opposite direction of the boundary edge a -> b
    return np.array(verts), np.concatenate([faces, np.array(new, dtype=np.int64)]), (None if uv is None else np.array(uv))


def drop_small_components(verts, faces, uv, fraction):
    """Removes the connected pieces (by shared vertices) with fewer than fraction x all triangles. Returns (verts, faces, uv, report)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    n = len(verts)
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    g = coo_matrix((np.ones(len(e), dtype=np.int8), (e[:, 0], e[:, 1])), shape=(n, n))
    ncomp, label = connected_components(g, directed=False)
    face_label = label[faces[:, 0]]
    sizes = np.bincount(face_label, minlength=ncomp)
    keep = sizes >= max(1, int(round(fraction * len(faces))))
    kept_faces = faces[keep[face_label]]
    report = "%d pieces, triangles %s; kept %d piece(s)" % (ncomp, sorted(sizes.tolist(), reverse=True)[:8], int(keep.sum()))
    return _compact(verts, kept_faces, uv) + (report,)


def _compact(verts, faces, uv):
    """Drops the vertices that no face uses."""
    used = np.zeros(len(verts), dtype=bool)
    used[faces.reshape(-1)] = True
    remap = np.cumsum(used) - 1
    return verts[used], remap[faces], (None if uv is None else uv[used])


_UP = {"y": np.eye(3), "-y": np.diag([1.0, -1.0, -1.0]),                          # a half turn about x: y down -> y up
       "z": np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]], dtype=np.float64),        # z up -> y up
       "-z": np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]], dtype=np.float64)}


def smooth_normals(verts, faces):
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    fn = np.cross(b - a, c - a)                  # the length is twice the area: an area weight
    n = np.zeros_like(verts)
    for k in range(3):
        np.add.at(n, faces[:, k], fn)
    length = np.linalg.norm(n, axis=1, keepdims=True)
    length[length == 0] = 1
    return n / length


def signed_volume(verts, faces):
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    return np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6


# ----------------------------------------------------------------------------------------------------------------------------------
# writing
# ----------------------------------------------------------------------------------------------------------------------------------

def write_obj(path, verts, normals, faces, uv, comment):
    with open(path, "w") as f:
        f.write("# %s\n" % comment)
        np.savetxt(f, verts, fmt="v %.6f %.6f %.6f")
        if uv is not None:
            np.savetxt(f, uv, fmt="vt %.6f %.6f")
        np.savetxt(f, normals, fmt="vn %.5f %.5f %.5f")
        i = faces + 1
        if uv is not None:
            np.savetxt(f, np.stack([i[:, 0]] * 3 + [i[:, 1]] * 3 + [i[:, 2]] * 3, axis=1), fmt="f %d/%d/%d %d/%d/%d %d/%d/%d")
        else:
            np.savetxt(f, np.stack([i[:, 0], i[:, 0], i[:, 1], i[:, 1], i[:, 2], i[:, 2]], axis=1), fmt="f %d//%d %d//%d %d//%d")


def write_ply(path, verts, normals, faces, uv, comment):
    """Binary little endian: x y z nx ny nz [s t] as float, faces as list uchar int."""
    props = ["x", "y", "z", "nx", "ny", "nz"] + (["s", "t"] if uv is not None else [])
    header = "ply\nformat binary_little_endian 1.0\ncomment %s\nelement vertex %d\n%s" % (comment, len(verts), "".join("property float %s\n" % p for p in props))
    header += "element face %d\nproperty list uchar int vertex_indices\nend_header\n" % len(faces)
    vdata = np.concatenate([verts, normals] + ([uv] if uv is not None else []), axis=1).astype("<f4")
    fdt = np.dtype([("n", "u1"), ("i", "<i4", (3,))])
    fdata = np.empty(len(faces), dtype=fdt)
    fdata["n"] = 3
    fdata["i"] = faces
    with open(path, "wb") as f:
        f.write(header.encode("ascii"))
        f.write(vdata.tobytes())
        f.write(fdata.tobytes())


# ----------------------------------------------------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="OBJ or PLY")
    ap.add_argument("output", help="OBJ or PLY (by the extension; PLY is binary and about 3x smaller)")
    ap.add_argument("--keep-uv", action="store_true", help="keep the texture coordinates (for a material with a \"texture\")")
    ap.add_argument("--up", choices=sorted(_UP), default="y", help="the input's up axis (COLMAP / 3DGS: write --up=-y); the output has y up")
    ap.add_argument("--flip", action="store_true", help="invert the winding (use when the signed volume comes out negative)")
    ap.add_argument("--drop-small", type=float, metavar="FRACTION", help="drop disconnected pieces with fewer than FRACTION of the triangles (floaters), e.g. 0.01")
    ap.add_argument("--height", type=float, metavar="METRES", help="scale the mesh so that its height (y extent) is this")
    ap.add_argument("--scale", type=float, help="multiply the coordinates by this (instead of --height)")
    ap.add_argument("--ground", action="store_true", help="put the lowest point at y = 0 and the center of the footprint at x = z = 0")
    ap.add_argument("--cap", action="store_true", help="close the open boundaries with a triangle fan")
    args = ap.parse_args()

    verts, faces, corner_uv = read_mesh(args.input)
    print("read: %d vertices, %d triangles%s" % (len(verts), len(faces), ", uv" if corner_uv is not None else ""))
    if not args.keep_uv:
        corner_uv = None
    verts = verts @ _UP[args.up].T
    if args.flip:
        faces = faces[:, [0, 2, 1]]
    verts, faces, uv = weld(verts, faces, corner_uv)
    print("welded: %d vertices, %d triangles" % (len(verts), len(faces)))
    if args.drop_small:
        verts, faces, uv, report = drop_small_components(verts, faces, uv, args.drop_small)
        print("pieces: %s -> %d vertices, %d triangles" % (report, len(verts), len(faces)))
    loops = boundary_loops(faces)
    print("open boundaries: %d loops %s" % (len(loops), sorted(map(len, loops), reverse=True)[:8]))
    if args.cap and loops:
        verts, faces, uv = cap(verts, faces, loops, uv)
        print("closed: %d vertices, %d triangles" % (len(verts), len(faces)))
    lo, hi = verts.min(axis=0), verts.max(axis=0)
    if args.height:
        s = args.height / (hi[1] - lo[1])
        verts = verts * s
        print("scaled by %.6g to a height of %g m" % (s, args.height))
    elif args.scale:
        verts = verts * args.scale
    if args.ground:
        lo, hi = verts.min(axis=0), verts.max(axis=0)
        verts = verts - np.array([(lo[0] + hi[0]) / 2, lo[1], (lo[2] + hi[2]) / 2])
    lo, hi = verts.min(axis=0), verts.max(axis=0)
    print("bounds: x %.4f..%.4f  y %.4f..%.4f  z %.4f..%.4f" % (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]))
    volume = signed_volume(verts, faces)
    print("signed volume %.6g (positive: outward normals)" % volume)
    if volume < 0:
        print("warning: negative volume - the winding is inside out (try --flip)", file=sys.stderr)
    normals = smooth_normals(verts, faces)
    steps = ["welded", "smooth normals recomputed"]
    if uv is not None:
        steps.append("uv kept")
    if args.up != "y":
        steps.append("up axis %s -> y" % args.up)
    if args.drop_small:
        steps.append("small pieces dropped")
    if args.height or args.scale:
        steps.append("scaled")
    if args.ground:
        steps.append("grounded")
    if args.cap:
        steps.append("boundaries closed")
    comment = "cleaned by tools/data/mesh_clean.py: " + ", ".join(steps)
    if os.path.splitext(args.output)[1].lower() == ".ply":
        write_ply(args.output, verts, normals, faces, uv, comment)
    else:
        write_obj(args.output, verts, normals, faces, uv, comment)
    print("wrote", args.output)


if __name__ == "__main__":
    main()
