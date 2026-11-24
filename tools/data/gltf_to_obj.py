"""
Writes the meshes of a glTF file (.gltf with embedded or external buffers, or .glb) as OBJ files the renderer reads
(positions, uv, triangles), one per mesh, and saves the images the meshes' materials use.

    python tools/data/gltf_to_obj.py scene.gltf out_folder [--world] [--list]

--list      only print the meshes, their materials and sizes
--world     apply the node transforms (translation, rotation, scale) of the default scene; otherwise each mesh stays in its own
            coordinates

glTF's uv origin is the top-left of the image; OBJ's is the bottom-left, so v is flipped.
"""
import argparse
import base64
import json
import os
import struct
import sys

import numpy as np

COMPONENT = {5120: np.int8, 5121: np.uint8, 5122: np.int16, 5123: np.uint16, 5125: np.uint32, 5126: np.float32}
COUNT = {"SCALAR": 1, "VEC2": 2, "VEC3": 3, "VEC4": 4, "MAT4": 16}


def load(path):
    if path.lower().endswith(".glb"):
        data = open(path, "rb").read()
        length = struct.unpack_from("<I", data, 8)[0]
        off, doc, binary = 12, None, None
        while off < length:
            clen, ctype = struct.unpack_from("<II", data, off)
            chunk = data[off + 8: off + 8 + clen]
            if ctype == 0x4E4F534A:
                doc = json.loads(chunk)
            elif ctype == 0x004E4942:
                binary = chunk
            off += 8 + clen
        buffers = [binary]
    else:
        doc = json.load(open(path, encoding="utf-8"))
        buffers = []
        for b in doc["buffers"]:
            uri = b["uri"]
            if uri.startswith("data:"):
                buffers.append(base64.b64decode(uri.split(",", 1)[1]))
            else:
                buffers.append(open(os.path.join(os.path.dirname(path), uri), "rb").read())
    return doc, buffers


def view_bytes(doc, buffers, i):
    v = doc["bufferViews"][i]
    start = v.get("byteOffset", 0)
    return buffers[v["buffer"]][start: start + v["byteLength"]], v.get("byteStride")


def accessor(doc, buffers, i):
    a = doc["accessors"][i]
    raw, stride = view_bytes(doc, buffers, a["bufferView"])
    dt = np.dtype(COMPONENT[a["componentType"]])
    n = COUNT[a["type"]]
    off = a.get("byteOffset", 0)
    if stride and stride != dt.itemsize * n:
        rows = [np.frombuffer(raw, dt, n, off + k * stride) for k in range(a["count"])]
        return np.array(rows)
    return np.frombuffer(raw, dt, a["count"] * n, off).reshape(a["count"], n) if n > 1 else np.frombuffer(raw, dt, a["count"], off)


def node_matrix(node):
    if "matrix" in node:
        return np.array(node["matrix"], dtype=np.float64).reshape(4, 4).T
    t = np.array(node.get("translation", [0, 0, 0]), dtype=np.float64)
    x, y, z, w = node.get("rotation", [0, 0, 0, 1])
    r = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                  [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                  [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
    s = np.array(node.get("scale", [1, 1, 1]), dtype=np.float64)
    m = np.eye(4)
    m[:3, :3] = r * s[None, :]
    m[:3, 3] = t
    return m


def mesh_transforms(doc):
    """mesh index -> world matrix of the first node of the default scene that draws it."""
    out = {}

    def walk(i, parent):
        node = doc["nodes"][i]
        m = parent @ node_matrix(node)
        if "mesh" in node and node["mesh"] not in out:
            out[node["mesh"]] = m
        for c in node.get("children", []):
            walk(c, m)

    scene = doc["scenes"][doc.get("scene", 0)]
    for i in scene["nodes"]:
        walk(i, np.eye(4))
    return out


def safe(name):
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gltf")
    ap.add_argument("out", nargs="?")
    ap.add_argument("--world", action="store_true")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    doc, buffers = load(args.gltf)
    xf = mesh_transforms(doc)
    if not args.list:
        if not args.out:
            sys.exit("give an output folder")
        os.makedirs(args.out, exist_ok=True)
    for mi, mesh in enumerate(doc["meshes"]):
        name = safe(mesh.get("name", f"mesh{mi}"))
        for pi, prim in enumerate(mesh["primitives"]):
            pos = accessor(doc, buffers, prim["attributes"]["POSITION"]).astype(np.float64)
            uv = accessor(doc, buffers, prim["attributes"]["TEXCOORD_0"]).astype(np.float64) if "TEXCOORD_0" in prim["attributes"] else None
            idx = accessor(doc, buffers, prim["indices"]).astype(np.int64).reshape(-1, 3) if "indices" in prim else np.arange(len(pos)).reshape(-1, 3)
            if args.world and mi in xf:
                m = xf[mi]
                pos = pos @ m[:3, :3].T + m[:3, 3]
            mat = doc["materials"][prim["material"]] if "material" in prim else {}
            tex = mat.get("pbrMetallicRoughness", {}).get("baseColorTexture", {}).get("index")
            img_name = None
            if tex is not None:
                img = doc["images"][doc["textures"][tex]["source"]]
                ext = ".png" if img.get("mimeType") == "image/png" else ".jpg"
                img_name = safe(img.get("name", f"image{tex}")) + ext
            lo, hi = pos.min(0), pos.max(0)
            label = name if len(mesh["primitives"]) == 1 else f"{name}_{pi}"
            print(f"{label}: {len(pos)} vertices, {len(idx)} triangles, size {np.round(hi - lo, 4)}, material {mat.get('name')}, image {img_name}")
            if args.list:
                continue
            path = os.path.join(args.out, label + ".obj")
            with open(path, "w") as f:
                f.write(f"# {label} from {os.path.basename(args.gltf)}" + (" (world coordinates)" if args.world else "") + "\n")
                np.savetxt(f, pos, fmt="v %.6f %.6f %.6f")
                if uv is not None:
                    np.savetxt(f, np.column_stack([uv[:, 0], 1 - uv[:, 1]]), fmt="vt %.6f %.6f")
                    np.savetxt(f, np.repeat(idx + 1, 2, axis=1), fmt="f %d/%d %d/%d %d/%d")
                else:
                    np.savetxt(f, idx + 1, fmt="f %d %d %d")
            if img_name:
                img = doc["images"][doc["textures"][tex]["source"]]
                if "bufferView" in img:
                    raw, _ = view_bytes(doc, buffers, img["bufferView"])
                    open(os.path.join(args.out, img_name), "wb").write(raw)


if __name__ == "__main__":
    main()
