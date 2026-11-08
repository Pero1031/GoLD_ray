"""Procedural meshes for the diamond ring scene: a round brilliant-cut diamond (flat facets) and a ring band with four prongs (smooth).

    python tools/data/make_ring.py [out_dir]          default: GoLD_rayt/assets/Meshes/ring

Writes   diamond.obj   57-facet-style round brilliant, girdle diameter 1, culet at y = -0.44 (flat facets: no normals in the file, the
                       loader flat-shades), closed, outward faces;
         ring.obj      band (a torus with a rounded rectangle section, tapered towards the back) and four prongs, with smooth normals.
Units: the diamond has a girdle diameter of 1 and its table up at y = +0.14; the ring is in metres for a 6.5 mm stone (girdle 6.5 mm) and is
placed so that the culet of the stone (scaled by 0.0065, moved up by the offset printed at the end) sits in the prongs.

Brilliant geometry (proportions of the "ideal" cut, in units of the girdle diameter): table 56 %, crown angle 34.5 degrees, pavilion angle 40.75
degrees, a thin girdle. Facets: table (8-gon), 8 star triangles on the table edges, 8 bezel kites with a corner at each table vertex, 16 upper-girdle
triangles, a thin girdle band, 16 lower-girdle triangles, 8 pavilion main kites meeting at the culet. The kites are planar by construction: the free
vertices (the star tips and the lower-girdle tips) are placed on the intersection of the two adjacent facet planes.
"""
import math
import os
import sys

import numpy as np


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def diamond():
    R = 0.5
    t = 0.008                                     # half the girdle thickness
    beta = math.radians(34.5)                     # crown angle
    alpha = math.radians(40.75)                   # pavilion angle
    rt = 0.28 / math.cos(math.radians(22.5))      # table: 56 % across the flats
    N = 8
    ang = lambda k: math.radians(45.0 * k)        # azimuth of table vertices and girdle "main" points

    # crown: plane k through the main girdle point G_k = (R, +t) with the crown angle; the table corner T_k lies in it
    hc = (R - rt) * math.tan(beta) + t
    T = [np.array([rt * math.cos(ang(k)), hc, rt * math.sin(ang(k))]) for k in range(N)]
    Gt_main = [np.array([R * math.cos(ang(k)), t, R * math.sin(ang(k))]) for k in range(N)]
    Gb_main = [np.array([R * math.cos(ang(k)), -t, R * math.sin(ang(k))]) for k in range(N)]
    half = lambda k: ang(k) + math.radians(22.5)
    Gt_half = [np.array([R * math.cos(half(k)), t, R * math.sin(half(k))]) for k in range(N)]
    Gb_half = [np.array([R * math.cos(half(k)), -t, R * math.sin(half(k))]) for k in range(N)]

    def crown_normal(k):
        phi = ang(k)
        return np.array([math.cos(phi) * math.sin(beta), math.cos(beta), math.sin(phi) * math.sin(beta)])

    # star tips S_k: on both bezel planes k and k+1, at azimuth half(k), radius rs
    rs = 0.405
    S = []
    for k in range(N):
        n = crown_normal(k)
        pos = np.array([rs * math.cos(half(k)), 0.0, rs * math.sin(half(k))])
        y = (np.dot(n, Gt_main[k]) - n[0] * pos[0] - n[2] * pos[2]) / n[1]
        S.append(np.array([pos[0], y, pos[2]]))

    # pavilion: main plane k through the culet C and G_k (below), containing the horizontal tangent at G_k
    pd = R * math.tan(alpha) + t
    C = np.array([0.0, -pd, 0.0])
    rl = 0.30
    Lm = []
    for k in range(N):
        phi = ang(k)
        tang = np.array([-math.sin(phi), 0, math.cos(phi)])
        n = unit(np.cross(Gb_main[k] - C, tang))
        if n[1] > 0:
            n = -n
        pos = np.array([rl * math.cos(half(k)), 0.0, rl * math.sin(half(k))])
        y = (np.dot(n, C) - n[0] * pos[0] - n[2] * pos[2]) / n[1]
        Lm.append(np.array([pos[0], y, pos[2]]))

    polys = []
    polys.append([T[k] for k in range(N)])                                                         # table (8-gon)
    for k in range(N):
        polys.append([T[k], T[(k + 1) % N], S[k]])                                                 # star (on the table edge)
        polys.append([T[k], S[(k - 1) % N], Gt_main[k], S[k]])                                     # bezel kite
        polys.append([S[k], Gt_main[k], Gt_half[k]])                                               # upper girdle
        polys.append([S[k], Gt_half[k], Gt_main[(k + 1) % N]])
    for k in range(N):                                                                             # girdle band (16 quads)
        polys.append([Gt_main[k], Gb_main[k], Gb_half[k], Gt_half[k]])
        polys.append([Gt_half[k], Gb_half[k], Gb_main[(k + 1) % N], Gt_main[(k + 1) % N]])
    for k in range(N):
        polys.append([Lm[k], Gb_main[k], Gb_half[k]])                                              # lower girdle
        polys.append([Lm[k], Gb_half[k], Gb_main[(k + 1) % N]])
        polys.append([C, Lm[(k - 1) % N], Gb_main[k], Lm[k]])                                      # pavilion main kite
    return polys


def orient_outward(polys, centre):
    out = []
    for p in polys:
        p = [np.asarray(v, float) for v in p]
        n = np.cross(p[1] - p[0], p[2] - p[0])
        cen = sum(p) / len(p)
        if np.dot(n, cen - centre) < 0:
            p = p[::-1]
        out.append(p)
    return out


def write_flat_obj(path, polys, comment):
    verts, faces = [], []
    for p in polys:
        idx = []
        for v in p:
            verts.append(v)
            idx.append(len(verts))
        for i in range(1, len(idx) - 1):               # fan triangulation (the facets are convex)
            faces.append((idx[0], idx[i], idx[i + 1]))
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# {comment}\n")
        for v in verts:
            f.write("v %.6f %.6f %.6f\n" % tuple(v))
        for a, b, c in faces:
            f.write(f"f {a} {b} {c}\n")
    return len(verts), len(faces)


def mesh_volume(polys):
    v = 0.0
    for p in polys:
        for i in range(1, len(p) - 1):
            v += np.dot(p[0], np.cross(p[i], p[i + 1])) / 6.0
    return v


def ring():
    """Band: a torus around the z axis (the finger axis), ring plane xy. Section: a rounded rectangle, thicker at the front (+y, where the stone is)."""
    Rc = 0.0089                       # centre radius of the band (inner diameter about 17 mm)
    nu, nv = 128, 22
    verts, norms, faces = [], [], []

    def section(u, s):
        """Point of the section at angle s around it: (radial offset, z offset) for the band at azimuth u (0 at the front)."""
        thick = 0.0013 + 0.0006 * math.cos(u) ** 2            # thicker at the front
        width = 0.0036 + 0.0006 * math.cos(u) ** 2
        # superellipse for a rounded rectangle
        e = 0.35
        c, sn = math.cos(s), math.sin(s)
        x = math.copysign(abs(c) ** e, c) * thick
        z = math.copysign(abs(sn) ** e, sn) * width / 2
        return x, z

    for i in range(nu):
        u = 2 * math.pi * i / nu
        for j in range(nv):
            s = 2 * math.pi * j / nv
            r, z = section(u, s)
            rad = Rc + r
            verts.append(np.array([rad * math.sin(u), rad * math.cos(u), z]))
    # smooth normals from the surface by finite differences
    P = lambda i, j: verts[(i % nu) * nv + (j % nv)]
    for i in range(nu):
        for j in range(nv):
            du = P(i + 1, j) - P(i - 1, j)
            dv = P(i, j + 1) - P(i, j - 1)
            n = unit(np.cross(du, dv))
            cen = np.array([0.0, 0.0, 0.0])
            norms.append(n)
    # orient normals outward: compare with the direction from the section centre
    for i in range(nu):
        u = 2 * math.pi * i / nu
        centre = np.array([Rc * math.sin(u), Rc * math.cos(u), 0.0])
        for j in range(nv):
            k = i * nv + j
            if np.dot(norms[k], verts[k] - centre) < 0:
                norms[k] = -norms[k]
    for i in range(nu):
        for j in range(nv):
            a, b = i * nv + j, ((i + 1) % nu) * nv + j
            c, d = ((i + 1) % nu) * nv + (j + 1) % nv, i * nv + (j + 1) % nv
            faces.append((a, b, c)); faces.append((a, c, d))

    # four prongs: slim tapered cylinders from the top of the band to the girdle of the stone
    stone_r = 0.00325
    base_y = Rc + 0.0020
    prongs = []
    for k in range(4):
        phi = math.radians(45 + 90 * k)
        p0 = np.array([0.0030 * math.cos(phi), base_y - 0.0004, 0.0030 * math.sin(phi)])
        p1 = np.array([(stone_r + 0.00045) * math.cos(phi), base_y + 0.0033, (stone_r + 0.00045) * math.sin(phi)])
        prongs.append((p0, p1))
    pv, pn, pf = [], [], []
    seg = 14
    for p0, p1 in prongs:
        axis = unit(p1 - p0)
        a = unit(np.cross(axis, [0, 0, 1]) if abs(axis[2]) < 0.9 else np.cross(axis, [1, 0, 0]))
        b = np.cross(axis, a)
        base = len(pv)
        for ring_i, (c, rad) in enumerate(((p0, 0.00042), (p1, 0.00030))):
            for s in range(seg):
                th = 2 * math.pi * s / seg
                d = a * math.cos(th) + b * math.sin(th)
                pv.append(c + d * rad)
                pn.append(d)
        for s in range(seg):
            q0, q1 = base + s, base + (s + 1) % seg
            q2, q3 = base + seg + (s + 1) % seg, base + seg + s
            pf.append((q0, q1, q2)); pf.append((q0, q2, q3))
        # caps
        for c, off, flip in ((p0, 0, True), (p1, seg, False)):
            ci = len(pv); pv.append(c); pn.append(-axis if flip else axis)
            for s in range(seg):
                q0, q1 = base + off + s, base + off + (s + 1) % seg
                pf.append((ci, q1, q0) if flip else (ci, q0, q1))
    return verts, norms, faces, pv, pn, pf, base_y


def write_smooth_obj(path, parts, comment):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# {comment}\n")
        off = 0
        for verts, norms, faces in parts:
            for v in verts:
                f.write("v %.6f %.6f %.6f\n" % tuple(v))
            for n in norms:
                f.write("vn %.5f %.5f %.5f\n" % tuple(n))
            for a, b, c in faces:
                f.write(f"f {a + 1 + off}//{a + 1 + off} {b + 1 + off}//{b + 1 + off} {c + 1 + off}//{c + 1 + off}\n")
            off += len(verts)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "GoLD_rayt", "assets", "Meshes", "ring")
    os.makedirs(out, exist_ok=True)
    polys = orient_outward(diamond(), np.array([0.0, -0.1, 0.0]))
    nv, nf = write_flat_obj(os.path.join(out, "diamond.obj"), polys, "round brilliant, girdle diameter 1 (tools/data/make_ring.py), flat facets")
    print(f"diamond.obj: {nv} vertices, {nf} triangles, volume {mesh_volume(polys):.4f} (a sphere of diameter 1 has 0.5236)")
    bv, bn, bf, pv, pn, pf, base_y = ring()
    write_smooth_obj(os.path.join(out, "ring.obj"), [(bv, bn, bf), (pv, pn, pf)], "ring band and four prongs, in metres, finger axis z (tools/data/make_ring.py)")
    print(f"ring.obj: band {len(bv)} vertices, prongs {len(pv)} vertices; prong base at y = {base_y:.4f}")
    print("stone placement for a 6.5 mm stone (scale 0.0065): the culet is at y = -0.44 * 0.0065 = -0.00286 in stone coordinates; translate the stone up by",
          "%.5f" % (base_y + 0.0004 + 0.0029 - 0.0009))


if __name__ == "__main__":
    main()
