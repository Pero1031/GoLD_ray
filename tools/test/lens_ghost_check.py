"""Ghost rays of the real lenses, traced by a program that is not the renderer: where a ray that the lens surfaces reflect twice meets
the film, and how much of it is left (Fresnel or a quarter-wave MgF2 coating at each crossing). The renderer's check (--scene camera,
check s) traces the same rays with LensSystem::traceGhost and compares them with the numbers this script prints.

Independent of the renderer's code: the lens files are read here, every surface (spheres too) is met by Newton's method from its vertex
plane in a frame whose z points toward the image, and the reflectances come from Airy's formula for a film on a substrate written out
here. The film lies at the paraxial focus (infinity), found with a ray a nanometre off the axis. Wavelength: the d line (587.56 nm), where
the renderer's glass has its n_d in both builds.

    python tools/test/lens_ghost_check.py
"""
import cmath
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
LENSES = os.path.join(os.path.dirname(os.path.dirname(HERE)), "GoLD_rayt", "assets", "Lenses")
NM = 587.56
MGF2, DESIGN_NM = 1.38, 550.0


def read_lens(name):
    """Surfaces front to back: dict(c, k, poly{i: A_i}, z (mm, toward the image), n (after), stop)."""
    surfs, z = [], 0.0
    for line in open(os.path.join(LENSES, name), encoding='utf-8'):
        tok = line.split('#')[0].split()
        if not tok or tok[0] in ('fnumber', 'image-height', 'half-field'):
            continue
        if tok[0] == 'asphere':
            kv = dict(zip(tok[1::2], map(float, tok[2::2])))
            surfs[-1]['k'] = kv.pop('k', 0.0)
            surfs[-1]['poly'] = {int(key[1:]): val for key, val in kv.items()}
            continue
        r = math.inf if tok[0] in ('inf', 'infinity') else float(tok[0])
        n = float(tok[2]) if float(tok[2]) > 1 else 1.0
        surfs.append(dict(c=0.0 if (r == 0 or math.isinf(r)) else 1.0 / r, k=0.0, poly={}, z=z, n=n, stop=(r == 0)))
        z += float(tok[1])
    return surfs


def sag(s, r):
    arg = 1 - (1 + s['k']) * s['c'] ** 2 * r * r
    q = math.sqrt(arg)
    z = s['c'] * r * r / (1 + q)
    dz = s['c'] * r / q
    for i, a in s['poly'].items():
        z += a * r ** i
        dz += i * a * r ** (i - 1)
    return z, dz


def meet(s, o, d, behind_ok=False):
    """Newton's method on the ray parameter from the vertex plane: the point where the ray meets the surface z = z_s + sag(r). It must
    lie ahead of the ray, unless the surface sits behind the previous one along the ray's direction (a negative spacing)."""
    t = (s['z'] - o[2]) / d[2]
    for _ in range(100):
        p = [o[i] + t * d[i] for i in range(3)]
        r = math.hypot(p[0], p[1])
        if s['stop']:
            if t <= 1e-12 and not behind_ok:
                raise ValueError('the stop lies behind the ray')
            return p, None
        z, dz = sag(s, r)
        g = p[2] - s['z'] - z
        dg = d[2] - (dz * (p[0] * d[0] + p[1] * d[1]) / r if r > 0 else 0)
        step = g / dg
        t -= step
        if abs(step) < 1e-13:
            break
    else:
        raise ValueError('Newton did not converge (a ray grazing the surface)')
    if t <= 1e-12 and not behind_ok:
        raise ValueError('the surface lies behind the ray')
    p = [o[i] + t * d[i] for i in range(3)]
    r = math.hypot(p[0], p[1])
    _, dz = sag(s, r)
    nrm = [-dz * p[0] / r, -dz * p[1] / r, 1.0] if r > 0 else [0.0, 0.0, 1.0]
    ln = math.sqrt(sum(x * x for x in nrm))
    return p, [x / ln for x in nrm]


def reflectance(cos_i, n_in, n_out, coated):
    """Unpolarized reflectance of ambient | film | substrate (Airy's formula); no film (uncoated or glass-glass): Fresnel."""
    nf = MGF2
    thick = DESIGN_NM / (4 * MGF2) if coated and ((n_in > 1) != (n_out > 1)) else 0.0
    sin_i = math.sqrt(max(0.0, 1 - cos_i * cos_i))
    n0, n1, n2 = complex(n_in), complex(nf), complex(n_out)
    c0 = complex(cos_i)
    c1 = cmath.sqrt(1 - (n0 * sin_i / n1) ** 2)
    c2 = cmath.sqrt(1 - (n0 * sin_i / n2) ** 2)
    if c2.imag < 0:
        c2 = -c2
    r01s = (n0 * c0 - n1 * c1) / (n0 * c0 + n1 * c1)
    r01p = (n1 * c0 - n0 * c1) / (n1 * c0 + n0 * c1)
    r12s = (n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2)
    r12p = (n2 * c1 - n1 * c2) / (n2 * c1 + n1 * c2)
    e2 = cmath.exp(2j * 2 * math.pi * thick / NM * n1 * c1)
    rs = (r01s + r12s * e2) / (1 + r01s * r12s * e2)
    rp = (r01p + r12p * e2) / (1 + r01p * r12p * e2)
    return min(1.0, 0.5 * (abs(rs) ** 2 + abs(rp) ** 2))


def ghost(surfs, i, j, o, d, coated, film_z):
    """Crosses 0..j, reflects at j, back to i, reflects there, on to the film. Returns (x, y, throughput, stop passes)."""
    n_before = lambda k: 1.0 if k == 0 else surfs[k - 1]['n']
    seq = [(k, k == j) for k in range(0, j + 1)] + [(k, k == i) for k in range(j - 1, i - 1, -1)] + [(k, False) for k in range(i + 1, len(surfs))]
    o, d, thr, passes, prev = list(o), list(d), 1.0, 0, None
    for k, refl in seq:
        s = surfs[k]
        behind = prev is not None and ((s['z'] < surfs[prev]['z']) if d[2] > 0 else (s['z'] > surfs[prev]['z']))
        prev = k
        p, nrm = meet(s, o, d, behind)
        o = p
        if s['stop']:
            passes += 1
            continue
        if sum(nrm[q] * d[q] for q in range(3)) > 0:
            nrm = [-x for x in nrm]
        cos_i = -sum(nrm[q] * d[q] for q in range(3))
        forward = d[2] > 0
        n_in, n_out = (n_before(k), s['n']) if forward else (s['n'], n_before(k))
        rr = reflectance(cos_i, n_in, n_out, coated)
        if refl:
            d = [d[q] + 2 * cos_i * nrm[q] for q in range(3)]
            thr *= rr
        else:
            eta = n_in / n_out
            s2 = eta * eta * (1 - cos_i * cos_i)
            if s2 >= 1:
                return None
            d = [eta * d[q] + (eta * cos_i - math.sqrt(1 - s2)) * nrm[q] for q in range(3)]
            thr *= 1 - rr
        ln = math.sqrt(sum(x * x for x in d))
        d = [x / ln for x in d]
    t = (film_z - o[2]) / d[2]
    return o[0] + t * d[0], o[1] + t * d[1], thr, passes


def paraxial_focus(surfs):
    """z of the paraxial focus: a ray parallel to the axis a nanometre off it."""
    o, d = [1e-6, 0.0, -10.0], [0.0, 0.0, 1.0]
    for k, s in enumerate(surfs):
        p, nrm = meet(s, o, d, k > 0 and s['z'] < surfs[k - 1]['z'])
        o = p
        if s['stop']:
            continue
        if sum(nrm[q] * d[q] for q in range(3)) > 0:
            nrm = [-x for x in nrm]
        cos_i = -sum(nrm[q] * d[q] for q in range(3))
        n_in, n_out = (1.0 if k == 0 else surfs[k - 1]['n']), s['n']
        eta = n_in / n_out
        d = [eta * d[q] + (eta * cos_i - math.sqrt(1 - eta * eta * (1 - cos_i * cos_i))) * nrm[q] for q in range(3)]
    return o[2] - o[0] / d[0] * d[2]


# (lens, i, j, field angle deg, x, y of the ray where it crosses the plane 1 mm in front of the lens's foremost vertex (mm), coated)
CASES = [
    ('dgauss_50mm.txt', 0, 1, 0.0, 3.0, 0.0, False), ('dgauss_50mm.txt', 2, 7, 5.0, -4.0, 2.0, True),
    ('dgauss_50mm.txt', 4, 9, 10.0, 2.0, -3.0, True), ('dgauss_50mm.txt', 4, 9, 4.0, 1.0, -2.0, True),
    ('dgauss_50mm.txt', 7, 10, 3.0, 0.0, 6.0, False),
    ('cellphone_5.6mm.txt', 1, 2, 0.0, 0.5, 0.0, True), ('cellphone_5.6mm.txt', 2, 4, 5.0, 0.3, -0.4, True),
    ('cellphone_5.6mm.txt', 5, 8, 5.0, -0.1, 0.2, False),
    ('cellphone_5.6mm.txt', 9, 10, 10.0, 0.7, 0.1, True),
    ('cellphone_5.6mm.txt', 7, 10, 20.0, -0.2, 0.6, False), ('phone_tele_12mm.txt', 1, 4, 3.0, 1.0, 0.5, True),
    ('phone_tele_12mm.txt', 6, 11, 8.0, -0.5, -1.0, True),
]


def main():
    print('ghost rays at %g nm: film x, y (mm, the renderer frame), throughput, stop crossings' % NM)
    rows = []
    for name, i, j, deg, x, y, coated in CASES:
        surfs = read_lens(name)
        z_top = min(s['z'] for s in surfs)   # the foremost vertex (z toward the image)
        film = paraxial_focus(surfs)
        th = math.radians(deg)
        d = [math.sin(th), 0.0, math.cos(th)]
        o = [x - 50 * d[0], y - 50 * d[1], z_top - 1.0 - 50 * d[2]]   # 50 mm back along the ray from (x, y) on the plane
        try:
            res = ghost(surfs, i, j, o, d, coated, film)
        except ValueError:   # the ray left an asphere's domain
            res = None
        print('  %-20s %2d %2d %5.1f deg (%5.2f, %5.2f) %-8s -> %s' % (name, i, j, deg, x, y, 'coated' if coated else 'bare', res))
        if res:
            rows.append((name, i, j, deg, x, y, coated) + res)
    print('C++ table for CameraChecks:')
    for r in rows:
        print('    { "%s", %d, %d, %g, %g, %g, %s, %.12f, %.12f, %.12e, %d },' % (r[0], r[1], r[2], r[3], r[4], r[5], 'true' if r[6] else 'false',
                                                                            r[7], r[8], r[9], r[10]))


if __name__ == '__main__':
    main()
