"""The aspheric phone lens of assets/Lenses/cellphone_5.6mm.txt (US patent 7,535,658, Example 1) in rayoptics, an optical-design library
that is not the renderer: its paraxial focal length and back focus, and where real rays land on the image plane. The renderer's check
(--scene camera, check p) compares its own trace of the same rays with the numbers this script prints.

Needs rayoptics (BSD). Its newest versions need Python 3.12; on Python 3.10 put rayoptics 0.9.5 in a separate environment that shares the
installed numpy and scipy, without its Qt user interface:
    python -m venv --system-site-packages D:\\venvs\\rayoptics
    D:\\venvs\\rayoptics\\Scripts\\python -m pip install --no-deps rayoptics==0.9.5 opticalglass==1.1.0 anytree transforms3d deprecation parsimonious regex xlrd
    D:\\venvs\\rayoptics\\Scripts\\python tools/test/lens_rayoptics_check.py

The prescription is built here from the patent's own numbers (its conic K = 1 + k, as rayoptics' "ec"; the coefficients A3 ... A10 as
rayoptics' documentation transcribes this patent), not from the renderer's file; the script then checks that the file says the same.
Rays: from the object side at the field angle theta in the x-z plane, through a point of the stop (the first surface, 0.20 mm behind the
first lens's vertex), at the d line (587.56 nm); the image plane lies 0.64 mm behind the cover glass, as in the patent.
"""
import math
import os
import sys

import numpy as np

from rayoptics.environment import OpticalModel, PupilSpec, FieldSpec, WvlSpec, RadialPolynomial
from rayoptics.raytr.trace import trace

HERE = os.path.dirname(os.path.abspath(__file__))
LENS_FILE = os.path.join(os.path.dirname(os.path.dirname(HERE)), "GoLD_rayt", "assets", "Lenses", "cellphone_5.6mm.txt")

# (radius, thickness, n_d, V_d, K, [A3 ... A10]); n_d = 1: air
SURFACES = [
    (1.962, 1.19, 1.471, 76.6, 2.153, [-1.895e-2, 2.426e-2, -5.123e-2, 8.371e-4, 7.850e-3, 4.091e-3, -7.732e-3, -4.265e-3]),
    (33.398, 0.93, 1, 0, 40.18, [-4.966e-3, -1.434e-2, -6.139e-3, -9.284e-5, 6.438e-3, -5.72e-3, -2.385e-2, 1.108e-2]),
    (-2.182, 0.75, 1.603, 27.5, 2.105, [-4.388e-2, -2.555e-2, 5.16e-2, -4.307e-2, -2.831e-2, 3.162e-2, 4.630e-2, -4.877e-2]),
    (-6.367, 0.10, 1, 0, 3.382, [-1.131e-1, -7.863e-2, 1.094e-1, 6.228e-3, -2.216e-2, -5.89e-3, 4.123e-3, 1.041e-3]),
    (5.794, 0.89, 1.510, 56.2, -221.1, [-7.876e-2, 7.02e-2, 5.175e-3, -9.958e-3, -7.322e-3, 6.914e-4, 2.54e-3, -7.65e-4]),
    (9.192, 0.16, 1, 0, 0.9331, [9.694e-3, -2.516e-3, -3.606e-3, -2.497e-4, -6.84e-4, -1.414e-4, 2.932e-4, -7.284e-5]),
    (1.674, 0.85, 1.510, 56.2, -7.617, [7.429e-2, -6.933e-2, -5.811e-3, 2.396e-3, 2.100e-3, -3.119e-4, -5.552e-5, 7.969e-6]),
    (1.509, 0.70, 1, 0, -2.707, [1.767e-3, -4.652e-2, 1.625e-2, -3.522e-3, -7.106e-4, 3.825e-4, 6.271e-5, -2.631e-5]),
]
STOP_GAP = -0.20          # stop to the first lens surface
COVER = (0.30, 1.516, 64.1)
IMAGE_GAP = 0.64          # cover glass to the image plane
OBJECT_GAP = 10.0         # where the rays start, in front of the stop (any distance: the ray is straight there)
FIELDS_DEG = [0, 10, 20, 30]
STOP_POINTS = [(0, 0), (0.8, 0), (-0.8, 0), (0, 0.8), (0.5, -0.5)]   # mm in the stop
LINE_D = 587.56


def build():
    opm = OpticalModel()
    sm, osp = opm['seq_model'], opm['optical_spec']
    osp['pupil'] = PupilSpec(osp, key=['image', 'f/#'], value=2.8)
    osp['fov'] = FieldSpec(osp, key=['object', 'angle'], value=30, flds=[0., 1.], is_relative=True)
    osp['wvls'] = WvlSpec([(LINE_D, 1.0)], ref_wl=0)
    opm.radius_mode = True
    sm.gaps[0].thi = OBJECT_GAP
    sm.add_surface([0., STOP_GAP])
    sm.set_stop()
    for r, t, n, v, K, a in SURFACES:
        sm.add_surface([r, t, n, v] if n > 1 else [r, t])
        sm.ifcs[sm.cur_surface].profile = RadialPolynomial(r=r, ec=K, coefs=[0.0, 0.0] + a)
    sm.add_surface([0., COVER[0], COVER[1], COVER[2]])
    sm.add_surface([0., IMAGE_GAP])
    opm.update_model()
    return opm


def check_file():
    """The renderer's file against the table above: radii, spacings, glasses, k = K - 1 and the coefficients."""
    rows, asph = [], []
    for line in open(LENS_FILE, encoding='utf-8'):
        line = line.split('#')[0].split()
        if not line or line[0] in ('fnumber', 'image-height'):
            continue
        if line[0] == 'asphere':
            kv = dict(zip(line[1::2], map(float, line[2::2])))
            asph.append(kv)
        else:
            rows.append(line)
    worst = 0.0
    assert float(rows[0][1]) == STOP_GAP
    for i, (r, t, n, v, K, a) in enumerate(SURFACES):
        row, kv = rows[i + 1], asph[i]
        assert float(row[0]) == r and float(row[1]) == t and float(row[2]) == n, row
        if n > 1:
            assert float(row[4]) == v, row
        worst = max(worst, abs(kv['k'] - (K - 1)) / max(1.0, abs(K)))
        for p, ap in enumerate(a, start=3):
            assert kv.get('A%d' % p, 0.0) == ap, (i, p)
    print('file: surfaces, glasses and coefficients as the patent; conic k = K - 1 within %.1e' % worst)


def main():
    check_file()
    opm = build()
    fod = opm['analysis_results']['parax_data'].fod
    print('rayoptics: efl %.6f mm, bfl (from the last surface, in air incl. the cover glass) %.6f mm' % (fod.efl, fod.bfl))
    sm = opm['seq_model']
    print('indices at %g nm:' % LINE_D, ['%.6f' % float(np.atleast_1d(n)[0]) for n in sm.rndx])
    rows = []
    for deg in FIELDS_DEG:
        th = math.radians(deg)
        d = np.array([math.sin(th), 0.0, math.cos(th)])   # rayoptics: +z toward the image
        for (px, py) in STOP_POINTS:
            p0 = np.array([px - OBJECT_GAP * d[0] / d[2], py - OBJECT_GAP * d[1] / d[2], 0.0])
            ray, _, _ = trace(sm, p0, d, LINE_D)
            pt, after = ray[-1][0], ray[-2][1]
            rows.append((deg, px, py, pt[0], pt[1], after[0], after[1]))
    print('rays (field deg, stop x, y mm) -> image x, y (mm), direction x, y:')
    for r in rows:
        print('  %4.1f  %5.2f %5.2f  ->  %.9f %.9f   %.9f %.9f' % r)
    print('C++ table for CameraChecks (field deg, stop x mm, stop y mm, image x mm, image y mm):')
    for r in rows:
        print('    { %g, %g, %g, %.9f, %.9f },' % r[:5])


if __name__ == '__main__':
    sys.exit(main())
