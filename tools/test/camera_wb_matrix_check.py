"""Does a camera's colour matrix need fitting under the light it is balanced to? A calculation separate from the renderer
(colour-science), on the lines of Color/CameraResponse: under D65 and under illuminant A, gains that make the light's white grey and a
3 x 3 matrix fitted on smooth reflectances (Meng 2015 spectra of 6 x 6 x 6 sRGB colours; the renderer uses its own smooth basis) whose
targets are the CIE 1931 colours adapted from the light's white to D65 with CAT16. Tested on the 24 ColorChecker patches (N Ohta),
which are not among the fitted reflectances: CIEDE2000 against the targets, with the matrix fitted under the light and, under A, with
D65's matrix and A's gains (white balance by gains alone).

    python tools/test/camera_wb_matrix_check.py GoLD_rayt/assets/CameraResponse/Canon_5DMarkII.csv     (about 6 minutes: Meng 2015)

Canon 5D Mark II: D65 mean 0.90 (max 2.89); A with its own matrix 1.29 (5.80); A with D65's matrix 2.48 (7.47).
"""
import sys, numpy as np, colour, warnings
warnings.filterwarnings('ignore')
shape = colour.SpectralShape(380, 780, 5)
wl = shape.wavelengths
cam = np.array([[float(v) for v in l.split(',')] for l in open(sys.argv[1]) if l.strip() and l[0].isdigit()])
curves = np.stack([np.interp(wl, cam[:, 0], cam[:, k], left=0, right=0) for k in (1, 2, 3)], 1)
cmf = colour.MSDS_CMFS['CIE 1931 2 Degree Standard Observer'].copy().align(shape).values
d65 = colour.SDS_ILLUMINANTS['D65'].copy().align(shape).values
illA = colour.SDS_ILLUMINANTS['A'].copy().align(shape).values
M_XYZ_SRGB = colour.models.RGB_COLOURSPACE_sRGB.matrix_XYZ_to_RGB
def xyz(E, R): return (cmf * (E * R)[:, None]).sum(0) / (cmf[:, 1] * d65).sum()
def raw(E, R): return (curves * (E * R)[:, None]).sum(0)
levels = [0.04, 0.2, 0.4, 0.6, 0.8, 0.96]
CMFS = colour.MSDS_CMFS['CIE 1931 2 Degree Standard Observer'].copy().align(shape)
ILL = colour.SDS_ILLUMINANTS['D65'].copy().align(shape)
train = []
for r in levels:
    for g in levels:
        for b in levels:
            sd = colour.recovery.XYZ_to_sd_Meng2015(colour.sRGB_to_XYZ([r, g, b], apply_cctf_decoding=False), cmfs=CMFS, illuminant=ILL)
            train.append(np.clip(sd.values, 0, 1))
test = [sd.copy().align(shape).values for sd in colour.SDS_COLOURCHECKERS['ColorChecker N Ohta'].values()]
def fit(E):
    W = xyz(E, np.ones_like(wl)); Wd = xyz(d65, np.ones_like(wl)) * W[1]
    cat = colour.adaptation.matrix_chromatic_adaptation_VonKries(W, Wd, transform='CAT16')
    rw = raw(E, np.ones_like(wl)) / W[1]
    X = np.array([raw(E, R) / rw for R in train]); T = np.array([M_XYZ_SRGB @ cat @ xyz(E, R) for R in train])
    A = X.T @ X; m = np.zeros((3, 3)); ones = np.ones(3); Ai = np.linalg.inv(A)
    for r in range(3):
        ls = Ai @ (X.T @ T[:, r]); m[r] = ls - Ai @ ones * ((ones @ ls - 1) / (ones @ Ai @ ones))
    return rw, m, cat
rwD, mD, _ = fit(d65)
rwA, mA, catA = fit(illA)
def de(E, rw, m, cat):
    out = []
    for R in test:
        rgb = m @ (raw(E, R) / rw); ref = M_XYZ_SRGB @ cat @ xyz(E, R)
        lab = lambda c: colour.XYZ_to_Lab(np.linalg.inv(M_XYZ_SRGB) @ c, colour.CCS_ILLUMINANTS['CIE 1931 2 Degree Standard Observer']['D65'])
        out.append(colour.delta_E(lab(rgb), lab(ref), method='CIE 2000'))
    return np.array(out)
eye = np.eye(3)
for name, E, rw, m, cat in (('D65, its own matrix', d65, rwD, mD, eye), ('A, its own matrix', illA, rwA, mA, catA), ('A, D65 matrix + A gains', illA, rwA, mD, catA)):
    d = de(E, rw, m, cat); print('%-26s ColorChecker CIEDE2000 mean %.2f, max %.2f' % (name, d.mean(), d.max()))
