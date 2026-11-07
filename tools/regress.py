"""Regression harness for refactoring: small renders (mean color) and self-checks of a build, compared with an earlier build.

  1. keep the old build:   copy out\\cli to out\\old-rgb (and out\\cli-spectral to out\\old-spectral)
  2. record it:            python tools/regress.py run old-rgb base_rgb 2          (a folder under out\\, a tag, the number of repeats)
  3. change the code, rebuild (tools\\build_cli.ps1), copy the new build to out\\new-rgb, record it:  python tools/regress.py run new-rgb new_rgb 2
  4. compare:              python tools/regress.py compare base_rgb new_rgb        (a tag may be a list: a,b)
  (run takes an optional fifth argument: the names of the cases to run, separated by commas)

A picture is "ok" when its mean color differs from the baseline by less than 3 x the spread of the baseline repeats (at least 1 %).
With only 2 repeats this is a weak test: a case that is slightly out is not evidence; repeat it 30-60 times with tools/ab_test.py.
The self-checks must give the same exit code and number of OK lines; some are flaky at the low spp used here (furnace, spectral build).
Results go to out/regress/.
"""
import io, json, os, re, subprocess, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJ = os.path.join(ROOT, "GoLD_rayt")
OUT = os.path.join(ROOT, "out")
S = os.path.join(OUT, "regress")
os.makedirs(S, exist_ok=True)
SC = os.path.join(PROJ, "scenes")
SMALL = ["--width", "192", "--spp", "16"]

# (name, args). Pictures: the mean color is compared. Checks: exit code and the number of OK / FAIL lines.
PICTURES = [
    ("gold", ["--scene", "gold"] + SMALL),
    ("goldrough", ["--scene", "goldrough"] + SMALL),
    ("glass", ["--scene", "glass"] + SMALL),
    ("glass_abbe", ["--scene", "glass", "--abbe", "30"] + SMALL),
    ("studio", ["--scene", "studio"] + SMALL),
    ("spd_tung", ["--scene", "spd", "--light", "tungsten"] + SMALL),
    ("goldbowl", ["--scene", "goldbowl", "--light", "tungsten"] + SMALL),
    ("blob_gray", ["--scene", "blob"] + SMALL),
    ("blob_chromatic", ["--scene", "blob", "--sigma-t-rgb", "250,150,250", "--albedo-rgb", "0.85,0.998,0.95", "--phase-g", "0.3"] + SMALL),
    ("blob_twolobe", ["--scene", "blob", "--phase-g", "0.9", "--phase-g2", "-0.4", "--phase-w", "0.7"] + SMALL),
    ("file_trio", [os.path.join(SC, "jade_wax_milk.json")] + SMALL),
    ("file_phase3b", [os.path.join(SC, "phase_backlit.json")] + SMALL),
    ("file_cloud", [os.path.join(SC, "cloud_noise.json")] + SMALL),
    ("water", ["--scene", "water"] + SMALL),
    ("sea", ["--scene", "sea"] + SMALL),
    ("skindemo_sss", ["--scene", "skindemo", "--sss"] + SMALL),
    ("manifold", ["--scene", "manifold"] + SMALL),
    ("furnace", ["--scene", "furnace"] + SMALL),
    ("mesh", ["--scene", "mesh"] + SMALL),
]
CHECKS = [
    ("chk_volume", ["--scene", "volume", "--spp", "128"]),
    ("chk_lights", ["--scene", "lights", "--spp", "256"]),
    ("chk_aov", ["--scene", "aov"]),
    ("chk_sss", ["--scene", "sss"]),
    ("chk_metals", ["--scene", "metals"]),
]


def pfm_mean(path):
    with open(path, 'rb') as f:
        t = f.readline().strip(); w, h = map(int, f.readline().split()); s = float(f.readline())
        a = np.frombuffer(f.read(), dtype='<f4' if s < 0 else '>f4').reshape(h, w, 3)
    return a.reshape(-1, 3).mean(0).tolist()


def run(build_dir, tag, repeats, only=None):
    exe = os.path.join(OUT, build_dir, "GoLD_rayt.exe")
    res = {}
    os.makedirs(os.path.join(S, "reg_" + tag), exist_ok=True)
    for rep in range(repeats):
        for name, args in [c for c in PICTURES + CHECKS if only is None or c[0] in only]:
            is_pic = (name, args) in PICTURES
            cmd = [exe] + args
            if is_pic:
                pf = os.path.join(S, "reg_" + tag, "%s_%d.pfm" % (name, rep))
                cmd += ["--out", pf]
            t0 = time.time()
            try:
                p = subprocess.run(cmd, cwd=PROJ, capture_output=True, timeout=900)
                text = p.stdout.decode('utf-8', errors='replace')
                code = p.returncode
            except subprocess.TimeoutExpired:
                text, code = "", -9
            dt = time.time() - t0
            r = res.setdefault(name, {"runs": []})
            entry = {"code": code, "sec": round(dt, 1)}
            if is_pic:
                entry["mean"] = pfm_mean(pf) if code == 0 and os.path.exists(pf) else None
            else:
                entry["ok"] = len(re.findall(r"\bOK\b", text)); entry["fail"] = len(re.findall(r"\bFAIL\b", text))
                entry["result"] = re.findall(r"RESULT: \w+", text)
            r["runs"].append(entry)
            print(tag, rep, name, entry, flush=True)
    json.dump(res, open(os.path.join(S, "results_%s.json" % tag), "w"), indent=1)


def compare(tagsA, tagsB):
    A, B = {}, {}
    for t in tagsA:
        for k, v in json.load(open(os.path.join(S, "results_%s.json" % t))).items(): A.setdefault(k, []).extend(v["runs"])
    for t in tagsB:
        for k, v in json.load(open(os.path.join(S, "results_%s.json" % t))).items(): B.setdefault(k, []).extend(v["runs"])
    bad = 0
    for name in A:
        a, b = A[name], B.get(name)
        if not b:
            print("MISSING", name); bad += 1; continue
        if "mean" in a[0]:
            if any(r["code"] != 0 or r["mean"] is None for r in a + b):
                print("FAILED RUN", name, [r["code"] for r in a], [r["code"] for r in b]); bad += 1; continue
            ma = np.array([r["mean"] for r in a]); mb = np.array([r["mean"] for r in b])
            spread = (ma.max(0) - ma.min(0)) if len(ma) > 1 else np.zeros(3)
            tol = np.maximum(3 * spread, 0.01 * np.abs(ma.mean(0))) + 1e-6
            d = np.abs(mb.mean(0) - ma.mean(0))
            ok = bool((d <= tol).all())
            bad += (not ok)
            print("%-4s %-16s base %s  new %s  diff/tol %s  secs %s -> %s" % ("ok" if ok else "DIFF", name, np.round(ma.mean(0), 4), np.round(mb.mean(0), 4),
                  np.round(d / tol, 2), round(np.mean([r['sec'] for r in a]), 1), round(np.mean([r['sec'] for r in b]), 1)))
        else:
            ca = [(r["code"], r["ok"], r["fail"], tuple(r["result"])) for r in a]
            cb = [(r["code"], r["ok"], r["fail"], tuple(r["result"])) for r in b]
            same = all(x[0] == 0 and x[2] == 0 for x in cb) and all(x[3] == cb[0][3] for x in cb)
            print("%-4s %-16s base %s  new %s" % ("ok" if same else "DIFF", name, ca[0], cb[0]))
            bad += (not same)
    print("TOTAL DIFFERENCES:", bad)


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5].split(",") if len(sys.argv) > 5 else None)
    else:
        compare(sys.argv[2].split(","), sys.argv[3].split(","))
