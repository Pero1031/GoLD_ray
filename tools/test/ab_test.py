"""A/B tests of two builds on the same scene: the mean color, the time, or how often a self-check fails.

    python tools/test/ab_test.py mean  <n> <buildA> <buildB> <scene args...>     n interleaved runs; mean and standard error per build, and a z score
    python tools/test/ab_test.py time  <n> <buildA> <buildB> <scene args...>     n interleaved runs; median and minimum time (this PC varies up to 1.4x between runs)
    python tools/test/ab_test.py flaky <n> <buildA,buildB,...> <scene args...>   how many of n runs fail (exit code != 0), and which lines

Builds are folders under out\\ (out\\cli, out\\cli-spectral, or copies you kept). Example:
    python tools/test/ab_test.py mean 30 old-rgb cli --scene blob --width 192 --spp 16

Lessons: 2-8 runs can show differences of 1-3 % that disappear with 30-60 runs; compare by z score, not by eye.
"""
import collections
import os
import re
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PROJ = os.path.join(ROOT, "GoLD_rayt")
OUT = os.path.join(ROOT, "out")
TMP = os.path.join(os.environ.get("TEMP", OUT), "ab_test.pfm")


def run(build, args):
    t0 = time.time()
    p = subprocess.run([os.path.join(OUT, build, "GoLD_rayt.exe")] + args + ["--out", TMP], cwd=PROJ, capture_output=True)
    return p, time.time() - t0


def mean_of(path):
    with open(path, "rb") as f:
        f.readline()
        w, h = map(int, f.readline().split())
        s = float(f.readline())
        a = np.frombuffer(f.read(), dtype="<f4" if s < 0 else ">f4").reshape(h, w, 3)
    return a.reshape(-1, 3).mean(0)


def main():
    mode, n = sys.argv[1], int(sys.argv[2])
    if mode == "flaky":
        builds, args = sys.argv[3].split(","), sys.argv[4:]
        for b in builds:
            fails, lines = 0, collections.Counter()
            for _ in range(n):
                p, _t = run(b, args)
                if p.returncode != 0:
                    fails += 1
                    for l in p.stdout.decode("utf-8", "replace").split("\n"):
                        if "FAIL" in l and "RESULT" not in l:
                            lines[re.sub(r"[-\d.e+]+", "#", l.strip())[:110]] += 1
            print("%-16s failed %2d of %d" % (b, fails, n))
            for l, c in lines.most_common(4):
                print("     %2dx  %s" % (c, l))
        return
    builds, args = sys.argv[3:5], sys.argv[5:]
    res = {b: [] for b in builds}
    for i in range(n):
        for b in (builds if i % 2 == 0 else builds[::-1]):      # interleaved: the speed of the PC cannot favor one build
            p, dt = run(b, args)
            res[b].append(mean_of(TMP) if mode == "mean" else dt)
    if mode == "time":
        for b in builds:
            print("%-16s median %.2f s  min %.2f s  (n=%d)" % (b, np.median(res[b]), np.min(res[b]), n))
        print("B / A: median %.3f   min %.3f" % (np.median(res[builds[1]]) / np.median(res[builds[0]]), np.min(res[builds[1]]) / np.min(res[builds[0]])))
        return
    m = {}
    for b in builds:
        a = np.array(res[b])
        m[b] = (a.mean(0), a.std(0, ddof=1) / np.sqrt(n))
        print("%-16s mean %s  se %s" % (b, np.round(m[b][0], 5), np.round(m[b][1], 5)))
    z = (m[builds[1]][0] - m[builds[0]][0]) / np.sqrt(m[builds[0]][1] ** 2 + m[builds[1]][1] ** 2)
    print("relative difference %s   z %s" % (np.round(m[builds[1]][0] / m[builds[0]][0] - 1, 4), np.round(z, 2)))


if __name__ == "__main__":
    main()
