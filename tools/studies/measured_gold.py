"""Measured gold (KAIST pBRDF 2020 sample 6, unpolarized) against the engine's rough conductor.

    python tools/studies/measured_gold.py     -> out/studies/measured_gold/{results.txt, measured_gold.png}

This is the gold entry point of measured_metal.py (see there for what is measured); the geometry and model helpers are
re-exported so that pbrdf_albedo.py and pbrdf_calibration.py can import them from here.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from measured_metal import directions, half_diff, fresnel_conductor, load_au, load_metal, main as _main   # noqa: F401,E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

if __name__ == "__main__":
    sys.argv = [sys.argv[0],
                os.path.join(ROOT, "references", "kaist_pbrdf", "gold_unpolarized.npz"),
                os.path.join(ROOT, "GoLD_rayt", "assets", "Metal_data", "Johnson_Au.csv"),
                os.path.join(ROOT, "out", "studies", "measured_gold"),
                "gold"]
    _main()
