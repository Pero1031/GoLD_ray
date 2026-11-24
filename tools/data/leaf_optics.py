"""
Data for the leaf material (docs/LEAF.md): the PROSPECT coefficient table and the ANGERS leaf database.

Neither is part of this repository; this script downloads them to git-ignored folders under GoLD_rayt/assets/Leaf/.

- PROSPECT: dataSpec_PRO_v2.txt of the authors' R package (J.-B. Feret, github.com/jbferet/prospect, GPL-3): refractive index and
  specific absorption coefficients of chlorophylls, carotenoids, anthocyanins, brown pigments, water and dry matter, 400-2500 nm, 1 nm.
  Cite Feret et al. 2017 (PROSPECT-D, Remote Sensing of Environment 193) and Feret et al. 2021 (PROSPECT-PRO, RSE 252).
- ANGERS (2003): directional-hemispherical reflectance and transmittance of 308 leaves with their measured pigments, as shared by
  J.-B. Feret (gitlab.com/jbferet/myshareddata, LOP/ANGERS), which adds 32 anthocyanin-rich Acer pseudoplatanus leaves to the
  276 of the official set. Cite Feret et al. 2008 (RSE 112) and 2017.

  python tools/data/leaf_optics.py fetch
"""
import os
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEST = os.path.join(ROOT, "GoLD_rayt", "assets", "Leaf")
FILES = {
    "prospect/dataSpec_PRO_v2.txt": "https://raw.githubusercontent.com/jbferet/prospect/master/data-raw/dataSpec_PRO_v2.txt",
    "angers/DataBioch.txt": "https://gitlab.com/jbferet/myshareddata/-/raw/master/LOP/ANGERS/DataBioch.txt",
    "angers/ReflectanceData.txt": "https://gitlab.com/jbferet/myshareddata/-/raw/master/LOP/ANGERS/ReflectanceData.txt",
    "angers/TransmittanceData.txt": "https://gitlab.com/jbferet/myshareddata/-/raw/master/LOP/ANGERS/TransmittanceData.txt",
}


def fetch():
    for rel, url in FILES.items():
        path = os.path.join(DEST, *rel.split("/"))
        if os.path.exists(path):
            print("have", rel)
            continue
        os.makedirs(os.path.dirname(path), exist_ok=True)
        print("downloading", url)
        urllib.request.urlretrieve(url, path + ".part")
        os.replace(path + ".part", path)
    print("done:", DEST)


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "fetch":
        fetch()
    else:
        print(__doc__)
