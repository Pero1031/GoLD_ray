# Erato (3D scan)

- **Title:** Erato
- **Author:** Geoffrey Marchal (© 2016)
- **License:** [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)
- **Source:** Morgan McGuire's Computer Graphics Archive, <https://casual-effects.com/data/> (412.7k triangles, 207.3k vertices, updated 2017-07-02)

The mesh (`erato.obj`, 107 MB) and its texture (`erato-101.jpg`, 20 MB) are **not part of this repository** (over GitHub's file size
limit): download the model from the archive above and put the files in this folder.

The scene files use a cleaned copy (the scan's own vertex normals disagree with the triangle winding on a quarter of the triangles, which
shades metal black in places). Make it, from the repository folder:

```
python tools/mesh_clean.py GoLD_rayt/assets/Meshes/erato/erato.obj GoLD_rayt/assets/Meshes/erato/erato_clean.obj
python tools/mesh_clean.py GoLD_rayt/assets/Meshes/erato/erato.obj GoLD_rayt/assets/Meshes/erato/erato_closed.obj --cap   # closed, for a volume
```

Changes made in this project's renderings: the model is scaled and placed by the scene file, its two open boundaries (the cut
at the base and at the top) are closed where a volume is needed, and the surface is given other materials (gold, wax, marble, jade);
the photographed texture is not used for those.
