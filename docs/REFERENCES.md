# References

The works that the code implements or follows, and the ones that the studies in `docs/` reproduce or test. Each entry says what it is used for.
Entries marked *(studied)* were implemented as a reference or evaluated, and are not part of the renderer.

**Rendering framework and sampling**

- M. Pharr, W. Jakob, G. Humphreys, *Physically Based Rendering: From Theory to Implementation*, 4th ed., MIT Press, 2023 (3rd ed., 2016). The overall design; the
  hero-wavelength walk through a medium (v4), light sampling, the dielectric BSDF, the probe-ray BSSRDF sampling (v3). [pbr-book.org](https://pbr-book.org)
- A. Wilkie, S. Nawaz, M. Droske, A. Weidlich, J. Hanika, "Hero Wavelength Spectral Sampling", *Computer Graphics Forum* 33(4), 2014. Four stratified wavelengths per path.
- E. Veach, L. J. Guibas, "Optimally Combining Sampling Techniques for Monte Carlo Rendering", SIGGRAPH 1995. Multiple importance sampling everywhere (next event estimation, strategies of the bidirectional estimator).
- L. Conty Estevez, C. Kulla, "Importance Sampling of Many Lights with Adaptive Tree Splitting", *Proc. ACM CGIT* 1(2), 2018. The spatial light sampler (the per-light importance of PBRT v4).
- T. Duff, J. Burgess, P. Christensen, C. Hery, A. Kensler, M. Liani, R. Villemin, "Building an Orthonormal Basis, Revisited", *JCGT* 6(1), 2017. Local frames.
- J. Ureña, M. Fajardo, A. King, "An Area-Preserving Parametrization for Spherical Rectangles", *Computer Graphics Forum* 32(4) (EGSR), 2013. Solid-angle sampling of rectangular lights.
- I. Wald, S. Woop, C. Benthin, G. Johnson, M. Ernst, "Embree: A Kernel Framework for Efficient CPU Ray Tracing", *ACM TOG* 33(4), 2014. The Embree backend.

**Metals and microfacets**

- P. B. Johnson, R. W. Christy, "Optical Constants of the Noble Metals", *Physical Review B* 6(12), 1972. The `n(λ)`, `k(λ)` of Au, Ag and Cu.
- O. Gulbrandsen, "Artist Friendly Metallic Fresnel", *JCGT* 3(4), 2014. A measured reflectance `R(λ)` to `(n, k)`, with the edge tint.
- B. Walter, S. R. Marschner, H. Li, K. E. Torrance, "Microfacet Models for Refraction through Rough Surfaces", EGSR 2007. Rough dielectrics.
- E. Heitz, "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs", *JCGT* 3(2), 2014; "Sampling the GGX Distribution of Visible Normals", *JCGT* 7(4), 2018. GGX and VNDF sampling.
- C. Kulla, A. Conty, "Revisiting Physically Based Shading at Imageworks", SIGGRAPH 2017 course. The multiple-scattering compensation of rough metals.
- C. Schlick, "An Inexpensive BRDF Model for Physically-based Rendering", *Computer Graphics Forum* 13(3), 1994. The Fresnel of the coated diffuse material.

**Spectra and color**

- C. Wyman, P.-P. Sloan, P. Shirley, "Simple Analytic Approximations to the CIE XYZ Color Matching Functions", *JCGT* 2(2), 2013. The film's color matching functions.
- CIE 15:2004, *Colorimetry* (superseded by CIE 015:2018). The 1931 observer, illuminants D65 and A.
- C. Li, Z. Li, Z. Wang, Y. Xu, M. R. Luo, G. Cui, M. Melgosa, M. H. Brill, M. Pointer, "Comprehensive color solutions: CAM16, CAT16, and CAM16-UCS", *Color Research & Application* 42(6), 2017. The chromatic adaptation of the white balance (`--wb`).
- CIE 142-2001, *Improvement to Industrial Colour-Difference Evaluation* (CIEDE2000); G. Sharma, W. Wu, E. N. Dalal, "The CIEDE2000 Color-Difference Formula: Implementation Notes, Supplementary Test Data, and Mathematical Observations", *Color Research & Application* 30(1), 2005. The color differences of the RGB-vs-spectral studies.
- X. Zhang, B. A. Wandell, "A Spatial Extension of CIELAB for Digital Color Image Reproduction", *Journal of the Society for Information Display* 5(1), 61-63, 1997; "Color image fidelity metrics evaluated using image distortion maps", *Signal Processing* 70(3), 1998. S-CIELAB: the difference after the eye's spatial filtering.
- W. Jakob, J. Hanika, "A Low-Dimensional Function Space for Efficient Spectral Upsampling", *Computer Graphics Forum* 38(2), 2019. *(studied)* An alternative to the RGB-to-spectrum conversion; not adopted.
- C. Peters, S. Merzbach, J. Hanika, C. Dachsbacher, "Using Moments to Represent Bounded Signals for Spectral Rendering", *ACM TOG* 38(4), 2019. *(studied)* The metamer-mismatch study ([SPECTRUM_FROM_RGB.md](SPECTRUM_FROM_RGB.md)).

**Dispersion and thin films**

- I. H. Malitson, "Refraction and Dispersion of Synthetic Sapphire", *JOSA* 52, 1962; "Interspecimen Comparison of the Refractive Index of Fused Silica", *JOSA* 55, 1965; the Schott catalog (BK7); F. Peter, *Z. Phys.* 15, 358 (1923) for diamond. The Sellmeier presets.
- M. Born, E. Wolf, *Principles of Optics*, Cambridge University Press. Airy's formula for a thin film.
- B. Smits, G. Meyer, "Newton's Colors: Simulating Interference Phenomena in Realistic Image Synthesis", Eurographics Workshop on Photosimulation, Realism and Physics in Computer Graphics, 1990, pp. 185-194 (in *Photorealism in Computer Graphics*, Springer, 1992). An early use of interference in graphics.

**Participating media and subsurface scattering**

- J. Novák, I. Georgiev, J. Hanika, W. Jarosz, "Monte Carlo Methods for Volumetric Light Transport Simulation", *Computer Graphics Forum* 37(2) (EG STAR), 2018. Delta tracking, null collisions, the bookkeeping of the walk.
- L. G. Henyey, J. L. Greenstein, "Diffuse Radiation in the Galaxy", *Astrophysical Journal* 93, 1941. The phase function.
- I. Gkioulekas, B. Xiao, S. Zhao, E. H. Adelson, T. Zickler, K. Bala, "Understanding the Role of Phase Function in Translucent Appearance", *ACM TOG* 32(5), 2013. The two-lobe phase function.
- P. H. Christensen, B. Burley, "Approximate Reflectance Profiles for Efficient Subsurface Scattering", Pixar technical memo, 2015. The normalized diffusion profile of the skin BSSRDF.
- S. L. Jacques, "Optical properties of biological tissues: a review", *Physics in Medicine & Biology* 58(11), 2013; S. Prahl, hemoglobin spectra ([omlc.org](https://omlc.org/spectra/hemoglobin/)). The chromophore model of the layered skin.

**Caustics and light transport**

- H. W. Jensen, "Global Illumination using Photon Maps", Eurographics Workshop on Rendering, 1996. The caustic photon map.
- J. Hanika, M. Droske, L. Fascione, "Manifold Next Event Estimation", *Computer Graphics Forum* 34(4) (EGSR), 2015; T. Zeltner, I. Georgiev, W. Jakob, "Specular Manifold Sampling for Rendering High-Frequency Caustics and Glints", *ACM TOG* 39(4), 2020. The specular manifold NEE of the prism caustic.
- E. Veach, *Robust Monte Carlo Methods for Light Transport Simulation*, PhD thesis, Stanford University, 1997; E. P. Lafortune, Y. D. Willems, "Bi-directional Path Tracing", Compugraphics 1993. Light tracing and bidirectional path tracing, with specular vertices.
- I. Georgiev, J. Křivánek, T. Davidovič, P. Slusallek, "Light Transport Simulation with Vertex Connection and Merging", *ACM TOG* 31(6) (SIGGRAPH Asia), 2012;
  T. Davidovič, J. Křivánek, M. Hašan, P. Slusallek, "Progressive Light Transport Simulation on the GPU: Survey and Improvements", *ACM TOG* 33(3), 2014. The vertex cache and the weights of the bidirectional estimators.
- shocker_0x15, "Renderer introduction: Ray Tracing Camp 10" (slides, in Japanese), [Speaker Deck](https://speakerdeck.com/shocker_0x15/renderer-introduction-ray-tracing-camp-10). Volumetric and dispersive caustics by light vertex cache bidirectional path tracing: the target of the fog scenes.

**3D Gaussian Splatting**

- J. Gao, C. Gu, Y. Lin, H. Zhu, X. Cao, L. Zhang, Y. Yao, "Relightable 3D Gaussian: Real-time Point Cloud Relighting with BRDF Decomposition and Ray Tracing", arXiv:2311.16043, 2023. *(studied)* The relightable reconstruction that the lighting study trains and relights.
- J.-N. Dihlmann, A. Majumdar, A. Engelhardt, R. Braun, H. P. A. Lensch, "Subsurface Scattering for 3D Gaussian Splatting", NeurIPS 2024 (arXiv:2408.12282). *(studied)* Reports that a relightable reconstruction bakes the subsurface scattering of translucent objects into the surface color.

**Software and background**

- W. Jakob, S. Speierer, N. Roussel, M. Nimier-David, D. Vicini, T. Zeltner, B. Nicolet, M. Crespo, V. Leroy, Z. Zhang, *Mitsuba 3 renderer*, 2022, [mitsuba-renderer.org](https://www.mitsuba-renderer.org). The RGB-vs-spectral study of the graduation research was done with it.
- J. Tessendorf, "Simulating Ocean Water", SIGGRAPH course notes, 2001 (the notes carry the copyright years 1999-2001). Background for the Gerstner-wave sea.
