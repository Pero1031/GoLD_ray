# References

The works that the code implements or follows, and the ones that the studies in `docs/` reproduce or test. Each entry says what it is used for.
Entries marked *(studied)* were implemented as a reference or evaluated, and are not part of the renderer.

**How far the bibliographic data were checked.** The entries that README and the documents on hair, skin, clouds, ocean, gloss and thin films carried on 2026-11-09 were checked by web search that day (authors, title, venue, volume, year, pages; 6 corrections were made). The entries of Haro et al., Graham et al., Barros and Walter were checked from search summaries only. The entries added later were taken from the title page of the PDF where one was at hand (Moon et al., Dihlmann et al., Lanza et al., Křivánek et al., Speierer et al., Yao et al., Kouros et al.) and otherwise written from memory; those are marked *(from memory)*. Not confirmed: the edition of Born & Wolf, the conference name (EGSR) of Hanika et al. 2015 (volume and issue are confirmed), and the source of the melanin absorption values (taken from PBRT). Remarks such as *(Text read.)* say whether the work was read in full or only cited second-hand. The entries added on 2026-11-19 (fluorescence, the camera, the gamut study, water drops and a few others that the documents cited but this list lacked) and the earlier ones marked *(from memory)* were checked by web search that day (authors, title, venue, volume, pages, DOI: DBLP, Crossref, the publishers' and authors' pages); one correction (the article number of Kerbl et al., which Crossref does not give, was removed).
The literature values behind the lotion and emulsion model are listed separately in [SKIN_MOISTURE_REFS.md](SKIN_MOISTURE_REFS.md) (search summaries; not cross-checked).

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
- E. Heitz, J. Hanika, E. d'Eon, C. Dachsbacher, "Multiple-Scattering Microfacet BSDFs with the Smith Model", *ACM TOG* 35(4), Article 58, 2016. The random walk on the microsurface, `"multipleScattering": "heitz"` and the reference of [MULTIBOUNCE.md](MULTIBOUNCE.md).
- S. Zhang, T. Funatomi, Y. Fujimura, Y. Mukaigawa, H. Kubo, "Successive Height Preintegration for a Height-Interdependent Path Formulation of Multiple Bounces in Smith Microfacet BRDFs", SIGGRAPH 2026 Conference Papers, doi:10.1145/3799902.3811194 (with its supplement). *(studied)* `"multipleScattering": "zhang"` ([MULTIBOUNCE.md](MULTIBOUNCE.md)).
- B. Wang, W. Jin, J. Fan, J. Yang, N. Holzschuch, L.-Q. Yan, "Position-Free Multiple-Bounce Computations for Smith Microfacet BSDFs", *ACM TOG* 41(4), Article 134, 2022, and Y. Cui, G. Pan, J. Yang, L. Zhang, L.-Q. Yan, B. Wang, "Multiple-bounce Smith Microfacet BRDFs using the Invariance Principle", SIGGRAPH Asia 2023. The vertex / segment term framework and the path estimator that Zhang et al. build on (read, not implemented).
- J. Dupuy, A. Benyoub, "Sampling Visible GGX Normals with Spherical Caps", *Computer Graphics Forum* 42(8), 2023. The visible-normal sampling from either hemisphere in the multiple-bounce walks.
- C. Schlick, "An Inexpensive BRDF Model for Physically-based Rendering", *Computer Graphics Forum* 13(3), 1994. The Fresnel of the coated diffuse material.
- A. D. Rakić, "Algorithm for the determination of intrinsic optical constants of metal films: application to aluminum", *Applied Optics* 34(22), 4755-4767, 1995, doi:10.1364/AO.34.004755; K. M. McPeak, S. V. Jayanti, S. J. P. Kress, S. Meyer, S. Iotti, A. Rossinelli, D. J. Norris, "Plasmonic Films Can Easily Be Better: Rules and Recipes", *ACS Photonics* 2(3), 326-333, 2015, doi:10.1021/ph5004237. The `n(λ)`, `k(λ)` of aluminum (`assets/Metal_data/Rakic_Al.csv`, `McPeak_Al.csv`), for the comparison with the measured aluminum.
- E. Heitz, J. Hanika, E. d'Eon, C. Dachsbacher, "Multiple-Scattering Microfacet BSDFs with the Smith Model", *ACM TOG* 35(4), Article 58, 2016 (SIGGRAPH), doi:10.1145/2897824.2925943. The exact random walk over microfacets that the multiple-scattering compensation approximates; named in [METAL_ACCURACY.md](METAL_ACCURACY.md) as a reference that could settle the remaining color error, not implemented.
- P. Kemppinen, L. Paulin, T. Thonat, J.-M. Thiery, J. Lehtinen, T. Boubekeur, "Evaluating and Sampling Glinty NDFs in Constant Time", *ACM TOG* 44(6), 2025 (SIGGRAPH Asia), doi:10.1145/3763282 (the authors' preprint read in full; its article number is printed as 255 there and listed as 256 by the publisher). The glinty GGX distribution (`Microfacet/GlintNDF`, [`GLINTS.md`](GLINTS.md)). Its lineage, from abstracts and from that paper's related work, not read: W. Jakob, M. Hašan, L.-Q. Yan, J. Lawrence, R. Ramamoorthi, S. Marschner, "Discrete Stochastic Microfacet Models", *ACM TOG* 33(4), 2014; T. Zirr, A. S. Kaplanyan, "Real-time Rendering of Procedural Multiscale Materials", I3D 2016. Read in full and not used (its authors rule it out for path tracing: no importance sampling, no energy conservation, point lights only): T. Deliot, L. Belcour, "Real-Time Rendering of Glinty Appearances using Distributed Binomial Laws on Anisotropic Grids", *Computer Graphics Forum* 42(8), 2023, doi:10.1111/cgf.14866.

**Measured reflectance and polarization**

- S.-H. Baek, T. Zeltner, H. J. Ku, I. Hwang, X. Tong, W. Jakob, M. H. Kim, "Image-Based Acquisition and Modeling of Polarimetric Reflectance", *ACM TOG* 39(4), 2020 (SIGGRAPH). *(studied)* The KAIST polarimetric BRDF dataset: its gold sample against the rough conductor ([MEASURED_GOLD.md](MEASURED_GOLD.md)).
- Y. Moon, R. Maeda, S. Shin, I. Hwang, Y. Kim, M. H. Kim, S.-H. Baek, "Hyperspectral Polarimetric BRDFs of Real-world Materials", 2025 (preprint in ACM format, September 2025; the PDF names no venue), [project page](https://yunseong0518.github.io/projects/hpBRDF/). The hpBRDF tables: the `measured` material reads one (PCA-compressed, with its measured Mueller matrices), and its aluminum is compared with the rough conductor (MEASURED_GOLD.md, section 6).
- M. Born, E. Wolf, *Principles of Optics* (see "Dispersion and thin films"): also the Fresnel amplitudes and the Stokes and Mueller formalism of the polarization code ([POLARIZATION.md](POLARIZATION.md)).

**Spectra and color**

- CIE 15:2004, *Colorimetry* (superseded by CIE 015:2018). The 1931 observer (its 1 nm table is the film's color matching functions), illuminants D65 and A.
- C. Wyman, P.-P. Sloan, P. Shirley, "Simple Analytic Approximations to the CIE XYZ Color Matching Functions", *JCGT* 2(2), 2013. The film's color matching functions until 2026-11-17 (CIEDE2000 0.4 median, 1.2 at most from the table on the `spd` scene); the studies in `tools/studies` that mirror the film still use it.
- C. Li, Z. Li, Z. Wang, Y. Xu, M. R. Luo, G. Cui, M. Melgosa, M. H. Brill, M. Pointer, "Comprehensive color solutions: CAM16, CAT16, and CAM16-UCS", *Color Research & Application* 42(6), 2017. The chromatic adaptation of the white balance (`--wb`).
- CIE 142-2001, *Improvement to Industrial Colour-Difference Evaluation* (CIEDE2000); G. Sharma, W. Wu, E. N. Dalal, "The CIEDE2000 Color-Difference Formula: Implementation Notes, Supplementary Test Data, and Mathematical Observations", *Color Research & Application* 30(1), 2005. The color differences of the RGB-vs-spectral studies.
- X. Zhang, B. A. Wandell, "A Spatial Extension of CIELAB for Digital Color Image Reproduction", *Journal of the Society for Information Display* 5(1), 61-63, 1997; "Color image fidelity metrics evaluated using image distortion maps", *Signal Processing* 70(3), 1998. S-CIELAB: the difference after the eye's spatial filtering.
- W. Jakob, J. Hanika, "A Low-Dimensional Function Space for Efficient Spectral Upsampling", *Computer Graphics Forum* 38(2), 2019. *(studied)* An alternative to the RGB-to-spectrum conversion; not adopted.
- C. Peters, S. Merzbach, J. Hanika, C. Dachsbacher, "Using Moments to Represent Bounded Signals for Spectral Rendering", *ACM TOG* 38(4), 2019. *(studied)* The metamer-mismatch study ([SPECTRUM_FROM_RGB.md](SPECTRUM_FROM_RGB.md)).
- G. Ward, E. Eydelberg-Vileshin, "Picture Perfect RGB Rendering Using Spectral Prefiltering and Sharp Color Primaries", Eurographics Workshop on Rendering 2002, 117-124, doi:10.2312/EGWR/EGWR02/117-124. *(the Sharp matrix used was written from memory)* The "Sharp" working space in the study of RGB working spaces ([METAL_ACCURACY.md](METAL_ACCURACY.md) 3.7).

**Dispersion and thin films**

- I. H. Malitson, "Refraction and Dispersion of Synthetic Sapphire", *JOSA* 52, 1962; "Interspecimen Comparison of the Refractive Index of Fused Silica", *JOSA* 55, 1965; the Schott catalog (BK7); F. Peter, *Z. Phys.* 15, 358 (1923) for diamond. The Sellmeier presets.
- M. Born, E. Wolf, *Principles of Optics*, Cambridge University Press. Airy's formula for a thin film.
- B. Smits, G. Meyer, "Newton's Colors: Simulating Interference Phenomena in Realistic Image Synthesis", Eurographics Workshop on Photosimulation, Realism and Physics in Computer Graphics, 1990, pp. 185-194 (in *Photorealism in Computer Graphics*, Springer, 1992). An early use of interference in graphics.
- L. Belcour, P. Barla, "A Practical Extension to Microfacet Theory for the Modeling of Varying Iridescence", *ACM TOG* 36(4), Article 65, 2017 (SIGGRAPH). Iridescence on rough surfaces; background for the thin-film material, which is smooth.

**Fluorescence**

- S. Gonzalez, M. D. Fairchild, "Evaluation of bispectral spectrophotometry for accurate colorimetry of printing materials", Color Imaging Conference 2000, pp. 14-23 (as cited by Hua et al.); the database: RIT Munsell Color Science Laboratory, "Fluorescent Database". The measured bispectral matrices of `Fluorescent` (downloaded by `tools/data/rit_bispectral.py`, not redistributed; [FLUORESCENCE.md](FLUORESCENCE.md)).
- Q. Hua, V. Tázlar, A. Fichet, A. Wilkie, "Efficient Storage and Importance Sampling for Fluorescent Reflectance", *Computer Graphics Forum* 42(1), 47-59, 2023, doi:10.1111/cgf.14716. Fluorescent reflectance in a spectral renderer; the source cited here for the RIT database.
- ASTM E2152, "Standard Practice for Computing the Colors of Fluorescent Objects from Bispectral Photometric Data". The definition of the matrix (the standard itself was not read).

**Layered materials and gloss**

- W. Jakob, E. d'Eon, O. Jakob, S. Marschner, "A Comprehensive Framework for Rendering Layered Materials", *ACM TOG* 33(4), Article 118, 2014 (SIGGRAPH). The framework of layered materials; `LayeredCoat` sums the round trips in the layers by a random walk, as the `LayeredBxDF` of PBRT v4 does ([GLOSS.md](GLOSS.md)).

**Participating media and subsurface scattering**

- J. Novák, I. Georgiev, J. Hanika, W. Jarosz, "Monte Carlo Methods for Volumetric Light Transport Simulation", *Computer Graphics Forum* 37(2) (EG STAR), 2018. Delta tracking, null collisions, the bookkeeping of the walk.
- L. G. Henyey, J. L. Greenstein, "Diffuse Radiation in the Galaxy", *Astrophysical Journal* 93, 1941. The phase function.
- I. Gkioulekas, B. Xiao, S. Zhao, E. H. Adelson, T. Zickler, K. Bala, "Understanding the Role of Phase Function in Translucent Appearance", *ACM TOG* 32(5), 2013. The two-lobe phase function.
- P. H. Christensen, B. Burley, "Approximate Reflectance Profiles for Efficient Subsurface Scattering", Pixar technical memo, 2015. The normalized diffusion profile of the skin BSSRDF.
- S. L. Jacques, "Optical properties of biological tissues: a review", *Physics in Medicine & Biology* 58(11), 2013; S. Prahl, hemoglobin spectra ([omlc.org](https://omlc.org/spectra/hemoglobin/)). The chromophore model of the layered skin.
- M. J.-Y. Chiang, P. Kutz, B. Burley, "Practical and Controllable Subsurface Scattering for Production Path Tracing", SIGGRAPH 2016 Talks. The volumetric random walk of the skin and its albedo inversion ([SSS.md](SSS.md)). (Read in full.)
- M. Wrenninge, R. Villemin, C. Hery, "Path Traced Subsurface Scattering using Anisotropic Phase Functions and Non-Exponential Free Flights", Pixar Technical Memo 17-07, 2017. Anisotropic random walks. (Read in full.)
- E. d'Eon, "A Hitchhiker's Guide to Multiple Scattering" (online book); S. Chandrasekhar, *Radiative Transfer*, Oxford, 1950 (Dover reprint, 1960). The H function and van de Hulst's formula used in the albedo inversion. The formula was taken from the source of Cycles (`subsurface_random_walk.h`), which cites d'Eon; the book was not read.

**Clouds, ocean and sky**

- C. F. Bohren, D. R. Huffman, *Absorption and Scattering of Light by Small Particles*, Wiley, 1983. Mie scattering of the cloud droplets (BHMIE; `tools/studies/mie_phase.py`, [CLOUD.md](CLOUD.md)).
- J. E. Hansen, L. D. Travis, "Light Scattering in Planetary Atmospheres", *Space Science Reviews* 16, 527-610, 1974. The size distribution of the droplets.
- T. Nishita, T. Sirai, K. Tadamura, E. Nakamae, "Display of the Earth Taking into Account Atmospheric Scattering", SIGGRAPH 1993, 175-182. The single-scattering sky (`tools/data/make_sky.py`).
- J. Tessendorf, "Simulating Ocean Water", SIGGRAPH course notes, 2001 (depending on the edition of the notes the year is 2001 or 2002; the notes carry the copyright years 1999-2001). Background for the wave-sum sea ([WATER.md](WATER.md)).
- W. J. Pierson, L. Moskowitz, "A proposed spectral form for fully developed wind seas based on the similarity theory of S. A. Kitaigorodskii", *J. Geophys. Res.* 69(24), 5181-5190, 1964. The wave spectrum.
- K. Hasselmann et al., "Measurements of wind-wave growth and swell decay during the Joint North Sea Wave Project (JONSWAP)", *Deutsche Hydrographische Zeitschrift*, Reihe A (8°), Nr. 12, 1973 (16 authors). The peak enhancement γ of the spectrum.
- C. Cox, W. Munk, "Measurement of the roughness of the sea surface from photographs of the sun's glitter", *J. Opt. Soc. Am.* 44(11), 838-850, 1954. The slope variance behind the roughness of the distant water.
- R. M. Pope, E. S. Fry, "Absorption spectrum (380-700 nm) of pure water. II. Integrating cavity measurements", *Applied Optics* 36(33), 8710-8723, 1997. The absorption of water (`assets/Water/`).
- K. Kaneda, T. Kagawa, H. Yamashita, "Animation of Water Droplets on a Glass Plate", Computer Animation '93, 177-189, 1993 (in *Models and Techniques in Computer Animation*, Springer, doi:10.1007/978-4-431-66911-1_17). The "part of a sphere" drop that the water-drop generator extends ([WATERDROPS.md](WATERDROPS.md)).
- C. M. Schmidt, B. Budge, "Simple Nested Dielectrics in Ray Traced Images", *Journal of Graphics Tools* 7(2), 1-8, 2002. The priority approach to overlapping media, which the water drops replace with explicit contact faces.

**Caustics and light transport**

- H. W. Jensen, "Global Illumination using Photon Maps", Eurographics Workshop on Rendering, 1996. The caustic photon map.
- J. Hanika, M. Droske, L. Fascione, "Manifold Next Event Estimation", *Computer Graphics Forum* 34(4) (EGSR), 2015; T. Zeltner, I. Georgiev, W. Jakob, "Specular Manifold Sampling for Rendering High-Frequency Caustics and Glints", *ACM TOG* 39(4), 2020. The specular manifold NEE of the prism caustic.
- E. Veach, *Robust Monte Carlo Methods for Light Transport Simulation*, PhD thesis, Stanford University, 1997; E. P. Lafortune, Y. D. Willems, "Bi-directional Path Tracing", Compugraphics 1993. Light tracing and bidirectional path tracing, with specular vertices.
- I. Georgiev, J. Křivánek, T. Davidovič, P. Slusallek, "Light Transport Simulation with Vertex Connection and Merging", *ACM TOG* 31(6) (SIGGRAPH Asia), 2012;
  T. Davidovič, J. Křivánek, M. Hašan, P. Slusallek, "Progressive Light Transport Simulation on the GPU: Survey and Improvements", *ACM TOG* 33(3), 2014. The vertex cache and the weights of the bidirectional estimators. (Read in part: in VCM, vertex merging, its weights and the radius schedule; in Davidovič et al., the light vertex cache and VCM on it.)
- J. Novák, V. Havran, C. Dachsbacher, "Path Regeneration for Interactive Path Tracing", Eurographics 2010 short papers, 61–64. (Cited second-hand, from Davidovič et al.) Bidirectional path tracing on the GPU that connects each camera vertex to all vertices of a randomly chosen light path: the connection scheme of the prototype.
- shocker_0x15, "Renderer introduction: Ray Tracing Camp 10" (slides, in Japanese), [Speaker Deck](https://speakerdeck.com/shocker_0x15/renderer-introduction-ray-tracing-camp-10). Volumetric and dispersive caustics by light vertex cache bidirectional path tracing: the target of the fog scenes.
- J. Křivánek, I. Georgiev, T. Hachisuka, P. Vévoda, M. Šik, D. Nowrouzezahrai, W. Jarosz, "Unifying Points, Beams, and Paths in Volumetric Light Transport Simulation", *ACM TOG* 33(4), 2014 (SIGGRAPH), doi:10.1145/2601097.2601219. (Read in part: the estimators and the results.) Related work of the bidirectional estimator in fog ([LVC_BPT.md](LVC_BPT.md)): point and beam estimators combined with paths by multiple importance sampling. The prototype has no beams.
- S. Speierer, C. Hery, R. Villemin, W. Jakob, "Caustic Connection Strategies for Bidirectional Path Tracing", Pixar Technical Memo #18-01, 2018. (Text read, from a Wayback Machine copy: the Pixar page now redirects.) Connection strategies through chains of refractive events, found by manifold exploration, for caustics that classical bidirectional path tracing finds poorly: the manifold connections of [LVC_BPT.md](LVC_BPT.md) section 8.
- D. van Antwerpen, "Improving SIMD Efficiency for Parallel Monte Carlo Light Transport on the GPU", High Performance Graphics 2011, 41-50, doi:10.1145/2018323.2018330. The recursive MIS weights that [LVC_BPT.md](LVC_BPT.md) names as a way to make its weights O(1); not implemented.

**Hair**

- S. R. Marschner, H. W. Jensen, M. Cammarano, S. Worley, P. Hanrahan, "Light Scattering from Human Hair Fibers", *ACM TOG* 22(3), 780-791, 2003 (SIGGRAPH). The three lobes R, TT, TRT ([HAIR.md](HAIR.md)).
- E. d'Eon, G. Francois, M. Hill, J. Letteri, J.-M. Aubry, "An Energy-Conserving Hair Reflectance Model", *Computer Graphics Forum* 30(4), 1181-1187, 2011 (EGSR). Energy conservation, the equation of M_p, the approximation that sets σ_a from the apparent color.
- M. J.-Y. Chiang, B. Bitterli, C. Tappan, B. Burley, "A Practical and Controllable Hair and Fur Model for Production Path Tracing", *Computer Graphics Forum* 35(2), 275-283, 2016 (Eurographics). The logistic distribution of N_p, the roughness equations, importance sampling.
- A. Zinke, C. Yuksel, A. Weber, J. Keyser, "Dual Scattering Approximation for Fast Multiple Scattering in Hair", *ACM TOG* 27(3), Article 32, 2008 (SIGGRAPH).

**Skin: micro-relief, pigment and wetness**

- S. Worley, "A Cellular Texture Basis Function", SIGGRAPH 1996, 291-294. The Worley noise of the micro-relief ([SKIN.md](SKIN.md)).
- A. Haro, B. Guenter, I. Essa, "Real-time, Photo-realistic, Physically Based Rendering of Fine Scale Human Skin Structure", Eurographics Workshop on Rendering, 53-62, 2001. Normal maps made from a mould, extended by texture synthesis.
- P. Graham, B. Tunwattanapong, J. Busch, X. Yu, A. Jones, P. Debevec, A. Ghosh, "Measurement-Based Synthesis of Facial Microgeometry", *Computer Graphics Forum* 32(2), 335-344, 2013 (Eurographics). Synthesis of measured close-up patches.
- S. Weiss, J. Moulin, P. Chandran, G. Zoss, P. Gotardo, D. Bradley, "Graph-Based Synthesis for Skin Micro Wrinkles", *Computer Graphics Forum* 42(5), 2023 (Symposium on Geometry Processing). Furrows from a graph with the pores as points; implemented in SKIN.md, section 2.4. (Read in full.)
- R. S. Barros, M. Walter, "Synthesis of Human Skin Pigmentation Disorders", *Computer Graphics Forum* 36(1), 330-344, 2017. Pigment patches by reaction-diffusion. Not used here.
- M. Sawayama, E. H. Adelson, S. Nishida, "Visual wetness perception based on image color statistics", *Journal of Vision* 17(5):7, 2017. The image statistics of the wet look ([WETNESS.md](WETNESS.md)). (Text read.)
- J. Lekner, M. C. Dorf, "Why some things are darker when wet", *Applied Optics* 27(7), 1988; S. A. Twomey, C. F. Bohren, J. L. Mergenthaler, "Reflectance and albedo differences between wet and dry surfaces", *Applied Optics* 25(3), 1986; H. W. Jensen, J. Legakis, J. Dorsey, "Rendering of wet materials", *Rendering Techniques '99*, 1999. (Cited second-hand from Sawayama et al. Not read.)

**Foundation make-up**

- D. Lanza, J. R. Padrón-Griffe, A. Pranovich, A. Muñoz, J. R. Frisvad, A. Jarabo, "Practical Appearance Model for Foundation Cosmetics", *Computer Graphics Forum* 43(4), 2024 (EGSR). The layer model of [COSMETICS.md](COSMETICS.md), implemented from the text of the paper (its public code is not used). (Text read.)
- Y. Guo, M. Hašan, S. Zhao, "Position-free Monte Carlo simulation for arbitrary layered BSDFs", *ACM TOG* 37(6), 2018. (Cited second-hand from Lanza et al. Not read.)
- E. Heitz, J. Dupuy, C. Crassin, C. Dachsbacher, "The SGGX microflake distribution", *ACM TOG* 34(4), 2015. The distribution of the platelet normals and its sampling; implemented from memory and checked by numerical integration.
- K. Yoshida, N. Okiyama, "Estimation of reflectance, transmittance, and absorbance of cosmetic foundation layer on skin using translucency of skin", *Optics Express* 29(24), 40038-40050, 2021, doi:10.1364/OE.442219. Measures r_F, t_F, a_F of a layer on skin; the formula R = r_F + t_F² r_S/(1 − r_F r_S). (Full text read.)
- M. Doi, R. Ohtsuki, S. Tominaga, "Spectral Estimation of Skin Color with Foundation Makeup", *SCIA 2005*, LNCS 3540, 95-104. Kubelka-Munk with two layers; K and S measured from a thick and a thin layer on a black background. (Full text read.)
- C.-G. Huang, T.-S. Huang, W.-C. Lin, J.-H. Chuang, "Physically based cosmetic rendering", *Computer Animation and Virtual Worlds* 24(3-4), 275-283, 2013, doi:10.1002/cav.1523. (Full text read.)
- F. Ikkai, C. Blaksley, R. Hayashi, A. Jaspers, "Cosmetic effects of red-light reflection by ultramarine blue pigment", *Int. J. Cosmet. Sci.* 45, 413-423, 2023, doi:10.1111/ics.12849. How to measure the intrinsic color, artificial skins on black and white, the definition of hiding power. (Full text read.)
- P. Kubelka, "New Contributions to the Optics of Intensely Light-Scattering Materials. Part II: Nonhomogeneous Layers", *J. Opt. Soc. Am.* 44(4), 330-335, 1954. (Cited second-hand from Ikkai et al. Not read.)

**Camera, lenses and sensor**

- C. Kolb, D. Mitchell, P. Hanrahan, "A Realistic Camera Model for Computer Graphics", SIGGRAPH 1995, 317-324. Tracing rays through a lens prescription surface by surface ([CAMERA.md](CAMERA.md)).
- M. B. Hullin, E. Eisemann, H.-P. Seidel, S. Lee, "Physically-Based Real-Time Lens Flare Rendering", *ACM TOG* 30(4), Article 108, 2011 (SIGGRAPH), doi:10.1145/2010324.1965003. The ghosts traced from the lights as grids of rays whose cells become triangles on the film, with the stop's edges from the interpolated stop coordinates (`Render/LensFlare`; traced exactly here, without their polynomial fits).
- W. J. Smith, *Modern Lens Design: A Resource Manual*, McGraw-Hill, 1992 (the edition Kolb et al. cite; 2nd ed. 2005). *(which edition pbrt's files follow is not confirmed)* The source of the bundled lenses' radii, through pbrt's lens files.
- A. W. Tronnier, US 2,673,491 (double-Gauss); A. Nakamura, US 4,211,472, embodiment 4 (wide angle); R. Muller, US 4,647,161, embodiment 1 (fisheye). The glasses (n_d, Abbe number) of the bundled lenses, read from the patents.
- M. Taniyama (Fujinon), US 7,535,658, *Imaging lens*, 2009 (priority 2005), Example 1. The aspheric phone-camera lens `cellphone_5.6mm`, read from the patent's figures.
- Y. Chiang (Huawei Technologies), US 11,435,561, *Lens system and image photographing apparatus and device*, 2022 (priority 2017), embodiment 1. The phone telephoto `phone_tele_12mm`, from the patent's text.
- J. Jiang, D. Liu, J. Gu, S. Süsstrunk, "What is the Space of Spectral Sensitivity Functions for Digital Color Cameras?", WACV 2013, 168-179, doi:10.1109/WACV.2013.6475015. The spectral sensitivities of 28 cameras (`--camera-response`; CC BY-NC-SA, not redistributed).
- H. S. Malvar, L.-W. He, R. Cutler, "High-Quality Linear Interpolation for Demosaicing of Bayer-Patterned Color Images", ICASSP 2004, doi:10.1109/ICASSP.2004.1326587. The default demosaicing of the sensor model.
- ISO 12232:2019, *Photography - Digital still cameras - Determination of exposure index, ISO speed ratings, standard output sensitivity, and recommended exposure index*. The ISO speed of the physical exposure (the standard itself was not read).
- Adobe, *Digital Negative (DNG) Specification*. The raw file the sensor model writes, with its color matrices for two illuminants.

**3D Gaussian Splatting**

- J. Gao, C. Gu, Y. Lin, H. Zhu, X. Cao, L. Zhang, Y. Yao, "Relightable 3D Gaussian: Real-time Point Cloud Relighting with BRDF Decomposition and Ray Tracing", arXiv:2311.16043, 2023. *(studied)* The relightable reconstruction that the lighting study trains and relights.
- J.-N. Dihlmann, A. Majumdar, A. Engelhardt, R. Braun, H. P. A. Lensch, "Subsurface Scattering for 3D Gaussian Splatting", NeurIPS 2024 (arXiv:2408.12282). *(studied)* Reports that a relightable reconstruction bakes the subsurface scattering of translucent objects into the surface color.
- B. Kerbl, G. Kopanas, T. Leimkühler, G. Drettakis, "3D Gaussian Splatting for Real-Time Radiance Field Rendering", *ACM TOG* 42(4), 2023 (SIGGRAPH), doi:10.1145/3592433. The representation that the `point_cloud.ply` reader expects ([CAPTURE.md](CAPTURE.md)).
- Y. Yao, Z. Zeng, C. Gu, X. Zhu, L. Zhang, "Reflective Gaussian Splatting", arXiv:2412.19282, 2024 (v2). (Abstract read.) Related work of the lighting study: relightable reconstruction of reflective objects. Not run here.
- G. Kouros, M. Wu, T. Tuytelaars, "Spec-Gloss Surfels and Normal-Diffuse Priors for Relightable Glossy Objects", arXiv:2510.02069, 2025 (v3). (Abstract read.) Related work of the lighting study: a microfacet BRDF in 2D Gaussian splatting for glossy objects. Not run here.
- N. Moenne-Loccoz, A. Mirzaei, O. Perel, R. de Lutio, J. Martinez Esturo, G. State, S. Fidler, N. Sharp, Z. Gojcic, "3D Gaussian Ray Tracing: Fast Tracing of Particle Scenes", *ACM TOG* 43(6), 2024 (SIGGRAPH Asia), doi:10.1145/3687934. The ray-Gaussian response with which a splat cloud is ray-traced ([CAPTURE.md](CAPTURE.md)).

**Software and background**

- W. Jakob, S. Speierer, N. Roussel, M. Nimier-David, D. Vicini, T. Zeltner, B. Nicolet, M. Crespo, V. Leroy, Z. Zhang, *Mitsuba 3 renderer*, 2022, [mitsuba-renderer.org](https://www.mitsuba-renderer.org). The RGB-vs-spectral study of the graduation research was done with it, and the renderer is cross-checked against it ([MITSUBA_CHECK.md](MITSUBA_CHECK.md)).
- M. Hayford, *ray-optics* (rayoptics), an optical-design library in Python (BSD), version 0.9.5, [github.com/mjhoptics/ray-optics](https://github.com/mjhoptics/ray-optics). The independent trace of the aspheric phone lens (`tools/test/lens_rayoptics_check.py`).
- Blender 4.5 and its OpenColorIO configuration (views Standard, Filmic, AgX), and the ACES 1.3 configuration built into OpenColorIO 2.4. The view transforms of the gamut study ([METAL_ACCURACY.md](METAL_ACCURACY.md) 3.8).
- rawpy (LibRaw) and colour-demosaicing. The independent reading of the sensor model's DNG files (`tools/test/sensor_dng_check.py`).
