# References

The works that the code implements or follows, and the ones that the studies in `docs/` reproduce or test. Each entry says what it is used for.
Entries marked *(studied)* were implemented as a reference or evaluated, and are not part of the renderer.

**How far the bibliographic data were checked.** The entries that README and the documents on hair, skin, clouds, ocean, gloss and thin films carried on 2026-11-09 were checked by web search that day (authors, title, venue, volume, year, pages; 6 corrections were made). The entries of Haro et al., Graham et al., Barros and Walter were checked from search summaries only. The entries added later were taken from the title page of the PDF where one was at hand (Moon et al., Dihlmann et al., Lanza et al., Křivánek et al., Speierer et al., Yao et al., Kouros et al.) and otherwise written from memory; those are marked *(from memory)*. Not confirmed: the edition of Born & Wolf, the conference name (EGSR) of Hanika et al. 2015 (volume and issue are confirmed), and the source of the melanin absorption values (taken from PBRT). Remarks such as *(Text read.)* say whether the work was read in full or only cited second-hand.
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
- C. Schlick, "An Inexpensive BRDF Model for Physically-based Rendering", *Computer Graphics Forum* 13(3), 1994. The Fresnel of the coated diffuse material.
- A. D. Rakić, "Algorithm for the determination of intrinsic optical constants of metal films: application to aluminum", *Applied Optics* 34(22), 4755-4767, 1995; K. M. McPeak, S. V. Jayanti, S. J. P. Kress, S. Meyer, S. Iotti, A. Rossinelli, D. J. Norris,"Plasmonic Films Can Easily Be Better: Rules and Recipes", *ACS Photonics* 2(3), 326-333, 2015. *(from memory)* The `n(λ)`, `k(λ)` of aluminum (`assets/Metal_data/Rakic_Al.csv`, `McPeak_Al.csv`), for the comparison with the measured aluminum.

**Measured reflectance and polarization**

- S.-H. Baek, T. Zeltner, H. J. Ku, I. Hwang, X. Tong, W. Jakob, M. H. Kim, "Image-Based Acquisition and Modeling of Polarimetric Reflectance", *ACM TOG* 39(4), 2020 (SIGGRAPH). *(studied)* The KAIST polarimetric BRDF dataset: its gold sample against the rough conductor ([MEASURED_GOLD.md](MEASURED_GOLD.md)).
- Y. Moon, R. Maeda, S. Shin, I. Hwang, Y. Kim, M. H. Kim, S.-H. Baek, "Hyperspectral Polarimetric BRDFs of Real-world Materials", 2025 (preprint in ACM format, September 2025; the PDF names no venue), [project page](https://yunseong0518.github.io/projects/hpBRDF/). The hpBRDF tables: the `measured` material reads one (PCA-compressed, with its measured Mueller matrices), and its aluminum is compared with the rough conductor (MEASURED_GOLD.md, section 6).
- M. Born, E. Wolf, *Principles of Optics* (see "Dispersion and thin films"): also the Fresnel amplitudes and the Stokes and Mueller formalism of the polarization code ([POLARIZATION.md](POLARIZATION.md)).

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
- L. Belcour, P. Barla, "A Practical Extension to Microfacet Theory for the Modeling of Varying Iridescence", *ACM TOG* 36(4), Article 65, 2017 (SIGGRAPH). Iridescence on rough surfaces; background for the thin-film material, which is smooth.

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

**Caustics and light transport**

- H. W. Jensen, "Global Illumination using Photon Maps", Eurographics Workshop on Rendering, 1996. The caustic photon map.
- J. Hanika, M. Droske, L. Fascione, "Manifold Next Event Estimation", *Computer Graphics Forum* 34(4) (EGSR), 2015; T. Zeltner, I. Georgiev, W. Jakob, "Specular Manifold Sampling for Rendering High-Frequency Caustics and Glints", *ACM TOG* 39(4), 2020. The specular manifold NEE of the prism caustic.
- E. Veach, *Robust Monte Carlo Methods for Light Transport Simulation*, PhD thesis, Stanford University, 1997; E. P. Lafortune, Y. D. Willems, "Bi-directional Path Tracing", Compugraphics 1993. Light tracing and bidirectional path tracing, with specular vertices.
- I. Georgiev, J. Křivánek, T. Davidovič, P. Slusallek, "Light Transport Simulation with Vertex Connection and Merging", *ACM TOG* 31(6) (SIGGRAPH Asia), 2012;
  T. Davidovič, J. Křivánek, M. Hašan, P. Slusallek, "Progressive Light Transport Simulation on the GPU: Survey and Improvements", *ACM TOG* 33(3), 2014. The vertex cache and the weights of the bidirectional estimators. (Read in part: in VCM, vertex merging, its weights and the radius schedule; in Davidovič et al., the light vertex cache and VCM on it.)
- J. Novák, V. Havran, C. Dachsbacher, "Path Regeneration for Interactive Path Tracing", Eurographics 2010 short papers, 61–64. (Cited second-hand, from Davidovič et al.) Bidirectional path tracing on the GPU that connects each camera vertex to all vertices of a randomly chosen light path: the connection scheme of the prototype.
- shocker_0x15, "Renderer introduction: Ray Tracing Camp 10" (slides, in Japanese), [Speaker Deck](https://speakerdeck.com/shocker_0x15/renderer-introduction-ray-tracing-camp-10). Volumetric and dispersive caustics by light vertex cache bidirectional path tracing: the target of the fog scenes.
- J. Křivánek, I. Georgiev, T. Hachisuka, P. Vévoda, M. Šik, D. Nowrouzezahrai, W. Jarosz, "Unifying Points, Beams, and Paths in Volumetric Light Transport Simulation", *ACM TOG* 33(4), 2014 (SIGGRAPH), doi:10.1145/2601097.2601219. (Read in part: the estimators and the results.) Related work of the bidirectional estimator in fog ([LVC_BPT.md](LVC_BPT.md)): point and beam estimators combined with paths by multiple importance sampling. The prototype has no beams.
- S. Speierer, C. Hery, R. Villemin, W. Jakob, "Caustic Connection Strategies for Bidirectional Path Tracing", Pixar Technical Memo #18-01, 2018. (Text read, from a Wayback Machine copy: the Pixar page now redirects.) Connection strategies through chains of refractive events, found by manifold exploration, for caustics that classical bidirectional path tracing finds poorly: the manifold connections of [LVC_BPT.md](LVC_BPT.md) section 8c.

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

**3D Gaussian Splatting**

- J. Gao, C. Gu, Y. Lin, H. Zhu, X. Cao, L. Zhang, Y. Yao, "Relightable 3D Gaussian: Real-time Point Cloud Relighting with BRDF Decomposition and Ray Tracing", arXiv:2311.16043, 2023. *(studied)* The relightable reconstruction that the lighting study trains and relights.
- J.-N. Dihlmann, A. Majumdar, A. Engelhardt, R. Braun, H. P. A. Lensch, "Subsurface Scattering for 3D Gaussian Splatting", NeurIPS 2024 (arXiv:2408.12282). *(studied)* Reports that a relightable reconstruction bakes the subsurface scattering of translucent objects into the surface color.
- B. Kerbl, G. Kopanas, T. Leimkühler, G. Drettakis, "3D Gaussian Splatting for Real-Time Radiance Field Rendering", *ACM TOG* 42(4), Article 139, 2023 (SIGGRAPH). *(from memory)* The representation that the `point_cloud.ply` reader expects ([CAPTURE.md](CAPTURE.md)).
- Y. Yao, Z. Zeng, C. Gu, X. Zhu, L. Zhang, "Reflective Gaussian Splatting", arXiv:2412.19282, 2024 (v2). (Abstract read.) Related work of the lighting study: relightable reconstruction of reflective objects. Not run here.
- G. Kouros, M. Wu, T. Tuytelaars, "Spec-Gloss Surfels and Normal-Diffuse Priors for Relightable Glossy Objects", arXiv:2510.02069, 2025 (v3). (Abstract read.) Related work of the lighting study: a microfacet BRDF in 2D Gaussian splatting for glossy objects. Not run here.

**Software and background**

- W. Jakob, S. Speierer, N. Roussel, M. Nimier-David, D. Vicini, T. Zeltner, B. Nicolet, M. Crespo, V. Leroy, Z. Zhang, *Mitsuba 3 renderer*, 2022, [mitsuba-renderer.org](https://www.mitsuba-renderer.org). The RGB-vs-spectral study of the graduation research was done with it.
